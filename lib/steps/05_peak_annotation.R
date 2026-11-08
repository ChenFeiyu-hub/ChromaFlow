#!/usr/bin/env Rscript
# ============================================================================
# 05_peak_annotation.R — ChIPseeker 注释 + 文献式图 + GO/KEGG
# ============================================================================
# 由 05_peak_annotation.sh 包装调用。参数（位置）：
#   1 GFF3   2 peaks.bed   3 outdir   4 promoter_up   5 promoter_down
#   6 annotation_info(可空)   7 target_dist_kb   8 pipeline   9 kegg_org(可空)   10 make_plots
# 输出契约（供 06/07/08/09 复用，保持与旧 Python 版一致）：
#   genes.bed  peak_annotation_full.tsv  promoter_target_genes.txt
#   target_genes_within_<N>kb_ranked.tsv  all_nearest_genes.txt
#   background_universe.txt  gene2go.tsv
#   + 图: anno_pie.png  dist_to_tss.png  anno_upset.png  go_dotplot.png  kegg_dotplot.png
# ============================================================================

# 全局禁用科学计数：基因组坐标一旦被写成 "8e+06"，下游按 BED/TSV 解析必然崩。
# R 的 write.table 对 double 用 as.character()，而 start()-1 之类运算会把
# integer 提升成 double，所以这是必须的（scipen 越大越倾向定点表示）。
options(scipen = 999)

a <- commandArgs(trailingOnly = TRUE)
gff        <- a[1]
bed        <- a[2]
outdir     <- a[3]
prom_up    <- as.integer(a[4])
prom_down  <- as.integer(a[5])
anno_info  <- if (length(a) >= 6) a[6] else ""
dist_kb    <- if (length(a) >= 7) as.numeric(a[7]) else 3
pipeline   <- if (length(a) >= 8) a[8] else "chipseq"
kegg_org   <- if (length(a) >= 9) a[9] else ""
plot_flag  <- if (length(a) >= 10) tolower(a[10]) else "true"
if (!plot_flag %in% c("true", "false")) stop("MAKE_PLOTS must be true or false")
make_plots <- plot_flag == "true"
# An optional existing table requests report-only rendering into a separate
# directory. Exact table equality is required before producing any reports.
report_table <- if (length(a) >= 11) a[11] else ""
chrom_sizes <- if (length(a) >= 12) a[12] else ""
peak_kind <- if (length(a) >= 13) a[13] else "unknown"
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(dirname(script_file), "tss_windows.R"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
say("== ChIPseeker annotation (%s) ==", pipeline)

suppressMessages({
    library(ChIPseeker); library(GenomicFeatures); library(rtracklayer)
})
makeTx <- if (requireNamespace("txdbmaker", quietly = TRUE))
    txdbmaker::makeTxDbFromGFF else GenomicFeatures::makeTxDbFromGFF
say("ChIPseeker version: %s; GenomicFeatures version: %s",
    as.character(packageVersion("ChIPseeker")), as.character(packageVersion("GenomicFeatures")))

# ── TxDb ──
say("Build TxDb: %s", gff)
# Keep construction warnings visible: they may identify discarded gene models.
txdb <- makeTx(gff, format = "gff3")

# ── 读 peaks；对齐 seqlevels，避免良性但刷屏的 seqinfo 警告 ──
peak <- readPeakFile(bed)
common_lv <- intersect(seqlevels(peak), seqlevels(txdb))
if (length(common_lv) == 0) stop("Peaks and TxDb share no sequence names; check reference naming consistency")
dropped <- setdiff(seqlevels(peak), common_lv)
if (length(dropped) > 0)
    say("  NOTE: %d peaks on %d sequence levels absent from annotation will be dropped: %s",
        sum(as.character(seqnames(peak)) %in% dropped), length(dropped), paste(head(dropped, 8), collapse = ","))
peak <- keepSeqlevels(peak, common_lv, pruning.mode = "coarse")

# ── 注释 ──
anno <- annotatePeak(peak, TxDb = txdb, tssRegion = c(-prom_up, prom_down),
                     level = "gene", overlap = "TSS", sameStrand = FALSE,
                     ignoreOverlap = FALSE, ignoreUpstream = FALSE, ignoreDownstream = FALSE,
                     genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"),
                     verbose = FALSE)
df <- as.data.frame(anno)
say("  Annotated %d peaks", nrow(df))

# peak id 列（readPeakFile 把 BED 第4列存为 V4）
df$peak_id <- if ("V4" %in% names(df)) as.character(df$V4) else paste0("peak_", seq_len(nrow(df)))
df$peak_score <- if ("V5" %in% names(df)) suppressWarnings(as.numeric(df$V5)) else NA_real_
# Column 7 is globalIDR only for the explicitly selected ChromaFlow IDR BED8.
# narrowPeak column 7 is MACS signalValue; custom files have unknown semantics.
df$peak_globalIDR <- peak_idr_values(df, peak_kind)

# ── 关键区分：annotation 里的基因 ≠ geneId ──────────────────────────────
# ChIPseeker 的两列回答的是两个不同问题：
#   annotation —— peak 物理上落在【哪个基因的哪个结构】里（可能是邻近基因）
#   geneId     —— 距离 peak 最近的【TSS 属于哪个基因】
# 典型情形：peak 位于 A 基因最后一个内含子内，最近 TSS 属于 B，但不在配置的启动子窗内，
# 于是 annotation="Intron (A...)" 而 geneId=B。这不是错误，但极易被误读，
# 且下游若「按 annotation 筛位置、按 geneId 记基因」就会张冠李戴。
# 这里把 annotation 里的基因单独解析成一列，并标记两者是否一致。
simplify_cat <- function(x) {
    x <- as.character(x)
    ifelse(grepl("^Promoter", x), "Promoter",
    ifelse(grepl("5' UTR", x), "5' UTR",
    ifelse(grepl("3' UTR", x), "3' UTR",
    ifelse(grepl("^Exon", x), "Exon",
    ifelse(grepl("^Intron", x), "Intron",
    ifelse(grepl("^Downstream", x), "Downstream", "Distal Intergenic"))))))
}
extract_feature_gene <- function(x) {
    x <- as.character(x); g <- rep(NA_character_, length(x))
    hit <- grepl("\\(.*/.*,", x)                 # 形如 "Intron (T1/GENE, intron 8 of 8)"
    g[hit] <- sub("^.*\\([^/]*/([^,]+),.*$", "\\1", x[hit])
    g
}
df$feature_simple  <- simplify_cat(df$annotation)
df$feature_geneId  <- extract_feature_gene(df$annotation)
df$gene_consistent <- is.na(df$feature_geneId) | df$feature_geneId == df$geneId

n_mis <- sum(!df$gene_consistent, na.rm = TRUE)
if (n_mis > 0) {
    say("  NOTE: %d/%d peaks overlap a feature from a gene other than the nearest-TSS gene", n_mis, nrow(df))
    say("        (The peak overlaps one gene body while another gene has a closer TSS.)")
    say("        See feature_geneId / gene_consistent in the annotation table")
    say("        Review these loci before changing PROMOTER_UP (currently %d bp)", prom_up)
}

# 完整注释表
verify_report_annotation <- function(df, saved, scratch_dir, up, down, policy_script, python = "python3") {
    reconstructed <- tempfile("reconstructed-annotation-", tmpdir = scratch_dir, fileext = ".tsv")
    on.exit(unlink(reconstructed), add = TRUE)
    write.table(df, reconstructed, sep = "\t", quote = FALSE, row.names = FALSE)
    status <- system2(python, c(shQuote(policy_script), "compare", shQuote(saved),
                     "--reconstructed", shQuote(reconstructed), "--up", as.character(up),
                     "--down", as.character(down)))
    if (status != 0L) stop("Report refresh refused: annotation comparison failed; see the detailed differences above")
    invisible(TRUE)
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
if (nzchar(report_table)) {
    verify_report_annotation(df, report_table, outdir, prom_up, prom_down,
                             file.path(dirname(script_file), "support", "annotation_policy.py"))
} else {
    write.table(df, file.path(outdir, "peak_annotation_full.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

# ── 基因组分布统计（细分类别，文献饼图口径）──
df$feature <- simplify_cat(df$annotation)
tab <- as.data.frame(table(df$feature))
names(tab) <- c("feature", "peaks"); tab$pct <- round(100 * tab$peaks / sum(tab$peaks), 2)
write.table(tab, file.path(outdir, "genomic_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
say("  Distribution: %s", paste(sprintf("%s=%d(%.1f%%)", tab$feature, tab$peaks, tab$pct), collapse = "  "))

# ── 图: 饼图 / 距 TSS(条形) / upset ──
source(file.path(dirname(script_file), "plot_tss.R"))
source(file.path(dirname(script_file), "annotation_plots.R"))
report <- new_plot_report(outdir)
if (make_plots) {
    report$attempt("anno_pie.png", plotAnnoPie(anno), 1000, 800)
    report$attempt("dist_to_tss.png", {
        if (!any(is.finite(df$distanceToTSS))) skip_plot("No finite TSS distances")
        plot_tss_distance(df$distanceToTSS, outdir)
    })
    report$attempt("anno_upset.png", plot_annotation_intersections(anno@detailGenomicAnnotation, outdir), 1400, 800)
    report$attempt("tss_avgprof.png", {
        prom <- ChIPseeker::getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
        tm <- reference_tss_matrix(peak, prom, chrom_sizes, outdir)
        plot_tss_profile(tm, outdir)
    }, 1200, 600)
} else {
    say("  MAKE_PLOTS=false: annotation plots are disabled; tables and audits remain enabled")
}
plot_status <- report$finish()
if (nzchar(report_table)) {
    say("[REPORT] Report-only rendering completed; scientific tables were not rewritten")
    quit(status = if (any(plot_status$status == "failed")) 1L else 0L)
}

# ── genes.bed（BED6, 供 06/07 复用）──
g <- genes(txdb)
# start(g) 是 integer，但 start(g)-1 会提升为 double，而 write.table 对 double
# 用 as.character()，8000000 这类整数会被写成 "8e+06" —— 下游按 BED 解析就崩。
# 所以显式转回 integer（顶部的 options(scipen=999) 是第二道保险）。
gb <- data.frame(chr = as.character(seqnames(g)),
                 start = as.integer(start(g) - 1L), end = as.integer(end(g)),
                 gid = names(g), score = 0L, strand = as.character(strand(g)))
write.table(gb[order(gb$chr, gb$start), ], file.path(outdir, "genes.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(sort(unique(gb$gid)), file.path(outdir, "background_universe.txt"))

# ── 靶基因三档（保持契约）──
strip <- function(x) sub("\\.[0-9]+$", "", x)
prom_genes <- unique(df$geneId[grepl("^Promoter", df$annotation) & !is.na(df$geneId)])
writeLines(sort(prom_genes), file.path(outdir, "promoter_target_genes.txt"))
writeLines(sort(unique(df$geneId[!is.na(df$geneId)])), file.path(outdir, "all_nearest_genes.txt"))

# 距离过滤 + 结合强度排序
sel <- !is.na(df$geneId) & abs(df$distanceToTSS) <= dist_kb * 1000
rk <- df[sel, c("geneId", "peak_score", "distanceToTSS")]
if (nrow(rk) > 0) {
    rk$peak_score[is.na(rk$peak_score)] <- 0
    # 每基因取结合最强的 peak
    best <- do.call(rbind, lapply(split(rk, rk$geneId), function(s) {
        i <- which.max(s$peak_score); s[i, ]
    }))
    best <- best[order(-best$peak_score), ]
    names(best) <- c("gene_id", "peak_score", "dist_to_TSS")
    write.table(best, file.path(outdir, sprintf("target_genes_within_%gkb_ranked.tsv", dist_kb)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(sort(unique(best$gene_id)),
               file.path(outdir, sprintf("target_genes_within_%gkb.txt", dist_kb)))
    say("  Target genes: promoter=%d  ≤%gkb=%d", length(prom_genes), dist_kb, nrow(best))
}

# ── gene2go（供 08 / GO 富集）──
t2g <- NULL
if (nzchar(anno_info) && file.exists(anno_info)) {
    say("Extract GO mappings from annotation_info: %s", basename(anno_info))
    ai <- tryCatch(read.delim(anno_info, stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) { warning(paste("Could not read annotation_info:", conditionMessage(e))); NULL })
    if (!is.null(ai)) {
        cn <- tolower(names(ai))
        gi <- which(cn %in% c("locusname", "gene", "genename", "locus"))[1]
        ti <- which(cn %in% c("transcriptname", "transcript", "pacid"))[1]
        goi <- which(cn == "go")[1]
        if (!is.na(goi)) {
            key <- if (!is.na(gi)) ai[[gi]] else ai[[ti]]
            pairs <- do.call(rbind, lapply(seq_len(nrow(ai)), function(i) {
                gos <- regmatches(ai[[goi]][i], gregexpr("GO:[0-9]+", ai[[goi]][i]))[[1]]
                if (length(gos) == 0) return(NULL)
                data.frame(GO = gos, gene = strip(key[i]), stringsAsFactors = FALSE)
            }))
            if (!is.null(pairs)) {
                pairs <- unique(pairs)
                write.table(pairs[, c("gene", "GO")], file.path(outdir, "gene2go.tsv"),
                            sep = "\t", quote = FALSE, row.names = FALSE)
                t2g <- pairs[, c("GO", "gene")]
                say("  gene2go: %d mappings (%d genes)", nrow(pairs), length(unique(pairs$gene)))
            }
        }
    }
}

# ── GO 富集 dotplot（clusterProfiler enricher + 自定义 TERM2GENE）──
if (!is.null(t2g) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    try({
        suppressMessages(library(clusterProfiler))
        t2n <- NULL
        if (requireNamespace("GO.db", quietly = TRUE)) {
            suppressMessages(library(GO.db))
            sel_go <- unique(t2g$GO)
            m <- suppressWarnings(AnnotationDbi::select(GO.db, keys = sel_go,
                    columns = "TERM", keytype = "GOID"))
            t2n <- data.frame(GO = m$GOID, name = m$TERM, stringsAsFactors = FALSE)
        }
        ego <- enricher(gene = strip(prom_genes), TERM2GENE = t2g, TERM2NAME = t2n,
                        universe = strip(unique(gb$gid)), pvalueCutoff = 0.05)
        if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
            write.table(as.data.frame(ego), file.path(outdir, "go_enrichment.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            suppressMessages(library(ggplot2))
            if (make_plots) ggsave(file.path(outdir, "go_dotplot.png"),
                   dotplot(ego, showCategory = 20) + ggtitle("GO enrichment"),
                   width = 8, height = 7, dpi = 150)
            say("  GO: %d significant terms → go_dotplot.png", nrow(as.data.frame(ego)))
        } else say("  GO: no significant terms")
    }, silent = FALSE)
}

# ── KEGG dotplot（best-effort, 需联网 + KEGG 支持该物种）──
if (nzchar(kegg_org) && requireNamespace("clusterProfiler", quietly = TRUE)) {
    try({
        suppressMessages(library(clusterProfiler))
        kk <- enrichKEGG(gene = strip(prom_genes), organism = kegg_org, pvalueCutoff = 0.05)
        if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
            write.table(as.data.frame(kk), file.path(outdir, "kegg_enrichment.tsv"),
                        sep = "\t", quote = FALSE, row.names = FALSE)
            suppressMessages(library(ggplot2))
            if (make_plots) ggsave(file.path(outdir, "kegg_dotplot.png"),
                   dotplot(kk, showCategory = 20) + ggtitle("KEGG enrichment"),
                   width = 8, height = 7, dpi = 150)
            say("  KEGG: %d pathways → kegg_dotplot.png", nrow(as.data.frame(kk)))
        } else say("  KEGG: No significant pathways (or IDs did not match KEGG)")
    }, silent = FALSE)
}

cat("\n[DONE] ChIPseeker annotation completed →", outdir, "\n")
