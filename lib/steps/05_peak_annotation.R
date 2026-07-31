#!/usr/bin/env Rscript
# ============================================================================
# 05_peak_annotation.R — ChIPseeker 注释 + 文献式图 + GO/KEGG
# ============================================================================
# 由 05_peak_annotation.sh 包装调用。参数（位置）：
#   1 GFF3   2 peaks.bed   3 outdir   4 promoter_up   5 promoter_down
#   6 annotation_info(可空)   7 target_dist_kb   8 pipeline   9 kegg_org(可空)
# 输出契约（供 06/07/08/09 复用，保持与旧 Python 版一致）：
#   genes.bed  peak_annotation_full.tsv  promoter_target_genes.txt
#   target_genes_within_<N>kb_ranked.tsv  all_nearest_genes.txt
#   background_universe.txt  gene2go.tsv
#   + 图: anno_pie.png  dist_to_tss.png  anno_upset.png  go_dotplot.png  kegg_dotplot.png
# ============================================================================

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

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
say("== ChIPseeker 注释 (%s) ==", pipeline)

suppressMessages({
    library(ChIPseeker); library(GenomicFeatures); library(rtracklayer)
})
makeTx <- if (requireNamespace("txdbmaker", quietly = TRUE))
    txdbmaker::makeTxDbFromGFF else GenomicFeatures::makeTxDbFromGFF

# ── TxDb ──
say("建 TxDb: %s", gff)
txdb <- suppressWarnings(makeTx(gff, format = "gff3"))

# ── 读 peaks；对齐 seqlevels，避免良性但刷屏的 seqinfo 警告 ──
peak <- readPeakFile(bed)
common_lv <- intersect(seqlevels(peak), seqlevels(txdb))
if (length(common_lv) == 0) stop("peak 与 TxDb 无共同染色体名，请检查命名一致性")
dropped <- setdiff(seqlevels(peak), common_lv)
if (length(dropped) > 0)
    say("  注意: %d 个 peak 所在 seqlevel 不在注释中，将忽略: %s",
        sum(as.character(seqnames(peak)) %in% dropped), paste(head(dropped, 8), collapse = ","))
peak <- keepSeqlevels(peak, common_lv, pruning.mode = "coarse")

# ── 注释 ──
anno <- suppressWarnings(
    annotatePeak(peak, TxDb = txdb, tssRegion = c(-prom_up, prom_down),
                 level = "gene", verbose = FALSE))
df <- as.data.frame(anno)
say("  注释 %d 个 peak", nrow(df))

# peak id 列（readPeakFile 把 BED 第4列存为 V4）
df$peak_id <- if ("V4" %in% names(df)) as.character(df$V4) else paste0("peak_", seq_len(nrow(df)))
df$peak_score <- if ("V5" %in% names(df)) suppressWarnings(as.numeric(df$V5)) else NA_real_
# globalIDR（03 输出的第7列，若存在）— 供反查报告可重复性
df$peak_globalIDR <- if ("V7" %in% names(df)) suppressWarnings(as.numeric(df$V7)) else NA_real_

# 完整注释表
write.table(df, file.path(outdir, "peak_annotation_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 基因组分布统计（细分类别，文献饼图口径）──
simplify_cat <- function(x) {
    x <- as.character(x)
    ifelse(grepl("^Promoter", x), "Promoter",
    ifelse(grepl("5' UTR", x), "5' UTR",
    ifelse(grepl("3' UTR", x), "3' UTR",
    ifelse(grepl("^Exon", x), "Exon",
    ifelse(grepl("^Intron", x), "Intron",
    ifelse(grepl("^Downstream", x), "Downstream", "Distal Intergenic"))))))
}
df$feature <- simplify_cat(df$annotation)
tab <- as.data.frame(table(df$feature))
names(tab) <- c("feature", "peaks"); tab$pct <- round(100 * tab$peaks / sum(tab$peaks), 2)
write.table(tab, file.path(outdir, "genomic_distribution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
say("  分布: %s", paste(sprintf("%s=%d(%.1f%%)", tab$feature, tab$peaks, tab$pct), collapse = "  "))

# ── 图: 饼图 / 距 TSS(条形) / upset ──
try({ png(file.path(outdir, "anno_pie.png"), 1000, 800, res = 130); plotAnnoPie(anno); dev.off() }, silent = TRUE)
try({ png(file.path(outdir, "dist_to_tss.png"), 1100, 500, res = 130); print(plotDistToTSS(anno)); dev.off() }, silent = TRUE)
try({ png(file.path(outdir, "anno_upset.png"), 1200, 700, res = 130); upsetplot(anno); dev.off() }, silent = TRUE)

# ── 折线形 metaplot: peak 在 TSS 附近的密度曲线（R 原生, 近似文献折线图）──
try({
    prom <- suppressWarnings(ChIPseeker::getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000))
    tm <- suppressWarnings(ChIPseeker::getTagMatrix(peak, windows = prom))
    png(file.path(outdir, "tss_avgprof.png"), 1100, 550, res = 130)
    print(ChIPseeker::plotAvgProf(tm, xlim = c(-3000, 3000),
          xlab = "Distance to TSS (bp)", ylab = "Peak frequency"))
    dev.off()
}, silent = TRUE)

# ── genes.bed（BED6, 供 06/07 复用）──
g <- genes(txdb)
gb <- data.frame(chr = as.character(seqnames(g)), start = start(g) - 1, end = end(g),
                 gid = names(g), score = 0, strand = as.character(strand(g)))
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
    say("  靶基因: promoter=%d  ≤%gkb=%d", length(prom_genes), dist_kb, nrow(best))
}

# ── gene2go（供 08 / GO 富集）──
t2g <- NULL
if (nzchar(anno_info) && file.exists(anno_info)) {
    say("解析 annotation_info 抽取 GO: %s", basename(anno_info))
    ai <- tryCatch(read.delim(anno_info, stringsAsFactors = FALSE, check.names = FALSE),
                   error = function(e) NULL)
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
                say("  gene2go: %d 条 (%d 基因)", nrow(pairs), length(unique(pairs$gene)))
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
            ggsave(file.path(outdir, "go_dotplot.png"),
                   dotplot(ego, showCategory = 20) + ggtitle("GO enrichment"),
                   width = 8, height = 7, dpi = 150)
            say("  GO: %d 显著 term → go_dotplot.png", nrow(as.data.frame(ego)))
        } else say("  GO: 无显著 term")
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
            ggsave(file.path(outdir, "kegg_dotplot.png"),
                   dotplot(kk, showCategory = 20) + ggtitle("KEGG enrichment"),
                   width = 8, height = 7, dpi = 150)
            say("  KEGG: %d 通路 → kegg_dotplot.png", nrow(as.data.frame(kk)))
        } else say("  KEGG: 无显著通路（或 ID 未匹配 KEGG）")
    }, silent = TRUE)
}

cat("\n[DONE] ChIPseeker 注释完成 →", outdir, "\n")
