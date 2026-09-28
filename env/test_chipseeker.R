#!/usr/bin/env Rscript
# ============================================================================
# test_chipseeker.R — ChIPseeker 环境安装自检（分阶段，定位失败点）
# ============================================================================
# 用法:
#   Rscript test_chipseeker.R                       # 用内置迷你数据自检
#   Rscript test_chipseeker.R genes.gff3 peaks.bed  # 用你的真实参考/peak
# 逐阶段测：① 载库 ② 从 GFF3 建 TxDb ③ annotatePeak ④ 出饼图
# 任一阶段失败会明确打印是哪一步，便于对症处理。
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
gff  <- if (length(args) >= 1) args[1] else ""
bed  <- if (length(args) >= 2) args[2] else ""
tmpd <- tempdir()
fail <- function(stage, e) {
    cat(sprintf("\n[FAIL @ %s]\n%s\n", stage, conditionMessage(e)))
    quit(status = 1)
}
say <- function(...) cat(sprintf(...), "\n")

# ── 阶段 1：载库 ──
say("[1/4] Load ChIPseeker and dependencies ...")
tryCatch({
    suppressMessages({
        library(ChIPseeker)
        library(GenomicFeatures)
        library(rtracklayer)
    })
    # makeTxDbFromGFF 在 Bioc>=3.19 迁到 txdbmaker
    makeTx <<- if (requireNamespace("txdbmaker", quietly = TRUE))
        txdbmaker::makeTxDbFromGFF else GenomicFeatures::makeTxDbFromGFF
    say("      OK — ChIPseeker %s", as.character(packageVersion("ChIPseeker")))
}, error = function(e) fail("package loading", e))

# ── 准备输入（无参数则用内置迷你数据）──
if (gff == "") {
    gff <- file.path(tmpd, "mini.gff3")
    writeLines(c(
        "##gff-version 3",
        "chr1\t.\tgene\t1000\t5000\t.\t+\t.\tID=gene1;Name=GENE1",
        "chr1\t.\tmRNA\t1000\t5000\t.\t+\t.\tID=mrna1;Parent=gene1",
        "chr1\t.\texon\t1000\t2000\t.\t+\t.\tID=ex1;Parent=mrna1",
        "chr1\t.\texon\t4000\t5000\t.\t+\t.\tID=ex2;Parent=mrna1",
        "chr1\t.\tCDS\t1000\t2000\t.\t+\t0\tID=cds1;Parent=mrna1",
        "chr1\t.\tCDS\t4000\t5000\t.\t+\t0\tID=cds2;Parent=mrna1",
        "chr1\t.\tgene\t8000\t9000\t.\t-\t.\tID=gene2;Name=GENE2",
        "chr1\t.\tmRNA\t8000\t9000\t.\t-\t.\tID=mrna2;Parent=gene2",
        "chr1\t.\texon\t8000\t9000\t.\t-\t.\tID=ex3;Parent=mrna2"
    ), gff)
    say("      (using the bundled tiny GFF3)")
}
if (bed == "") {
    bed <- file.path(tmpd, "mini_peaks.bed")
    writeLines(c(
        "chr1\t900\t950\tpeak1\t100",
        "chr1\t3000\t3050\tpeak2\t80",
        "chr1\t9100\t9200\tpeak3\t60"
    ), bed)
}

# ── 阶段 2：从 GFF3 建 TxDb ──（Chlamydomonas 等非模式物种最容易出问题的一步）
say("[2/4] Build TxDb from GFF3: %s", gff)
txdb <- tryCatch(
    suppressWarnings(makeTx(gff, format = "gff3")),
    error = function(e) fail("Build TxDb (makeTxDbFromGFF)", e))
say("      OK — %d transcripts", length(transcripts(txdb)))

# ── 阶段 3：annotatePeak ──（会触发 GenomeInfoDb 数据，genomeinfodbdata 坏则在此报错）
say("[3/4] annotatePeak ...")
anno <- tryCatch({
    peak <- readPeakFile(bed)
    annotatePeak(peak, TxDb = txdb, tssRegion = c(-1000, 200), verbose = FALSE)
}, error = function(e) fail("annotatePeak", e))
say("      OK — Annotated %d peaks", length(anno@anno))

# ── 阶段 4：出饼图（验证 R 绘图链）──
say("[4/4] Render plotAnnoPie ...")
out <- file.path(getwd(), "chipseeker_selftest_pie.png")
tryCatch({
    png(out, width = 900, height = 700, res = 130)
    plotAnnoPie(anno)
    dev.off()
}, error = function(e) fail("plotAnnoPie", e))
say("      OK — Wrote %s", out)

cat("\n[ALL PASS] ChIPseeker TxDb construction, annotation and plotting passed\n")
