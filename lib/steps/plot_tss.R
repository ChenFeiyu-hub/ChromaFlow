#!/usr/bin/env Rscript
# Base R only: can also redraw an existing annotation table without rebuilding TxDb.
tss_distance_summary <- function(distance) {
    x <- suppressWarnings(as.numeric(distance))
    valid <- x[is.finite(x)]
    breaks <- c(-Inf, 1000, 3000, 5000, 10000, 100000, Inf)
    band <- as.integer(cut(abs(valid), breaks, right = TRUE))
    index <- ifelse(valid < 0, 7L - band, ifelse(valid > 0, 7L + band, 7L))
    labels <- c("Upstream >100 kb", "Upstream 10-100 kb", "Upstream 5-10 kb",
                "Upstream 3-5 kb", "Upstream 1-3 kb", "Upstream 0-1 kb",
                "At TSS (0 bp)", "Downstream 0-1 kb", "Downstream 1-3 kb",
                "Downstream 3-5 kb", "Downstream 5-10 kb", "Downstream 10-100 kb",
                "Downstream >100 kb")
    counts <- tabulate(as.integer(index), nbins = 13L)
    data.frame(distance_bin = labels,
               direction = c(rep("Upstream", 6), "At TSS", rep("Downstream", 6)),
               count = counts,
               min_abs_bp = c(100001,10001,5001,3001,1001,1,0,1,1001,3001,5001,10001,100001),
               max_abs_bp = c(Inf,100000,10000,5000,3000,1000,0,1000,3000,5000,10000,100000,Inf),
               percent = if (length(valid)) counts / length(valid) * 100 else rep(0, 13),
               valid_peaks = length(valid), excluded_missing_distance = length(x) - length(valid),
               stringsAsFactors = FALSE)
}

plot_tss_distance <- function(distance, outdir) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    tab <- tss_distance_summary(distance)
    write.table(tab, file.path(outdir, "dist_to_tss.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    png(file.path(outdir, "dist_to_tss.png"), width = 1500, height = 950, res = 150)
    on.exit(dev.off(), add = TRUE)
    par(mar = c(4.8, 12.8, 4.6, 2.5), mgp = c(2.8, 0.8, 0), las = 1)
    palette <- c("Upstream" = "#347CA5", "At TSS" = "#777777", "Downstream" = "#CF8540")
    order <- rev(seq_len(nrow(tab)))
    xmax <- max(10, tab$percent) * 1.28
    bars <- barplot(tab$percent[order], names.arg = tab$distance_bin[order], horiz = TRUE,
                    col = palette[tab$direction[order]], border = NA,
                    xlim = c(0, xmax), xlab = "Peaks with a valid TSS distance (%)",
                    cex.names = 0.9, axes = FALSE)
    axis(1)
    abline(v = axTicks(1), col = "#EEEEEE", lwd = 0.7)
    # Redraw bars above the grid, using the identical geometry.
    barplot(tab$percent[order], horiz = TRUE, col = palette[tab$direction[order]],
            border = NA, add = TRUE, axes = FALSE, axisnames = FALSE)
    text(tab$percent[order] + xmax * 0.015, bars,
         sprintf("%.1f%% (n=%s)", tab$percent[order], format(tab$count[order], trim = TRUE)),
         adj = 0, cex = 0.8)
    title("Peak distance to annotated TSS", adj = 0, line = 2.4, cex.main = 1.15)
    mtext(sprintf("Strand-relative direction | valid peaks: %d | missing distance: %d",
                  tab$valid_peaks[1], tab$excluded_missing_distance[1]),
          side = 3, line = 0.8, adj = 0, cex = 0.78)
    invisible(tab)
}

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2L) stop("Usage: Rscript plot_tss.R peak_annotation_full.tsv output_directory")
    data <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
    if (!"distanceToTSS" %in% names(data)) stop("Annotation table has no distanceToTSS column")
    plot_tss_distance(data$distanceToTSS, args[2])
    cat("[OK] Wrote dist_to_tss.png and dist_to_tss.tsv\n")
}
