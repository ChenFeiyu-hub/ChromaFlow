# Dependency-free annotation diagnostics. Assignment policy is never changed here.
skip_plot <- function(reason) {
    stop(structure(list(message = reason, call = NULL), class = c("plot_skip", "error", "condition")))
}

new_plot_report <- function(outdir) {
    rows <- list()
    attempt <- function(name, draw, width = NULL, height = NULL) {
        opened <- FALSE
        temporary <- tempfile(pattern = "plot-", tmpdir = outdir, fileext = ".png")
        status <- "written"; reason <- ""
        tryCatch({
            if (!is.null(width)) {
                png(temporary, width = width, height = height, res = 130)
                opened <- TRUE
            }
            value <- force(draw)
            if (inherits(value, "ggplot")) print(value)
            if (opened) { dev.off(); opened <- FALSE }
            if (!is.null(width) && !file.copy(temporary, file.path(outdir, name), overwrite = TRUE))
                stop("Could not publish plot")
        }, plot_skip = function(e) {
            status <<- "skipped"; reason <<- conditionMessage(e)
        }, error = function(e) {
            status <<- "failed"; reason <<- conditionMessage(e)
        }, finally = {
            if (opened) dev.off()
            unlink(temporary)
        })
        if (status != "written") {
            unlink(file.path(outdir, name))
            tables <- c("anno_upset.png" = "annotation_intersections.tsv", "tss_avgprof.png" = "tss_avgprof.tsv",
                        "dist_to_tss.png" = "dist_to_tss.tsv")
            if (name %in% names(tables)) unlink(file.path(outdir, tables[[name]]))
        }
        reason <- gsub("[\t\r\n]+", " ", reason)
        rows[[length(rows) + 1L]] <<- data.frame(plot = name, status = status, reason = reason)
        cat(sprintf("[PLOT] %s: %s%s\n", name, status, if (nzchar(reason)) paste0("; ", reason) else ""))
        invisible(status)
    }
    finish <- function() {
        tab <- if (length(rows)) do.call(rbind, rows) else
            data.frame(plot = character(), status = character(), reason = character())
        write.table(tab, file.path(outdir, "plot_status.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
        invisible(tab)
    }
    list(attempt = attempt, finish = finish)
}

annotation_intersections <- function(detail) {
    m <- as.matrix(detail)
    if (!nrow(m) || !ncol(m)) skip_plot("No detailed annotation memberships")
    if (!is.logical(m) && (!is.numeric(m) || any(!m %in% c(0, 1))))
        stop("Detailed annotation memberships must be logical or 0/1")
    if (anyNA(m) || is.null(colnames(m)) || anyDuplicated(colnames(m)))
        stop("Missing or ambiguous detailed annotation memberships")
    m <- m != 0
    labels <- apply(m, 1, function(row) {
        if (any(row)) paste(colnames(m)[row], collapse = " & ") else "Unassigned"
    })
    counts <- table(labels)
    tab <- data.frame(intersection = names(counts), peaks = as.integer(counts), stringsAsFactors = FALSE)
    tab <- tab[order(-tab$peaks, tab$intersection), , drop = FALSE]
    rownames(tab) <- NULL
    membership <- m[match(tab$intersection, labels), , drop = FALSE]
    list(table = tab, membership = membership, total_peaks = nrow(m))
}

plot_annotation_intersections <- function(detail, outdir, top = 20L) {
    data <- annotation_intersections(detail)
    write.table(data$table, file.path(outdir, "annotation_intersections.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    n <- min(nrow(data$table), top)
    tab <- head(data$table, n)
    m <- data$membership[seq_len(n), , drop = FALSE]
    oldpar <- par(no.readonly = TRUE)
    on.exit({ layout(1); par(oldpar) }, add = TRUE)
    layout(matrix(c(1, 2), ncol = 1), heights = c(2, 1.5))
    par(mar = c(0.5, 11, 4, 2), mgp = c(2.8, 0.7, 0), las = 1)
    limits <- c(0, n * 1.2 + 0.2)
    xs <- barplot(tab$peaks, col = "#347CA5", border = NA, xlim = limits, xaxs = "i",
                  ylim = c(0, max(tab$peaks) * 1.2), ylab = "Peaks", axes = FALSE)
    axis(2); text(xs, tab$peaks, labels = tab$peaks, pos = 3, cex = 0.8)
    title("Overlapping genomic annotations", adj = 0)
    mtext(sprintf("%d peaks | %d of %d intersections shown | memberships may overlap",
                  data$total_peaks, n, nrow(data$table)), side = 3, line = 0.4, cex = 0.75, adj = 0)
    par(mar = c(2, 11, 0.5, 2), mgp = c(2.8, 0.7, 0), las = 1)
    ys <- rev(seq_len(ncol(m)))
    plot(NA, xlim = limits, ylim = c(0.5, ncol(m) + 0.5),
         xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    axis(2, at = ys, labels = gsub("_", " ", colnames(m)), tick = FALSE, cex.axis = 0.85)
    for (i in seq_len(n)) {
        points(rep(xs[i], length(ys)), ys, pch = 16, col = "#E1E5E8", cex = 0.9)
        present <- ys[m[i, ]]
        if (length(present) > 1) segments(xs[i], min(present), xs[i], max(present), col = "#347CA5", lwd = 2)
        if (length(present)) points(rep(xs[i], length(present)), present, pch = 16, col = "#347CA5", cex = 1.1)
    }
    mtext("Each column is one exact intersection; blank columns are unassigned", side = 1, line = 0.5, cex = 0.75)
    invisible(data$table)
}

tss_profile <- function(tm, upstream = 3000L, downstream = 3000L) {
    if (!is.matrix(tm) || !is.numeric(tm) || ncol(tm) != upstream + downstream + 1L)
        stop("Unexpected TSS tag-matrix dimensions or type")
    if (!nrow(tm)) skip_plot("No TSS windows with retained peak overlap; no curve generated")
    if (any(!is.finite(tm)) || any(tm < 0))
        stop("Non-finite or negative TSS tag-matrix values; no curve generated")
    if (!any(tm > 0)) skip_plot("All TSS window overlaps are zero; no curve generated")
    data.frame(distance_bp = seq.int(-upstream, downstream),
               mean_peak_overlap = colMeans(tm), windows = nrow(tm))
}

plot_tss_profile <- function(tm, outdir) {
    tab <- tss_profile(tm)
    write.table(tab, file.path(outdir, "tss_avgprof.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    draw_tss_profile(tab)
    invisible(tab)
}

draw_tss_profile <- function(tab) {
    # Leave separate space for horizontal decimal tick labels and the rotated
    # title. This changes presentation only; the exported profile is untouched.
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mar = c(4.5, 6.5, 4, 2), mgp = c(2.8, 0.8, 0), las = 1)
    plot(tab$distance_bp, tab$mean_peak_overlap, type = "l", col = "#347CA5", lwd = 2,
         xlab = "Strand-relative distance to TSS (bp)", ylab = "",
         main = "Peak overlap around annotated TSS")
    mtext("Mean peak overlap per retained window", side = 2, line = 4.5, las = 0)
    abline(v = 0, col = "#777777", lty = 2)
    mtext(sprintf("%d retained TSS windows | display window: -3000 to +3000 bp | no smoothing", tab$windows[1]),
          side = 3, line = 0.4, cex = 0.75, adj = 0)
    invisible(tab)
}
