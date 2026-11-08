# Full-width, strand-oriented TSS peak overlap with explicit reference bounds.
# Coordinates in this helper are 1-based, inclusive. Peaks are not weighted.
read_tss_sizes <- function(path) {
    if (!nzchar(path) || !file.exists(path))
        stop("TSS profiles require CHROM_SIZES from the prepared reference")
    tab <- read.table(path, header = FALSE, sep = "\t", quote = "", comment.char = "",
                      colClasses = "character", stringsAsFactors = FALSE)
    if (ncol(tab) != 2L || !nrow(tab) || anyDuplicated(tab[[1]]) || any(!nzchar(tab[[1]])))
        stop("Invalid or duplicate sequence names in CHROM_SIZES")
    values <- suppressWarnings(as.numeric(tab[[2]]))
    if (any(!is.finite(values)) || any(values < 1 | values != floor(values)))
        stop("CHROM_SIZES must contain positive integer lengths")
    setNames(values, tab[[1]])
}

tss_matrix_from_hits <- function(peaks, windows, sizes, query, subject,
                                 upstream = 3000L, downstream = 3000L) {
    columns <- c("seqnames", "start", "end", "strand")
    if (!all(columns %in% names(windows)) || !all(columns %in% names(peaks)))
        stop("Missing interval columns for TSS profile")
    peaks$seqnames <- as.character(peaks$seqnames); windows$seqnames <- as.character(windows$seqnames)
    peaks$strand <- as.character(peaks$strand); windows$strand <- as.character(windows$strand)
    if (anyNA(sizes) || any(!is.finite(sizes)) || any(sizes < 1) || anyDuplicated(names(sizes)))
        stop("Invalid chromosome lengths")
    for (x in list(peaks, windows)) {
        if (anyNA(x[, columns]) || any(!is.finite(x$start) | !is.finite(x$end)) ||
            any(x$start != floor(x$start) | x$end != floor(x$end) | x$end < x$start))
            stop("Invalid TSS input coordinates")
        if (any(!x$seqnames %in% names(sizes)))
            stop("TSS input sequence is absent from CHROM_SIZES")
    }
    if (any(peaks$start < 1 | peaks$end > sizes[peaks$seqnames]))
        stop("Peak extends outside the reference sequence")
    width <- upstream + downstream + 1L
    if (any(windows$end - windows$start + 1 != width) || any(!windows$strand %in% c("+", "-")))
        stop("Unexpected TSS window width or strand")
    if (length(query) != length(subject) || anyNA(query) || anyNA(subject) ||
        any(query < 1 | query > nrow(windows) | query != floor(query)) ||
        any(subject < 1 | subject > nrow(peaks) | subject != floor(subject)) ||
        anyDuplicated(paste(query, subject))) stop("Invalid TSS overlap pairs")
    if (length(query) && any(windows$seqnames[query] != peaks$seqnames[subject] |
        windows$start[query] > peaks$end[subject] | windows$end[query] < peaks$start[subject]))
        stop("Non-overlapping intervals supplied as TSS overlap pairs")
    audit <- windows[, columns, drop = FALSE]
    audit$window_id <- seq_len(nrow(windows))
    audit$tss <- ifelse(audit$strand == "+", audit$start + upstream, audit$start + downstream)
    audit$chrom_length <- unname(sizes[audit$seqnames])
    audit$overlapping_peaks <- tabulate(query, nbins = nrow(windows))
    audit$status <- ifelse(audit$start < 1 | audit$end > audit$chrom_length,
                           "outside_reference", ifelse(audit$overlapping_peaks == 0L, "no_peak_overlap", "retained"))
    last_peak <- tapply(peaks$end, peaks$seqnames, max)
    audit$legacy_end_bound_exclusion <- audit$status == "retained" &
        audit$end > unname(last_peak[audit$seqnames])
    audit$legacy_end_bound_exclusion[is.na(audit$legacy_end_bound_exclusion)] <- FALSE
    keep <- which(audit$status == "retained")
    matrix <- matrix(0, nrow = length(keep), ncol = width,
                     dimnames = list(as.character(keep), NULL))
    hits <- split(subject, query)
    for (i in seq_along(keep)) {
        w <- keep[i]; p <- hits[[as.character(w)]]
        left <- pmax(peaks$start[p], windows$start[w]) - windows$start[w] + 1L
        right <- pmin(peaks$end[p], windows$end[w]) - windows$start[w] + 1L
        changes <- tabulate(left, nbins = width + 1L) - tabulate(right + 1L, nbins = width + 1L)
        values <- cumsum(changes)[seq_len(width)]
        matrix[i, ] <- if (windows$strand[w] == "-") rev(values) else values
    }
    list(matrix = matrix, audit = audit)
}

reference_tss_matrix <- function(peaks, windows, chrom_sizes, outdir) {
    sizes <- read_tss_sizes(chrom_sizes)
    # findOverlaps supplies interval pairs only; coverage limits cannot truncate
    # a window. The numeric helper preserves one-row matrices and row orientation.
    hits <- GenomicRanges::findOverlaps(windows, peaks, ignore.strand = TRUE)
    result <- tss_matrix_from_hits(as.data.frame(peaks), as.data.frame(windows), sizes,
                                  S4Vectors::queryHits(hits), S4Vectors::subjectHits(hits))
    write.table(result$audit, file.path(outdir, "tss_windows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    counts <- table(factor(result$audit$status, levels = c("retained", "no_peak_overlap", "outside_reference")))
    summary <- data.frame(metric = c(names(counts), "legacy_end_bound_exclusion"),
        windows = c(as.integer(counts), sum(result$audit$legacy_end_bound_exclusion)))
    write.table(summary, file.path(outdir, "tss_window_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("[TSS] %d retained windows; %d valid windows would be lost with the old last-peak boundary\n",
                nrow(result$matrix), sum(result$audit$legacy_end_bound_exclusion)))
    result$matrix
}

peak_idr_values <- function(df, kind) {
    if (!kind %in% c("idr", "non_idr", "unknown")) stop("Invalid explicit peak-set kind")
    if (kind != "idr") return(rep(NA_real_, nrow(df)))
    if (!all(c("V7", "V8") %in% names(df)) || "V9" %in% names(df))
        stop("IDR peak sets require the ChromaFlow BED8 schema")
    values <- suppressWarnings(as.numeric(df$V7))
    pairs <- suppressWarnings(as.numeric(df$V8))
    if (any(!is.finite(values) | values < 0) || any(!is.finite(pairs) | pairs < 1 | pairs != floor(pairs)))
        stop("Invalid IDR values or supporting-pair counts")
    values
}
