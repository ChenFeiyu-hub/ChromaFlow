args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[1], "lib/steps/tss_windows.R"))
bad <- function(expr) tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
sizes <- c(Chr1 = 100L, Chr2 = 90L, Chr3 = 40L)
# Deliberately use factor sequence columns in non-alphabetic level order, as
# GRanges -> data.frame does. Positive and negative windows have the same TSS.
p <- data.frame(seqnames = factor(c("Chr2", "Chr2", "Chr1"), levels = c("Chr3", "Chr2", "Chr1")),
                start = c(10L, 12L, 7L), end = c(12L, 13L, 7L), strand = "*")
w <- data.frame(seqnames = factor(c("Chr2", "Chr2", "Chr1", "Chr3", "Chr1", "Chr2")),
                start = c(8L, 8L, 5L, 5L, -1L, 88L), end = c(14L, 14L, 11L, 11L, 5L, 94L),
                strand = c("+", "-", "+", "+", "+", "-"))
result <- tss_matrix_from_hits(p, w, sizes, c(1L, 1L, 2L, 2L, 3L), c(1L, 2L, 1L, 2L, 3L), 3L, 3L)
stopifnot(identical(dim(result$matrix), c(3L, 7L)),
          identical(as.numeric(result$matrix[1, ]), c(0, 0, 1, 1, 2, 1, 0)),
          identical(as.numeric(result$matrix[2, ]), c(0, 1, 2, 1, 1, 0, 0)),
          identical(result$audit$status, c(rep("retained", 3), "no_peak_overlap", rep("outside_reference", 2))),
          sum(result$audit$legacy_end_bound_exclusion) == 3L)
one <- tss_matrix_from_hits(p[3, , drop=FALSE], w[3, , drop=FALSE], sizes, 1L, 1L, 3L, 3L)
stopifnot(identical(dim(one$matrix), c(1L, 7L)))
zero <- tss_matrix_from_hits(p, w[4, , drop=FALSE], sizes, integer(), integer(), 3L, 3L)
stopifnot(identical(dim(zero$matrix), c(0L, 7L)))
stopifnot(bad(tss_matrix_from_hits(p, w, sizes[-2], integer(), integer(), 3L, 3L)),
          bad(tss_matrix_from_hits(p, w, sizes, 4L, 1L, 3L, 3L)),
          bad(tss_matrix_from_hits(p, w, sizes, c(1L, 1L), c(1L, 1L), 3L, 3L)))
# Independent per-base enumeration across randomized overlapping peaks/windows.
set.seed(26)
for (iteration in seq_len(40L)) {
    pp <- data.frame(seqnames="Chr1", start=sample(1:90, 20, replace=TRUE), strand="*")
    pp$end <- pmin(100L, pp$start + sample(0:10, 20, replace=TRUE))
    ww <- data.frame(seqnames="Chr1", start=sample(-5:100, 15), strand=sample(c("+", "-"), 15, replace=TRUE))
    ww$end <- ww$start+6L
    hit <- which(outer(ww$start, pp$end, "<=") & outer(ww$end, pp$start, ">="), arr.ind=TRUE)
    actual <- tss_matrix_from_hits(pp, ww, sizes, hit[,1], hit[,2], 3L, 3L)
    for (j in seq_len(nrow(actual$matrix))) {
        i <- as.integer(rownames(actual$matrix)[j])
        coordinates <- seq.int(ww$start[i], ww$end[i])
        if (ww$strand[i] == "-") coordinates <- rev(coordinates)
        expected <- vapply(coordinates, function(x) sum(pp$start <= x & pp$end >= x), 0L)
        stopifnot(all(actual$matrix[j, ] == expected))
    }
}
# Explicit IDR BED8 only; narrowPeak and custom column 7 must stay NA.
df <- data.frame(V7=c(4, 7), V8=c(1, 3))
stopifnot(all(is.na(peak_idr_values(df, "non_idr"))), all(is.na(peak_idr_values(df, "unknown"))),
          identical(peak_idr_values(df, "idr"), c(4, 7)),
          bad(peak_idr_values(cbind(df, V9=8), "idr")),
          bad(peak_idr_values(data.frame(V7=3, V8=0), "idr")))
cat("[OK] Reference bounds, factor sequences, strand orientation, single/zero windows, overlapping peaks, randomized per-base comparison and IDR semantics passed\n")
