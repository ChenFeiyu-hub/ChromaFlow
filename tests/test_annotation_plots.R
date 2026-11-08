args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
out <- if (length(args) > 1) args[2] else tempfile("annotation-plots-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "lib/steps/annotation_plots.R"))
for (file in c("annotation_plots.R", "05_peak_annotation.R")) parse(file.path(root, "lib/steps", file))
detail <- data.frame(Promoter = c(TRUE, TRUE, FALSE, FALSE, TRUE),
                     Exon = c(FALSE, TRUE, TRUE, FALSE, FALSE),
                     Intron = c(FALSE, FALSE, FALSE, FALSE, TRUE))
detail <- rbind(detail, detail[c(1, 1), ])
tab <- annotation_intersections(detail)
stopifnot(sum(tab$table$peaks) == 7L, nrow(tab$table) == 5L, tab$table$peaks[1] == 3L,
          tab$table$peaks[tab$table$intersection == "Promoter & Exon"] == 1L)
skip <- function(expr) tryCatch({ force(expr); FALSE }, plot_skip = function(e) TRUE)
bad <- function(expr) tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
stopifnot(skip(annotation_intersections(detail[FALSE, ])))
invalid <- detail; invalid$Exon[1] <- NA
stopifnot(bad(annotation_intersections(invalid)))
tm <- matrix(0, nrow = 3, ncol = 6001)
tm[1, 2000:3000] <- 1; tm[2, 2500:3500] <- 1; tm[3, 3000:4000] <- 1
profile <- tss_profile(tm)
stopifnot(profile$distance_bp[3001] == 0, all(profile$mean_peak_overlap == colMeans(tm)),
          skip(tss_profile(tm[FALSE, , drop = FALSE])), skip(tss_profile(tm * 0)),
          bad(tss_profile(tm[, 1:100])), bad(tss_profile(tm * NA)), bad(tss_profile(-tm)))
report <- new_plot_report(out)
report$attempt("anno_upset.png", plot_annotation_intersections(detail, out), 1400, 800)
report$attempt("tss_avgprof.png", plot_tss_profile(tm, out), 1200, 600)
writeLines("stale plot", file.path(out, "empty.png"))
report$attempt("empty.png", plot_tss_profile(tm[FALSE, , drop = FALSE], out), 1200, 600)
report$attempt("invalid.png", plot_tss_profile(tm * NA, out), 1200, 600)
status <- report$finish()
stopifnot(identical(status$status, c("written", "written", "skipped", "failed")),
          !file.exists(file.path(out, "empty.png")), !file.exists(file.path(out, "invalid.png")),
          file.size(file.path(out, "anno_upset.png")) > 1000,
          file.size(file.path(out, "tss_avgprof.png")) > 1000,
          sum(read.delim(file.path(out, "annotation_intersections.tsv"))$peaks) == 7)
cat("[OK] Intersection counts, profile values, empty/invalid matrix handling and rendered plots verified\n")
