# Exercise the production R guard and actual Python postprocessing without
# substituting a renderer stub or requiring Bioconductor to construct a TxDb.
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
python <- if (length(args) > 1) args[2] else "python3"
out <- if (length(args) > 2) args[3] else tempfile("report-guard-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
expressions <- parse(file.path(root, "lib/steps/05_peak_annotation.R"))
guard <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                    identical(x[[2]], as.name("verify_report_annotation")), as.list(expressions))
stopifnot(length(guard) == 1L)
eval(guard[[1]])
policy <- file.path(root, "lib/steps/support/annotation_policy.py")
df <- data.frame(seqnames = c("Chr1", "Chr2"), start = c(101L, 501L), end = c(200L, 700L),
                 peak_id = c("P1", "P2"), geneId = c("gene_1", 'gene "two"'),
                 annotation = c("Promoter (<=1kb)", "Intron (T2/gene_two, intron 1 of 3)"),
                 distanceToTSS = c(-10L, 3000L), peak_score = c(2.123456, 4.5678),
                 gene_consistent = c(TRUE, FALSE), missing_value = c(NA_real_, NA_real_))
saved <- file.path(out, "saved.tsv")
write.table(df, saved, sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(system2(python, c(shQuote(policy), "validate", shQuote(saved), "--up", "1500", "--down", "0")) == 0)
original <- readLines(saved)
# Reproduce the former failure even when the biological content is unchanged.
stopifnot(!identical(original, capture.output(write.table(df, sep = "\t", quote = FALSE, row.names = FALSE))))
stopifnot(verify_report_annotation(df, saved, out, 1500L, 0L, policy, python))
stopifnot(identical(original, readLines(saved)))
bad <- function(expr) tryCatch({ force(expr); FALSE }, error = function(e) {
    cat("[EXPECTED REJECTION]", conditionMessage(e), "\n"); TRUE
})
changed <- df; changed$geneId[1] <- "changed_gene"
stopifnot(bad(verify_report_annotation(changed, saved, out, 1500L, 0L, policy, python)))
changed <- df; changed$peak_score[1] <- 2.123457
stopifnot(bad(verify_report_annotation(changed, saved, out, 1500L, 0L, policy, python)))
stopifnot(bad(verify_report_annotation(df, saved, out, 0L, 0L, policy, python)))
stopifnot(identical(original, readLines(saved)),
          length(list.files(out, pattern = "^reconstructed-annotation-")) == 0L)
cat("[OK] Real R serialization + Python feature policy + production R guard passed; changed evidence rejected\n")
