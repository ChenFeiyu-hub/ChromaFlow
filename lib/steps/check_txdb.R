#!/usr/bin/env Rscript
# Export the actual installed importer's identities. Python performs strict
# comparison using explicit GFF ID/Name/gene_id/transcript_id relationships.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || args[4] != "export-for-identity-v2")
    stop("Use the patched validate_annotation.py with check_txdb.R; mixed old/new files are not supported")
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
options(warn = 1, scipen = 999)
suppressPackageStartupMessages(library(GenomicFeatures))
makeTx <- if (requireNamespace("txdbmaker", quietly = TRUE))
    txdbmaker::makeTxDbFromGFF else GenomicFeatures::makeTxDbFromGFF
writeLines(capture.output(sessionInfo()), file.path(args[3], "R-sessionInfo.txt"))
warnings <- character()
txdb <- withCallingHandlers(makeTx(args[1], format = "gff3"), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
})
writeLines(warnings, file.path(args[3], "TxDb-warnings.txt"))
tx <- transcripts(txdb, columns = c("tx_id", "tx_name", "gene_id"))
observed <- data.frame(tx_id = as.character(mcols(tx)$tx_name),
                       seqid = as.character(seqnames(tx)), strand = as.character(strand(tx)),
                       tss = ifelse(as.character(strand(tx)) == "-", end(tx), start(tx)),
                       txdb_id = as.character(mcols(tx)$tx_id), tx_start = start(tx), tx_end = end(tx),
                       stringsAsFactors = FALSE)
gene_lists <- lapply(mcols(tx)$gene_id, as.character)
observed$gene_ids <- vapply(gene_lists, paste, collapse = ",", FUN.VALUE = character(1))
write.table(observed, file.path(args[3], "TxDb-transcripts.tsv"), sep = "\t", quote = TRUE,
            qmethod = "double", na = "", row.names = FALSE)
# One row per edge avoids ambiguity if a gene identifier itself contains commas.
edges <- data.frame(txdb_id = rep(observed$txdb_id, lengths(gene_lists)),
                    gene_id = as.character(unlist(gene_lists, use.names = FALSE)))
write.table(edges, file.path(args[3], "TxDb-transcript-genes.tsv"), sep = "\t", quote = TRUE,
            qmethod = "double", na = "", row.names = FALSE)
cat(sprintf("[EXPORT] %d TxDb transcripts; identity comparison follows in Python\n", nrow(observed)))
