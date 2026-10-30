#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(SeuratObject); library(Matrix); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
tmp <- tempfile("sclra-dominance-"); dir.create(tmp)
cells <- c("A_c1", "A_c2", "B_c1", "B_tied", "A_empty")
counts <- Matrix(matrix(c(2,0,0,0, 0,3,0,0, 0,4,0,0, 0,0,1,1, 0,0,0,0), nrow = 4), sparse = TRUE)
dimnames(counts) <- list(paste0("tx", 1:4), cells)
options(Seurat.object.assay.version = "v4")
gtf <- paste0('chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g"; gene_name "G"; transcript_id "tx', 1:4, '";')
writeLines(gtf, file.path(tmp, "transcripts.gtf"))
fwrite(data.table(sample = c("A", "B"), stage = c("case", "control")), file.path(tmp, "mapping.tsv"), sep = "\t")
fwrite(data.table(group1 = "case", group2 = "control", name = "case_vs_control"), file.path(tmp, "comparisons.tsv"), sep = "\t")
run <- function(matrix, label) {
  saveRDS(CreateSeuratObject(matrix, min.cells = 0, min.features = 0), file.path(tmp, "object.rds"))
  output <- file.path(tmp, label)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
    file.path(root, "bin/AlternativeSplicing/dominant_isoform_detection.R"),
    "--seurat_rds", file.path(tmp, "object.rds"), "--gtf_file", file.path(tmp, "transcripts.gtf"),
    "--stage_mapping", file.path(tmp, "mapping.tsv"), "--comparisons", file.path(tmp, "comparisons.tsv"),
    "--output_dir", output, "--n_cores", "1")))
  stopifnot(status == 0)
  output
}
out <- run(counts, "mixed")
calls <- fread(file.path(out, "dominant_transcripts_per_cell.txt"))
stopifnot(nrow(calls) == 3L, !"B_tied" %in% calls$Cellid)
audit <- fread(file.path(out, "dominance_calls.tsv"))
stopifnot(audit[Cellid == "B_tied", dominance_status] == "tie",
          is.na(audit[Cellid == "B_tied", DominantTranscript]))
stats <- fread(file.path(out, "case_vs_control/statistics_case_vs_control.txt"))
stopifnot(stats[transcript_id == "tx2", control_pct] == .5,
          stats[transcript_id == "tx2", case_pct] == .5)
usage <- fread(file.path(out, "transcript_usage.txt"))
stopifnot(all(is.na(usage$A_empty)))
all_ties <- counts
all_ties[,] <- 1
empty <- run(all_ties, "all_ties")
stopifnot(nrow(fread(file.path(empty, "dominant_transcripts_per_cell.txt"))) == 0L,
          nrow(fread(file.path(empty, "case_vs_control/statistics_case_vs_control.txt"))) == 0L)
cat("Dominance integration passed: unique calls, ties, empty genes and complete denominators.\n")
