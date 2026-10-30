#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Matrix); library(SeuratObject); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".")
source(file.path(root, "bin/AlternativeSplicing/count_helpers.R"))
source(file.path(root, "bin/AlternativeSplicing/ase_statistics.R"))
tmp <- tempfile("sclra-splicing-")
dir.create(tmp)
counts <- Matrix(matrix(c(1, 3, 0, 3, 1, 4, 0, 0, 0, 1, 1, 0), nrow = 3), sparse = TRUE)
dimnames(counts) <- list(c("tx1", "tx2", "tx3"), c("A_c1", "B_c2", "A_empty", "B_c4"))
linear <- linear_abundance(counts)
stopifnot(identical(dimnames(linear), dimnames(counts)), linear[1, 1] == 250000,
          linear[2, 1] == 750000, linear[1, 2] == 375000, all(linear[, 3] == 0))
# Transcript-wise log transforms change the relative abundance; library scaling does not.
stopifnot(abs(linear[1, 1] / sum(linear[1:2, 1]) - 0.25) < 1e-12,
          abs(log1p(2500) / sum(log1p(c(2500, 7500))) - 0.25) > 0.1)
for (version in c("v4", "v5")) {
  options(Seurat.object.assay.version = version)
  object <- CreateSeuratObject(counts, min.cells = 0, min.features = 0)
  if (version == "v5") object[["RNA"]] <- split(object[["RNA"]], f = c("A", "B", "A", "B"))
  actual <- rna_counts(object)
  stopifnot(isTRUE(all.equal(as.matrix(actual)[rownames(counts), colnames(counts)], as.matrix(counts))))
  saveRDS(object, file.path(tmp, paste0(version, ".rds")))
}
stopifnot(dominant_call(c(1, 1), c("t1", "t2"))$dominance_status == "tie",
          is.na(dominant_call(c(1, 1), c("t2", "t1"))$DominantTranscript),
          dominant_call(c(1, 3), c("t1", "t2"))$DominantTranscript == "t2",
          dominant_call(c(0, 0), c("t1", "t2"))$dominance_status == "unexpressed")
set.seed(24)
cases <- list(list(a = c(0, .1, .4), b = c(0, 0, .3)),
              list(a = rep(0, 60), b = rep(0, 70)),
              list(a = c(rep(0, 40), rep(.1, 15), rep(.8, 10)),
                   b = c(rep(0, 45), rep(.1, 15), rep(.9, 10))))
for (i in 1:20) cases[[length(cases) + 1L]] <- list(
  a = sample(c(0, .1, .2, .5, 1), sample(50:200, 1), replace = TRUE),
  b = sample(c(0, .1, .2, .5, 1), sample(50:200, 1), replace = TRUE))
for (case in cases) {
  a <- case$a; b <- case$b
  expected <- suppressWarnings(wilcox.test(a, b)$p.value)
  actual <- zero_wilcox(a[a > 0], b[b > 0], length(a), length(b))
  stopifnot(identical(is.na(expected), is.na(actual)),
            is.na(expected) || abs(expected - actual) < 1e-12)
}
mapping <- data.table(sample = c("A", "B"), stage = c("case", "control"))
fwrite(mapping, file.path(tmp, "mapping.tsv"), sep = "\t")
psi <- data.table(ASE = c("event1", "event2", "all_zero"),
                  A_c1 = c(.9, .1, 0), A_c2 = c(.8, .2, 0),
                  B_c3 = c(.1, .9, 0), B_c4 = c(.2, .8, 0))
fwrite(psi, file.path(tmp, "merged.psi"), sep = "\t")
opt <- list(psi_file = file.path(tmp, "merged.psi"), stage_mapping = file.path(tmp, "mapping.tsv"),
            group1 = "case", group2 = "control", output_dir = file.path(tmp, "results"),
            fdr = .05, fold_change = 1.5, pseudocount = 1e-6)
run_statistics(opt)
result <- fread(file.path(opt$output_dir, "all_events.tsv"))
stopifnot(nrow(result) == 3L, result$log2FC[1] > 0, result$log2FC[2] < 0,
          is.na(result$p_value[3]), !result$selected[3],
          isTRUE(all.equal(result$FDR, p.adjust(result$p_value, "BH"))))
check_merge <- function() {
  previous <- getwd(); on.exit(setwd(previous)); setwd(tmp)
  fwrite(data.table(event = c("e1", "e2"), A_c1 = c(.5, .1), A_empty = c(NA_real_, NA_real_)),
         "part1.psi", sep = "\t", na = "NA")
  fwrite(data.table(event = c("e1", "e2"), B_c2 = c(.7, .4)), "part2.psi", sep = "\t")
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
    file.path(root, "bin/AlternativeSplicing/merge_psi.R"), "part1.psi", "part2.psi")))
  stopifnot(status == 0)
  merged <- fread("merged_project_events.psi")
  missing <- fread("undefined_psi_by_sample.tsv")
  stopifnot(all(merged$A_empty == 0), all(missing$A == 1), all(missing$B == 0))
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
    file.path(root, "bin/AlternativeSplicing/merge_psi.R"), "part1.psi", "part1.psi")),
    stdout = TRUE, stderr = TRUE))
  stopifnot(!is.null(attr(output, "status")), attr(output, "status") != 0,
            any(grepl("Duplicate PSI cell IDs", output)))
}
check_merge()
cat("Splicing numerical checks passed (v4/v5, scale, ties, Wilcoxon, all-event BH).\n")
