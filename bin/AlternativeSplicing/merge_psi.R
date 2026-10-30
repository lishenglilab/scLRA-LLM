#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: merge_psi.R part.psi [...]")
merged <- NULL
missing <- NULL
seen_cells <- character()
for (path in sort(args)) {
  dt <- fread(path, na.strings = c("NA", "nan", "NaN"))
  setnames(dt, 1L, "ASE")
  if (ncol(dt) < 2L || anyDuplicated(dt$ASE)) stop("Invalid PSI event table: ", path)
  cells <- names(dt)[-1L]
  if (anyDuplicated(cells) || any(cells %in% seen_cells)) stop("Duplicate PSI cell IDs")
  seen_cells <- c(seen_cells, cells)
  if (!is.null(merged) && !identical(dt$ASE, merged$ASE)) stop("PSI event order mismatch")
  if (is.null(missing)) missing <- data.table(ASE = dt$ASE)
  samples <- sub("_[^_]+$", "", cells)
  for (sample in unique(samples)) {
    columns <- cells[samples == sample]
    count <- rowSums(is.na(dt[, ..columns]))
    if (sample %in% names(missing)) count <- count + missing[[sample]]
    set(missing, j = sample, value = count)
  }
  for (column in cells) {
    values <- dt[[column]]
    if (!is.numeric(values) && !is.logical(values)) stop("Nonnumeric PSI column")
    if (any(!is.na(values) & (!is.finite(values) | values < 0 | values > 1)))
      stop("PSI values must be in [0, 1] or undefined")
    set(dt, j = column, value = as.numeric(values))
    set(dt, which(is.na(values)), column, 0)
  }
  if (is.null(merged)) merged <- dt else merged <- cbind(merged, dt[, -1L, with = FALSE])
}
fwrite(merged, "merged_project_events.psi", sep = "\t", quote = FALSE)
fwrite(missing, "undefined_psi_by_sample.tsv", sep = "\t", quote = FALSE)
