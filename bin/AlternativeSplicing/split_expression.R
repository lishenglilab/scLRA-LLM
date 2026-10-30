#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(optparse)
  library(jsonlite)
})
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
source(file.path(dirname(normalizePath(script)), "count_helpers.R"))
parser <- OptionParser(option_list = list(
  make_option("--seurat_file", type = "character"),
  make_option("--split_size", type = "integer", default = 100),
  make_option("--output_dir", type = "character", default = "split_exp")
))
opt <- parse_args(parser)
if (is.null(opt$seurat_file)) stop("--seurat_file is required")
if (is.na(opt$split_size) || opt$split_size < 1L) stop("--split_size must be positive")
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(opt$output_dir, pattern = "^part_[0-9]+\\.txt$")))
  stop("Output directory already contains expression parts; use an empty directory")
counts <- rna_counts(readRDS(opt$seurat_file))
abundance <- linear_abundance(counts)
totals <- Matrix::colSums(counts)
fwrite(data.table(cell_id = colnames(counts), transcript_library_count = totals),
       file.path(opt$output_dir, "cell_libraries.tsv"), sep = "\t")
parts <- split(seq_len(ncol(abundance)), ceiling(seq_len(ncol(abundance)) / opt$split_size))
for (i in seq_along(parts)) {
  block <- as.matrix(abundance[, parts[[i]], drop = FALSE])
  output <- file.path(opt$output_dir, paste0("part_", i, ".txt"))
  # SUPPA expects a header of cell IDs, then rows beginning with transcript ID.
  writeLines(paste(colnames(block), collapse = "\t"), output)
  fwrite(as.data.frame(block), output, append = TRUE, row.names = TRUE,
         col.names = FALSE, quote = FALSE, sep = "\t")
}
provenance <- list(assay = "RNA", layer = "counts", log = FALSE,
  transformation = "1e6 * count / full retained transcript library count per cell",
  length_adjusted = FALSE, zero_library_cells = sum(totals == 0),
  zero_library_behavior = "all transcript abundances remain zero; PSI may be undefined",
  transcripts = nrow(counts), cells = ncol(counts), parts = length(parts))
write_json(provenance, file.path(opt$output_dir, "input_provenance.json"),
           auto_unbox = TRUE, pretty = TRUE)
writeLines(c("SUPPA input: linear count-based per-million abundance; no log or length adjustment",
             paste("Cells:", ncol(counts)), paste("Transcripts:", nrow(counts)),
             paste("Parts:", length(parts))),
           file.path(opt$output_dir, "..", "split_info.txt"))
