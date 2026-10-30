#!/usr/bin/env Rscript
# Descriptive QC from the actual post-filter RNA count matrices.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(jsonlite)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2L) stop("Arguments must be --name value pairs")
opt <- setNames(as.list(args[seq(2, length(args), 2)]),
                sub("^--", "", args[seq(1, length(args), 2)]))
for (key in c("gene", "transcript", "isoquant_dir", "reference_gtf", "output_dir"))
  if (is.null(opt[[key]])) stop("Missing --", key)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- function(obj) {
  if (!"RNA" %in% names(obj@assays)) stop("An RNA assay with linear counts is required")
  assay <- obj@assays$RNA
  if ("counts" %in% slotNames(assay)) return(assay@counts)
  layers <- Layers(assay, search = "^counts($|\\.)")
  if (!length(layers)) stop("RNA counts layer is absent")
  if (length(layers) > 1L) assay <- JoinLayers(assay, layers = "counts", new = "counts")
  LayerData(assay, layer = if (length(layers) > 1L) "counts" else layers[1])
}
check_counts <- function(x) {
  if (is.null(colnames(x)) || is.null(rownames(x)) ||
      anyDuplicated(colnames(x)) || anyDuplicated(rownames(x)))
    stop("Count matrices require unique feature and cell IDs")
  values <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
  if (any(!is.finite(values)) || any(values < 0)) stop("Invalid linear RNA counts")
}
flag_value <- function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(!is.na(x) & x %in% c("true", "false"), x == "true", NA)
}
median_ok <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
cor_ok <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L)
    return(NA_real_)
  cor(x[keep], y[keep], method = method)
}

gene <- readRDS(opt$gene)
g <- counts(gene)
check_counts(g)
md <- gene@meta.data[colnames(g), , drop = FALSE]
n_source <- ncol(g)
doublet_source <- NULL
doublet <- rep(NA, ncol(g))
if ("double.cell" %in% names(md)) {
  doublet <- flag_value(md$double.cell)
  keep <- !is.na(doublet)
  g <- g[, keep, drop = FALSE]
  md <- md[keep, , drop = FALSE]
  doublet <- doublet[keep]
  doublet_source <- "double.cell"
} else if ("predicted_doublets" %in% names(md)) {
  performed <- if ("doublet_detection_performed" %in% names(md)) {
    flag_value(md$doublet_detection_performed) %in% TRUE
  } else if ("doublet_rate" %in% names(md)) {
    is.finite(suppressWarnings(as.numeric(md$doublet_rate)))
  } else rep(FALSE, nrow(md))
  doublet[performed] <- flag_value(md$predicted_doublets[performed])
  if (any(!is.na(doublet))) doublet_source <- "predicted_doublets"
}
if (!ncol(g)) stop("No eligible filtered cells remain for the QC report")
cells <- colnames(g)
# MAKE_READ_GROUPS uses sample_barcode; sample IDs cannot contain underscores.
if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9.-]*_.+$", cells)))
  stop("Expected pipeline cell IDs in sample_barcode format")
sample_id <- sub("_.*", "", cells)
sample_order <- sort(unique(sample_id))
tr <- readRDS(opt$transcript)
t <- counts(tr)
check_counts(t)
if (!all(cells %in% colnames(t))) stop("Filtered transcript object does not cover the report cell cohort")
t <- t[, match(cells, colnames(t)), drop = FALSE]
stopifnot(identical(colnames(g), colnames(t)))

mt_features <- grepl("^MT-", rownames(g), ignore.case = TRUE)
mt_source <- "RNA feature names"
if (!any(mt_features)) {
  # Resolve Ensembl IDs through the supplied reference; do not reuse inherited percent.mt.
  ref <- if (grepl("\\.gz$", opt$reference_gtf)) gzfile(opt$reference_gtf, "rt") else file(opt$reference_gtf, "rt")
  mt_ids <- character()
  repeat {
    lines <- readLines(ref, n = 100000L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines, "#") & grepl("\tgene\t", lines, fixed = TRUE)]
    for (line in lines) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(parts) < 9L) next
      attr <- parts[9]
      name <- if (grepl('gene_name "', attr, fixed = TRUE))
        sub('.*gene_name "([^"]+)".*', '\\1', attr) else ""
      if (grepl("^MT-", name, ignore.case = TRUE) || tolower(parts[1]) %in% c("chrm", "mt", "m")) {
        id <- sub('.*gene_id "([^"]+)".*', '\\1', attr)
        mt_ids <- c(mt_ids, id, name)
      }
    }
  }
  close(ref)
  strip_version <- function(x) sub("\\.[0-9]+$", "", x)
  mt_features <- strip_version(rownames(g)) %in% strip_version(unique(mt_ids))
  mt_source <- "Reference GTF gene IDs/names and mitochondrial contigs"
}
gc <- as.numeric(Matrix::colSums(g))
tc <- as.numeric(Matrix::colSums(t))
gd <- as.numeric(Matrix::colSums(g > 0))
td <- as.numeric(Matrix::colSums(t > 0))
per_cell <- data.frame(cell_id = cells, sample = sample_id, gene_counts = gc,
  transcript_counts = tc, genes_detected = gd, transcripts_detected = td,
  detection_ratio = ifelse(gd > 0, td/gd, NA_real_))
if (any(mt_features)) per_cell$mitochondrial_percent <- ifelse(gc > 0,
  100 * as.numeric(Matrix::colSums(g[mt_features, , drop = FALSE])) / gc, NA_real_)
if (!is.null(doublet_source)) per_cell$doublet_label <- doublet

classification_path <- opt$classification
if (is.null(classification_path) || !nzchar(classification_path)) {
  candidates <- file.path(opt$isoquant_dir, "OUT", c("OUT.classification.txt",
    "OUT.SQANTI-like.tsv", "OUT.novel_vs_known.SQANTI-like.tsv"))
  classification_path <- candidates[file.exists(candidates)][1]
}
cl <- NULL
if (length(classification_path) && !is.na(classification_path) && file.exists(classification_path)) {
  first <- strsplit(readLines(classification_path, n = 1L), "\t", fixed = TRUE)[[1]]
  header <- "isoform" %in% sub("^#", "", first)
  raw <- fread(classification_path, header = header, na.strings = c("NA", ""))
  if (header) {
    setnames(raw, sub("^#", "", names(raw)))
    required <- c("isoform", "length", "exons", "structural_category")
    if (!all(required %in% names(raw))) stop("Classification table lacks required SQANTI columns")
    cl <- as.data.frame(raw[, ..required])
  } else {
    if (ncol(raw) < 6L) stop("Headerless SQANTI table must have at least six columns")
    cl <- as.data.frame(raw[, c(1L, 4L, 5L, 6L), with = FALSE])
    names(cl) <- c("isoform", "length", "exons", "structural_category")
  }
  if (anyDuplicated(cl$isoform)) stop("Duplicate transcript IDs in classification table")
  for (key in c("length", "exons")) {
    cl[[key]] <- suppressWarnings(as.numeric(cl[[key]]))
    if (any(is.finite(cl[[key]]) & cl[[key]] <= 0)) stop("Invalid model length/exon count")
  }
}
categories <- c(fsm = "full-splice_match", ism = "incomplete-splice_match",
  nic = "novel_in_catalog", nnc = "novel_not_in_catalog", antisense = "antisense",
  intergenic = "intergenic", fusion = "fusion", genic = "genic")
records <- lapply(sample_order, function(s) {
  keep <- sample_id == s
  d <- per_cell[keep, , drop = FALSE]
  pos <- Matrix::rowSums(t[, keep, drop = FALSE]) > 0
  row <- list(sample = s, cells = nrow(d),
    gene_features = sum(Matrix::rowSums(g[, keep, drop = FALSE]) > 0), transcript_features = sum(pos),
    gene_counts_median = median_ok(d$gene_counts), transcript_counts_median = median_ok(d$transcript_counts),
    genes_detected_median = median_ok(d$genes_detected), transcripts_detected_median = median_ok(d$transcripts_detected),
    detection_ratio_median = median_ok(d$detection_ratio),
    r_gene_transcript_detection = cor_ok(d$genes_detected, d$transcripts_detected),
    rho_gene_transcript_detection = cor_ok(d$genes_detected, d$transcripts_detected, "spearman"),
    r_gene_transcript_counts = cor_ok(d$gene_counts, d$transcript_counts))
  if ("mitochondrial_percent" %in% names(d)) row$mt_median_pct <- median_ok(d$mitochondrial_percent)
  if ("doublet_label" %in% names(d) && any(!is.na(d$doublet_label))) {
    row$double_cell_available <- sum(!is.na(d$doublet_label))
    row$double_cell_true <- sum(d$doublet_label, na.rm = TRUE)
    row$double_cell_false <- sum(!d$doublet_label, na.rm = TRUE)
    row$double_cell_true_pct <- 100 * mean(d$doublet_label, na.rm = TRUE)
  }
  if (!is.null(cl) && sum(pos)) {
    m <- match(rownames(t)[pos], cl$isoform)
    row$annotation_coverage_pct <- 100 * mean(!is.na(m))
    if (any(!is.na(m))) {
      labels <- cl$structural_category[m]
      # Denominator is all positively detected transcript features, including unmatched IDs.
      for (key in names(categories)) row[[paste0(key, "_pct")]] <-
        100 * sum(labels == categories[[key]], na.rm = TRUE) / sum(pos)
      for (key in c("length", "exons")) {
        x <- cl[[key]][m]
        x <- x[is.finite(x)]
        if (length(x)) {
          values <- c(median = median(x), q25 = unname(quantile(x, .25)),
            q75 = unname(quantile(x, .75)), min = min(x), max = max(x))
          for (nm in names(values)) row[[paste0("model_", key, "_", nm)]] <- values[[nm]]
        }
      }
    }
  }
  row
})
fwrite(per_cell, file.path(opt$output_dir, "cells.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")
write_json(records, file.path(opt$output_dir, "sample_metrics.json"), auto_unbox = TRUE,
  pretty = TRUE, digits = 16, na = "null")
manifest <- list(source_gene_cells = n_source, report_cells = nrow(per_cell),
  excluded_without_valid_double_cell = n_source - nrow(per_cell),
  cohort_policy = if (identical(doublet_source, "double.cell")) "Filtered cells with valid existing double.cell records" else "Filtered gene/transcript cell intersection",
  doublet_source = doublet_source, mitochondrial_source = if (any(mt_features)) mt_source else NULL,
  mitochondrial_features = rownames(g)[mt_features], counts_scale = "linear RNA counts",
  classification_file = if (!is.null(cl)) basename(classification_path) else NULL,
  classification_records = if (!is.null(cl)) nrow(cl) else 0L,
  exact_transcript_cell_alignment = identical(colnames(g), colnames(t)),
  software = list(R = R.version.string, Matrix = as.character(packageVersion("Matrix")),
    SeuratObject = as.character(packageVersion("SeuratObject"))))
write_json(manifest, file.path(opt$output_dir, "extraction.json"), auto_unbox = TRUE,
  pretty = TRUE, digits = 16, na = "null", null = "null")
cat("Extracted", nrow(per_cell), "filtered cells in", length(records), "samples\n")
