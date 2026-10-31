#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(SeuratObject); library(Matrix); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript demo/scripts/check_objects.R RESULTS")
folder <- file.path(args[1], "PreliminarySingleCellAnalysis")
objects <- lapply(c("sce.gene.rds", "sce.tr.rds"), function(name) readRDS(file.path(folder, name)))
stopifnot(identical(colnames(objects[[1]]), colnames(objects[[2]])), ncol(objects[[1]]) > 30)
summary <- lapply(objects, function(object) {
  counts <- GetAssayData(object, assay = "RNA", layer = "counts")
  stopifnot(all(is.finite(counts@x)), all(counts@x >= 0), all(colSums(counts) > 0))
  stopifnot(setequal(unique(sub("_[^_]+$", "", colnames(object))), c("DemoA", "DemoB")))
  stopifnot(identical(as.character(object$orig.ident), sub("_[^_]+$", "", colnames(object))))
  list(features = nrow(counts), cells = ncol(counts),
       cells_by_sample = as.list(table(object$orig.ident)),
       min_detected_features = min(colSums(counts > 0)), max_detected_features = max(colSums(counts > 0)),
       doublet_detection_performed = if ("doublet_detection_performed" %in% colnames(object[[]]))
         all(object$doublet_detection_performed) else NULL)
})
names(summary) <- c("gene", "transcript")
write_json(summary, file.path(args[1], "demo_object_checks.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
cat("Demo Seurat objects verified: matching ordered cells, two technical groups, positive raw counts.\n")
