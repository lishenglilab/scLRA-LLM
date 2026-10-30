#!/usr/bin/env Rscript
# Tiny count matrices with hand-checkable medians, fractional counts and labels.
suppressPackageStartupMessages({library(Matrix); library(SeuratObject)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: make_qc_fixture.R output_directory")
out <- args[1]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells <- c("alpha_AAAA", "alpha_AAAC", "alpha_AAAG", "beta_CCCC", "beta_CCCA", "beta_CCCG")
g <- Matrix(matrix(c(1,3,0, 2,4,1, 3,5,2, 4,6,3, 5,7,4, 6,8,5), nrow = 3), sparse = TRUE)
rownames(g) <- c("MT-CO1", "G1", "G2")
colnames(g) <- cells
t <- Matrix(matrix(c(1,2,0,0, 2,3,0,1, 3,4,1,1, 4,5,2,1, 5,6,3,2, 6,7,4,3), nrow = 4), sparse = TRUE) * 1.5
rownames(t) <- c("tx1", "tx2", "tx3", "unannotated")
colnames(t) <- cells
options(Seurat.object.assay.version = "v4")
gene <- CreateSeuratObject(g, min.cells = 0, min.features = 0)
gene$double.cell <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
gene$percent.mt <- 99 # Must not be reused by the extractor.
tr <- CreateSeuratObject(t[, rev(seq_along(cells))], min.cells = 0, min.features = 0)
saveRDS(gene, file.path(out, "gene.rds"))
saveRDS(tr, file.path(out, "transcript.rds"))
missing_flag <- gene
missing_flag$double.cell[2] <- NA
saveRDS(missing_flag, file.path(out, "gene_missing_flag.rds"))
skipped <- gene
skipped$double.cell <- NULL
skipped$predicted_doublets <- FALSE
skipped$doublet_detection_performed <- FALSE
skipped$doublet_rate <- NA_real_
saveRDS(skipped, file.path(out, "gene_skipped_doublets.rds"))
performed <- skipped
performed$doublet_detection_performed <- TRUE
performed$doublet_rate <- 0
saveRDS(performed, file.path(out, "gene_performed_doublets.rds"))
no_mt <- skipped
no_mt[["RNA"]] <- CreateAssayObject(counts = g[-1, ])
saveRDS(no_mt, file.path(out, "gene_no_mt.rds"))
options(Seurat.object.assay.version = "v5")
ensembl <- g
rownames(ensembl)[1] <- "ENSGMT.2"
v5 <- CreateSeuratObject(ensembl, min.cells = 0, min.features = 0)
v5[["RNA"]] <- split(v5[["RNA"]], f = sub("_.*", "", colnames(v5)))
saveRDS(v5, file.path(out, "gene_v5.rds"))
saveRDS(tr[, -1], file.path(out, "transcript_missing_cell.rds"))
dir.create(file.path(out, "isoquant", "OUT"), recursive = TRUE, showWarnings = FALSE)
classification <- data.frame(isoform = c("tx1", "tx2", "tx3"), chr = "chr1", strand = "+",
  length = c(100,200,300), exons = c(1,2,3),
  structural_category = c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog"))
write.table(classification, file.path(out, "classification.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(classification, file.path(out, "isoquant", "OUT", "OUT.novel_vs_known.SQANTI-like.tsv"),
  sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
writeLines('chrM\tfixture\tgene\t1\t100\t.\t+\t.\tgene_id "ENSGMT.2"; gene_name "MT-CO1";', file.path(out, "reference.gtf"))
read_summary <- data.frame(CellBarcode = cells, unique_reads = seq(10,60,10), unique_UMIs = seq(8,48,8))
write.table(read_summary, file.path(out, "bam_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
groups <- unlist(lapply(cells, function(cell) {
  barcode <- sub("^[^_]+_", "", cell)
  ids <- paste0(barcode, "_", c("AAAA", "AAAC", "AAAC"), "#fixture", 1:3, "_+")
  c(paste(ids, cell, sep = "\t"), paste(ids[1], cell, sep = "\t"))
}))
writeLines(groups, file.path(out, "read_groups.tsv"))
cat("QC fixture created in", out, "\n")
