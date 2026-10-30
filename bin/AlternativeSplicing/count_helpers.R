rna_counts <- function(object) {
  if (!"RNA" %in% names(object@assays)) stop("RNA counts are required")
  assay <- object@assays$RNA
  if ("counts" %in% methods::slotNames(assay)) return(assay@counts)
  layers <- SeuratObject::Layers(assay, search = "^counts($|\\.)")
  if (!length(layers)) stop("RNA counts layer is absent; normalized data cannot replace counts")
  if (length(layers) > 1L) {
    assay <- SeuratObject::JoinLayers(assay, layers = "counts", new = "counts")
    return(SeuratObject::LayerData(assay, layer = "counts"))
  }
  SeuratObject::LayerData(assay, layer = layers[1])
}

validate_counts <- function(x) {
  if (!nrow(x) || !ncol(x) || is.null(rownames(x)) || is.null(colnames(x)) ||
      anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
    stop("Counts require nonempty, unique feature and cell IDs")
  values <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
  if (any(!is.finite(values)) || any(values < 0)) stop("Counts must be finite and nonnegative")
  invisible(x)
}

linear_abundance <- function(counts) {
  validate_counts(counts)
  totals <- Matrix::colSums(counts)
  # Normalize before splitting to keep each cell's full transcript-library denominator.
  scale <- numeric(length(totals))
  scale[totals > 0] <- 1e6 / totals[totals > 0]
  result <- counts %*% Matrix::Diagonal(x = scale)
  dimnames(result) <- dimnames(counts)
  result
}

dominant_call <- function(values, transcripts) {
  best <- max(values)
  if (!is.finite(best) || best <= 0)
    return(list(DominantTranscript = NA_character_, dominance_status = "unexpressed", n_tied = 0L))
  winners <- which(values == best)
  list(DominantTranscript = if (length(winners) == 1L) transcripts[winners] else NA_character_,
       dominance_status = if (length(winners) == 1L) "unique" else "tie",
       n_tied = length(winners))
}
