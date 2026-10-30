#!/usr/bin/env Rscript
# Scientific distributions are rendered in R, using the full report cohort.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: plot_qc.R cells.tsv.gz output.svg")
d <- as.data.frame(fread(args[1]))
d$sample <- factor(d$sample, levels = sort(unique(d$sample)))
metrics <- c(unique_read_ids = "Read IDs / cell", gene_counts = "Gene counts / cell",
  transcript_counts = "Transcript counts / cell", genes_detected = "Detected genes / cell",
  transcripts_detected = "Detected transcripts / cell", mitochondrial_percent = "MT (%) / cell")
metrics <- metrics[names(metrics) %in% names(d)]
metrics <- metrics[vapply(names(metrics), function(k) any(is.finite(d[[k]])), logical(1))]
if (!length(metrics)) stop("No observed cell-level metrics to plot")
palette <- grDevices::hcl.colors(nlevels(d$sample), "Dark 3")
plots <- lapply(names(metrics), function(key) {
  z <- d[, c("sample", key)]
  names(z)[2] <- "value"
  is_mt <- key == "mitochondrial_percent"
  if (!is_mt) z$value <- log10(1 + z$value)
  ggplot(z, aes(sample, value, fill = sample)) +
    geom_violin(scale = "width", linewidth = .25, color = "#687782", na.rm = TRUE) +
    geom_boxplot(width = .15, outlier.shape = NA, fill = "white", linewidth = .3, na.rm = TRUE) +
    scale_fill_manual(values = palette, guide = "none") +
    labs(title = unname(metrics[[key]]), x = NULL, y = if (is_mt) "Percent" else "log10(1 + value)") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 11),
      axis.text.x = element_text(angle = if (nlevels(d$sample) > 6) 45 else 0,
                                 hjust = if (nlevels(d$sample) > 6) 1 else .5))
})
combined <- wrap_plots(plots, ncol = min(3L, length(plots))) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold")))
grDevices::svg(args[2], width = max(10.5, nlevels(d$sample) * .8),
  height = 3.2 * ceiling(length(plots)/3), bg = "white")
print(combined)
invisible(dev.off())
