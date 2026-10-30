#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(optparse); library(jsonlite)})

zero_wilcox <- function(a, b, na, nb) {
  na <- as.double(na); nb <- as.double(nb)
  if (na < 50L && nb < 50L) return(suppressWarnings(
    wilcox.test(c(rep(0, na - length(a)), a), c(rep(0, nb - length(b)), b))$p.value))
  nzero <- na + nb - length(a) - length(b)
  values <- c(a, b)
  n <- as.double(na + nb)
  statistic <- (na - length(a)) * (nzero + 1) / 2 +
    sum(rank(values)[seq_along(a)]) + length(a) * nzero - na * (na + 1) / 2
  ties <- c(nzero, tabulate(match(values, unique(values))))
  variance <- (as.double(na) * nb / 12) * (n + 1 - sum(ties^3 - ties) / (n * (n - 1)))
  if (!is.finite(variance) || variance <= 0) return(NA_real_)
  z <- statistic - as.double(na) * nb / 2
  z <- (z - sign(z) * 0.5) / sqrt(variance)
  min(1, 2 * pnorm(abs(z), lower.tail = FALSE))
}

run_statistics <- function(opt) {
  if (!is.finite(opt$fdr) || opt$fdr <= 0 || opt$fdr > 1 ||
      !is.finite(opt$fold_change) || opt$fold_change <= 1 ||
      !is.finite(opt$pseudocount) || opt$pseudocount <= 0) stop("Invalid AS thresholds")
  mapping <- fread(opt$stage_mapping)
  if (!all(c("sample", "stage") %in% names(mapping)) || anyDuplicated(mapping$sample))
    stop("Stage mapping requires unique sample IDs and sample/stage columns")
  input <- file(opt$psi_file, "rt")
  on.exit(close(input))
  header <- strsplit(readLines(input, n = 1L), "\t", fixed = TRUE)[[1]]
  if (length(header) < 3L || header[1] != "ASE" || anyDuplicated(header)) stop("Invalid PSI header")
  samples <- sub("_[^_]+$", "", header[-1L])
  if (any(!samples %in% mapping$sample)) stop("PSI cells have unmapped samples")
  groups <- mapping$stage[match(samples, mapping$sample)]
  ix1 <- which(groups == opt$group1); ix2 <- which(groups == opt$group2)
  if (opt$group1 == opt$group2 || !length(ix1) || !length(ix2)) stop("Both distinct groups need cells")
  result <- list()
  seen <- new.env(hash = TRUE, parent = emptyenv())
  repeat {
    lines <- readLines(input, n = 128L)
    if (!length(lines)) break
    for (line in lines) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) != length(header)) stop("PSI row/header mismatch")
      event <- fields[1]
      if (!nzchar(event) || exists(event, envir = seen, inherits = FALSE)) stop("Duplicate/empty event ID")
      assign(event, TRUE, envir = seen)
      values <- suppressWarnings(as.numeric(fields[-1L]))
      if (any(!is.finite(values)) || any(values < 0 | values > 1))
        stop("Statistics require merged, zero-filled PSI in [0, 1]")
      a <- values[ix1]; b <- values[ix2]
      ma <- mean(a); mb <- mean(b)
      p <- zero_wilcox(a[a > 0], b[b > 0], length(a), length(b))
      result[[length(result) + 1L]] <- data.table(ASE = event,
        mean_group1 = ma, mean_group2 = mb, delta_PSI = ma - mb,
        log2FC = log2((ma + opt$pseudocount) / (mb + opt$pseudocount)), p_value = p)
    }
  }
  if (!length(result)) stop("PSI table contains no events")
  table <- rbindlist(result)
  table[, FDR := p.adjust(p_value, method = "BH")]
  table[, selected := !is.na(FDR) & FDR < opt$fdr & abs(log2FC) > log2(opt$fold_change)]
  table[, direction := fifelse(!selected, "not_selected", fifelse(log2FC > 0, opt$group1, opt$group2))]
  dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(table, file.path(opt$output_dir, "all_events.tsv"), sep = "\t", na = "NA")
  fwrite(table[selected == TRUE], file.path(opt$output_dir, "candidates.tsv"), sep = "\t", na = "NA")
  write_json(list(group1 = opt$group1, group2 = opt$group2,
    cells_group1 = length(ix1), cells_group2 = length(ix2), tested_events = nrow(table),
    selected_events = sum(table$selected), test = "two-sided unpaired Wilcoxon rank-sum",
    exact = "R default for small samples; tie-corrected normal approximation otherwise",
    correction = "BH across all input events; no detection-proportion prefilter",
    undefined_PSI = "zero-filled by merge; undefined_psi_by_sample.tsv preserves missing counts",
    limitation = "pooled-cell exploration; cells are not independent donors; zero filling conflates missing PSI and exon exclusion",
    fdr = opt$fdr, fold_change = opt$fold_change, pseudocount = opt$pseudocount,
    threshold_comparators = "strict < FDR and > absolute log2FC"),
    file.path(opt$output_dir, "statistics.json"), auto_unbox = TRUE, pretty = TRUE)
}

if (sys.nframe() == 0L) {
  parser <- OptionParser(option_list = list(
    make_option("--psi_file", type = "character"), make_option("--stage_mapping", type = "character"),
    make_option("--group1", type = "character"), make_option("--group2", type = "character"),
    make_option("--output_dir", type = "character", default = "ase_statistics"),
    make_option("--fdr", type = "double", default = 0.05),
    make_option("--fold_change", type = "double", default = 1.5),
    make_option("--pseudocount", type = "double", default = 1e-6)))
  opt <- parse_args(parser)
  for (key in c("psi_file", "stage_mapping", "group1", "group2"))
    if (is.null(opt[[key]]) || !nzchar(opt[[key]])) stop("Missing --", key)
  run_statistics(opt)
}
