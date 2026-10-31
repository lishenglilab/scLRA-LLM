# Changes

## Compact demos, 2026-10-31

- Add a matrix entry point that starts from matched, unfiltered IsoQuant gene/transcript
  counts and runs the existing downstream analysis processes.
- Provide versioned matrix and cell-focused FASTQ downloads, SHA256 verification,
  source provenance, reconstruction utilities and output validation commands.
- Add demo profiles with normal Seurat filtering and clustering settings, disabled
  LLM calls and explicitly artificial technical groups from one public specimen.
  Replicate-based isoform-switch modelling is disabled in these profiles.
- Keep sample IDs paired with their FASTQs when upstream tasks finish in different orders.
- Avoid Python bytecode writes in the staged source directory to keep cache inputs stable.
- Report a clear error when sparse input has no cells passing the configured filters.
- Allow matrix inputs to provide matched ORF sequences from the FASTQ run, while retaining
  ORF prediction as a fallback when `ORF.fa` is absent.
- Pin NanoPlot to 1.44.1 and allow interactive-only output with `nanoplot_static=false`
  to avoid static-export browser requirements in demo runs.

## Pipeline update, 2026-10-30

- Build SUPPA inputs from linear RNA count-based per-million abundance, before splitting cells;
  retain Seurat v4/v5 compatibility and record library denominators and input provenance.
- Add configurable pooled-cell AS comparisons across all events, with two-sided Wilcoxon tests,
  BH adjustment and explicit thresholds. Retain zero-filled PSI and publish original undefined
  counts by sample. See [splicing conventions](docs/alternative_splicing.md) for interpretation.
- Exclude tied maxima from unique dominant calls and publish tie/unexpressed summaries.
  Use all quantified gene transcripts in conditional detection denominators. Keep the legacy
  effect codes with explicit selection metadata; these codes are not P values.
- Generate a standalone filtered-cell QC report with HTML/TSV/JSON outputs, optional aggregate-only
  LLM interpretation, preserved numerical metrics and failure fallback. Add a QC dependency
  environment and synthetic numerical/integration checks.
- Annotate clusters in independent requests using a versioned prompt, retain raw answers and
  qualitative confidence grades, and include the whole gene intersection in functional reports.
  Remove unsupported significance/causality claims from interpretation instructions.
- Reserve the underscore for the sample/barcode delimiter at input validation. Sample IDs may
  contain letters, numbers, periods and hyphens; rename underscored sample IDs and update their
  mappings before running this version.

These changes affect PSI, dominance tables and optional LLM results. Use a new output directory
for comparison with an earlier run. Nextflow invalidates tasks whose staged scripts/inputs change;
review the new provenance and do not mix old analysis tables with new outputs.

Validation uses synthetic data, native SUPPA and local mock LLM responses. It does not establish
full-data biological equivalence or an end-to-end FASTQ production run. The new environment YAML
specifies dependencies; a fresh Conda solve is separate from these installed-environment checks.
