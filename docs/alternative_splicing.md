# Alternative-splicing input and output conventions

SUPPA expression parts are built exclusively from RNA counts, supporting Seurat v4 counts and
Seurat v5 joined or split counts layers. Each cell is scaled by its full retained transcript count
library to one million **before** partitioning. Values are linear, without logarithms or additional
transcript-length adjustment; this is count-based per-million abundance, not length-adjusted TPM.
Cells with an empty library retain zero abundance. Normalized RNA `data` is never a substitute
when counts are absent. Transcript IDs are preserved for matching to the generated IOE.

The splitter publishes `split_files/input_provenance.json`, `cell_libraries.tsv` and expression
parts. Per-cell undefined PSI is zero-filled when parts are merged, preserving the previous
workflow convention. `undefined_psi_by_sample.tsv` separately records the original undefined
count for each event and sample. Missing PSI and measured exon exclusion are different states;
the zero-filled test must be interpreted with that limitation.

## Differential results

| Parameter | Default | Effect |
| --- | --- | --- |
| `ase_differential_enabled` | `true` | Run each configured comparison on merged PSI |
| `ase_fdr` | `0.05` | Strict BH-adjusted P-value cutoff |
| `ase_fold_change` | `1.5` | Strict absolute log2 mean-ratio cutoff, expressed as a fold change |
| `ase_pseudocount` | `1e-6` | Added to both group means for the ratio |
| `ase_min_cell_prop` | `0.01` | Detection threshold for descriptive plots only |

Results appear in `AlternativeSplicing/Differential/<comparison_name>/`: `all_events.tsv`,
`candidates.tsv` and `statistics.json`. Positive `log2FC` and `delta_PSI` indicate group1 higher
than group2. All input events enter a two-sided unpaired Wilcoxon rank-sum test with R's
small-sample default and tie/continuity-corrected normal approximation otherwise; BH adjustment
is applied across the full input event family. Constant, untestable events retain an undefined
P value and cannot be selected. Plot detection thresholds do not prefilter this test family.

The statistics reader streams event rows rather than loading a second complete PSI matrix.
Zero multiplicities are used to calculate ranks without sorting redundant zeros; the numerical
tests compare this implementation with `stats::wilcox.test`, including ties and degenerate data.

This is exploratory pooled-cell analysis. It does not model donor pairing or make cells
independent biological replicates. No donor-specific examples or fixed study results are embedded
in the pipeline. Disable this additional step with `--ase_differential_enabled false`.

## Dominance and conditional detection fractions

Only a unique positive maximum produces a dominant-transcript call. Tied maxima have no
assigned winner and are retained in `dominance_calls.tsv`; `dominance_by_cell.tsv` counts unique,
tied and unexpressed multi-transcript genes per cell. Zero gene denominators remain `NA` in the
usage table. The legacy `dominant_transcripts_per_cell.txt` contains unique calls only.

The comparison table's existing `<group>_pct` columns describe transcript detection among cells
expressing **any quantified transcript of the gene**, including transcripts that never win a
dominant call. They are conditional detection fractions, not counts usage or the fraction of cells
in which that transcript dominates. Legacy `type`/`notSig` codes remain for compatibility; they
represent the configured effect threshold and are not statistical significance. `evaluable` and
`selection_basis` make that distinction explicit. Ties and the corrected denominator can change
results from previous versions.

## Rechecking the implementation

```bash
Rscript tests/test_splicing.R .
Rscript tests/test_dominance_integration.R .
python3 -m unittest discover -s tests -p test_pipeline_updates.py -v
# Requires native SUPPA's suppa.py and Nextflow on PATH; uses six synthetic cells.
python3 tests/run_splicing_smoke.py --output /tmp/new-splicing-smoke --rscript /path/to/Rscript
```
