# scLRA-LLM

scLRA-LLM is a Nextflow DSL2 pipeline for single-cell long-read RNA-seq analysis. It connects read preprocessing, IsoQuant transcript reconstruction, gene- and transcript-level single-cell analysis, alternative-splicing analysis, ORF analysis, and optional LLM-assisted reporting in one workflow.

LLM steps use an OpenAI-compatible Chat Completions endpoint and are not tied to a specific provider. All biological analysis steps can run with LLM support disabled.

## Workflow

| Module | Main tasks |
| --- | --- |
| `InputAndPreprocess` | Metadata validation, BLAZE barcode processing, NanoFilt filtering, NanoPlot QC, optional LLM read-QC report |
| `AlignmentAndTranscriptReconstruction` | IsoQuant transcript reconstruction and quantification, transcript sequence extraction, ORF prediction, optional LLM assembly-QC report |
| `PreliminarySingleCellAnalysis` | Gene/transcript matrices, Seurat QC and clustering, marker detection, English HTML QC report from filtered cells, optional LLM annotation and QC interpretation |
| `AlternativeSplicing` | SUPPA2 event quantification, dominant-isoform analysis, isoform-switch preparation and consequences, optional LLM interpretation |
| `IntegrativeAnalysis` | DDTU, DTE, DGE and DGU integration, ORF clustering, optional LLM interpretation |

The five modules are orchestrated from [`main.nf`](main.nf), with process definitions under [`modules/`](modules/).

The current update adds linear count-based SUPPA inputs, optional pooled-cell AS differential
analysis, explicit dominance ties, and filtered-cell QC reports. See [CHANGELOG.md](CHANGELOG.md)
for changes that affect cached results and input conventions.

## Requirements

- Linux
- Nextflow with DSL2 support; the current workflow was validated with Nextflow 25.10.2
- A Java version supported by the installed Nextflow release (Java 17 or newer for current Nextflow releases)
- Conda or Mamba for the process environments
- Sufficient CPU, memory and storage for IsoQuant and single-cell analysis

The pipeline provides Conda YAML files for its Python and command-line dependencies. The main environment currently does not fully lock the R/Bioconductor stack; see [R environment](#r-environment).

## Installation

```bash
git clone https://github.com/lishenglilab/scLRA-LLM.git
cd scLRA-LLM
unzip assets/3M-february-2018.zip -d assets
```

The archive contains the 10x 3M barcode whitelist referenced by the default configuration. If a different chemistry or whitelist is used, update `params.whitelist` instead. No API credentials are required when LLM support is disabled.

## Input configuration

Before running, edit [`nextflow.config`](nextflow.config) and check at least:

- `params.metadata`
- `params.whitelist`
- `params.reference_fa`
- `params.genedb_gtf`
- `params.sample_stage_mapping`
- `params.comparisons`
- Thread counts and analysis thresholds appropriate for the server

### Metadata

Metadata may be TSV or CSV and must contain `sample_id` and `fastq` columns. Absolute FASTQ paths are recommended.

```text
sample_id<TAB>fastq
AD1<TAB>/absolute/path/AD1.fastq
AD2<TAB>/absolute/path/AD2.fastq
N1<TAB>/absolute/path/N1.fastq
```

Sample IDs must be unique and may contain letters, numbers, `.` and `-`. The workflow uses `_` to separate sample IDs from cell barcodes.

### Groups and comparisons

Each sample must map to a biological group:

```groovy
sample_stage_mapping = [
  ["AD1", "AD"],
  ["AD2", "AD"],
  ["N1",  "Normal"]
]
```

Comparisons use `[group1, group2, comparison_name]`:

```groovy
comparisons = [
  ["AD", "Normal", "AD_vs_Normal"]
]
```

Comparison names must be unique and may contain letters, numbers, `.`, `_` and `-`. Both groups must occur in `sample_stage_mapping`.

## Provider-neutral LLM configuration

The optional LLM processes use the OpenAI Python SDK against an OpenAI-compatible Chat Completions API. Configure credentials with environment variables:

```bash
export LLM_BASE_URL='https://provider.example/v1'
export LLM_API_KEY='your-key'
export LLM_MODEL='provider-model-name'
```

- `LLM_BASE_URL` and `LLM_API_KEY` identify and authenticate the endpoint.
- `LLM_MODEL` is optional for backward compatibility, but should normally be set for a non-DeepSeek endpoint. OpenAI-compatible APIs require a model string and cannot reliably infer one from a multi-model URL.
- If `LLM_API_KEY` is absent, LLM steps are disabled automatically and all non-LLM analyses continue.
- Use `--llm_enabled false` to disable LLM calls explicitly.
- Enabling LLM without a key or a valid absolute HTTP(S) base URL fails during workflow validation.

The API key is deliberately not stored in `nextflow.config`, Nextflow parameters, task command lines, metadata or committed `.env` files. Do not pass real credentials on the command line.

Some output directories and JSON filenames retain historical `deepseek_*` names so existing downstream consumers do not break. These names do not restrict the configured provider; provider-neutral CLI entry points live in [`bin/llm/`](bin/llm/).

Cell annotation uses a separate request for each cluster with the ranked marker list and a
versioned [prompt](bin/llm/prompts/cell_annotation.txt). The CSV retains a qualitative confidence
grade from 0 to 3; grade 0 is a forced best estimate, and must not be treated as a resolved
annotation or a calibrated probability. Requests, raw answers and status are retained in
`cell_ann_by_llm.requests.json`; the gene object also retains `cell_anno_llm_confidence`.
The maximum marker count remains 30 by default. Changing the prompt or request grouping can
change labels; outputs from the preceding batch prompt are not equivalent runs.

Integrative interpretation includes all intersection genes. DGU/DDTU effect labels are described
separately from expression-test P values, and the prompt allows abstention when a unifying
function is not supported. Large intersections require a model context window that can hold the
full request; no hidden top-gene truncation is applied. LLM annotations and interpretations remain
optional and are not independent validation of the biological results.

## Running the pipeline

Select the compatible R interpreter when more than one R installation is available:

```bash
export RSCRIPT_BIN='/usr/bin/Rscript'
"$RSCRIPT_BIN" --version
```

Run the full workflow:

```bash
nextflow run main.nf -resume
```

Run without any LLM requests:

```bash
nextflow run main.nf -resume --llm_enabled false
```

Useful structure-only checks:

```bash
NXF_OFFLINE=true nextflow lint main.nf nextflow.config modules/*.nf
NXF_OFFLINE=true nextflow run main.nf -preview --llm_enabled false
```

Nextflow caching is enabled. Keep `work/` and use `-resume` to reuse successfully completed tasks.

## Output layout

The default output root is `results/`:

```text
results/
├── InputAndPreprocess/
├── AlignmentAndTranscriptReconstruction/
├── PreliminarySingleCellAnalysis/
├── AlternativeSplicing/
└── IntegrativeAnalysis/
```

Each module publishes its reports, tables, plots and reusable objects below its corresponding directory. Large intermediate files remain in the Nextflow `work/` directory.

### Filtered-cell QC report

The concise English report is generated automatically from the final filtered `sce.gene.rds` and `sce.tr.rds` objects in this run:

```text
results/PreliminarySingleCellAnalysis/qc_report/
  scLRA_LLM_QC_report.html      # Standalone HTML; open locally without a web server
  scLRA_LLM_QC_report.tsv       # All sample metrics, with Assessment last
  sc_qc_report.txt             # Short text summary
  data/                       # Count-derived metrics, cell table, JSON and R-rendered SVG
  llm/                        # Aggregate request, raw response and status, when enabled
```

The report retains the approved compact layout: overall assessment, up to two review priorities, a sample-by-metric table, a distribution figure, and tabs for read identifiers/labels, explicit correlation pairs, and transcript models. Definitions are collapsed. No study-specific sample names, cell counts or conclusions are embedded in the template.

Metrics use linear RNA counts from the same ordered cell IDs. MT is recomputed from mitochondrial features, resolving reference GTF IDs when necessary. Existing valid `double.cell` records define the report cohort when that field exists; TRUE records are reported without another filtering step. Otherwise, `predicted_doublets` is shown only when detection was actually performed. Skipped detection is not reported as a negative result. Unavailable measurements and undefined correlations are omitted rather than replaced by zero.

Read identifiers default to distinct query names in this run's IsoQuant BAMs, matched to retained cells through `all_read_groups.tsv`. All stored alignment records, including unmapped records, contribute identifiers; repeated alignments of the same ReadID are counted once per cell. If no BAM is retained, the report explicitly states that identifiers come from filtered FASTQs **before alignment**. UMI-string counts require recognizable BLAZE headers and remain distinct from gene-aware molecule counts. Disk sorting bounds memory use; large runs need temporary storage in the Nextflow work directory.

Structural annotation uses the available IsoQuant SQANTI-like table. An IsoQuant novel-model table can provide partial coverage; the report shows that coverage and uses **all detected transcript features** as the class-percentage denominator. Supply a complete classification table to use a wider annotation snapshot. Neither structural classes nor within-platform correlations establish isoform accuracy.

| Parameter | Default | Effect |
| --- | --- | --- |
| `qc_report_enabled` | `true` | Generate the filtered-cell QC report |
| `qc_project` | `scLRA-LLM` | Report project label |
| `qc_read_source` | `auto` | Prefer BAMs; `bam` requires BAMs; `fastq` uses pre-alignment read groups |
| `qc_read_summary` | empty | Override with a BAM-derived per-cell TSV: `CellBarcode`, `unique_reads`, optional `unique_UMIs`; exact retained-cell coverage is required |
| `qc_classification` | empty | Override with a SQANTI classification TSV, headered or standard headerless format |
| `qc_llm_strict` | `false` | Fail on QC interpretation errors instead of keeping the metric-based report |

With LLM enabled, only aggregate metrics and their definitions are sent to the configured endpoint. The model supplies short English assessments; it cannot overwrite numerical metrics. Requests, raw responses and status are retained for review, without API credentials. With LLM disabled, the same report is built with metric-based descriptions. Provider failures or malformed responses also retain the numerical report and clearly label interpretation as unavailable, unless `qc_llm_strict` is enabled. The historical `deepseek_qc/sc_qc_report.txt` path remains available when LLM interpretation runs.

```bash
# Standard run, including the report, with all LLM calls disabled
nextflow run main.nf -resume --llm_enabled false

# Match an existing BAM-summary / full-classification snapshot
nextflow run main.nf -resume \
  --qc_read_summary /absolute/path/cell_read_summary.tsv \
  --qc_classification /absolute/path/classification.tsv
```

## Alternative splicing

See [docs/alternative_splicing.md](docs/alternative_splicing.md) for input scaling, differential
parameters, output definitions and migration notes.

## R environment

The new QC processes use the portable [`envs/qc.yml`](envs/qc.yml) environment, including Python 3,
pysam, the LLM client, SeuratObject, Matrix, jsonlite, data.table and the R plotting dependencies.
When `RSCRIPT_BIN` is unset, they use that process environment's `Rscript`. An explicit override
must provide the same R dependencies. This environment is a dependency specification, not a
fully resolved lock file. The standalone smoke tests disable Conda and require installed R/Python
dependencies; they do not validate a fresh environment solve.

The R scripts use packages including Seurat, tidyverse, data.table, harmony, rtracklayer, IsoformSwitchAnalyzeR, DESeq2, muscat, BSgenome.Hsapiens.UCSC.hg38, Biostrings, ggVennDiagram and UpSetR.

The current [`envs/sclong.yml`](envs/sclong.yml) mainly locks Python and command-line packages. For production deployment, create a dedicated R/Bioconductor environment matched to the server's R version and set `RSCRIPT_BIN` to that environment's `Rscript`. Avoid installing packages dynamically during a production run.

## Validation and compatibility

The current update is checked using:

- Nextflow lint for the main workflow, configuration and all five modules
- Seurat v4/v5 linear-count, PSI-ratio, tied-dominance and all-event BH numerical checks
- Python regression checks for QC metrics, independent annotation requests and interpretation input scope
- QC-only Nextflow runs with LLM disabled, a successful local mock response and a malformed-response fallback
- Native SUPPA split/merge and two-comparison Nextflow integration on six synthetic cells

These checks use installed dependencies and local mock LLM responses. They do not rerun a complete
FASTQ-to-results study or a fresh Conda solve, and do not measure biological or model accuracy.

QC-specific checks include Seurat v4/v5 count layers, exact cell alignment, valid label cohorts, skipped doublet detection, recomputed MT, partial structural annotation, distinct BAM/read-group identifiers, sample totals and medians, immutable metrics across interpretation, and local mock-API success/failure. Run the QC regression and the actual Nextflow QC-only tasks with:

```bash
RSCRIPT_BIN=/path/to/Rscript python -m unittest discover -s tests -p test_qc_report.py -v
python tests/run_qc_smoke.py --output /tmp/new-qc-smoke --rscript /path/to/Rscript
```

The QC-only smoke runs use six synthetic cells and a local mock endpoint, covering LLM-disabled, successful interpretation, and malformed-response fallback. They do not rerun alignment or the full biological workflow.

Isoform-switch preparation still uses the first entry in `params.comparisons`; the new AS differential step, dominant-isoform analysis and integrative differential analysis iterate over all comparisons.

## Security

- Never commit API keys, tokens or private endpoint credentials.
- Keep local secret files outside the repository; `.env*` files are ignored by default.
- If a credential has ever been committed, deleting it in a later commit is not sufficient. Revoke and rotate the credential immediately, then clean the Git history if required.

## License

This project is distributed under the terms in [`LICENSE`](LICENSE).
