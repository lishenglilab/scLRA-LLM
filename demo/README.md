# Compact, traceable demos

Two downloadable inputs exercise the same analysis processes as normal runs. Data are
versioned GitHub Release assets, with SHA256 checksums in `manifest.json`; they are not
committed as large binary files in the code repository. Python 3.7+ can download them.

| Input | Compressed download | Entry point | Whole-genome reference |
| --- | --- | --- | --- |
| Matrix | 31.4 MiB | Unfiltered counts and matched upstream models/ORFs | Not required |
| FASTQ | 270 MiB | Original reads through the full preprocessing/alignment path | GRCh38 and GENCODE v48 supplied separately |

## Quick matrix demo

From the repository root, in the configured Linux analysis environment:

```bash
python3 demo/download.py matrix
nextflow run main.nf -profile demo_matrix --llm_enabled false -with-trace trace_matrix.tsv
python3 demo/validate.py results_demo_matrix --trace trace_matrix.tsv
"${RSCRIPT_BIN:-Rscript}" demo/scripts/check_objects.R results_demo_matrix
```

This entry starts with paired, unfiltered raw-count matrices from the FASTQ demo's
IsoQuant run. It runs the usual Seurat filtering and clustering, filtered-cell QC,
SUPPA, AS comparisons, dominance, ORF clustering and integrative analysis. Matched ORF
sequences from the FASTQ run are bundled; matrix mode predicts ORFs if `ORF.fa` is absent.
No genome FASTA or alignment index is required. It does require the normal R analysis
packages, Python/ORF dependencies and SUPPA; downloading data does not install software.
See the repository's [R environment guidance](../README.md#r-environment).

The package includes:

- Gene and transcript count matrices in the pipeline's `IsoQuant_sc/OUT` layout.
- Matched transcript models and sequences, and a reference GTF subset retaining their
  original genomic coordinates. Count features and library denominators are preserved.
- Original partial IsoQuant structural-classification output and filtered FASTQ read
  groups. BAMs are omitted. QC explicitly labels read support as pre-alignment FASTQ
  identifiers and retains the original partial classification coverage.
- A provenance record identifying the source and exact subsetting protocol.

The transcript-model GTF omits the machine-local command-line comment; all annotation
records are retained, with the original and bundled file hashes recorded in provenance.

## Cell-focused FASTQ demo

```bash
python3 demo/download.py fastq
unzip -n assets/3M-february-2018.zip -d assets
nextflow run main.nf -profile demo_fastq --llm_enabled false \
  --reference_fa /absolute/path/GRCh38.primary_assembly.genome.fa \
  --genedb_gtf /absolute/path/gencode.v48.annotation.gtf \
  -with-trace trace_fastq.tsv
python3 demo/validate.py results_demo_fastq --trace trace_fastq.tsv
"${RSCRIPT_BIN:-Rscript}" demo/scripts/check_objects.R results_demo_fastq
```

The downloader writes absolute paths into the local `metadata.tsv`. Both compressed
FASTQs contain unmodified original read records. They exercise BLAZE, read filtering,
NanoPlot, IsoQuant and the downstream analysis entry described above. BLAZE must identify
the barcodes independently; source barcode assignments are used only to select reads.

The data contain 120 selected cells and 240,000 raw reads (2,000 reads/cell), split into
two technical groups with 60 cells each. Reads are downsampled, not complete original
cell libraries. Both demos retain the normal `min_cells=10` and `min_features=200`
settings and the original clustering settings.
Seurat applies its native filtering order; detected features in a final object can fall
below the initial cell threshold after feature filtering. Demo NanoPlot outputs interactive HTML (`nanoplot_static=false`),
so this mode does not require a headless browser for static-image export.

Use the genome sequence and comprehensive annotation from [GENCODE human release 48](https://www.gencodegenes.org/human/release_48.html).
The full human reference is separate from the demo download. Reconstruction uses
substantially more disk/RAM and time than matrix mode.
In the validation run, IsoQuant alone took about 82 minutes with 8 threads and reached
18.5 GB peak RSS as reported by Nextflow for that task. These are server-specific measurements, not runtime or
whole-workflow memory guarantees. Matrix mode avoids this reference-dependent work.
An existing IsoQuant database may be supplied with `--isoquant_annotation_db /path/file.db`
only if built from exactly the same GTF and a compatible IsoQuant version. Continue using
`--genedb_gtf` for the actual reference annotation needed by QC.

## Source and interpretation

The public source is [GSE178175](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178175),
[FCtx1_ONT / GSM5381317](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5381317),
run [SRR14808743](https://www.ebi.ac.uk/ena/browser/view/SRR14808743).
Published matching 10x cell barcodes come from
[GSM5381311](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5381311).
Please acknowledge the original study when reusing the data. The repository's software
license does not replace the original data source's terms.

`DemoA` and `DemoB` are artificial technical groups from **one specimen**: selected barcodes
are alternated after SHA256 ordering. They are neither disease/control conditions nor
independent biological replicates. No biological finding, significant event or specific
cluster assignment is required for success. All measured outputs must still be generated
and internally consistent. LLM calls are disabled by default.

Replicate-based isoform-switch preparation/consequence analysis is explicitly disabled in
both demo profiles (`isoform_switch_enabled=false`). The regular pipeline default remains
enabled. Full-data biological validation and model-provider testing are separate activities.

## Rebuilding and checking

`scripts/select_fastq.py` selects the top 120 cells by scisorseqr read assignments among the
published 10x cell barcodes, breaking ties lexicographically, then retains the first 2,000
raw records per selected barcode in source FASTQ order. It preserves sequences, quality
strings and read headers. Selection details and retained cells accompany the FASTQ asset.
The `DemoA.read_ids.txt.gz` and `DemoB.read_ids.txt.gz` lists also allow exact extraction
from the public source without rerunning the original barcode-assignment software:

```bash
gzip -dc demo/data/fastq/DemoA.read_ids.txt.gz > DemoA.read_ids.txt
seqkit grep -f DemoA.read_ids.txt SRR14808743_1.fastq.gz | gzip > reproduced_DemoA.fastq.gz
```

Check reproduced reads with:

```bash
python3 demo/scripts/check_fastq.py reproduced_DemoA.fastq.gz \
  --selection demo/data/fastq/selection.json --sample DemoA
```

This compares read IDs, sequences, qualities and ordering after removing header comments
and normalizing the `+` separator. ENA and SRA FASTQ dumps use different header comments;
gzip container bytes also differ between encoders. The demo retains the local SRA-style
records unchanged, and records their original uncompressed SHA256 separately.

After running the FASTQ demo, `scripts/package_matrix.py` extracts the exact raw matrices,
models, transcript/ORF sequences and support tables for matrix mode. `--isoquant-dir` may
point directly to the successful IsoQuant process output. It does not reuse a preclustered RDS
or select genes by a differential result. `validate.py` checks required outputs, AS event
uniqueness, PSI ranges, exact candidate thresholds and completed trace statuses;
`check_objects.R` checks ordered matching cells and raw-count integrity.

Actual validation evidence is recorded in [validation.json](validation.json). Full demos
use installed analysis environments. Preprocessing dependencies were also installed fresh
and tested separately; a fresh installation of the remaining R/alignment/splicing stack
is not implied. Older small
`test1.fastq`/`test2.fastq` files remain only as sparse upstream fixtures and are insufficient
for the normal downstream cell filters.
