# Study filtered-cell QC report

This is the original English QC report prepared on **22 October 2026** from the
study's filtered single-cell long-read RNA data: **AD1–AD3 and N1–N3**, representing
three paired donors and **53,489 retained cells**. It is preserved as a completed
study-report snapshot. The compact public-data demos use a different dataset.

## Open the report

1. [Download the standalone HTML](https://github.com/lishenglilab/scLRA-LLM/raw/refs/heads/main/docs/reports/filtered-cell-qc-20261022/scLRA_LLM_QC_report.html).
2. Open the downloaded file in a browser. It works offline; no server, account or
   API key is required.

GitHub's normal [HTML file view](scLRA_LLM_QC_report.html) displays source code.
Use the download link above to view the interactive report. In a local checkout,
open `docs/reports/filtered-cell-qc-20261022/scLRA_LLM_QC_report.html` directly.

The report embeds its sample data, styles, interaction code and distribution figure.
It supports sample filtering, metric tabs, definitions, figure enlargement, export
of all 42 metric columns for all six samples, and browser printing.
The same full table is also available as [a TSV file](scLRA_LLM_QC_report.tsv).

## Measurement scope

- Read identifiers are distinct BAM-derived read IDs per retained cell; the sample
  totals sum to **492,206,930**.
- Gene/transcript counts use the corresponding linear RNA matrices. UMI strings
  are pooled across genes within a cell and are distinct from gene-aware molecule counts.
- `double.cell` labels were inherited from matched filtered Illumina metadata.
  All labels in this retained cohort are False; this does not measure the residual
  doublet rate.
- Assessments describe this run. Structural classes and aggregate correlations do
  not establish isoform accuracy. Detailed definitions accompany the report.

## Preservation and reproducibility

The HTML and TSV are byte-for-byte copies of the original approved report files;
the embedded numbers and assessments have not been regenerated for publication.
[provenance.json](provenance.json) records their SHA256 hashes, cohort and publication
checks. [SHA256SUMS](SHA256SUMS) can be checked from this directory with
`sha256sum -c SHA256SUMS`.

This snapshot is distinct from a new execution of the current pipeline. To generate
a report for another dataset, follow the main
[filtered-cell QC instructions](../../../README.md#filtered-cell-qc-report), or run
the [compact demos](../../../demo/README.md).
