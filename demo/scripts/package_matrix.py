#!/usr/bin/env python3
"""Create the matrix demo from a completed FASTQ-demo run, retaining count libraries."""
import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--isoquant-dir", type=Path, help="Optional successful IsoQuant_sc work output instead of the published copy")
    parser.add_argument("--reference-gtf", required=True, type=Path)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    alignment = args.results / "AlignmentAndTranscriptReconstruction"
    source = (args.isoquant_dir or alignment / "isoquant" / "IsoQuant_sc") / "OUT"
    destination = args.output / "IsoQuant_sc" / "OUT"
    destination.mkdir(parents=True)
    for kind in ("gene", "transcript"):
        prefix = "OUT.discovered_" + kind + "_grouped_counts."
        for extension in ("barcodes.tsv", "features.tsv", "matrix.mtx"):
            shutil.copyfile(source / (prefix + extension), destination / (prefix + extension))
    classification = source / "OUT.novel_vs_known.SQANTI-like.tsv"
    shutil.copyfile(classification, destination / classification.name)
    shutil.copyfile(alignment / "all_read_groups.tsv", args.output / "all_read_groups.tsv")
    shutil.copyfile(args.selection.with_name("selected_cells.tsv"), args.output / "selected_cells.tsv")
    models = (alignment / "orf/OUT.transcript_models.gtf").read_bytes()
    public_models = b"".join(line for line in models.splitlines(keepends=True)
                              if not line.startswith(b"# Command line:"))
    (args.output / "transcripts.gtf").write_bytes(public_models)
    shutil.copyfile(alignment / "orf/OUT.transcript_models.1line.fa", args.output / "transcripts.fa")
    shutil.copyfile(alignment / "orf/ORF.fa", args.output / "ORF.fa")
    genes = set(re.findall(r'gene_id "([^"]+)"', (args.output / "transcripts.gtf").read_text()))
    with args.reference_gtf.open() as reference, (args.output / "reference.gtf").open("w") as output:
        for line in reference:
            match = re.search(r'gene_id "([^"]+)"', line)
            if match and match.group(1) in genes:
                output.write(line)
    provenance = json.loads(args.selection.read_text())
    provenance["matrix_origin"] = "Unfiltered discovered raw-count matrices from this FASTQ demo's IsoQuant run"
    provenance["feature_selection"] = "No further count-feature filtering or rescaling; entire retained count libraries"
    provenance["reference_gtf"] = "GENCODE v48, subset to genes present in transcript models; genome coordinates preserved"
    provenance["read_support"] = "Filtered FASTQ read groups before alignment; BAMs are not bundled"
    provenance["transcript_models"] = {
        "original_sha256": hashlib.sha256(models).hexdigest(),
        "bundled_sha256": hashlib.sha256(public_models).hexdigest(),
        "comments": "Machine-local command-line comment omitted; all annotation records preserved"}
    provenance["orf_sequences"] = {
        "origin": "ORF prediction from these same transcript model sequences in the FASTQ run",
        "sha256": hashlib.sha256((args.output / "ORF.fa").read_bytes()).hexdigest()}
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Matrix demo packaged:", args.output)


if __name__ == "__main__":
    main()
