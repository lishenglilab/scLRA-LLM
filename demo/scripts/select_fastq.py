#!/usr/bin/env python3
"""Freeze a small cell-focused subset of public raw FASTQ, without editing reads.

Uses scisorseqr barcode assignments only for subset selection; the released reads
must still pass the pipeline's independent BLAZE identification. Group A/B are
deterministic technical splits of one public specimen, not biological conditions.
"""
import argparse
import collections
import gzip
import hashlib
import json
from pathlib import Path


def assignments(path):
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5:
                yield fields[0].lstrip("@"), fields[4]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", required=True)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cells", type=int, default=120)
    parser.add_argument("--reads-per-cell", type=int, default=2000)
    parser.add_argument("--source-compressed-md5", help="Optional independently measured MD5 of the local source archive")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with open(args.barcodes) as handle:
        barcodes = [line.strip() for line in handle]
    valid = {barcode[:-2] if barcode.endswith("-1") else barcode for barcode in barcodes}
    counts = collections.Counter(bc for _, bc in assignments(args.assignments) if bc in valid)
    selected = sorted(counts, key=lambda bc: (-counts[bc], bc))[:args.cells]
    if len(selected) != args.cells or min(counts[bc] for bc in selected) < args.reads_per_cell:
        raise SystemExit("Insufficient cells/depth in source barcode assignments")
    # Alternate a hash order, avoiding ranking-derived technical groups.
    ordered = sorted(selected, key=lambda bc: hashlib.sha256(bc.encode()).hexdigest())
    groups = {bc: "DemoA" if i % 2 == 0 else "DemoB" for i, bc in enumerate(ordered)}
    lookup = {}
    for qname, bc in assignments(args.assignments):
        if bc in groups:
            lookup[qname] = bc
    print("Selected", len(selected), "cells;", len(lookup), "source read assignments", flush=True)
    retained = collections.Counter()
    outputs = {name: gzip.open(args.output / (name + ".fastq.gz"), "wt") for name in ("DemoA", "DemoB")}
    id_outputs = {name: gzip.open(args.output / (name + ".read_ids.txt.gz"), "wt") for name in outputs}
    raw_hashes = {name: hashlib.sha256() for name in outputs}
    canonical_hashes = {name: hashlib.sha256() for name in outputs}
    sample_counts = collections.Counter()
    seen = 0
    try:
        with gzip.open(args.fastq, "rt") as handle:
            while True:
                record = [handle.readline() for _ in range(4)]
                if not record[0]:
                    break
                if not record[3]:
                    raise ValueError("Truncated source FASTQ")
                seen += 1
                qname = record[0][1:].split()[0]
                bc = lookup.get(qname)
                if bc is not None and retained[bc] < args.reads_per_cell:
                    sample = groups[bc]
                    outputs[sample].writelines(record)
                    id_outputs[sample].write(qname + "\n")
                    raw_hashes[sample].update("".join(record).encode())
                    normalized = "@" + qname + "\n" + record[1].rstrip("\r\n") + "\n+\n" + record[3].rstrip("\r\n") + "\n"
                    canonical_hashes[sample].update(normalized.encode())
                    sample_counts[sample] += 1
                    retained[bc] += 1
                if seen % 1000000 == 0:
                    print("Scanned", seen, "retained", sum(retained.values()), flush=True)
                if len(retained) == len(selected) and min(retained.values()) >= args.reads_per_cell:
                    break
    finally:
        for output in outputs.values():
            output.close()
        for output in id_outputs.values():
            output.close()
    if len(retained) != len(selected) or min(retained.values()) < args.reads_per_cell:
        raise SystemExit("Raw FASTQ did not supply the requested reads for every selected cell")
    (args.output / "selected_cells.tsv").write_text("sample\tbarcode\tsource_assigned_reads\tretained_reads\n" +
        "".join(f"{groups[bc]}\t{bc}\t{counts[bc]}\t{retained[bc]}\n" for bc in ordered))
    (args.output / "selection.json").write_text(json.dumps({
        "accession": "SRR14808743", "series": "GSE178175", "specimen": "FCtx1_ONT",
        "selection": "top source barcode counts among published 10x cell barcodes; tie-break lexicographic",
        "grouping": "alternating SHA256 order; two artificial technical groups from one specimen",
        "read_sampling": "first requested number of raw reads per selected cell in source FASTQ order",
        "cells": len(selected), "reads_per_cell": args.reads_per_cell,
        "scanned_reads": seen, "retained_reads": sum(retained.values()),
        "read_integrity": "source FASTQ record text preserved, including original sequence and qualities",
        "source_ena_fastq_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR148/043/SRR14808743/SRR14808743_1.fastq.gz",
        "source_ena_fastq_md5": "8f56e4d9af0a2447c47c8a64b0ae9165",
        "local_sra_dump_compressed_md5": args.source_compressed_md5,
        "header_format": "Local SRA FASTQ-dump headers; ENA header comments differ; canonical hash normalizes comments and + lines",
        "subsets": {name: {"records": sample_counts[name], "uncompressed_record_sha256": raw_hashes[name].hexdigest(),
                           "canonical_record_sha256": canonical_hashes[name].hexdigest()} for name in outputs}
    }, indent=2) + "\n")
    print("Selection complete", flush=True)


if __name__ == "__main__":
    main()
