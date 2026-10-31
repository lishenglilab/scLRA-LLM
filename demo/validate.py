#!/usr/bin/env python3
"""Check actual demo outputs and (optionally) its Nextflow trace."""
import argparse
import csv
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--trace", type=Path)
    args = parser.parse_args()
    required = ["PreliminarySingleCellAnalysis/sce.gene.rds",
                "PreliminarySingleCellAnalysis/sce.tr.rds",
                "PreliminarySingleCellAnalysis/qc_report/scLRA_LLM_QC_report.html",
                "PreliminarySingleCellAnalysis/qc_report/scLRA_LLM_QC_report.tsv",
                "AlternativeSplicing/merged_project_events.psi",
                "AlternativeSplicing/undefined_psi_by_sample.tsv",
                "AlternativeSplicing/split_files/input_provenance.json",
                "AlternativeSplicing/DominantIsoform/dominant_isoform_output/dominance_calls.tsv",
                "IntegrativeAnalysis/ORF_Cluster/ORF_Cluster/sce.orf.rds"]
    comparison = "TechnicalA_vs_TechnicalB"
    required += [f"IntegrativeAnalysis/{comparison}/{comparison}/Differential_{kind}_{comparison}.txt"
                 for kind in ("Gene_Expression", "Transcript_Expression", "Gene_Usage", "Dominant_Transcript_Usage")]
    for relative in required:
        path = args.results / relative
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit("Missing/empty required output: " + str(path))
    diff = args.results / "AlternativeSplicing/Differential/TechnicalA_vs_TechnicalB"
    statistics = json.loads((diff / "statistics.json").read_text())
    with (diff / "all_events.tsv").open() as handle:
        events = list(csv.DictReader(handle, delimiter="\t"))
    ids = [row["ASE"] for row in events]
    if not events or len(ids) != len(set(ids)):
        raise SystemExit("AS events absent or duplicated")
    selected = []
    for row in events:
        if not all(0 <= float(row[key]) <= 1 for key in ("mean_group1", "mean_group2")):
            raise SystemExit("Invalid PSI means")
        expected = row["FDR"] != "NA" and float(row["FDR"]) < statistics["fdr"] and abs(float(row["log2FC"])) > math.log2(statistics["fold_change"])
        if (row["selected"] == "TRUE") != expected:
            raise SystemExit("AS threshold/selection disagreement")
        if expected:
            selected.append(row["ASE"])
    with (diff / "candidates.tsv").open() as handle:
        actual = [row["ASE"] for row in csv.DictReader(handle, delimiter="\t")]
    if selected != actual:
        raise SystemExit("Candidate list disagrees with all-event table")
    statuses = []
    if args.trace:
        with args.trace.open() as handle:
            tasks = list(csv.DictReader(handle, delimiter="\t"))
        statuses = [row["status"] for row in tasks]
        if not tasks or any(status not in ("COMPLETED", "CACHED") for status in statuses):
            raise SystemExit("Trace contains incomplete or failed tasks")
    summary = {"required_outputs": len(required), "AS_events": len(events),
               "selected_AS_events": len(selected), "trace_tasks": len(statuses),
               "scope": "Software demo; artificial technical groups, not biological inference"}
    (args.results / "demo_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
