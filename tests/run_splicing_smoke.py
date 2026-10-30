"""Exercise the actual split, native SUPPA, merge and multi-comparison Nextflow tasks."""
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--rscript", default=os.getenv("RSCRIPT_BIN", "Rscript"))
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    fixture = output / "fixture"
    subprocess.run([args.rscript, str(ROOT / "tests/make_qc_fixture.R"), str(fixture)], check=True)
    event = "gene;SE:chr1:1-2:3-4:+"
    (fixture / "events.ioe").write_text(
        "seqname\tgene_id\tevent_id\talternative_transcripts\ttotal_transcripts\n"
        "chr1\tgene\t{}\ttx1\ttx1,tx2\n".format(event), encoding="utf-8")
    command = ["nextflow", "-c", str(ROOT / "tests/qc_smoke.config"), "run",
               str(ROOT / "tests/splicing_smoke.nf"), "-ansi-log", "false",
               "-work-dir", str(output / "work"), "--fixture", str(fixture),
               "--pipeline_bin", str(ROOT / "bin"), "--rscript_bin", args.rscript,
               "--outdir", str(output / "results"), "-with-trace", str(output / "trace.tsv")]
    env = dict(os.environ, NXF_OFFLINE="true")
    process = subprocess.run(command, cwd=str(output), env=env, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, universal_newlines=True)
    (output / "nextflow.log").write_text(process.stdout, encoding="utf-8")
    if process.returncode:
        raise RuntimeError("Splicing smoke failed: " + str(output / "nextflow.log"))
    results = output / "results/AlternativeSplicing"
    with (results / "merged_project_events.psi").open() as f:
        row = next(csv.DictReader(f, delimiter="\t"))
    expected = {"alpha_AAAA": 1 / 3, "alpha_AAAC": 2 / 5, "alpha_AAAG": 3 / 7,
                "beta_CCCC": 4 / 9, "beta_CCCA": 5 / 11, "beta_CCCG": 6 / 13}
    assert row["ASE"] == event
    for cell, value in expected.items():
        assert abs(float(row[cell]) - value) < 1e-10, (cell, row[cell], value)
    statistics = []
    for name in ("case_vs_control", "control_vs_case"):
        directory = results / "Differential" / name
        meta = json.loads((directory / "statistics.json").read_text())
        assert meta["tested_events"] == 1
        with (directory / "all_events.tsv").open() as f:
            statistics.append(next(csv.DictReader(f, delimiter="\t")))
    assert abs(float(statistics[0]["log2FC"]) + float(statistics[1]["log2FC"])) < 1e-12
    assert abs(float(statistics[0]["p_value"]) - float(statistics[1]["p_value"])) < 1e-12
    assert (results / "split_files/input_provenance.json").exists()
    assert (results / "undefined_psi_by_sample.tsv").exists()
    (output / "summary.json").write_text(json.dumps({"status": "passed", "cells": 6,
        "comparisons": 2, "native_suppa": True, "ratios_match_raw_counts": True}, indent=2))
    print("Splicing Nextflow smoke passed: native SUPPA ratios and two comparisons", flush=True)


if __name__ == "__main__":
    main()
