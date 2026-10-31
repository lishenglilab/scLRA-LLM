#!/usr/bin/env python3
"""Freeze deterministic release archives and update their public checksum manifest."""
import argparse
import gzip
import hashlib
import json
import tarfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", required=True, choices=["matrix", "fastq"])
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    version = "v1.0.0"
    filename = f"scLRA-demo-{args.kind}-{version}.tar.gz"
    args.output.mkdir(parents=True, exist_ok=True)
    archive = args.output / filename
    with archive.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode="w|") as tar:
            for path in sorted(args.source.rglob("*")):
                if not path.is_file() or path.name == "metadata.tsv":
                    continue
                name = str(Path(args.kind) / path.relative_to(args.source)).replace("\\", "/")
                info = tar.gettarinfo(str(path), arcname=name)
                info.uid = info.gid = info.mtime = 0
                info.uname = info.gname = ""
                info.mode = 0o644
                with path.open("rb") as source:
                    tar.addfile(info, source)
    digest = hashlib.sha256()
    with archive.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    manifest = json.loads(args.manifest.read_text()) if args.manifest.exists() else {
        "schema_version": 1, "dataset_version": version, "source_run": "SRR14808743",
        "source_series": "GSE178175", "assets": {}}
    manifest["assets"][args.kind] = {"filename": filename, "bytes": archive.stat().st_size,
        "sha256": digest.hexdigest(),
        "url": f"https://github.com/lishenglilab/scLRA-LLM/releases/download/demo-{version}/{filename}"}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["assets"][args.kind], indent=2))


if __name__ == "__main__":
    main()
