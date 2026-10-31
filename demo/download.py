#!/usr/bin/env python3
"""Download and verify a versioned demo asset using Python 3's standard library."""
import argparse
import hashlib
import json
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=["matrix", "fastq"])
    parser.add_argument("--archive", type=Path, help="Use an already downloaded archive")
    parser.add_argument("--dest", type=Path, default=ROOT / "data")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "manifest.json").read_text())
    asset = manifest["assets"][args.kind]
    destination = args.dest.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    target = destination / args.kind
    if target.exists():
        raise SystemExit(f"Destination already exists: {target}. Use a fresh --dest directory.")
    with tempfile.TemporaryDirectory(prefix=".download-", dir=destination) as tmp:
        tmp = Path(tmp)
        archive = tmp / "asset.tar.gz"
        if args.archive:
            shutil.copyfile(args.archive, archive)
        else:
            print("Downloading", asset["url"], flush=True)
            request = urllib.request.Request(asset["url"], headers={"User-Agent": "scLRA-LLM-demo"})
            with urllib.request.urlopen(request, timeout=120) as response, archive.open("wb") as output:
                shutil.copyfileobj(response, output)
        hasher = hashlib.sha256()
        with archive.open("rb") as downloaded:
            for block in iter(lambda: downloaded.read(1024 * 1024), b""):
                hasher.update(block)
        digest = hasher.hexdigest()
        if digest != asset["sha256"]:
            raise SystemExit("Archive SHA256 mismatch; download discarded")
        unpacked = tmp / "unpacked"
        unpacked.mkdir()
        with tarfile.open(archive) as handle:
            for member in handle.getmembers():
                path = Path(member.name)
                if path.is_absolute() or ".." in path.parts or path.parts[0] != args.kind:
                    raise SystemExit("Archive contains an unexpected path")
                if not (member.isfile() or member.isdir()):
                    raise SystemExit("Archive links and special files are not permitted")
                output = unpacked / path
                if member.isdir():
                    output.mkdir(parents=True, exist_ok=True)
                else:
                    output.parent.mkdir(parents=True, exist_ok=True)
                    with handle.extractfile(member) as source, output.open("wb") as dest:
                        shutil.copyfileobj(source, dest)
        shutil.move(str(unpacked / args.kind), str(target))
    if args.kind == "fastq":
        (target / "metadata.tsv").write_text("sample_id\tfastq\n" + "".join(
            f"{sample}\t{target / (sample + '.fastq.gz')}\n" for sample in ("DemoA", "DemoB")))
    print("Verified and prepared:", target)
    flag = "--matrix_input_dir" if args.kind == "matrix" else "--metadata"
    input_path = target if args.kind == "matrix" else target / "metadata.tsv"
    print(f'Run: nextflow run main.nf -profile demo_{args.kind} {flag} "{input_path}" --llm_enabled false')


if __name__ == "__main__":
    main()
