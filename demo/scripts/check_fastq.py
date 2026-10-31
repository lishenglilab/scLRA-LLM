#!/usr/bin/env python3
"""Verify read IDs, sequences and quality strings, ignoring FASTQ header comments."""
import argparse
import gzip
import hashlib
import json


def canonical_digest(path):
    digest = hashlib.sha256()
    count = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            separator = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            if not header.startswith("@") or not separator.startswith("+") or len(sequence) != len(quality):
                raise ValueError("Invalid/truncated FASTQ record")
            record = "@" + header[1:].split()[0] + "\n" + sequence + "\n+\n" + quality + "\n"
            digest.update(record.encode())
            count += 1
    return count, digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq")
    parser.add_argument("--selection", required=True)
    parser.add_argument("--sample", required=True, choices=["DemoA", "DemoB"])
    args = parser.parse_args()
    count, digest = canonical_digest(args.fastq)
    expected = json.load(open(args.selection))["subsets"][args.sample]
    if count != expected["records"] or digest != expected["canonical_record_sha256"]:
        raise SystemExit("Read content/order mismatch")
    print(args.sample, count, "records: identifiers, sequences, qualities and ordering verified")


if __name__ == "__main__":
    main()
