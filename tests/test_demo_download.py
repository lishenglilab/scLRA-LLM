import hashlib
import io
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class DemoDownloadTest(unittest.TestCase):
    def exercise(self, member, digest_ok=True):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            shutil.copyfile(ROOT / "demo/download.py", tmp / "download.py")
            archive = tmp / "input.tar.gz"
            with tarfile.open(archive, "w:gz") as handle:
                info = tarfile.TarInfo(member)
                info.size = 4
                handle.addfile(info, io.BytesIO(b"data"))
            digest = hashlib.sha256(archive.read_bytes()).hexdigest() if digest_ok else "0" * 64
            (tmp / "manifest.json").write_text(json.dumps({"assets": {"fastq": {"sha256": digest}}}))
            result = subprocess.run([sys.executable, str(tmp / "download.py"), "fastq", "--archive", str(archive),
                                     "--dest", str(tmp / "output")], capture_output=True, text=True)
            if result.returncode == 0:
                self.assertIn(str((tmp / "output/fastq/DemoA.fastq.gz").resolve()),
                              (tmp / "output/fastq/metadata.tsv").read_text())
            else:
                self.assertFalse((tmp / "output/fastq").exists())
            return result

    def test_verified_archive_and_absolute_metadata(self):
        self.assertEqual(self.exercise("fastq/DemoA.fastq.gz").returncode, 0)

    def test_hash_mismatch_is_rejected(self):
        result = self.exercise("fastq/DemoA.fastq.gz", digest_ok=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("SHA256 mismatch", result.stderr)

    def test_archive_traversal_is_rejected(self):
        result = self.exercise("fastq/../escape.txt")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unexpected path", result.stderr)


if __name__ == "__main__":
    unittest.main()
