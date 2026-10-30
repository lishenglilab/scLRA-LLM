import argparse
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / "deepseek" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class PipelineUpdateTest(unittest.TestCase):
    def test_annotation_requests_are_independent_and_low_confidence_is_explicit(self):
        module = load("cell_ann_deekseek")
        client = MagicMock()
        client.chat.completions.create.side_effect = [
            SimpleNamespace(model="mock", choices=[SimpleNamespace(message=SimpleNamespace(content=text))])
            for text in ("Cluster0:T_cells[3]", "Cluster1:Keratinocytes[0]")]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            markers = tmp / "markers.csv"
            markers.write_text("0,CD3D\n0,CD3E\n1,KRT14\n", encoding="utf-8")
            args = argparse.Namespace(csv_file_path=str(markers), output_csv=str(tmp / "annotations.csv"),
                                      max_genes_to_display=30, species_type="human", tissue_type="skin", model="mock")
            with patch.object(module, "client_from_args", return_value=client):
                module.cell_ann_llm(args)
            calls = client.chat.completions.create.call_args_list
            self.assertEqual(len(calls), 2)
            self.assertEqual(len(calls[1].kwargs["messages"]), 2)
            self.assertNotIn("CD3D", calls[1].kwargs["messages"][1]["content"])
            table = module.pd.read_csv(args.output_csv)
            self.assertEqual(table["Confidence"].tolist(), [3, 0])
            self.assertEqual(table["Forced estimate"].tolist(), [False, True])
            audit = json.loads((tmp / "annotations.requests.json").read_text())
            self.assertEqual(audit["status"], "completed")
            self.assertEqual(audit["runs"][0]["response"], "Cluster0:T_cells[3]")

    def test_annotation_rejects_wrong_cluster_multiple_lines_and_invalid_labels(self):
        module = load("cell_ann_deekseek")
        for raw in ("Cluster1:T_cells[3]", "Cluster0:T_cells[3]\nCluster1:B_cells[2]",
                    "Cluster0:T_cells/B_cells[1]", "Cluster0:Unknown[0]", "Cluster0:T_cells[9]"):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                module.parse_annotation(raw, "0")

    def test_functional_prompt_keeps_all_genes_and_measurement_limits(self):
        module = load("differential_interpretation")
        genes = ["GENE{}".format(i) for i in range(35)]
        prompt = module.generate_prompt(genes, {}, "study", "case_vs_control")
        for gene in genes:
            self.assertIn(gene, prompt)
        self.assertIn("not P values", prompt)
        self.assertIn("abstain", prompt)

    def test_underscored_sample_id_is_rejected_at_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            fastq = tmp / "read.fastq"
            fastq.write_text("@read\nA\n+\nI\n")
            metadata = tmp / "metadata.tsv"
            metadata.write_text("sample_id\tfastq\ninvalid_sample\t{}\n".format(fastq))
            result = subprocess.run([sys.executable, str(ROOT / "bin/parse_metadata.py"),
                                     "--metadata", str(metadata), "--out", str(tmp / "parsed.tsv")],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Invalid sample_id", result.stderr)


if __name__ == "__main__":
    unittest.main()
