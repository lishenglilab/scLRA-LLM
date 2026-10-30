import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from llm.client import add_llm_arguments, client_from_args, response_text


def normalize_cluster_id(value):
    return re.sub(r"^cluster\s*", "", str(value).strip(), flags=re.IGNORECASE).strip()


def parse_annotation(raw, cluster):
    match = re.fullmatch(r"\s*Cluster\s*([^:]+?)\s*:\s*([^\[\]\n]+)\[([0-3])\]\s*", raw)
    if not match or normalize_cluster_id(match.group(1)) != cluster:
        raise ValueError("Expected exactly Cluster<ID>:CellType[0-3] for the requested cluster")
    label = match.group(2).strip()
    if not label or any(x in label for x in ("/", "*", "<", ">")) or re.search(
        r"unknown|unassigned|insufficient|unresolved", label, re.IGNORECASE
    ):
        raise ValueError("Annotation must contain one concrete major cell type")
    return label, int(match.group(3))


def cell_ann_llm(args):
    if args.max_genes_to_display < 1:
        raise ValueError("--max_genes_to_display must be positive")
    df = pd.read_csv(args.csv_file_path, header=None, dtype=str, keep_default_na=False)
    if df.shape[1] < 2:
        raise ValueError("Marker CSV must contain cluster and gene columns")
    clusters = {}
    for raw_cluster, raw_genes in df.iloc[:, :2].itertuples(index=False, name=None):
        cluster = normalize_cluster_id(raw_cluster)
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", cluster):
            raise ValueError("Invalid cluster identifier")
        genes = clusters.setdefault(cluster, [])
        for gene in str(raw_genes).split(","):
            gene = gene.strip()
            if gene and gene not in genes:
                genes.append(gene)
    if not clusters or any(not genes for genes in clusters.values()):
        raise ValueError("Each cluster must have nonempty markers")
    client = client_from_args(args)
    prompt = (Path(__file__).resolve().parents[1] / "llm" / "prompts" / "cell_annotation.txt").read_text(encoding="utf-8").format(
        species_type=args.species_type, tissue_type=args.tissue_type)
    audit_path = Path(args.output_csv).with_suffix(".requests.json")
    history, annotations = [], []
    audit = {"protocol": "independent-cluster-v1", "model": args.model,
             "confidence": "qualitative evidence grade; 0 is a forced estimate, not a resolved label",
             "status": "running", "runs": history}
    def save_audit():
        audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    try:
        for cluster, genes in clusters.items():
            messages = [{"role": "system", "content": prompt},
                        {"role": "user", "content": "Cluster{}: {}".format(cluster, ", ".join(genes[:args.max_genes_to_display]))}]
            record = {"cluster": cluster, "messages": messages, "status": "requested"}
            history.append(record)
            save_audit()
            response = client.chat.completions.create(model=args.model, messages=messages)
            raw = response_text(response)
            record.update(response=raw, returned_model=getattr(response, "model", None))
            label, confidence = parse_annotation(raw, cluster)
            record["status"] = "completed"
            annotations.append({"Cluster": "Cluster" + cluster,
                                "Cell Type": "{}[{}]".format(label, confidence),
                                "Confidence": confidence, "Forced estimate": confidence == 0})
            save_audit()
        pd.DataFrame(annotations).to_csv(args.output_csv, index=False)
        audit["status"] = "completed"
    except Exception as exc:
        audit["status"] = "failed"
        audit["error_type"] = type(exc).__name__
        raise
    finally:
        save_audit()
    return "\n".join(record["response"] for record in history), history


def main():
    parser = argparse.ArgumentParser(description="Annotate each cluster in an independent LLM request")
    parser.add_argument("--csv_file_path", required=True)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--tissue_type", default="Skin")
    parser.add_argument("--species_type", default="human")
    parser.add_argument("--max_genes_to_display", type=int, default=30)
    add_llm_arguments(parser)
    args = parser.parse_args()
    cell_ann_llm(args)


if __name__ == "__main__":
    main()
