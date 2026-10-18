#!/usr/bin/env python3
"""All annotated isoforms and undirected homology components, before fitting models."""
import argparse
import json
import subprocess
from collections import Counter
from pathlib import Path

from diffuprimer.acquisition import load_reference_sequences
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.splits import assign_cluster_splits


parser = argparse.ArgumentParser()
parser.add_argument("--mmseqs", default="tools/mmseqs/bin/mmseqs")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
destination = root / "data/splits/pilot.jsonl"
manifest_path = destination.with_suffix(".manifest.json")
selection_path = root / "data/processed/pilot_gene_selection.jsonl"
if destination.exists():
    manifest = json.loads(manifest_path.read_text())
    if manifest["selection_sha256"] != sha256(selection_path) or manifest["split_sha256"] != sha256(destination):
        raise ValueError("Frozen split or input selection was modified")
    print("Frozen split intact; no repartitioning performed")
    raise SystemExit(0)
genes = read_jsonl(selection_path)
acc_to_gene = {acc: gene["gene_id"] for gene in genes for acc in gene["reference_transcripts"]}
seqs = load_reference_sequences(root / "data/raw/refseq/GCF_000001635.27_GRCm39_rna.fna.gz", set(acc_to_gene))
fasta = root / "data/processed/pilot_all_isoforms.fasta"
with fasta.open("w") as handle:
    for acc, sequence in sorted(seqs.items()):
        handle.write(f">{acc}\n{sequence}\n")
hits_path = root / "data/processed/isoform_homology.tsv"
command = [str((root / args.mmseqs).resolve()), "easy-search", str(fasta), str(fasta), str(hits_path),
           str(root / "data/processed/mmseqs_tmp"), "--search-type", "3", "--min-seq-id", "0.8",
           "-c", "0.8", "--cov-mode", "0", "--max-seqs", "100000", "--threads", "8",
           "--format-output", "query,target,fident,alnlen,qcov,tcov"]
subprocess.run(command, check=True, stdout=(root / "results/mmseqs_split.log").open("w"), stderr=subprocess.STDOUT)
edges = set()
for line in hits_path.read_text().splitlines():
    q, t, identity, length, qcov, tcov = line.split("\t")
    if float(identity) >= 0.8 and float(qcov) >= 0.8 and float(tcov) >= 0.8:
        a, b = acc_to_gene[q], acc_to_gene[t]
        if a != b:
            edges.add(tuple(sorted((a, b))))
config = json.loads((root / "configs/pilot.json").read_text())
splits = assign_cluster_splits([gene["gene_id"] for gene in genes], edges, config["split_fractions"], config["seed"])
write_jsonl(destination, splits)
version = subprocess.check_output([str((root / args.mmseqs).resolve()), "version"], text=True).strip()
write_json(manifest_path, {"frozen_at": utc_now(), "selection_sha256": sha256(selection_path),
                         "all_isoforms_sha256": sha256(fasta), "homology_hits_sha256": sha256(hits_path),
                         "split_sha256": sha256(destination), "seed": config["seed"],
                         "homology_rule": "connected genes if any annotated NM_ isoform aligns at >=80% identity and >=80% query AND target coverage",
                         "mmseqs_version": version, "command": command,
                         "similarity_search_limitation": "Heuristic search; absence of a hit is not proof of unrelatedness. Additional family and domain audit required for final study.",
                         "cross_gene_edges": len(edges), "counts": dict(Counter(x["split"] for x in splits)),
                         "clusters": len({x["cluster_id"] for x in splits})})
print(manifest_path.read_text())
