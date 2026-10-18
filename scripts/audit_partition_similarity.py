#!/usr/bin/env python3
"""Exploratory similarity sensitivity audit; never alter the frozen partition."""
import argparse
import json
import subprocess
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    genes_path = root / "data/processed/pilot_gene_selection.jsonl"
    split_path = root / "data/splits/pilot.jsonl"
    fasta = root / "data/processed/pilot_all_isoforms.fasta"
    output_dir = root / "results/partition_audit"
    output_dir.mkdir(parents=True, exist_ok=True)
    split_hash = sha256(split_path)
    genes = read_jsonl(genes_path)
    accession_gene = {acc: row["gene_id"] for row in genes for acc in row["reference_transcripts"]}
    assignments = {row["gene_id"]: row["split"] for row in read_jsonl(split_path)}
    hits = output_dir / "nucleotide_70.tsv"
    command = [str(root / "tools/mmseqs/bin/mmseqs"), "easy-search", str(fasta), str(fasta),
        str(hits), str(output_dir / "tmp70"), "--search-type", "3", "--min-seq-id", "0.7",
        "-c", "0.7", "--cov-mode", "0", "--max-seqs", "100000", "--threads", str(args.threads),
        "--format-output", "query,target,fident,alnlen,qcov,tcov"]
    with (output_dir / "search.log").open("w") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)
    rules = [(0.7, 0.7), (0.8, 0.8), (0.9, 0.9)]
    edges = {rule: {} for rule in rules}
    for line in hits.read_text().splitlines():
        q, t, identity, length, qcov, tcov = line.split("\t")
        left, right = accession_gene[q], accession_gene[t]
        if left == right:
            continue
        pair = tuple(sorted((left, right), key=int))
        for rule in rules:
            if float(identity) >= rule[0] and min(float(qcov), float(tcov)) >= rule[1]:
                if pair not in edges[rule] or float(identity) > edges[rule][pair]["identity"]:
                    edges[rule][pair] = {"genes": list(pair), "query": q, "target": t,
                        "identity": float(identity), "qcov": float(qcov), "tcov": float(tcov),
                        "splits": [assignments[gene] for gene in pair]}
    reports = []
    for rule in rules:
        all_edges = list(edges[rule].values())
        cross = [row for row in all_edges if len(set(row["splits"])) > 1]
        train_test = [row for row in cross if set(row["splits"]) == {"train", "test"}]
        reports.append({"identity": rule[0], "bidirectional_coverage": rule[1],
            "cross_gene_edges": len(all_edges), "cross_partition_edges": len(cross),
            "train_test_edges": len(train_test), "cross_partition_details": cross,
            "test_genes_with_train_edge": sorted({gene for row in train_test for gene in row["genes"]
                if assignments[gene] == "test"}, key=int)})
    if sha256(split_path) != split_hash:
        raise ValueError("Frozen partition changed during read-only audit")
    report = {"generated_at": utc_now(), "classification": "exploratory nucleotide-threshold sensitivity audit",
        "partition_modified": False, "split_sha256": split_hash, "selection_sha256": sha256(genes_path),
        "fasta_sha256": sha256(fasta), "hits_sha256": sha256(hits), "command": command,
        "rules": reports, "limitations": ["Heuristic alignment search is not exhaustive",
            "RNA similarity thresholds do not establish absence of shared gene families or domains",
            "Primary evaluation keeps the originally frozen partition and denominator"]}
    write_json(output_dir / "summary.json", report)
    text = "An exploratory nucleotide-similarity sensitivity audit retained the frozen partition. "
    for row in reports:
        text += (f"At {100*row['identity']:.0f}\\% identity and {100*row['bidirectional_coverage']:.0f}\\% bidirectional coverage, "
            f"the search found {row['cross_gene_edges']} cross-gene edges, including {row['cross_partition_edges']} across partitions "
            f"and {row['train_test_edges']} between training and test genes. ")
    text += "These heuristic nucleotide checks do not replace a gene-family or domain audit.\n"
    (root / "paper/tables/partition_audit.tex").write_text(text)
    print(json.dumps(reports, indent=2))


if __name__ == "__main__":
    main()
