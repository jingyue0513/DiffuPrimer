#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.specificity import classify_transcript_products, parse_blast_tsv, potential_amplicons


parser = argparse.ArgumentParser()
parser.add_argument("evaluation_directory")
parser.add_argument("--genome", action="store_true")
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
directory = root / args.evaluation_directory
source = directory / "candidates.jsonl"
candidates = read_jsonl(source)
sequence_ids, sequence_map = {}, {}
for row in candidates:
    if row["assessment"]["computational_feasible"]:
        for key in ["forward_primer", "reverse_primer"]:
            sequence = row[key]
            if sequence not in sequence_ids:
                qid = f"P{len(sequence_ids):06d}"
                sequence_ids[sequence] = qid
                sequence_map[qid] = sequence
if not sequence_ids:
    write_json(directory / "specificity_status.json", {"status": "no_computationally_feasible_candidates", "candidate_sha256": sha256(source)})
    raise SystemExit(0)
queries = directory / "specificity_queries.fasta"
queries.write_text("".join(f">{qid}\n{sequence}\n" for qid, sequence in sequence_map.items()))
db_dir = root / "data/processed/blastdb"
transcript_map = json.loads((db_dir / "transcript_gene_map.json").read_text())
screens, commands = {}, []
names = ["transcriptome"] + (["genome"] if args.genome else [])
for name in names:
    manifest = db_dir / f"{name}.manifest.json"
    if not manifest.exists():
        raise ValueError(f"Specificity database not built: {name}")
    output = directory / f"blast_{name}.tsv"
    command = [str(root / "tools/ncbi-blast-2.17.0+/bin/blastn"), "-task", "blastn-short",
               "-query", str(queries), "-db", str(db_dir / name), "-out", str(output),
               "-dust", "no", "-soft_masking", "false", "-word_size", "7", "-evalue", "1000",
               "-max_target_seqs", "1000000", "-max_hsps", "1000000", "-num_threads", "8",
               "-outfmt", "6 qseqid sseqid qstart qend sstart send length mismatch gaps qlen qseq sseq bitscore evalue"]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL,
                   stderr=(directory / f"blast_{name}.stderr").open("w"))
    commands.append(command)
    screens[name] = parse_blast_tsv(output)
output = []
for row in candidates:
    if not row["assessment"]["computational_feasible"]:
        output.append({"candidate_id": row["candidate_id"], "status": "not_screened_computational_failure"})
        continue
    forward = sequence_ids[row["forward_primer"]]
    reverse = sequence_ids[row["reverse_primer"]]
    hits = screens["transcriptome"]
    products = potential_amplicons(hits.get(forward, []), hits.get(reverse, []))
    screen = classify_transcript_products(products, row["gene_id"], transcript_map)
    if "genome" in screens:
        genome = screens["genome"]
        screen["genomic_potential_amplicons"] = potential_amplicons(genome.get(forward, []), genome.get(reverse, []))
        screen["genomic_DNA_product_detected"] = bool(screen["genomic_potential_amplicons"])
        screen["genome_status"] = "limited_alignment_rule_assessed"
    else:
        screen["genome_status"] = "not_assessed"
    output.append({"candidate_id": row["candidate_id"], "status": "limited_alignment_rule_assessed", **screen})
write_jsonl(directory / "specificity.jsonl", output)
write_json(directory / "specificity_status.json", {"status": "limited_alignment_rule_assessed",
    "completed_at": utc_now(), "candidate_sha256": sha256(source), "queries_sha256": sha256(queries),
    "commands": commands, "databases": {name: sha256(db_dir / f"{name}.manifest.json") for name in names},
    "limitations": ["Heuristic local BLAST alignments; not exhaustive mismatch enumeration", "Full-length ungapped hit rule excludes partial and gapped priming", "Genome search, when requested, reports possible genomic-DNA amplification separately", "In-silico screening is not experimental specificity validation"]})
print(json.dumps({"screened_candidates": sum(x["status"] == "limited_alignment_rule_assessed" for x in output),
                  "gene_screen_pass": sum(x.get("gene_specificity_screen_pass", False) for x in output)}), flush=True)
