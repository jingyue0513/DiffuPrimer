#!/usr/bin/env python3
"""Apply the frozen specificity rule, reusing identical-sequence BLAST searches."""
import argparse
import fcntl
import json
import subprocess
from pathlib import Path

from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.specificity import classify_transcript_products, parse_blast_tsv
from amplicon_pairing_fast import potential_amplicons
from blast_query_batches import run_query_batches


def fasta_sequences(path):
    records, name, parts = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directories", nargs="+")
    parser.add_argument("--reuse-from", nargs="*", default=[])
    parser.add_argument("--batch-size", type=int, default=1000)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    pairing_check = json.loads((root / "results/checks/fast_pairing/check.json").read_text())
    if not pairing_check["exact_ordered_output_equality"] or pairing_check["indexed_sha256"] != sha256(root / "scripts/amplicon_pairing_fast.py"):
        raise ValueError("Indexed pairing equivalence check missing or stale")
    db = root / "data/processed/blastdb"
    cache_dir = root / "results/specificity_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    identity = {"specificity_code_sha256": sha256(root / "src/diffuprimer/specificity.py"),
        "blastn_sha256": sha256(root / "tools/ncbi-blast-2.17.0+/bin/blastn"),
        "databases": {name: sha256(db / (name + ".manifest.json")) for name in ["transcriptome", "genome"]},
        "search_parameters": ["blastn-short", "dust=no", "soft_masking=false", "word_size=7", "evalue=1000", "max_target_seqs=1000000", "max_hsps=1000000"]}
    identity_path = cache_dir / "identity.json"
    if identity_path.exists() and json.loads(identity_path.read_text()) != identity:
        raise ValueError("Specificity rule or reference changed; cached hits cannot be reused")
    write_json(identity_path, identity)
    cache_path = cache_dir / "sequence_hits.jsonl"
    cache = {row["sequence"]: row for row in read_jsonl(cache_path)} if cache_path.exists() else {}
    sources = []
    for raw in args.reuse_from:
        directory = root / raw
        status = json.loads((directory / "specificity_status.json").read_text())
        if status["databases"] != identity["databases"]:
            raise ValueError("Source screen uses a different reference")
        if status["candidate_sha256"] != sha256(directory / "candidates.jsonl"):
            raise ValueError("Source screen candidate file changed")
        queries = directory / "specificity_queries.fasta"
        if status["queries_sha256"] != sha256(queries):
            raise ValueError("Source query file changed")
        mapping = fasta_sequences(queries)
        hits = {name: parse_blast_tsv(directory / ("blast_" + name + ".tsv")) for name in identity["databases"]}
        for qid, sequence in mapping.items():
            cache.setdefault(sequence, {"sequence": sequence, **{name: hits[name].get(qid, []) for name in hits}})
        sources.append({"directory": str(directory), "status_sha256": sha256(directory / "specificity_status.json"),
            "query_sha256": sha256(queries), "blast_sha256": {name: sha256(directory / ("blast_" + name + ".tsv")) for name in hits}})
    directories = [root / raw for raw in args.directories]
    candidates = {}
    sequences = set()
    for directory in directories:
        completion = json.loads((directory / "completion.json").read_text())
        if not completion["all_expected_targets_evaluated"] or completion["candidates_sha256"] != sha256(directory / "candidates.jsonl"):
            raise ValueError("Evaluation incomplete or changed")
        own = read_jsonl(directory / "candidates.jsonl")
        candidates[directory] = own
        for row in own:
            if row["assessment"]["computational_feasible"]:
                sequences.update((row["forward_primer"], row["reverse_primer"]))
    missing = sorted(sequences - set(cache))
    commands = []
    for offset in range(0, len(missing), args.batch_size):
        batch = missing[offset:offset + args.batch_size]
        mapping = {f"Q{index:07d}": sequence for index, sequence in enumerate(batch)}
        query = cache_dir / "batch_query.fasta"
        query.write_text("".join(f">{qid}\n{sequence}\n" for qid, sequence in mapping.items()))
        batch_hash = sha256(query)[:16]
        query_record = cache_dir / (batch_hash + ".fasta")
        query_record.write_bytes(query.read_bytes())
        hits = {}
        for name in identity["databases"]:
            output = cache_dir / (batch_hash + "_" + name + ".tsv")
            command = [str(root / "tools/ncbi-blast-2.17.0+/bin/blastn"), "-task", "blastn-short",
                "-query", str(query_record), "-db", str(db / name), "-out", str(output), "-dust", "no",
                "-soft_masking", "false", "-word_size", "7", "-evalue", "1000", "-max_target_seqs", "1000000",
                "-max_hsps", "1000000", "-num_threads", "8", "-outfmt",
                "6 qseqid sseqid qstart qend sstart send length mismatch gaps qlen qseq sseq bitscore evalue"]
            with output.with_suffix(".stderr").open("w") as handle:
                if name == "genome":
                    run_query_batches(command)
                else:
                    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=handle)
            hits[name] = parse_blast_tsv(output)
            commands.append({"command": command, "output_sha256": sha256(output), "queries_sha256": sha256(query_record)})
        for qid, sequence in mapping.items():
            cache[sequence] = {"sequence": sequence, **{name: hits[name].get(qid, []) for name in hits}}
        write_jsonl(cache_path, [cache[sequence] for sequence in sorted(cache)])
        print(json.dumps({"new_sequences_completed": min(offset + len(batch), len(missing)), "new_sequences_total": len(missing)}), flush=True)
    write_jsonl(cache_path, [cache[sequence] for sequence in sorted(cache)])
    map_by_transcript = json.loads((db / "transcript_gene_map.json").read_text())
    for directory in directories:
        output = []
        for row in candidates[directory]:
            if not row["assessment"]["computational_feasible"]:
                output.append({"candidate_id": row["candidate_id"], "status": "not_screened_computational_failure"})
                continue
            forward, reverse = cache[row["forward_primer"]], cache[row["reverse_primer"]]
            products = potential_amplicons(forward["transcriptome"], reverse["transcriptome"])
            screen = classify_transcript_products(products, row["gene_id"], map_by_transcript)
            genome_products = potential_amplicons(forward["genome"], reverse["genome"])
            screen.update(genomic_potential_amplicons=genome_products, genomic_DNA_product_detected=bool(genome_products),
                genome_status="limited_alignment_rule_assessed")
            output.append({"candidate_id": row["candidate_id"], "status": "limited_alignment_rule_assessed", **screen})
        write_jsonl(directory / "specificity.jsonl", output)
        write_json(directory / "specificity_status.json", {"status": "limited_alignment_rule_assessed", "completed_at": utc_now(),
            "candidate_sha256": sha256(directory / "candidates.jsonl"), "databases": identity["databases"],
            "cache_sha256": sha256(cache_path), "cache_identity_sha256": sha256(identity_path),
            "indexed_pairing_sha256": sha256(root / "scripts/amplicon_pairing_fast.py"),
            "pairing_equivalence_check_sha256": sha256(root / "results/checks/fast_pairing/check.json"),
            "reused_source_screens": sources, "new_search_batches": commands,
            "interpretation": "Identical synthesis sequences reuse qualifying BLAST hits under unchanged reference and screening rules"})
        print(json.dumps({"directory": directory.name, "gene_screen_pass": sum(row.get("gene_specificity_screen_pass", False) for row in output)}), flush=True)
    write_json(cache_dir / "last_run.json", {"completed_at": utc_now(), "reused_sources": sources, "commands": commands,
        "new_sequences": len(missing), "total_cached_sequences": len(cache), "cache_sha256": sha256(cache_path)})


if __name__ == "__main__":
    lock_dir = Path(__file__).resolve().parents[1] / "results/specificity_cache"
    lock_dir.mkdir(parents=True, exist_ok=True)
    with (lock_dir / "writer.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        main()
