#!/usr/bin/env python3
"""Independent all-strata ungapped mismatch audit, cached by oligo sequence."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

from bowtie_audit import alignment_command
from bowtie_sqlite_cache import AlignmentCache
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.specificity import classify_transcript_products


parser = argparse.ArgumentParser()
parser.add_argument("directories", nargs="+")
parser.add_argument("--threads", type=int, default=4)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
cache_dir = root / "results/bowtie_audit"
cache_dir.mkdir(parents=True, exist_ok=True)
binary = root / "tools/bowtie-1.3.1-linux-x86_64/bowtie"
db = root / "data/processed/bowtiedb"
identity = {"parser_sha256": sha256(root / "scripts/bowtie_audit.py"),
    "sqlite_cache_sha256": sha256(root / "scripts/bowtie_sqlite_cache.py"),
    "pair_classifier_sha256": sha256(root / "src/diffuprimer/specificity.py"),
    "tool_manifest_sha256": sha256(root / "tools/bowtie_manifest.json"),
    "databases": {name: sha256(db / f"{name}.manifest.json") for name in ["transcriptome", "genome"]},
    "alignment_rule": "All full-length ungapped ACGT-reference alignments with <=2 substitutions and exact original 3prime terminal 5 bases",
    "product_rule": "Inward F/R,R/F,F/F,R/R products 36..1000 bases; same-gene isoforms accepted",
    "product_evidence_format": "one_coordinate_witness_per_subject_and_primer_role"}
check = json.loads((root / "results/checks/bowtie_exhaustive/check.json").read_text())
if not check["exact_set_equality"] or check["parser_sha256"] != identity["parser_sha256"]:
    raise ValueError("Actual aligner/parser exhaustive reference check missing or stale")
pairing_check = json.loads((root / "results/checks/sqlite_specificity/check.json").read_text())
if not all(pairing_check[key] for key in ["exact_subject_role_coverage", "exact_gene_specificity_decisions", "exact_genomic_existence_decisions"]):
    raise ValueError("SQLite specificity equivalence check failed")
if pairing_check["cache_sha256"] != identity["sqlite_cache_sha256"] or pairing_check["parser_sha256"] != identity["parser_sha256"]:
    raise ValueError("SQLite specificity equivalence check is stale")
identity_path = cache_dir / "identity.json"
if identity_path.exists() and json.loads(identity_path.read_text()) != identity:
    previous = json.loads(identity_path.read_text())
    unchanged = ["pair_classifier_sha256", "tool_manifest_sha256", "databases", "alignment_rule", "product_rule"]
    if not (cache_dir / "raw_alignment_resume.json").exists() or any(previous[k] != identity[k] for k in unchanged):
        raise ValueError("Existing audit cache uses different algorithm/reference")
    if (cache_dir / "sequence_hits.jsonl").exists():
        raise ValueError("Unexpected prior materialized cache; explicit migration required")
    write_json(cache_dir / "identity_before_sqlite.json", previous)
    write_json(cache_dir / "sqlite_migration.json", {"at": utc_now(), "previous": previous,
        "current": identity, "equivalence_check_sha256": sha256(root / "results/checks/sqlite_specificity/check.json")})
write_json(identity_path, identity)
cache_path = cache_dir / "alignments.sqlite"
cache = AlignmentCache(cache_path)


def query_map(path):
    lines = path.read_text().splitlines()
    if len(lines) % 2 or any(not lines[i].startswith(">") for i in range(0, len(lines), 2)):
        raise ValueError("Expected two-line query FASTA")
    return {lines[i][1:]: lines[i+1] for i in range(0, len(lines), 2)}


# Preserve and ingest the already completed large batch before writing new queries.
resume_path = cache_dir / "raw_alignment_resume.json"
resume_done = cache_dir / "raw_alignment_resume_ingested.json"
if resume_path.exists() and not resume_done.exists():
    resume = json.loads(resume_path.read_text())
    queries = cache_dir / "new_queries.fasta"
    if sha256(queries) != resume["queries_sha256"]:
        raise ValueError("Preserved query file changed")
    qmap = query_map(queries)
    for reference, name in enumerate(["transcriptome", "genome"]):
        output = cache_dir / f"new_{name}.sam"
        evidence = resume["alignments"][name]
        if evidence["exit_code"] != 0 or sha256(output) != evidence["sam_sha256"]:
            raise ValueError("Preserved completed alignment changed")
        print(json.dumps({"ingesting_preserved_alignment": name, "at": utc_now()}), flush=True)
        cache.ingest(output, qmap, reference)
    write_json(resume_done, {"completed_at": utc_now(), "resume_sha256": sha256(resume_path)})
inputs, sequences = [], set()
for raw in args.directories:
    directory = root / raw
    source = directory / "candidates.jsonl"
    completion = json.loads((directory / "completion.json").read_text())
    if not completion["all_expected_targets_evaluated"] or sha256(source) != completion["candidates_sha256"]:
        raise ValueError("Incomplete or changed evaluation")
    candidates = read_jsonl(source)
    for row in candidates:
        if row["assessment"]["computational_feasible"]:
            sequences.update([row["forward_primer"], row["reverse_primer"]])
    inputs.append((directory, candidates))
registered = cache.register(sequences)
missing = sorted(s for s in sequences if not all(registered[s][1:]))
commands = []
if missing:
    qmap = {"Q" + hashlib.sha256(sequence.encode()).hexdigest(): sequence for sequence in missing}
    batch_dir = cache_dir / "batches" / hashlib.sha256("\n".join(missing).encode()).hexdigest()
    batch_dir.mkdir(parents=True, exist_ok=True)
    queries = batch_dir / "queries.fasta"
    queries.write_text("".join(f">{qid}\n{sequence}\n" for qid, sequence in qmap.items()))
    for reference, name in enumerate(["transcriptome", "genome"]):
        output = batch_dir / f"{name}.sam"
        command = alignment_command(binary, db / name, queries, output, args.threads)
        commands.append(command)
        print(json.dumps({"alignment": name, "new_oligos": len(missing), "at": utc_now()}), flush=True)
        marker = batch_dir / f"{name}.completed.json"
        valid = False
        if marker.exists():
            old = json.loads(marker.read_text())
            valid = old["command"] == command and old["queries_sha256"] == sha256(queries) and old["sam_sha256"] == sha256(output) and old["identity_sha256"] == sha256(identity_path)
        if not valid:
            with (batch_dir / f"{name}.log").open("w") as handle:
                subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)
            write_json(marker, {"completed_at": utc_now(), "command": command,
                "queries_sha256": sha256(queries), "sam_sha256": sha256(output), "identity_sha256": sha256(identity_path)})
        cache.ingest(output, qmap, reference)
    registered = cache.register(sequences)
    write_json(cache_dir / "last_alignment.json", {"completed_at": utc_now(), "new_oligos": len(missing),
        "queries_sha256": sha256(queries), "commands": commands, "cache_path": str(cache_path)})
transcript_map = json.loads((root / "data/processed/blastdb/transcript_gene_map.json").read_text())
for directory, candidates in inputs:
    records = []
    for row in candidates:
        if not row["assessment"]["computational_feasible"]:
            records.append({"candidate_id": row["candidate_id"], "status": "not_screened_computational_failure"})
            continue
        forward, reverse = registered[row["forward_primer"]], registered[row["reverse_primer"]]
        if not all(forward[1:]) or not all(reverse[1:]):
            raise ValueError("Pair contains an incompletely aligned oligo")
        products = cache.products(forward[0], reverse[0], 0)
        result = classify_transcript_products(products, row["gene_id"], transcript_map)
        genomic = cache.products(forward[0], reverse[0], 1)
        result.update(genomic_potential_amplicons=genomic, genomic_DNA_product_detected=bool(genomic),
                      genome_status="full_length_ungapped_mismatch_rule_assessed")
        records.append({"candidate_id": row["candidate_id"], "status": "full_length_ungapped_mismatch_rule_assessed", **result})
    output = directory / "specificity_bowtie.jsonl"
    write_jsonl(output, records)
    write_json(directory / "specificity_bowtie_status.json", {"status": "completed", "completed_at": utc_now(),
        "candidate_sha256": sha256(directory / "candidates.jsonl"), "output_sha256": sha256(output),
        "product_evidence_format": identity["product_evidence_format"],
        "identity_sha256": sha256(identity_path), "unique_oligos": len({s for row in candidates
             if row["assessment"]["computational_feasible"] for s in [row["forward_primer"], row["reverse_primer"]]})})
    print(json.dumps({"directory": str(directory), "RNA_pass_proposals": sum(r.get("gene_specificity_screen_pass", False) for r in records)}), flush=True)
