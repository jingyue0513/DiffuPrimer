#!/usr/bin/env python3
"""Run frozen BLAST screening with equivalent indexed pairing and raw-file reuse."""
import json
import runpy
import subprocess
import sys
from pathlib import Path

from amplicon_pairing_fast import potential_amplicons
from blast_query_batches import fasta, run_query_batches
from diffuprimer import specificity
from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
use_witnesses = "--product-witnesses" in sys.argv
if use_witnesses:
    sys.argv.remove("--product-witnesses")
    from cached_amplicon_witnesses import potential_amplicons
    witness_check = json.loads((root / "results/checks/cached_witnesses/check.json").read_text())
    if not all(witness_check[key] for key in ["exact_subject_role_coverage", "witness_coordinates_belong_to_full_catalog", "exact_gene_and_genome_existence_decisions"]):
        raise ValueError("Product-existence equivalence check failed")
    if witness_check["source_sha256"] != sha256(root / "scripts/cached_amplicon_witnesses.py") or witness_check["frozen_classifier_sha256"] != sha256(root / "src/diffuprimer/specificity.py"):
        raise ValueError("Product-existence equivalence check is stale")
check = json.loads((root / "results/checks/fast_pairing/check.json").read_text())
if not check["exact_ordered_output_equality"] or check["original_sha256"] != sha256(root / "src/diffuprimer/specificity.py") or check["indexed_sha256"] != sha256(root / "scripts/amplicon_pairing_fast.py"):
    raise ValueError("Pairing equivalence check missing or stale")
original_run = subprocess.run
commands_used = []


def verify_alignment_tail(path):
    if not path.stat().st_size:
        return
    with path.open("rb") as handle:
        handle.seek(max(0,path.stat().st_size-4096))
        tail=handle.read()
    if not tail.endswith(b"\n") or len(tail.splitlines()[-1].split(b"\t")) != 14:
        raise ValueError("BLAST output ends in an incomplete alignment record")


def run(command, *args, **kwargs):
    if not isinstance(command, list) or "-task" not in command or "blastn-short" not in command:
        return original_run(command, *args, **kwargs)
    output = Path(command[command.index("-out") + 1])
    query = Path(command[command.index("-query") + 1])
    database = Path(command[command.index("-db") + 1])
    marker = output.with_suffix(".completed.json")
    identity = {"command": command, "query_sha256": sha256(query),
                "database_manifest_sha256": sha256(database.with_suffix(".manifest.json"))}
    if marker.exists():
        verify_alignment_tail(output)
        recorded = json.loads(marker.read_text())
        if any(recorded[key] != value for key, value in identity.items()) or recorded["output_sha256"] != sha256(output):
            raise ValueError("Completed raw alignment identity differs")
        commands_used.append({"command": command, "action": "reuse_verified_complete_alignment"})
        return subprocess.CompletedProcess(command, 0)
    batched = database.name == "genome" and len(fasta(query)) > 1
    result = run_query_batches(command, run=original_run) if batched else original_run(command, *args, **kwargs)
    if result.returncode == 0:
        verify_alignment_tail(output)
        write_json(marker, {**identity, "output_sha256": sha256(output), "completed_at": utc_now(),
            "completion_evidence": "Successful subprocess exit checked before marker creation"})
    commands_used.append({"command": command, "action": "execute_equivalent_independent_queries" if batched else "execute_frozen_alignment_command"})
    return result


specificity.potential_amplicons = potential_amplicons
subprocess.run = run
directory = root / sys.argv[1]
sys.argv[0] = str(root / "scripts/screen_specificity.py")
runpy.run_path(sys.argv[0], run_name="__main__")
if use_witnesses:
    status_path = directory / "specificity_status.json"
    status = json.loads(status_path.read_text())
    status.update(product_evidence_format="one_coordinate_witness_per_subject_and_primer_role",
                  product_existence_equivalence_check_sha256=sha256(root / "results/checks/cached_witnesses/check.json"),
                  output_sha256=sha256(directory / "specificity.jsonl"))
    write_json(status_path, status)
write_json(directory / "specificity_execution.json", {"completed_at": utc_now(),
    "frozen_screen_script_sha256": sha256(root / "scripts/screen_specificity.py"),
    "frozen_classifier_sha256": sha256(root / "src/diffuprimer/specificity.py"),
    "indexed_pairing_sha256": sha256(root / "scripts/amplicon_pairing_fast.py"),
    "equivalence_check_sha256": sha256(root / "results/checks/fast_pairing/check.json"),
    "product_evidence_format": "one_coordinate_witness_per_subject_and_primer_role" if use_witnesses else "full_ordered_product_catalog",
    "product_witnesses_sha256": sha256(root / "scripts/cached_amplicon_witnesses.py") if use_witnesses else None,
    "commands": commands_used, "semantics": "Identical qualifying hits and every subject/role existence decision" if use_witnesses else "Identical qualifying hits and exact ordered products; distance indexing only"})
