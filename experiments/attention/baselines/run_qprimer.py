#!/usr/bin/env python3
"""Run the upstream qPrimer 1.0.6 design module with shared qPCR bounds.

Retains qPrimer's additional design rules. Candidate evaluation and specificity
assessment remain shared with other methods, so the tool's own screening is not
substituted for the common benchmark assessment.
"""
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import summarize_candidates
from diffuprimer.sequence import PrimerPair

UPSTREAM_COMMIT = "34f047e5d77ec0353846a9a18b0b94f74013db3b"
UPSTREAM_DESIGN_SHA256 = "43fdb0329db15d1188c266cb5291614acc8fb1cbe309aaed9b972b89ec6eb1e5"


def load_upstream():
    vendor = Path(__file__).parent / "vendor/qPrimer"
    if (vendor / ".git").exists():
        commit = subprocess.check_output(["git", "-C", str(vendor), "rev-parse", "HEAD"], text=True).strip()
    else:
        commit = json.loads((vendor / "source_snapshot.json").read_text())["commit"]
    if commit != UPSTREAM_COMMIT:
        raise ValueError(f"Expected qPrimer {UPSTREAM_COMMIT}, got {commit}")
    if sha256(vendor / "qPrimer/design.py") != UPSTREAM_DESIGN_SHA256:
        raise ValueError("Pinned upstream qPrimer design source changed")
    spec = importlib.util.spec_from_file_location("qprimer_upstream_design", vendor / "qPrimer/design.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure(module, conditions, budget):
    # All methods share biological bounds and thermodynamic concentrations.
    # Remaining upstream defaults (including end stability and alignment scores)
    # retain qPrimer's distinct native design profile.
    c = conditions
    module.GLOBAL_ARGS.update({
        "PRIMER_TASK": "generic", "PRIMER_NUM_RETURN": budget,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_MIN_SIZE": c.primer_min_length, "PRIMER_MAX_SIZE": c.primer_max_length,
        "PRIMER_MIN_TM": c.tm_min, "PRIMER_OPT_TM": c.tm_opt, "PRIMER_MAX_TM": c.tm_max,
        "PRIMER_PAIR_MAX_DIFF_TM": c.max_tm_difference,
        "PRIMER_MIN_GC": c.gc_min, "PRIMER_MAX_GC": c.gc_max,
        "PRIMER_PRODUCT_SIZE_RANGE": [[c.amplicon_min, c.amplicon_max]],
        "PRIMER_MAX_POLY_X": c.max_homopolymer,
        "PRIMER_SALT_MONOVALENT": c.monovalent_mM,
        "PRIMER_SALT_DIVALENT": c.divalent_mM,
        "PRIMER_DNTP_CONC": c.dntp_mM, "PRIMER_DNA_CONC": c.primer_nM,
        "PRIMER_MAX_HAIRPIN_TH": c.max_structure_tm,
        "PRIMER_MAX_SELF_ANY_TH": c.max_structure_tm,
        "PRIMER_MAX_SELF_END_TH": c.max_structure_tm,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": c.max_structure_tm,
        "PRIMER_PAIR_MAX_COMPL_END_TH": c.max_structure_tm,
    })
    module.SEQ_ARGS.clear()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--task", choices=["design", "completion", "redesign"], default="design")
    p.add_argument("--fixed-side", choices=["forward", "reverse"])
    p.add_argument("--split", choices=["validation", "test", "external_test"], default="test")
    p.add_argument("--targets-path",default="data/processed/targets_pilot.jsonl")
    p.add_argument("--output-root",default="experiments/attention/results/evaluations")
    p.add_argument("--amplicon-min", type=int, default=80)
    p.add_argument("--amplicon-max", type=int, default=200)
    p.add_argument("--budget", type=int, default=20)
    p.add_argument("--seed", type=int, default=20261013)
    a = p.parse_args()
    if a.task != "design" and a.fixed_side is None:
        p.error("Completion and redesign require --fixed-side")
    if a.task != "design" and a.split != "test":
        p.error("Frozen anchor tasks are defined on the test partition")
    module = load_upstream()
    config = json.loads((ROOT / "configs/pilot.json").read_text())
    target_path=ROOT / a.targets_path
    targets = [t for t in read_jsonl(target_path) if t["split"] == a.split]
    targets.sort(key=lambda t: int(hashlib.sha256(f"{t['gene_id']}:{a.seed}".encode()).hexdigest()[:8], 16))
    anchors = {t["gene_id"]: t for t in read_jsonl(ROOT / "experiments/final_followup/redesign_tasks.jsonl")}
    name = f"qprimer_{a.task}_{a.fixed_side or str(a.amplicon_min)+'-'+str(a.amplicon_max)}_{a.split}_seed{a.seed}"
    out = ROOT / a.output_root / name
    out.mkdir(parents=True, exist_ok=True)
    base = dict(config["conditions"], amplicon_min=a.amplicon_min, amplicon_max=a.amplicon_max)
    configure(module, Conditions(**base), a.budget)
    manifest = {"method": "qprimer", "stage": "attention_external_comparison" if a.split=="external_test" else "attention_comparison", "arguments": vars(a),
        "task": a.task, "split":a.split,"conditions": base, "upstream_repository": "https://github.com/swu1019lab/qPrimer",
        "upstream_commit": UPSTREAM_COMMIT, "upstream_design_sha256": sha256(Path(module.__file__)),
        "adapter_sha256": sha256(Path(__file__)), "upstream_global_arguments_after_shared_overrides": module.GLOBAL_ARGS,
        "target_sha256": sha256(target_path),
        "expected_targets": [t["gene_id"] for t in targets],
        "candidate_accounting": "At most K returned candidates; shortages and duplicates count against K",
        "screening": "Upstream design module; shared thermodynamic assessment and full-reference Bowtie downstream",
        "pair_only_adjustment": "Internal probe generation disabled for the common F/R pair task"}
    write_json(out / "manifest.json", manifest)
    rows, all_candidates = [], []
    for target in targets:
        started = time.monotonic()
        gene, template = target["gene_id"], target["template_sequence"]
        data = dict(base)
        task = anchors.get(gene, {})
        if a.task == "redesign":
            data.update(amplicon_min=task["new_length_band"][0], amplicon_max=task["new_length_band"][1])
        conditions = Conditions(**data)
        configure(module, conditions, a.budget)
        anchor = task.get("anchor") if a.task != "design" else None
        available = a.task == "design" or anchor is not None
        fixed = None
        if anchor is not None:
            if a.fixed_side == "forward":
                module.SEQ_ARGS["SEQUENCE_PRIMER"] = anchor["forward_primer"]
                fixed = (anchor["forward_start"], anchor["forward_end"])
            else:
                module.SEQ_ARGS["SEQUENCE_PRIMER_REVCOMP"] = anchor["reverse_primer"]
                fixed = (anchor["reverse_start"], anchor["reverse_end"])
        raw = module.design_primers(SeqRecord(Seq(template), id=gene)) if available else {}
        candidates = []
        for i in range(raw.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            left, llen = raw[f"PRIMER_LEFT_{i}"]
            right, rlen = raw[f"PRIMER_RIGHT_{i}"]
            pair = PrimerPair(left, left + llen, right - rlen + 1, right + 1)
            assert pair.sequences(template) == (raw[f"PRIMER_LEFT_{i}_SEQUENCE"], raw[f"PRIMER_RIGHT_{i}_SEQUENCE"])
            candidate = {**pair.as_dict(template), "structural_valid": True, "gene_id": gene,
                "transcript_accession_version": target["transcript_accession_version"],
                "method": "qprimer", "split":a.split,"proposal_index": i, "candidate_id": f"{name}:{gene}:{i}",
                "primer3_penalty": float(raw[f"PRIMER_PAIR_{i}_PENALTY"]),
                "requested_amplicon_min": conditions.amplicon_min, "requested_amplicon_max": conditions.amplicon_max,
                "assessment": assess_pair(template, pair, conditions)}
            actual = (pair.forward_start, pair.forward_end) if a.fixed_side == "forward" else (pair.reverse_start, pair.reverse_end)
            candidate["anchor_preserved"] = fixed is None or actual == fixed
            if not candidate["anchor_preserved"]:
                candidate["assessment"]["computational_feasible"] = False
                candidate["assessment"]["failures"].append("fixed_binding_interval_changed")
            candidates.append(candidate)
        row = {"gene_id": gene, "method": "qprimer", "split":a.split,"task": a.task, "fixed_side": a.fixed_side,
            "conditions": data, "anchor_available": available,
            **summarize_candidates(candidates, a.budget), "elapsed_seconds": time.monotonic() - started,
            "primer3_explain": {k: v for k, v in raw.items() if k.endswith("EXPLAIN") or k in ("PRIMER_ERROR", "PRIMER_WARNING")}}
        rows.append(row)
        all_candidates.extend(candidates)
        write_jsonl(out / "per_gene.jsonl", rows)
        write_jsonl(out / "candidates.jsonl", all_candidates)
        print(json.dumps({"run": name, "completed_genes": len(rows), "unique_feasible_pairs": row["unique_feasible_pairs"]}), flush=True)
    write_json(out / "completion.json", {"completed_at": utc_now(), "targets": len(rows),
        "all_expected_targets_evaluated": len(rows) == len(targets),
        "per_gene_sha256": sha256(out / "per_gene.jsonl"), "candidates_sha256": sha256(out / "candidates.jsonl")})


if __name__ == "__main__":
    main()
