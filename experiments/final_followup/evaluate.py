#!/usr/bin/env python3
"""Fixed-seed component and condition-change experiments using existing models."""
import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "experiments/joint_tm/src"), str(ROOT / "scripts")]
import torch
from diffuprimer import guidance
from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.metrics import summarize_candidates
from diffuprimer.model import ConditionalDenoiser
from diffuprimer.sampling import generate
from diffuprimer.sequence import PrimerPair
from span_cache import install

install()
spec = importlib.util.spec_from_file_location("diffuprimer.uniform_compatible_pairs", ROOT / "experiments/uniform_pairs/src/diffuprimer/guidance.py")
uniform = importlib.util.module_from_spec(spec)
spec.loader.exec_module(uniform)


def stable_seed(gene, seed):
    return int(hashlib.sha256(f"{gene}:{seed}".encode()).hexdigest()[:8], 16)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", choices=["ablation", "redesign", "completion"], required=True)
    parser.add_argument("--method", choices=["diffusion", "direct", "primer3", "rule_sampler", "uniform_pair_sampler"], required=True)
    parser.add_argument("--variant", choices=["full", "single_oligo", "geometry_only", "raw"], default="full")
    parser.add_argument("--fixed-side", choices=["forward", "reverse"])
    parser.add_argument("--seed", type=int, default=20261013)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--anchor-prefilter", action="store_true")
    args = parser.parse_args()
    if args.anchor_prefilter:
        decoder_spec = importlib.util.spec_from_file_location("diffuprimer.anchored_decode", ROOT / "experiments/final_followup/anchored_decode.py")
        decoder = importlib.util.module_from_spec(decoder_spec)
        decoder_spec.loader.exec_module(decoder)
        from diffuprimer import sampling
        sampling.sample_grammar = decoder.sample_grammar
    protocol = json.loads((ROOT / "experiments/final_followup/protocol.json").read_text())
    config = json.loads((ROOT / "configs/pilot.json").read_text())
    torch.set_num_threads(args.threads)
    targets = [r for r in read_jsonl(ROOT / "data/processed/targets_pilot.jsonl") if r["split"] == "test"]
    targets.sort(key=lambda r: stable_seed(r["gene_id"], args.seed))
    anchors = {r["gene_id"]: r for r in read_jsonl(ROOT / "experiments/final_followup/redesign_tasks.jsonl")}
    name = f"{args.method}_{args.task}_{args.fixed_side or args.variant}_seed{args.seed}"
    if args.anchor_prefilter:
        name += "_anchor_prefilter"
    out = ROOT / "experiments/final_followup/results/evaluations" / name
    out.mkdir(parents=True, exist_ok=True)
    checkpoint_path = ROOT / f"results/checkpoints/{args.method}_weak_partial_pilot_seed{args.seed}.pt"
    model = checkpoint = None
    if args.method in ["diffusion", "direct"]:
        checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
        model = ConditionalDenoiser(**checkpoint["model_spec"])
        model.load_state_dict(checkpoint["state_dict"])
    write_json(out / "manifest.json", {"method": args.method, "stage": "fixed_seed_final_followup",
        "arguments": {**vars(args), "split": "test", "budget": 20, "temperature": 4.0},
        "task": args.task, "variant": args.variant, "fixed_side": args.fixed_side,
        "protocol_sha256": sha256(ROOT / "experiments/final_followup/protocol.json"),
        "tasks_sha256": sha256(ROOT / "experiments/final_followup/redesign_tasks.jsonl"),
        "target_sha256": sha256(ROOT / "data/processed/targets_pilot.jsonl"),
        "checkpoint_sha256": sha256(checkpoint_path) if model is not None else None,
        "conditions": config["conditions"], "expected_targets": [r["gene_id"] for r in targets]})
    ledger_path, candidates_path = out / "per_gene.jsonl", out / "candidates.jsonl"
    ledger = read_jsonl(ledger_path) if ledger_path.exists() else []
    all_candidates = read_jsonl(candidates_path) if candidates_path.exists() else []
    completed = {r["gene_id"] for r in ledger}
    for target in targets:
        gene = target["gene_id"]
        if gene in completed:
            continue
        started = time.monotonic()
        template = target["template_sequence"]
        task = anchors[gene]
        data = dict(config["conditions"])
        fixed_forward = fixed_reverse = None
        forward_sequence = reverse_sequence = None
        if args.task in ["redesign", "completion"]:
            if args.task == "redesign":
                data.update(amplicon_min=task["new_length_band"][0], amplicon_max=task["new_length_band"][1])
            if task["anchor"] is not None:
                anchor = task["anchor"]
                if args.fixed_side == "forward":
                    fixed_forward = (anchor["forward_start"], anchor["forward_end"])
                    forward_sequence = anchor["forward_primer"]
                else:
                    fixed_reverse = (anchor["reverse_start"], anchor["reverse_end"])
                    reverse_sequence = anchor["reverse_primer"]
        conditions = Conditions(**data)
        span_data = temperatures = None
        if args.variant in ["full", "single_oligo"]:
            span_data = guidance.feasible_oligo_spans(template, conditions)
            if args.variant == "full":
                temperatures = guidance.compatible_span_temperatures(template, conditions, span_data)
        candidates, explain = [], {}
        available = args.task == "ablation" or task["anchor"] is not None
        if available:
            if model is not None:
                proposals = generate(model, template, conditions, 20, checkpoint["steps"], stable_seed(gene, args.seed),
                    args.method, fixed_forward, fixed_reverse, 4.0, constrained=args.variant != "raw",
                    span_allowed=span_data[1] if span_data else None, span_temperatures=temperatures)
                candidates = [{"structural_valid": p["structural_valid"], **(p["pair"].as_dict(template) if p["pair"] else {})} for p in proposals]
            elif args.method == "primer3":
                candidates, explain = primer3_candidates(template, conditions, 20, forward_sequence, reverse_sequence)
                for c in candidates:
                    c["structural_valid"] = True
            else:
                sampler = uniform.rule_span_sampler if args.method == "uniform_pair_sampler" else guidance.rule_span_sampler
                pairs = sampler(template, conditions, 20, stable_seed(gene, args.seed), span_data,
                                fixed_forward, fixed_reverse, temperatures)
                candidates = [{"structural_valid": p is not None, **(p.as_dict(template) if p else {})} for p in pairs]
        for index, candidate in enumerate(candidates):
            candidate.update(gene_id=gene, method=args.method, proposal_index=index,
                transcript_accession_version=target["transcript_accession_version"], candidate_id=f"{name}:{gene}:{index}",
                requested_amplicon_min=conditions.amplicon_min, requested_amplicon_max=conditions.amplicon_max)
            if candidate["structural_valid"]:
                pair = PrimerPair(*(candidate[k] for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"]))
                candidate["assessment"] = assess_pair(template, pair, conditions)
                preserved = ((fixed_forward is None or (pair.forward_start, pair.forward_end) == fixed_forward)
                             and (fixed_reverse is None or (pair.reverse_start, pair.reverse_end) == fixed_reverse))
                candidate["anchor_preserved"] = preserved
                if not preserved:
                    candidate["assessment"]["computational_feasible"] = False
                    candidate["assessment"]["failures"].append("fixed_binding_interval_changed")
            else:
                candidate["anchor_preserved"] = False if args.task == "redesign" else None
                candidate["assessment"] = {"computational_feasible": False, "failures": ["invalid_structure"]}
        row = {"gene_id": gene, "method": args.method, "task": args.task, "variant": args.variant,
            "fixed_side": args.fixed_side, "anchor_available": available,
            "conditions": data, "original_amplicon_length": task.get("original_amplicon_length"),
            "old_pair_violates_new_length": task.get("old_pair_violates_new_length"),
            **summarize_candidates(candidates, 20), "elapsed_seconds": time.monotonic() - started,
            "primer3_explain": explain}
        ledger.append(row)
        all_candidates = [c for c in all_candidates if c["gene_id"] != gene] + candidates
        write_jsonl(candidates_path, all_candidates)
        write_jsonl(ledger_path, ledger)
        print(json.dumps({"run": name, "gene": gene, "completed_genes": len(ledger),
            "unique_feasible_pairs": row["unique_feasible_pairs"], "seconds": row["elapsed_seconds"]}), flush=True)
    write_json(out / "completion.json", {"completed_at": utc_now(), "targets": len(ledger),
        "all_expected_targets_evaluated": len(ledger) == len(targets),
        "per_gene_sha256": sha256(ledger_path), "candidates_sha256": sha256(candidates_path)})


if __name__ == "__main__":
    main()
