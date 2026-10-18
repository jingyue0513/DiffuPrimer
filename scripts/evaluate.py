#!/usr/bin/env python3
import argparse
import hashlib
import json
import time
from pathlib import Path

import torch

from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.guidance import feasible_oligo_spans, rule_span_sampler
from diffuprimer.metrics import summarize_candidates
from diffuprimer.model import ConditionalDenoiser
from diffuprimer.sampling import generate
from diffuprimer.sequence import PrimerPair


def stable_seed(gene_id, seed):
    return int(hashlib.sha256(f"{gene_id}:{seed}".encode()).hexdigest()[:8], 16)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", help="Omit for Primer3-only evaluation")
    parser.add_argument("--rule-sampler", action="store_true", help="Uniform forward then reverse span sampler with shared single-oligo constraints")
    parser.add_argument("--oligo-guidance", action="store_true")
    parser.add_argument("--split", choices=["train", "validation", "test"], default="validation")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--budget", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261013)
    parser.add_argument("--raw", action="store_true", help="Disable grammar decoder ablation")
    parser.add_argument("--repair", choices=["forward", "reverse"])
    parser.add_argument("--amplicon-min", type=int, default=80)
    parser.add_argument("--amplicon-max", type=int, default=200)
    parser.add_argument("--temperature", type=float, default=1.0)
    parser.add_argument("--steps", type=int, help="Sampling-step ablation")
    parser.add_argument("--calibration", action="store_true")
    args = parser.parse_args()
    if args.rule_sampler and args.checkpoint:
        parser.error("Rule sampler and learned checkpoint are mutually exclusive")
    if args.raw and args.oligo_guidance:
        parser.error("Oligo guidance requires the grammar decoder")
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "configs/pilot.json").read_text())
    torch.set_num_threads(config["model"]["cpu_threads"])
    target_path = root / "data/processed/targets_pilot.jsonl"
    targets = [row for row in read_jsonl(target_path) if row["split"] == args.split]
    targets.sort(key=lambda x: stable_seed(x["gene_id"], config["seed"]))
    if args.limit:
        if not args.calibration and args.split == "test":
            raise ValueError("Final test evaluations must include all frozen targets")
        targets = targets[:args.limit]
    conditions_data = dict(config["conditions"], amplicon_min=args.amplicon_min, amplicon_max=args.amplicon_max)
    conditions = Conditions(**conditions_data)
    model = checkpoint = None
    if args.checkpoint:
        checkpoint = torch.load(root / args.checkpoint, weights_only=True, map_location="cpu")
        if checkpoint["stage"] == "calibration" and not args.calibration:
            raise ValueError("Calibration checkpoints cannot provide final results")
        split_hash = sha256(root / "data/splits/pilot.jsonl")
        if checkpoint["split_sha256"] != split_hash:
            raise ValueError("Checkpoint split does not match evaluation split")
        model = ConditionalDenoiser(**checkpoint["model_spec"])
        model.load_state_dict(checkpoint["state_dict"])
        method = checkpoint["mode"]
        stage = checkpoint["stage"]
    else:
        method, stage = "rule_sampler" if args.rule_sampler else "primer3", "calibration" if args.calibration else "pilot"
    descriptor = f"{method}_{stage}_{args.split}_seed{args.seed}_amp{args.amplicon_min}-{args.amplicon_max}"
    if args.raw:
        descriptor += "_raw"
    if args.repair:
        descriptor += "_fix_" + args.repair
    if args.oligo_guidance or args.rule_sampler:
        descriptor += "_oligo_guided"
    if args.steps:
        descriptor += f"_steps{args.steps}"
    out_dir = root / "results/evaluations" / descriptor
    out_dir.mkdir(parents=True, exist_ok=True)
    ledger, all_candidates = [], []
    ledger_path, candidates_path = out_dir / "per_gene.jsonl", out_dir / "candidates.jsonl"
    provenance_path = out_dir / "manifest.json"
    provenance = {"arguments": vars(args), "target_sha256": sha256(target_path),
                  "checkpoint_sha256": sha256(root / args.checkpoint) if args.checkpoint else None,
                  "conditions": conditions_data, "method": method, "stage": stage,
                  "expected_targets": [row["gene_id"] for row in targets],
                  "candidate_accounting": "K proposals per target, no replacement for invalid or duplicate model outputs; Primer3 shortage counts against K",
                  "interpretation": "Computational pilot, not qPCR validation. Raw and post-decoder validity are reported separately by ablation."}
    if provenance_path.exists():
        existing = json.loads(provenance_path.read_text())
        if existing != provenance:
            raise ValueError("Output directory already exists with different experiment provenance")
        ledger = read_jsonl(ledger_path) if ledger_path.exists() else []
        all_candidates = read_jsonl(candidates_path) if candidates_path.exists() else []
    else:
        write_json(provenance_path, provenance)
    completed = {row["gene_id"] for row in ledger}
    for target in targets:
        if target["gene_id"] in completed:
            continue
        template = target["template_sequence"]
        started = time.monotonic()
        fixed_forward = fixed_reverse = None
        forward_sequence = reverse_sequence = None
        anchor_status = "not_requested"
        if args.repair:
            anchor_rows, _ = primer3_candidates(template, conditions, budget=1)
            if not anchor_rows:
                anchor_status = "no_primer3_anchor"
            else:
                anchor = anchor_rows[0]
                anchor_status = "primer3_top_pair_anchor"
                if args.repair == "forward":
                    fixed_forward = (anchor["forward_start"], anchor["forward_end"])
                    forward_sequence = anchor["forward_primer"]
                else:
                    fixed_reverse = (anchor["reverse_start"], anchor["reverse_end"])
                    reverse_sequence = anchor["reverse_primer"]
        candidates, explain = [], {}
        span_data = None
        if args.oligo_guidance or args.rule_sampler:
            span_data = feasible_oligo_spans(template, conditions)
        if anchor_status != "no_primer3_anchor":
            if args.rule_sampler:
                pairs = rule_span_sampler(template, conditions, args.budget, stable_seed(target["gene_id"], args.seed),
                                          span_data, fixed_forward, fixed_reverse)
                candidates = [{"structural_valid": pair is not None, **(pair.as_dict(template) if pair else {})} for pair in pairs]
            elif model is None:
                candidates, explain = primer3_candidates(template, conditions, args.budget, forward_sequence, reverse_sequence)
                for candidate in candidates:
                    candidate["structural_valid"] = True
            else:
                generated = generate(model, template, conditions, args.budget, args.steps or checkpoint["steps"],
                    stable_seed(target["gene_id"], args.seed), checkpoint["mode"],
                    fixed_forward, fixed_reverse, args.temperature, constrained=not args.raw,
                    span_allowed=span_data[1] if span_data is not None else None)
                for proposal in generated:
                    candidate = {"proposal_index": proposal["proposal_index"], "structural_valid": proposal["structural_valid"]}
                    if proposal["pair"] is not None:
                        candidate.update(proposal["pair"].as_dict(template))
                    candidates.append(candidate)
        for index, candidate in enumerate(candidates):
            candidate.update({"gene_id": target["gene_id"], "transcript_accession_version": target["transcript_accession_version"],
                              "method": method, "proposal_index": index,
                              "candidate_id": f"{descriptor}:{target['gene_id']}:{index}"})
            if candidate["structural_valid"]:
                pair = PrimerPair(*(candidate[key] for key in ["forward_start", "forward_end", "reverse_start", "reverse_end"]))
                candidate["assessment"] = assess_pair(template, pair, conditions)
                if fixed_forward is not None:
                    assert (pair.forward_start, pair.forward_end) == fixed_forward or model is None
                if fixed_reverse is not None:
                    assert (pair.reverse_start, pair.reverse_end) == fixed_reverse or model is None
            else:
                candidate["assessment"] = {"computational_feasible": False, "failures": ["invalid_structure"],
                                            "specificity_status": "not_assessed", "experimental_status": "not_tested"}
        summary = {"gene_id": target["gene_id"], "method": method, "stage": stage,
                   "anchor_status": anchor_status, **summarize_candidates(candidates, args.budget),
                   "elapsed_seconds": time.monotonic() - started, "primer3_explain": explain}
        ledger.append(summary)
        all_candidates.extend(candidates)
        # Candidate file committed first. A crash before ledger commit may cause a
        # rerun, so remove existing candidates for this gene before re-adding below.
        all_candidates = [row for row in all_candidates if row["gene_id"] != target["gene_id"]] + candidates
        write_jsonl(candidates_path, all_candidates)
        write_jsonl(ledger_path, ledger)
        print(json.dumps(summary), flush=True)
    write_json(out_dir / "completion.json", {"completed_at": utc_now(), "targets": len(ledger),
        "all_expected_targets_evaluated": len(ledger) == len(targets),
        "per_gene_sha256": sha256(ledger_path), "candidates_sha256": sha256(candidates_path)})


if __name__ == "__main__":
    main()
