#!/usr/bin/env python3
"""Measure identical recorded designs under a common, explicit timing policy."""
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import torch

from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions, assess_pair
from diffuprimer.guidance import compatible_span_temperatures, feasible_oligo_spans, rule_span_sampler
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl
from diffuprimer.model import ConditionalDenoiser
from diffuprimer.sampling import generate
from diffuprimer.sequence import PrimerPair


parser = argparse.ArgumentParser()
parser.add_argument("--repetitions", type=int, default=3)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
selection_path = root / "results/temperature_tuning/selection.json"
selection = json.loads(selection_path.read_text())
config = json.loads((root / "configs/pilot.json").read_text())
conditions = Conditions(**config["conditions"])
torch.set_num_threads(config["model"]["cpu_threads"])
seed, budget = 20261013, 20


def stable_seed(gene):
    return int(hashlib.sha256(f"{gene}:{seed}".encode()).hexdigest()[:8], 16)


targets = [row for row in read_jsonl(root / "data/processed/targets_pilot.jsonl") if row["split"] == "test"]
targets.sort(key=lambda row: stable_seed(row["gene_id"]))
output = root / "results/runtime_benchmark"
output.mkdir(parents=True, exist_ok=True)
sources, models, temperatures, loads = {}, {}, {}, {}
uniform_path = root / "experiments/uniform_pairs/src/diffuprimer/guidance.py"
spec = importlib.util.spec_from_file_location("diffuprimer.uniform_guidance", uniform_path)
uniform_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(uniform_module)
for method in ["direct", "diffusion", "rule_sampler", "uniform_pair_sampler", "primer3"]:
    if method in ["direct", "diffusion"]:
        temperatures[method] = selection["selected"][method]["temperature"]
        label = {0.5: "05", 2.0: "2", 4.0: "4"}.get(temperatures[method])
        base = root / ("experiments/temperature_" + label if label else "experiments/joint_tm")
        source = base / f"results/evaluations/{method}_weak_partial_pilot_test_seed{seed}_amp80-200_oligo_guided_joint_tm"
        checkpoint_path = root / f"results/checkpoints/{method}_weak_partial_pilot_seed{seed}.pt"
        started = time.perf_counter()
        checkpoint = torch.load(checkpoint_path, weights_only=True, map_location="cpu")
        model = ConditionalDenoiser(**checkpoint["model_spec"])
        model.load_state_dict(checkpoint["state_dict"])
        model.eval()
        models[method] = (model, checkpoint["steps"])
        loads[method] = {"seconds": time.perf_counter()-started, "checkpoint_sha256": sha256(checkpoint_path)}
    elif method == "rule_sampler":
        source = root / f"experiments/joint_tm/results/evaluations/rule_sampler_pilot_test_seed{seed}_amp80-200_oligo_guided_joint_tm"
    elif method == "uniform_pair_sampler":
        source = root / f"experiments/uniform_pairs/results/evaluations/uniform_pair_sampler_pilot_test_seed{seed}_amp80-200_oligo_guided_joint_tm"
    else:
        source = root / f"results/evaluations/primer3_pilot_test_seed{seed}_amp80-200"
    completion = json.loads((source / "completion.json").read_text())
    if not completion["all_expected_targets_evaluated"] or completion["candidates_sha256"] != sha256(source / "candidates.jsonl"):
        raise ValueError("Incomplete or changed reference designs")
    manifest = json.loads((source / "manifest.json").read_text())
    if manifest["conditions"] != config["conditions"]:
        raise ValueError("Runtime constraints differ from recorded experiment")
    own = read_jsonl(source / "candidates.jsonl")
    sources[method] = {"directory": str(source), "candidates_sha256": sha256(source / "candidates.jsonl"),
        "per_gene": {t["gene_id"]: [r for r in own if r["gene_id"] == t["gene_id"]] for t in targets}}

protocol = {"created_at": utc_now(), "targets": [t["gene_id"] for t in targets], "budget": budget,
    "seed": seed, "repetitions": args.repetitions, "cpu_threads": torch.get_num_threads(),
    "temperatures": temperatures, "selection_sha256": sha256(selection_path),
    "target_sha256": sha256(root / "data/processed/targets_pilot.jsonl"),
    "sources": {m: {k:v for k,v in s.items() if k != "per_gene"} for m,s in sources.items()},
    "timing_policy": "One shared cold span/Tm preparation per target; timed generation and pair assessment separately; span preparation charged to guided methods; no model loading or file I/O in per-target times; all methods repeated in rotating order; no specificity search in generation runtime",
    "script_sha256": sha256(Path(__file__)),
    "joint_source_sha256": {p.name: sha256(p) for p in (root / "experiments/joint_tm/src/diffuprimer").glob("*.py")}}
write_json(output / "protocol.json", protocol)
write_json(output / "environment.json", {"at": utc_now(), "platform": platform.platform(),
    "torch": torch.__version__, "torch_threads": torch.get_num_threads(), "cpu_count": os.cpu_count(),
    "affinity": sorted(os.sched_getaffinity(0)), "cpu_max": Path('/sys/fs/cgroup/cpu.max').read_text().strip(),
    "memory_max": Path('/sys/fs/cgroup/memory.max').read_text().strip(), "initial_load_average": os.getloadavg(),
    "model_loading": loads})
records = []
methods = list(sources)
for index, target in enumerate(targets):
    template = target["template_sequence"]
    started = time.perf_counter()
    spans = feasible_oligo_spans(template, conditions)
    tm = compatible_span_temperatures(template, conditions, spans)
    prepare_seconds = time.perf_counter()-started
    for repetition in range(args.repetitions):
        order = (index+repetition) % len(methods)
        for method in methods[order:]+methods[:order]:
            started = time.perf_counter()
            if method in models:
                model, steps = models[method]
                proposals = generate(model, template, conditions, budget, steps, stable_seed(target["gene_id"]),
                    method, temperature=temperatures[method], span_allowed=spans[1], span_temperatures=tm)
                pairs = [p["pair"] for p in proposals]
            elif method == "rule_sampler":
                pairs = rule_span_sampler(template, conditions, budget, stable_seed(target["gene_id"]), spans, span_temperatures=tm)
            elif method == "uniform_pair_sampler":
                pairs = uniform_module.rule_span_sampler(template, conditions, budget, stable_seed(target["gene_id"]), spans, span_temperatures=tm)
            else:
                rows, _ = primer3_candidates(template, conditions, budget)
                pairs = [PrimerPair(*(r[k] for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"])) for r in rows]
            generation_seconds = time.perf_counter()-started
            expected = sources[method]["per_gene"][target["gene_id"]]
            coordinates = lambda r: tuple(r[k] for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"]) if r["structural_valid"] else None
            observed = [tuple(getattr(p,k) for k in ["forward_start", "forward_end", "reverse_start", "reverse_end"]) if p else None for p in pairs]
            if observed != [coordinates(r) for r in expected]:
                raise ValueError(f"Benchmark candidates differ: {method}/{target['gene_id']}")
            started = time.perf_counter()
            assessments = [assess_pair(template,p,conditions) for p in pairs if p is not None]
            assessment_seconds = time.perf_counter()-started
            records.append({"gene_id":target["gene_id"],"template_length":len(template),"method":method,
                "repetition":repetition,"generation_seconds":generation_seconds,"assessment_seconds":assessment_seconds,
                "shared_preparation_seconds":prepare_seconds if method != "primer3" else 0.0,
                "total_seconds_with_preparation":generation_seconds+assessment_seconds+(prepare_seconds if method != "primer3" else 0.0),
                "candidate_coordinates_exactly_match_recorded_run":True,"load_average":os.getloadavg()})
    write_jsonl(output / "per_gene.jsonl",records)
    write_json(output / "status.json", {"stage":"benchmarking","completed_targets":index+1,"expected_targets":len(targets),"at":utc_now()})
    print(json.dumps({"benchmark_completed_targets":index+1,"at":utc_now()}),flush=True)
summary = {}
for method in methods:
    own = [r for r in records if r["method"] == method]
    summary[method] = {field: {"mean":float(np.mean([r[field] for r in own])),"median":float(np.median([r[field] for r in own])),
        "p95":float(np.quantile([r[field] for r in own],.95))} for field in ["generation_seconds","assessment_seconds","shared_preparation_seconds","total_seconds_with_preparation"]}
write_json(output / "summary.json", {"completed_at":utc_now(),"targets":len(targets),"repetitions":args.repetitions,
    "methods":summary,"all_candidate_coordinates_match":True,"protocol_sha256":sha256(output / "protocol.json")})
write_json(output / "status.json", {"stage":"completed","at":utc_now(),"targets":len(targets)})
