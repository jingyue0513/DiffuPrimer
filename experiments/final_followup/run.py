#!/usr/bin/env python3
"""Run the user-authorized fixed-seed additions and common reference assessment."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from diffuprimer.baselines import primer3_candidates
from diffuprimer.constraints import Conditions
from diffuprimer.io import read_jsonl, sha256, utc_now, write_json, write_jsonl

base = ROOT / "experiments/final_followup"
results = base / "results"
results.mkdir(exist_ok=True)
config = json.loads((ROOT / "configs/pilot.json").read_text())
targets = [r for r in read_jsonl(ROOT / "data/processed/targets_pilot.jsonl") if r["split"] == "test"]
python = str(ROOT / ".venv/bin/python")
env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
seed = 20261013
bands = [(80, 120), (121, 160), (161, 200)]


def status(stage, **kw):
    write_json(results / "status.json", {"at": utc_now(), "stage": stage, "seed": seed, **kw})
    print(json.dumps({"stage": stage, **kw}), flush=True)


try:
    tasks_path = base / "redesign_tasks.jsonl"
    if not tasks_path.exists():
        tasks = []
        for target in targets:
            anchors, _ = primer3_candidates(target["template_sequence"], Conditions(**config["conditions"]), budget=1)
            anchor = anchors[0] if anchors else None
            length = anchor["amplicon_length"] if anchor else None
            old_band = next((i for i, (lo, hi) in enumerate(bands) if lo <= length <= hi), 0) if anchor else 0
            new_band = bands[(old_band + 1) % 3]
            tasks.append({"gene_id": target["gene_id"], "anchor": anchor, "original_amplicon_length": length,
                "new_length_band": list(new_band), "old_pair_violates_new_length": bool(anchor and not new_band[0] <= length <= new_band[1])})
        write_jsonl(tasks_path, tasks)
    commands = []
    for method, variant in [("diffusion", "raw"), ("diffusion", "geometry_only"),
                            ("diffusion", "single_oligo"), ("direct", "single_oligo")]:
        commands.append([python, str(base / "evaluate.py"), "--task", "ablation", "--method", method,
                         "--variant", variant, "--seed", str(seed)])
    for side in ["forward", "reverse"]:
        for method in ["primer3", "direct", "diffusion", "rule_sampler", "uniform_pair_sampler"]:
            commands.append([python, str(base / "evaluate.py"), "--task", "redesign", "--method", method,
                             "--fixed-side", side, "--seed", str(seed)])
    protocol_path = base / "protocol.json"
    if not protocol_path.exists():
        write_json(protocol_path, {"declared_at": utc_now(), "seed": seed, "seed_basis": "User explicitly requests a fixed seed for additions",
            "training": "Reuse existing master-seed checkpoints; no new training", "targets": 75, "budget": 20,
            "temperature": 4.0, "temperature_basis": "Existing frozen validation-only selection",
            "redesign_rule": "Retain one side of the broad-band Primer3 top pair, then cyclically change the product-length band: 80-120 -> 121-160 -> 161-200 -> 80-120",
            "anchor_tasks_sha256": sha256(tasks_path), "targets_sha256": sha256(ROOT / "data/processed/targets_pilot.jsonl"),
            "commands": commands, "generation_script_sha256": sha256(base / "evaluate.py"),
            "analysis": "Compare fixed-seed candidates at the gene level; original multi-seed results remain separate",
            "checkpoint_sha256": {m: sha256(ROOT / f"results/checkpoints/{m}_weak_partial_pilot_seed{seed}.pt") for m in ["direct", "diffusion"]}})
    status("generating", expected_new_runs=len(commands))
    # Different result directories; at most two CPU workers share the deterministic span cache.
    active, next_index, finished = [], 0, 0
    while finished < len(commands):
        while len(active) < 2 and next_index < len(commands):
            index = next_index
            handle = (results / f"generation_{index:02d}.log").open("a")
            process = subprocess.Popen(commands[index], cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
            active.append((index, process, handle))
            next_index += 1
        for entry in list(active):
            index, process, handle = entry
            code = process.poll()
            if code is not None:
                handle.close()
                active.remove(entry)
                if code:
                    raise RuntimeError(f"Generation {index} failed; see generation_{index:02d}.log")
                finished += 1
                status("generating", finished_runs=finished, expected_new_runs=len(commands), active_jobs=[p.pid for _, p, _ in active])
        time.sleep(2)
    directories = sorted((results / "evaluations").glob("*"))
    status("screening", new_runs=len(directories))
    with (results / "screening.log").open("a") as handle:
        subprocess.run([python, str(ROOT / "scripts/screen_bowtie.py"), *map(str, directories), "--threads", "4"],
                       cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
    status("followup_experiments_complete", new_runs=len(directories))
except Exception as exc:
    status("error", error=str(exc))
    raise
