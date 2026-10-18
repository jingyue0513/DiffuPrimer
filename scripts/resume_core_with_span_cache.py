#!/usr/bin/env python3
"""Switch execution caching after all ongoing original BLAST screens finish."""
import json
import os
import signal
import subprocess
import time
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
state_path = root / "results/study_status.json"
while True:
    state = json.loads(state_path.read_text())
    stage = state["stage"]
    if stage == "pipeline_error":
        raise ValueError("Core driver failed before cache transition")
    if state["execution_history"][-1].get("deterministic_span_cache"):
        print("Current driver already uses verified span cache", flush=True)
        raise SystemExit(0)
    if stage.startswith(("condition_bands", "fixed_primer_completion", "decoder_ablations")):
        break
    if stage == "computational_pilot_complete_independent_qpcr_pending":
        print("Core sampling already completed", flush=True)
        raise SystemExit(0)
    time.sleep(3)
rows = []
for line in subprocess.check_output(["ps", "-eo", "pid,ppid,args"], text=True).splitlines()[1:]:
    pid, parent, *command = line.split()
    rows.append((int(pid), int(parent), " ".join(command)))
parents = [pid for pid, _, command in rows if command.startswith(".venv/bin/python -u scripts/run_study.py")
           and Path(f"/proc/{pid}/cwd").resolve() == root]
if len(parents) != 1:
    raise ValueError("Cannot identify one core study driver")
children = set(parents)
for _ in range(3):
    children.update(pid for pid, parent, _ in rows if parent in children)
if any("blastn " in command for pid, _, command in rows if pid in children):
    raise ValueError("Refuse cache transition while an original BLAST search is running")
write_json(root / "results/core_span_cache_transition.json", {"at": utc_now(), "old_stage": stage,
    "old_state_sha256": sha256(state_path), "terminated_processes": sorted(children),
    "scope": "Resume unchanged scientific evaluations with verified deterministic span caching; preserve completed screens and per-gene results"})
os.kill(parents[0], signal.SIGTERM)
for pid in children - set(parents):
    try:
        os.kill(pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
time.sleep(1)
python = str(root / ".venv/bin/python")
with (root / "results/study.log").open("a") as handle:
    subprocess.run([python, "-u", "scripts/run_study.py", "--extended", "--resume", "--workers", "6"],
        cwd=root, env=dict(os.environ, PYTHONPATH=str(root / "src")), check=True, stdout=handle, stderr=subprocess.STDOUT)
