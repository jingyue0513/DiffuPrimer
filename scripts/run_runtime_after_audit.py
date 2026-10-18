#!/usr/bin/env python3
"""Wait for model selection and our other jobs, then measure generation."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import utc_now,write_json


root=Path(__file__).resolve().parents[1]
output=root/"results/runtime_benchmark"
output.mkdir(parents=True,exist_ok=True)
write_json(output/"driver_status.json",{'stage':'waiting_for_completed_experiments','at':utc_now()})
paths=[root/'results/bowtie_audit/driver_status.json',root/'results/result_checks_status.json',root/'results/step_sweep/analysis_status.json']
expected=['all_completed_experiments_audited','core_experiments_analyzed','analyzed']
if (root/'experiments/epoch100/three_seed_evaluation_protocol.json').exists():
    paths.append(root/'experiments/epoch100/results/followup_status.json')
    expected.append('duration_followup_analyzed')
while True:
    states=[json.loads(p.read_text())['stage'] if p.exists() else 'not_started' for p in paths]
    if states == expected:
        break
    if any(s in ['audit_error','analysis_error','duration_followup_error'] for s in states):
        write_json(output/"driver_status.json",{'stage':'upstream_error','at':utc_now(),'upstream_states':states})
        raise ValueError('Experiment analysis requires repair before timing')
    time.sleep(10)
write_json(output/"driver_status.json",{'stage':'benchmarking','at':utc_now()})
with (output/'benchmark.log').open('a') as handle:
    subprocess.run([str(root/'.venv/bin/python'),'scripts/check_temperature_selection_consistency.py'],cwd=root,
        env=dict(os.environ,PYTHONPATH=str(root/'src')),check=True,stdout=handle,stderr=subprocess.STDOUT)
    subprocess.run([str(root/'.venv/bin/python'),'-u','scripts/benchmark_generation_runtime.py'],cwd=root,
        env=dict(os.environ,PYTHONPATH=str(root/'experiments/joint_tm/src')),check=True,stdout=handle,stderr=subprocess.STDOUT)
write_json(output/"driver_status.json",{'stage':'completed','at':utc_now()})
