#!/usr/bin/env python3
"""Finalize diagnostics and scientific figures after recorded experiments finish."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import utc_now,write_json


root=Path(__file__).resolve().parents[1]
state_path=root/'results/computational_artifact_status.json'
write_json(state_path,{'stage':'waiting_for_complete_selected_matrix','at':utc_now(),'manuscript_work_paused':True})
expected={root/'results/temperature_tuning/selected_conditions_analysis_status.json':'analyzed',
    root/'results/temperature_tuning/selection_status.json':'selected_test_analyzed',
    root/'results/runtime_benchmark/driver_status.json':'completed'}
if (root/'experiments/epoch100/three_seed_protocol.json').exists():
    expected[root/'experiments/epoch100/results/three_seed_validation_status.json']='validation_comparison_complete'
if (root/'experiments/epoch100/three_seed_evaluation_protocol.json').exists():
    expected[root/'experiments/epoch100/results/followup_status.json']='duration_followup_analyzed'
while True:
    states={str(p):json.loads(p.read_text())['stage'] if p.exists() else 'not_started' for p in expected}
    if all(states[str(p)]==stage for p,stage in expected.items()):break
    if any(s in ['analysis_error','selection_error','upstream_error','duration_validation_error','duration_followup_error'] for s in states.values()):
        write_json(state_path,{'stage':'upstream_error','at':utc_now(),'upstream_states':states})
        raise ValueError('Required experiment analysis failed')
    time.sleep(10)
inputs=json.loads((root/'results/temperature_tuning/selected_conditions_analysis_inputs.json').read_text())['directories']
commands=[['scripts/audit_candidate_integrity.py'], ['scripts/compare_specificity_screens.py'],
    ['scripts/inspect_design_results.py',*inputs,'--screen','bowtie','--output','results/selected_design_diagnostics.json'],
    ['scripts/audit_exon_context.py',*inputs,'--output','results/selected_exon_context_audit.json'],
    ['scripts/plot_completed_study.py'], ['scripts/plot_training_curves.py'], ['scripts/experiment_snapshot.py']]
logs=root/'results/study_logs'
for command in commands:
    label=Path(command[0]).stem
    write_json(state_path,{'stage':label,'at':utc_now(),'manuscript_work_paused':True})
    with (logs/('final_'+label+'.log')).open('w') as handle:
        subprocess.run([str(root/'.venv/bin/python'),*command],cwd=root,env=dict(os.environ,PYTHONPATH=str(root/'src')),
            check=True,stdout=handle,stderr=subprocess.STDOUT)
write_json(state_path,{'stage':'computational_artifacts_complete_independent_assays_pending','at':utc_now(),
    'selected_matrix_runs':len(inputs),'manuscript_work_paused':True})
