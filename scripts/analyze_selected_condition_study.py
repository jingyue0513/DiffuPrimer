#!/usr/bin/env python3
"""Analyze the complete matched matrix after every required audit finishes."""
import json
import os
import subprocess
import time
from pathlib import Path

from diffuprimer.io import utc_now,write_json


root=Path(__file__).resolve().parents[1]
output=root/'results/temperature_tuning'
state_path=output/'selected_conditions_analysis_status.json'
write_json(state_path,{'stage':'waiting_for_selected_conditions','at':utc_now()})
protocol_path=output/'selected_conditions_protocol.json'
while not protocol_path.exists():time.sleep(5)
directories=json.loads(protocol_path.read_text())['directories']
directories += json.loads((output/'selected_test_protocol.json').read_text())['directories']
for method in ['rule_sampler','uniform_pair_sampler']:
    exp=root/('experiments/uniform_pairs' if method=='uniform_pair_sampler' else 'experiments/joint_tm')
    for seed in [20261013,20261014,20261015]:
        for lower,upper,role in [(80,200,None),(80,120,None),(121,160,None),(161,200,None),(80,200,'forward'),(80,200,'reverse')]:
            suffix=f'{method}_pilot_test_seed{seed}_amp{lower}-{upper}'
            if role:suffix += '_fix_'+role
            suffix += '_oligo_guided_joint_tm'
            directories.append(str(exp/'results/evaluations'/suffix))
for lower,upper,role in [(80,200,None),(80,120,None),(121,160,None),(161,200,None),(80,200,'forward'),(80,200,'reverse')]:
    suffix=f'primer3_pilot_test_seed20261013_amp{lower}-{upper}'
    if role:suffix += '_fix_'+role
    directories.append(str(root/'results/evaluations'/suffix))
assert len(directories)==87 and len(set(directories))==87
write_json(output/'selected_conditions_analysis_inputs.json',{'at':utc_now(),'directories':directories,'runs':87})
while not all((Path(p)/'specificity_bowtie_status.json').exists() for p in directories):
    if json.loads((root/'results/bowtie_audit/driver_status.json').read_text())['stage']=='audit_error':
        write_json(state_path,{'stage':'analysis_error','at':utc_now(),'error':'Upstream reference audit failed'})
        raise ValueError('Upstream reference audit failed')
    time.sleep(10)
write_json(state_path,{'stage':'analyzing','at':utc_now()})
with (output/'selected_conditions_statistics.log').open('w') as handle:
    subprocess.run([str(root/'.venv/bin/python'),'scripts/summarize_specificity_audit.py',*directories,
        '--selection',str(output/'selection.json'),'--selected-conditions','--output','results/temperature_tuning/selected_conditions_summary.json'],
        cwd=root,env=dict(os.environ,PYTHONPATH=str(root/'src')),check=True,stdout=handle,stderr=subprocess.STDOUT)
write_json(state_path,{'stage':'analyzed','at':utc_now(),'runs':87})
