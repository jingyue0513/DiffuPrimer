#!/usr/bin/env python3
"""Evaluate condition generalization with validation-selected temperatures."""
import json
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from diffuprimer.io import sha256, utc_now, write_json


root=Path(__file__).resolve().parents[1]
output=root/'results/temperature_tuning'
selection_path=output/'selection.json'
status_path=output/'selected_conditions_status.json'
python=str(root/'.venv/bin/python')
write_json(status_path,{'stage':'waiting_for_validation_selection','at':utc_now()})
while not selection_path.exists():
    if json.loads((output/'selection_status.json').read_text())['stage']=='selection_error':
        raise ValueError('Temperature selection requires repair')
    time.sleep(5)
selection=json.loads(selection_path.read_text())
jobs=[];directories=[]
for method in ['direct','diffusion']:
    temperature=selection['selected'][method]['temperature']
    label={0.5:'05',2.0:'2',4.0:'4'}.get(temperature)
    exp=root/('experiments/temperature_'+label if label else 'experiments/joint_tm')
    configurations=[(a,b,None,None) for a,b in [(80,120),(121,160),(161,200)]]
    configurations += [(80,200,role,None) for role in ['forward','reverse']]
    if method=='diffusion':
        configurations += [(80,200,None,steps) for steps in [1,4,8]]
    for seed in [20261013,20261014,20261015]:
        for lower,upper,role,steps in configurations:
            descriptor=f'{method}_weak_partial_pilot_test_seed{seed}_amp{lower}-{upper}'
            command=[python,str(root/'scripts/evaluate_cached_spans.py'),'--evaluation-script',str(exp/'scripts/evaluate.py'),
                '--split','test','--seed',str(seed),'--joint-tm-guidance','--oligo-guidance','--temperature',str(temperature),
                '--checkpoint',str(root/f'results/checkpoints/{method}_weak_partial_pilot_seed{seed}.pt'),
                '--amplicon-min',str(lower),'--amplicon-max',str(upper)]
            if role:
                command += ['--repair',role];descriptor += '_fix_'+role
            descriptor += '_oligo_guided'
            if steps:
                command += ['--steps',str(steps)];descriptor += '_steps'+str(steps)
            descriptor += '_joint_tm'
            directory=exp/'results/evaluations'/descriptor
            jobs.append((descriptor,command,directory,exp,temperature))
            directories.append(str(directory))
protocol={'frozen_at':utc_now(),'selection_sha256':sha256(selection_path),'commands':[j[1] for j in jobs],
    'directories':directories,'criterion':'Use broad-band validation-selected temperatures unchanged for all held-out condition/completion/step tasks. No task-specific retuning.',
    'temperature_1_policy':'Reuse completed owned joint-guidance outputs with identical manifests; no duplicate evaluation',
    'reused_runs':sum(j[-1]==1.0 for j in jobs),'new_runs':sum(j[-1]!=1.0 for j in jobs)}
protocol_path=output/'selected_conditions_protocol.json'
if protocol_path.exists():
    old=json.loads(protocol_path.read_text())
    if old['commands'] != protocol['commands'] or old['selection_sha256'] != protocol['selection_sha256']:
        raise ValueError('Frozen condition protocol changed')
else:
    write_json(protocol_path,protocol)
state={'stage':'selected_condition_sampling','at':utc_now(),'completed':[],'directories':directories}
write_json(status_path,state)


def run(job):
    descriptor,command,directory,exp,temperature=job
    if temperature==1.0:
        while not (directory/'completion.json').exists():
            if json.loads((root/'experiments/joint_tm/results/test_status.json').read_text())['stage']=='test_sampling_error':
                raise ValueError('Owned joint condition run failed')
            time.sleep(5)
    else:
        with (output/('condition_'+descriptor+'.log')).open('a') as handle:
            subprocess.run(command,cwd=root,env=dict(os.environ,PYTHONPATH=str(exp/'src')),check=True,
                stdout=handle,stderr=subprocess.STDOUT)
    manifest=json.loads((directory/'manifest.json').read_text())
    if manifest['arguments']['temperature'] != temperature:
        raise ValueError('Selected condition output uses a different temperature')
    return descriptor


try:
    with ThreadPoolExecutor(max_workers=3) as executor:
        for future in as_completed([executor.submit(run,j) for j in jobs]):
            state['completed'].append(future.result());state['at']=utc_now();write_json(status_path,state)
            print(json.dumps({'completed':state['completed'][-1],'at':state['at']}),flush=True)
    state['stage']='selected_condition_sampling_complete_specificity_pending';write_json(status_path,state)
except Exception as exc:
    state.update(stage='selected_condition_error',error=str(exc));write_json(status_path,state)
    raise
