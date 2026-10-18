#!/usr/bin/env python3
"""Generate constrained primer pairs from one FASTA; seed defaults to 42."""
import argparse
import dataclasses
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments/hybrid_attention/src'),str(ROOT/'scripts')]
import torch
from diffuprimer.cnn_model import ConditionalDenoiser
from diffuprimer.constraints import Conditions,assess_pair
from diffuprimer import guidance
from diffuprimer.sampling import generate
from diffuprimer.io import write_json,sha256
from span_cache import install


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fasta',required=True,help='Single-record DNA FASTA, binding template')
    p.add_argument('--output',required=True,help='Output JSON with oligos, coordinates and assessments')
    p.add_argument('--checkpoint',type=Path,default=ROOT/'results/checkpoints/diffusion_weak_partial_pilot_seed20261013.pt',help='CNN checkpoint (historical training seed is file provenance)')
    p.add_argument('--seed',type=int,default=42,help='Constant RNG seed, never derived from time (default 42)')
    p.add_argument('--budget',type=int,default=20,help='Requested proposal slots; shortages remain failures (default 20)')
    p.add_argument('--steps',type=int,default=4,help='Recovery steps: 1 for one-step; 4 is a useful tested short setting; historical main benchmark used 16')
    p.add_argument('--temperature',type=float,default=4.,help='Compatible-pair sampling temperature (default 4)')
    p.add_argument('--amplicon-min',type=int,default=80,help='Minimum product bp (default 80)')
    p.add_argument('--amplicon-max',type=int,default=200,help='Maximum product bp (default 200)')
    p.add_argument('--fixed-forward',help='Observed forward interval START:END, zero-based half-open')
    p.add_argument('--fixed-reverse',help='Observed reverse binding interval START:END on forward template')
    p.add_argument('--threads',type=int,default=4,help='CPU intraop threads (default 4)')
    a=p.parse_args()
    if a.seed<0 or a.budget<1 or a.steps<1 or a.temperature<=0 or a.threads<1:p.error('Invalid seed, budget, steps, temperature or threads')
    lines=Path(a.fasta).read_text().splitlines()
    headers=[line[1:] for line in lines if line.startswith('>')]
    if len(headers)!=1:p.error('Expected exactly one FASTA record')
    sequence=''.join(line.strip() for line in lines if line and not line.startswith('>')).upper()
    if not sequence or set(sequence)-set('ACGTN'):p.error('FASTA must contain A/C/G/T/N')
    def interval(value):
        if value is None:return None
        try:start,end=map(int,value.split(':'))
        except ValueError:p.error('Anchor must be START:END')
        if not 0<=start<end<=len(sequence):p.error('Anchor outside sequence')
        return start,end
    f,r=interval(a.fixed_forward),interval(a.fixed_reverse)
    config=json.loads((ROOT/'configs/pilot.json').read_text())
    conditions=Conditions(**dict(config['conditions'],amplicon_min=a.amplicon_min,amplicon_max=a.amplicon_max))
    if a.amplicon_min<1 or a.amplicon_min>a.amplicon_max:p.error('Invalid amplicon bounds')
    torch.set_num_threads(a.threads);torch.set_num_interop_threads(1)
    torch.manual_seed(a.seed);torch.use_deterministic_algorithms(True)
    ck=torch.load(a.checkpoint,map_location='cpu',weights_only=True)
    model=ConditionalDenoiser(**ck['model_spec']);model.load_state_dict(ck['state_dict']);model.eval()
    install();spans=guidance.feasible_oligo_spans(sequence,conditions)
    temperatures=guidance.compatible_span_temperatures(sequence,conditions,spans)
    mode=ck['mode']
    proposals=generate(model,sequence,conditions,a.budget,a.steps,a.seed,mode,f,r,
        temperature=a.temperature,span_allowed=spans[1],span_temperatures=temperatures)
    rows=[]
    for proposal in proposals:
        pair=proposal['pair'];row={'proposal_index':proposal['proposal_index'],'structural_valid':pair is not None}
        if pair is not None:row.update(pair.as_dict(sequence));row['assessment']=assess_pair(sequence,pair,conditions)
        rows.append(row)
    write_json(Path(a.output),{'reference_id':headers[0],'reference_length':len(sequence),
        'seed':a.seed,'steps':a.steps if mode=='diffusion' else 1,'mode':mode,'checkpoint_sha256':sha256(a.checkpoint),
        'conditions':dataclasses.asdict(conditions),'fixed_forward':f,'fixed_reverse':r,
        'reference_specificity_screen_performed':False,'candidates':rows})
    print(f'Wrote {len(rows)} proposal slots to {a.output}')


if __name__=='__main__':main()
