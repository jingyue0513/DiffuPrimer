#!/usr/bin/env python3
"""Design primer pairs with the local-encoder attention diffusion checkpoint."""
import argparse,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'experiments/hybrid_attention/src'),str(ROOT/'scripts')]
import torch
from diffuprimer.model import ConditionalDenoiser
from diffuprimer.constraints import Conditions,assess_pair
from diffuprimer import guidance
from diffuprimer.sampling import generate
from span_cache import install
from diffuprimer.io import write_json,sha256

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--fasta',required=True);p.add_argument('--output',required=True)
p.add_argument('--checkpoint',default='experiments/hybrid_attention/results/checkpoints/hybrid_attention_diffusion_seed20261013.pt')
p.add_argument('--budget',type=int,default=20);p.add_argument('--seed',type=int,default=20261013)
p.add_argument('--amplicon-min',type=int,default=80);p.add_argument('--amplicon-max',type=int,default=200)
p.add_argument('--fixed-forward',help='0-based binding interval START:END');p.add_argument('--fixed-reverse',help='0-based binding interval START:END')
p.add_argument('--device',default='cuda' if torch.cuda.is_available() else 'cpu')
a=p.parse_args();torch.set_num_threads(4)
lines=Path(a.fasta).read_text().splitlines();headers=[x[1:] for x in lines if x.startswith('>')]
if len(headers)!=1:raise ValueError('Provide one reference sequence per FASTA')
sequence=''.join(x.strip() for x in lines if not x.startswith('>')).upper()
config=json.loads((ROOT/'configs/pilot.json').read_text());conditions=Conditions(**dict(config['conditions'],amplicon_min=a.amplicon_min,amplicon_max=a.amplicon_max))
def interval(value):
 if value is None:return None
 start,end=map(int,value.split(':'))
 if not 0<=start<end<=len(sequence):raise ValueError('Fixed interval is outside the reference')
 return start,end
f,r=interval(a.fixed_forward),interval(a.fixed_reverse)
ckpath=ROOT/a.checkpoint;checkpoint=torch.load(ckpath,map_location='cpu',weights_only=True);model=ConditionalDenoiser(**checkpoint['model_spec']);model.load_state_dict(checkpoint['state_dict']);model.to(a.device)
install();spans=guidance.feasible_oligo_spans(sequence,conditions);temperatures=guidance.compatible_span_temperatures(sequence,conditions,spans)
proposals=generate(model,sequence,conditions,a.budget,checkpoint['steps'],a.seed,checkpoint['mode'],f,r,temperature=4,span_allowed=spans[1],span_temperatures=temperatures)
candidates=[]
for i,proposal in enumerate(proposals):
 row=dict(proposal_index=i,structural_valid=proposal['structural_valid'])
 if proposal['pair'] is not None:
  row.update(proposal['pair'].as_dict(sequence));row['assessment']=assess_pair(sequence,proposal['pair'],conditions)
 candidates.append(row)
write_json(Path(a.output),dict(reference_id=headers[0],reference_length=len(sequence),architecture=checkpoint['architecture'],checkpoint_sha256=sha256(ckpath),seed=a.seed,requested_budget=a.budget,conditions=dict(config['conditions'],amplicon_min=a.amplicon_min,amplicon_max=a.amplicon_max),fixed_forward=f,fixed_reverse=r,candidates=candidates))
print(f'Wrote {len(candidates)} primer-pair proposals to {a.output}')
