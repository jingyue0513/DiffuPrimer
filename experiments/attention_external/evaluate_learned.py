#!/usr/bin/env python3
"""Fixed-weight learned-model evaluation on the independent RefSeq gene cohort."""
import argparse,hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT/'experiments/hybrid_attention/src'),str(ROOT/'scripts')]
import torch
from diffuprimer import guidance
from diffuprimer.constraints import Conditions,assess_pair
from diffuprimer.io import read_jsonl,sha256,write_json,write_jsonl,utc_now
from diffuprimer.metrics import summarize_candidates
from diffuprimer.sequence import PrimerPair
from diffuprimer.sampling import generate
from span_cache import install
install()
p=argparse.ArgumentParser();p.add_argument('--method',choices=['hybrid_diffusion','hybrid_direct','cnn_diffusion','cnn_direct','rope_diffusion'],required=True);p.add_argument('--seed',type=int,default=20261013);p.add_argument('--threads',type=int,default=4);a=p.parse_args();torch.set_num_threads(a.threads)
mode='direct' if a.method.endswith('direct') else 'diffusion';config=json.loads((ROOT/'configs/pilot.json').read_text());conditions=Conditions(**config['conditions'])
if a.method.startswith('hybrid'):
 from diffuprimer.model import ConditionalDenoiser
 ck=ROOT/f'experiments/hybrid_attention/results/checkpoints/hybrid_attention_{mode}_seed{a.seed}.pt';canonical='attention_hybrid_'+mode;architecture='conditional_hybrid_transformer_rope'
elif a.method.startswith('cnn'):
 from diffuprimer.cnn_model import ConditionalDenoiser
 ck=ROOT/f'results/checkpoints/{mode}_weak_partial_pilot_seed{a.seed}.pt';canonical='cnn_'+mode;architecture='conditional_dilated_cnn'
else:
 import importlib.util
 spec=importlib.util.spec_from_file_location('rope_model',ROOT/'experiments/attention/src/diffuprimer/model.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);ConditionalDenoiser=module.ConditionalDenoiser
 ck=ROOT/f'experiments/attention/results/checkpoints/attention_{mode}_seed{a.seed}.pt';canonical='attention_rope_'+mode;architecture='conditional_transformer_rope'
checkpoint=torch.load(ck,map_location='cpu',weights_only=True);model=ConditionalDenoiser(**checkpoint['model_spec']);model.load_state_dict(checkpoint['state_dict']);model.to('cuda')
target_path=ROOT/'data/processed/attention_external_targets.jsonl';targets=read_jsonl(target_path)
def stable(g):return int(hashlib.sha256(f'{g}:{a.seed}'.encode()).hexdigest()[:8],16)
targets.sort(key=lambda r:stable(r['gene_id']));name=f'{canonical}_external_test_amp80-200_seed{a.seed}';out=ROOT/'experiments/attention_external/results/evaluations'/name;out.mkdir(parents=True,exist_ok=True)
write_json(out/'manifest.json',dict(method=canonical,architecture=architecture,mode=mode,stage='independent_external_gene_test',task='design',variant='full',fixed_side=None,arguments=dict(split='external_test',seed=a.seed,budget=20,amplicon_min=80,amplicon_max=200,steps=16,temperature=4.0),conditions=config['conditions'],target_sha256=sha256(target_path),checkpoint_sha256=sha256(ck),expected_targets=[t['gene_id'] for t in targets],runtime_versions=dict(torch=torch.__version__,numpy=__import__('numpy').__version__,primer3=__import__('primer3').__version__)))
ledger_path=out/'per_gene.jsonl';candidates_path=out/'candidates.jsonl';ledger=read_jsonl(ledger_path) if ledger_path.exists() else [];allc=read_jsonl(candidates_path) if candidates_path.exists() else [];done={r['gene_id'] for r in ledger}
for target in targets:
 gene=target['gene_id']
 if gene in done:continue
 started=time.monotonic();template=target['template_sequence'];span=guidance.feasible_oligo_spans(template,conditions);temps=guidance.compatible_span_temperatures(template,conditions,span)
 proposals=generate(model,template,conditions,20,16,stable(gene),mode,temperature=4,span_allowed=span[1],span_temperatures=temps)
 candidates=[]
 for i,proposal in enumerate(proposals):
  c=dict(gene_id=gene,method=canonical,proposal_index=i,candidate_id=f'{name}:{gene}:{i}',transcript_accession_version=target['transcript_accession_version'],structural_valid=proposal['structural_valid'],requested_amplicon_min=80,requested_amplicon_max=200)
  if proposal['pair'] is not None:c.update(proposal['pair'].as_dict(template));c['assessment']=assess_pair(template,proposal['pair'],conditions)
  else:c['assessment']=dict(computational_feasible=False,failures=['invalid_structure'])
  candidates.append(c)
 row=dict(gene_id=gene,method=canonical,task='design',variant='full',fixed_side=None,conditions=config['conditions'],**summarize_candidates(candidates,20),elapsed_seconds=time.monotonic()-started)
 ledger.append(row);allc=[c for c in allc if c['gene_id']!=gene]+candidates;write_jsonl(candidates_path,allc);write_jsonl(ledger_path,ledger)
 print(json.dumps(dict(run=name,gene=gene,completed_genes=len(ledger),unique_feasible_pairs=row['unique_feasible_pairs'],seconds=row['elapsed_seconds'])),flush=True)
write_json(out/'completion.json',dict(completed_at=utc_now(),targets=len(ledger),all_expected_targets_evaluated=len(ledger)==len(targets),candidates_sha256=sha256(candidates_path),per_gene_sha256=sha256(ledger_path)))
