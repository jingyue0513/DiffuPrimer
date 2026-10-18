#!/usr/bin/env python3
"""Train conditional Transformer diffusion/direct models on the frozen gene split."""
import argparse,sys,json,random,time,math,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/attention/src'))
import numpy as np
import torch
from torch.nn import functional as F
from torch.utils.data import DataLoader
from diffuprimer.model import ConditionalDenoiser,MASK,corrupt,diffusion_loss
from diffuprimer.training import PrimerDataset
from diffuprimer.io import read_jsonl,sha256,write_json,utc_now
p=argparse.ArgumentParser();p.add_argument('--mode',choices=['diffusion','direct'],required=True);p.add_argument('--epochs',type=int,default=30);p.add_argument('--seed',type=int,default=20261013);p.add_argument('--lr',type=float,default=.001);p.add_argument('--device',default='cuda');a=p.parse_args()
torch.set_num_threads(4);random.seed(a.seed);np.random.seed(a.seed);torch.manual_seed(a.seed)
if torch.cuda.is_available():torch.cuda.manual_seed_all(a.seed)
device=torch.device(a.device);path=ROOT/'data/processed/dataset_partial_pilot.jsonl';manifest=json.loads(path.with_suffix('.manifest.json').read_text());assert manifest['dataset_sha256']==sha256(path)
rows=[r for r in read_jsonl(path) if r['evidence_type'] in {'experimental_success','unverified','computational_candidate'}]
train=PrimerDataset([r for r in rows if r['split']=='train'],512,a.seed,True);val=PrimerDataset([r for r in rows if r['split']=='validation'],512,a.seed+1,False)
loader=DataLoader(train,batch_size=16,shuffle=True,generator=torch.Generator().manual_seed(a.seed));vl=DataLoader(val,batch_size=16)
model=ConditionalDenoiser().to(device);opt=torch.optim.AdamW(model.parameters(),lr=a.lr,weight_decay=.01);steps=16
out=ROOT/'experiments/attention/results';ck=out/'checkpoints'/f'attention_{a.mode}_seed{a.seed}.pt';history=[];best=math.inf;start=time.monotonic()
def objective(raw,generator=None):
 batch={k:v.to(device) for k,v in raw.items()};clean,valid=batch['clean'],batch['valid']
 if a.mode=='diffusion':
  t=torch.randint(1,steps+1,(len(clean),),generator=generator,device=device);noisy,masked=corrupt(clean,t,steps,valid,generator,batch['fixed'])
  logits=model(batch['dna'],noisy,t.float()/steps,batch['conditions'],valid)
  return diffusion_loss(logits.float(),clean,masked,t,steps,valid)
 noisy=torch.full_like(clean,MASK);noisy[batch['fixed']]=clean[batch['fixed']];logits=model(batch['dna'],noisy,torch.ones(len(clean),device=device),batch['conditions'],valid)
 ce=F.cross_entropy(logits.float().transpose(1,2),clean,reduction='none');supervised=valid&~batch['fixed'];return (ce*supervised).sum()/supervised.sum().clamp(min=1)
for epoch in range(a.epochs):
 train.epoch=epoch;model.train();losses=[]
 for batch in loader:
  opt.zero_grad(set_to_none=True)
  with torch.autocast(device_type=device.type,dtype=torch.bfloat16,enabled=device.type=='cuda'):loss=objective(batch)
  if not torch.isfinite(loss):raise ValueError('Nonfinite training loss')
  loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),1);opt.step();losses.append(loss.item())
 model.eval();vals=[];g=torch.Generator(device=device).manual_seed(a.seed+1000)
 with torch.no_grad():
  for batch in vl:
   with torch.autocast(device_type=device.type,dtype=torch.bfloat16,enabled=device.type=='cuda'):vals.append(objective(batch,g).item())
 row=dict(epoch=epoch+1,training_loss=float(np.mean(losses)),validation_loss=float(np.mean(vals)),elapsed_seconds=time.monotonic()-start);history.append(row)
 state=dict(state_dict=model.state_dict(),model_spec=model.spec,architecture='conditional_transformer_rope',mode=a.mode,epoch=epoch+1,seed=a.seed,steps=steps,dataset_sha256=sha256(path),split_sha256=manifest['split_sha256'],stage='attention_fixed_seed',allowed_evidence=['experimental_success','unverified','computational_candidate'])
 if row['validation_loss']<best:best=row['validation_loss'];torch.save(state,ck)
 torch.save(state,ck.with_name(ck.stem+'_last.pt'))
 write_json(out/f'attention_{a.mode}_training.json',dict(run=ck.stem,history=history,train_genes=len(train.genes),validation_genes=len(val.genes),parameters=sum(p.numel() for p in model.parameters()),arguments=vars(a),checkpoint=str(ck),best_validation_loss=best,model_spec=model.spec,device=str(device),torch_version=torch.__version__,dataset_sha256=sha256(path),split_sha256=manifest['split_sha256'],updated_at=utc_now()))
 print(json.dumps(dict(mode=a.mode,**row)),flush=True)
