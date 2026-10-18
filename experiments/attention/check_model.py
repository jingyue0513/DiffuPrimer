"""Check global label coupling, padding isolation and checkpoint roundtrip."""
from pathlib import Path
import sys,torch,json
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/attention/src'))
from diffuprimer.model import ConditionalDenoiser
from diffuprimer.io import write_json,utc_now

torch.manual_seed(20261013);torch.set_num_threads(2);m=ConditionalDenoiser().eval();dna=torch.randint(0,4,(1,128));labels=torch.full_like(dna,3);c=torch.ones(1,9);t=torch.ones(1)
with torch.no_grad():
 base=m(dna,labels,t,c);changed=labels.clone();changed[:,110:114]=1;out=m(dna,changed,t,c)
 distance_effect=(base[:,5]-out[:,5]).abs().max().item();assert distance_effect>1e-6
 padded_dna=torch.cat([dna,torch.zeros((1,20),dtype=torch.long)],1);padded_labels=torch.cat([labels,torch.full((1,20),3,dtype=torch.long)],1);valid=torch.arange(148)[None,:]<128
 padded=m(padded_dna,padded_labels,t,c,valid);padding_error=(base-padded[:,:128]).abs().max().item();assert padding_error<1e-5
 copy=ConditionalDenoiser(**m.spec);copy.load_state_dict(m.state_dict());copy.eval();assert torch.equal(copy(dna,labels,t,c),base)
loss=m(dna,labels,t,c).square().mean();loss.backward();assert all(p.grad is None or torch.isfinite(p.grad).all() for p in m.parameters())
result=dict(at=utc_now(),passed=True,parameters=sum(p.numel() for p in m.parameters()),global_anchor_label_effect=distance_effect,padding_isolation_max_error=padding_error,checkpoint_roundtrip_exact=True,finite_backward=True)
write_json(ROOT/'experiments/attention/results/model_check.json',result);print(json.dumps(result,indent=2))
