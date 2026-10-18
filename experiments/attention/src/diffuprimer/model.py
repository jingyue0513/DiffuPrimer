"""Conditional Transformer clean-position predictor with rotary attention."""
from __future__ import annotations
import math
import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

MASK=3
DNA_CODES={"A":0,"C":1,"G":2,"T":3,"N":4}
CONDITION_FIELDS=["primer_min_length","primer_max_length","amplicon_min","amplicon_max","tm_min","tm_max","max_tm_difference","gc_min","gc_max"]
CONDITION_SCALES=[25,25,350,350,65,65,5,100,100]
def encode_dna(sequence):return np.array([DNA_CODES.get(b,4) for b in sequence.upper()],dtype=np.int64)
def encode_conditions(c):return np.array([getattr(c,k)/v for k,v in zip(CONDITION_FIELDS,CONDITION_SCALES)],dtype=np.float32)
def rotary(x):
    d=x.shape[-1];pos=torch.arange(x.shape[-2],device=x.device,dtype=torch.float32)
    inv=10000.0**(-torch.arange(0,d,2,device=x.device,dtype=torch.float32)/d)
    angles=pos[:,None]*inv[None,:];cos=angles.cos().to(x.dtype);sin=angles.sin().to(x.dtype)
    even,odd=x[...,0::2],x[...,1::2]
    return torch.stack((even*cos-odd*sin,even*sin+odd*cos),dim=-1).flatten(-2)
class AttentionBlock(nn.Module):
    def __init__(self,width,heads,ff_mult):
        super().__init__();self.heads=heads;self.head_width=width//heads
        self.norm_attention=nn.LayerNorm(width);self.norm_ff=nn.LayerNorm(width)
        self.modulation=nn.Linear(width,2*width)
        self.qkv=nn.Linear(width,3*width);self.output=nn.Linear(width,width)
        self.ff=nn.Sequential(nn.Linear(width,ff_mult*width),nn.GELU(),nn.Linear(ff_mult*width,width))
    def forward(self,x,condition,valid):
        scale,shift=self.modulation(F.silu(condition)).chunk(2,dim=-1)
        h=self.norm_attention(x)*(1+scale[:,None,:])+shift[:,None,:]
        b,l,w=h.shape;q,k,v=self.qkv(h).reshape(b,l,3,self.heads,self.head_width).permute(2,0,3,1,4).unbind(0)
        q,k=rotary(q),rotary(k)
        mask=None if valid is None else valid[:,None,None,:]
        h=F.scaled_dot_product_attention(q,k,v,attn_mask=mask,dropout_p=0.0).transpose(1,2).reshape(b,l,w)
        x=x+self.output(h);x=x+self.ff(self.norm_ff(x))
        return x if valid is None else x*valid[:,:,None]
class ConditionalDenoiser(nn.Module):
    def __init__(self,width=96,layers=4,condition_dim=len(CONDITION_FIELDS),heads=4,ff_mult=2):
        super().__init__()
        if width%heads or (width//heads)%2:raise ValueError('Rotary attention requires an even head width')
        self.spec=dict(width=width,layers=layers,condition_dim=condition_dim,heads=heads,ff_mult=ff_mult)
        self.dna=nn.Embedding(5,width);self.labels=nn.Embedding(4,width)
        self.condition=nn.Sequential(nn.Linear(condition_dim+4,width),nn.SiLU(),nn.Linear(width,width))
        self.blocks=nn.ModuleList([AttentionBlock(width,heads,ff_mult) for _ in range(layers)])
        self.head=nn.Sequential(nn.LayerNorm(width),nn.Linear(width,3))
    def forward(self,dna,noisy_labels,time_fraction,conditions,valid=None):
        t=time_fraction.reshape(-1,1);features=torch.cat([t,t*t,torch.sin(math.pi*t),torch.cos(math.pi*t)],dim=1)
        condition=self.condition(torch.cat([conditions,features],dim=1))
        x=self.dna(dna)+self.labels(noisy_labels)+condition[:,None,:]
        if valid is not None:x=x*valid[:,:,None]
        for block in self.blocks:x=block(x,condition,valid)
        return self.head(x)
def corrupt(labels,t,steps,valid=None,generator=None,fixed=None):
    mask=torch.rand(labels.shape,generator=generator,device=labels.device)<t[:,None]/steps
    if valid is not None:mask &= valid
    if fixed is not None:mask &= ~fixed
    return torch.where(mask,MASK,labels),mask
def diffusion_loss(logits,clean,masked,t,steps,valid):
    ce=F.cross_entropy(logits.transpose(1,2),clean,reduction='none')
    return (((ce*masked*valid).sum(1)/valid.sum(1).clamp(min=1))*steps/t).mean()
def reveal_probability(t):return 1.0/t
