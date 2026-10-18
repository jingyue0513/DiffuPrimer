#!/usr/bin/env python3
"""Compare the sampler's pool with independent exhaustive span enumeration."""
import itertools
from pathlib import Path

import numpy as np

from diffuprimer.constraints import Conditions
from diffuprimer.guidance import rule_span_sampler
from diffuprimer.io import sha256, utc_now, write_json


root = Path(__file__).resolve().parents[1]
intervals = np.array([(0,20),(5,25),(60,80),(70,90),(80,100),(1000,1020)])
allowed = [np.array([1,1,0,0,0,1],bool),np.array([0,0,1,1,1,0],bool)]
tm = np.array([[60,61,np.nan,np.nan,np.nan,60],[np.nan,np.nan,58,60,63,np.nan]])
conditions = Conditions(amplicon_min=80,amplicon_max=100,max_tm_difference=2)
comparisons = 0
for f,r in itertools.product([None,(0,20),(5,25),(1000,1020)],[None,(60,80),(70,90),(80,100)]):
    expected = set()
    for i,(fs,fe) in enumerate(intervals):
        for j,(rs,re) in enumerate(intervals):
            if allowed[0][i] and allowed[1][j] and rs >= fe and 80 <= re-fs <= 100 and abs(tm[0,i]-tm[1,j]) <= 2:
                if (f is None or tuple((fs,fe)) == f) and (r is None or tuple((rs,re)) == r):
                    expected.add(tuple(map(int,(fs,fe,rs,re))))
    actual = rule_span_sampler('A'*1020,conditions,1000,20261013,(intervals,allowed),f,r,tm)
    observed = {None if p is None else (p.forward_start,p.forward_end,p.reverse_start,p.reverse_end) for p in actual}
    assert observed == (expected or {None})
    again = rule_span_sampler('A'*1020,conditions,1000,20261013,(intervals,allowed),f,r,tm)
    assert actual == again
    comparisons += 1
write_json(root/'experiments/uniform_pairs/results/check.json',{'at':utc_now(),'independent_pool_comparisons':comparisons,
    'cases':['broad range','fixed forward','fixed reverse','no compatible pair','delta-Tm boundary','inclusive product-length boundaries'],
    'exact_pool_equality':True,'seed_reproducibility':True,
    'source_sha256':sha256(root/'experiments/uniform_pairs/src/diffuprimer/guidance.py')})
print('Uniform sampler agrees with exhaustive compatible-pair pool',flush=True)
