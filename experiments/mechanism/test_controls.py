"""Scientific-control invariants, checked before the frozen study starts."""
import sys
from pathlib import Path
import unittest
import numpy as np
import torch

ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'experiments/hybrid_attention/src'),str(Path(__file__).parent)]
from diffuprimer.constraints import Conditions
from diffuprimer.decode import interval_table, grammar_candidates
from diffuprimer.sequence import PrimerPair
from diffuprimer.sampling import generate
from controls import all_support, iterative, weighted_order, diverse_indices, timed_direct


class CountingModel(torch.nn.Module):
    def __init__(self):
        super().__init__();self.weight=torch.nn.Parameter(torch.tensor(0.3));self.calls=0
    def forward(self,dna,labels,t,conditions):
        self.calls+=1
        pos=torch.arange(dna.shape[1]).float()[None].expand(len(dna),-1)
        return torch.stack([pos*0,pos*-.02+self.weight+(labels==1)*.2,
                            pos*.02+self.weight+(labels==2)*.2],-1)


class ControlsTest(unittest.TestCase):
    def setUp(self):
        self.c=Conditions(primer_min_length=2,primer_max_length=3,amplicon_min=6,amplicon_max=12)
        self.seq='ACGT'*4
        self.intervals=interval_table(len(self.seq),self.c,1)
        n=len(self.intervals)
        self.allowed=[np.arange(n)%3!=0,np.arange(n)%4!=0]
        self.temps=np.vstack([58+np.arange(n)%5,58+np.arange(n)%4])

    def test_exhaustive_support_matches_independent_cartesian_enumeration(self):
        logits=np.random.default_rng(42).normal(size=(len(self.seq),3))
        coords,scores=all_support(logits,self.c,self.allowed,self.temps)
        actual={tuple(map(int,c)):s for c,s in zip(coords,scores)}
        expected={}
        for i,(fs,fe) in enumerate(self.intervals):
            for j,(rs,re) in enumerate(self.intervals):
                if self.allowed[0][i] and self.allowed[1][j] and fe<=rs and 6<=re-fs<=12 and abs(self.temps[0,i]-self.temps[1,j])<=2:
                    expected[(fs,fe,rs,re)]=sum(logits[fs:fe,1]-logits[fs:fe,0])+sum(logits[rs:re,2]-logits[rs:re,0])
        self.assertEqual(set(actual),set(expected))
        np.testing.assert_allclose([actual[k] for k in expected],list(expected.values()),atol=1e-12)

    def test_full_process_matches_existing_sampler_and_frozen_calls_once(self):
        # Broad support avoids artificial failures on the tiny reference.
        n=len(self.intervals);allowed=[np.ones(n,bool),np.ones(n,bool)];temps=np.ones((2,n))*60
        for steps in [1,4,16]:
            model=CountingModel()
            actual,stats=iterative(model,self.seq,self.c,allowed,temps,42,steps,budget=3)
            self.assertEqual(model.calls,steps)
            baseline=generate(CountingModel(),self.seq,self.c,budget=3,steps=steps,seed=42,
                              temperature=4,span_allowed=allowed,span_temperatures=temps)
            self.assertEqual(actual,[r['pair'] for r in baseline])
            frozen=CountingModel()
            result,stats=iterative(frozen,self.seq,self.c,allowed,temps,42,steps,frozen=True,budget=3)
            self.assertEqual(frozen.calls,1)
            self.assertTrue(all(p is not None for p in result))
            if steps==1:self.assertEqual(actual,result)

    def test_without_replacement_is_repeatable_and_full_support(self):
        a=weighted_order(np.array([1.,2.,3.,-10000.]),42,4)
        np.testing.assert_array_equal(a,weighted_order(np.array([1.,2.,3.,-10000.]),42,4))
        self.assertEqual(set(a),{0,1,2,3})
        self.assertEqual(len(weighted_order([],42,20)),0)

    def test_diversity_never_uses_infeasible_and_keeps_distant_candidate(self):
        pairs=[PrimerPair(0,2,8,10),PrimerPair(1,3,9,11),PrimerPair(20,22,28,30),PrimerPair(40,42,48,50)]
        assessments=[{'computational_feasible':i!=3,'shared_temperature_penalty':float(i)} for i in range(4)]
        self.assertEqual(diverse_indices(pairs,assessments,2),[0,2])

    def test_deadline_does_not_accept_late_output(self):
        n=len(self.intervals);allowed=[np.ones(n,bool),np.ones(n,bool)];temps=np.ones((2,n))*60
        # Avoid thermodynamics on the toy 2bp oligos; zero deadline has no valid
        # output even if the minimum non-preemptible operation has completed.
        import controls
        original=controls.assess_pair
        controls.assess_pair=lambda *a:{'computational_feasible':True,'shared_temperature_penalty':0.}
        try:
            pairs,assessments,details=timed_direct(CountingModel(),self.seq,self.c,allowed,temps,42,0)
            self.assertEqual(pairs,[]);self.assertEqual(assessments,[])
            self.assertTrue(all(not r['accepted_before_deadline'] for r in details['trace']))
        finally:controls.assess_pair=original


if __name__=='__main__':unittest.main()
