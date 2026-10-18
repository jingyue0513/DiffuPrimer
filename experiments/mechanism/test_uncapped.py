import unittest
from unittest.mock import patch
import numpy as np
import run_uncapped as experiment
from test_controls import CountingModel
from diffuprimer.constraints import Conditions
from diffuprimer.decode import interval_table


class DeadlineTest(unittest.TestCase):
    def setUp(self):
        experiment.MODEL=CountingModel()
        experiment.CONDITIONS=Conditions(primer_min_length=2,primer_max_length=3,amplicon_min=6,amplicon_max=12)
        self.sequence='ACGT'*4
        n=len(interval_table(len(self.sequence),experiment.CONDITIONS,1))
        self.allowed=[np.ones(n,bool),np.ones(n,bool)]
        self.temps=np.ones((2,n))*60
        self.assessment=patch.object(experiment,'assess_pair',lambda *args:{'computational_feasible':True,'shared_temperature_penalty':0.})
        self.assessment.start()

    def tearDown(self):self.assessment.stop()

    def test_zero_budget_never_commits_output(self):
        pairs,_,details=experiment.search(self.sequence,self.allowed,self.temps,42,0.000000001)
        self.assertEqual(pairs,[])
        self.assertEqual(details['accepted_checkpoint_cpu_seconds'],0.)

    def test_finite_support_stops_without_wasting_budget(self):
        pairs,_,details=experiment.search(self.sequence,self.allowed,self.temps,42,100.)
        self.assertEqual(details['stop_reason'],'entire_support_exhausted')
        self.assertIsNone(details['candidate_count_cap'])
        self.assertGreater(len(pairs),0)
        self.assertEqual(len(pairs),len(set(p.sequences(self.sequence) for p in pairs)))
        self.assertLess(details['generation_cpu_seconds'],100.)

    def test_virtual_deadline_only_accepts_on_time_results(self):
        value=[0.]
        def clock():value[0]+=.001;return value[0]
        with patch.object(experiment.time,'process_time',clock):
            pairs,_,details=experiment.search(self.sequence,self.allowed,self.temps,42,.08)
        self.assertLessEqual(details['accepted_checkpoint_cpu_seconds'],.08)
        self.assertTrue(all(r['cpu_seconds']<=.08 for r in details['trace'] if r['accepted_before_deadline']))
        self.assertGreater(len(pairs),0)


if __name__=='__main__':unittest.main()
