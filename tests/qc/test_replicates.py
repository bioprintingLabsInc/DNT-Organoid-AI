import unittest
from src.qc.replicates import spearman,condition_results
from src.qc.result import SampleResult

def sample(s,r,t=None,c="A"):
 return SampleResult(s,c,r,t,1,1,0,0,0,0,1,(0,0,1),1,1,())

class ReplicateTests(unittest.TestCase):
 def test_spearman(self): self.assertAlmostEqual(1,spearman((1,2,3),(2,4,6))); self.assertAlmostEqual(-1,spearman((1,2,3),(6,4,2)))
 def test_constant_not_evaluable(self): self.assertIsNone(spearman((1,1,1),(1,2,3)))
 def test_all_zero_not_evaluable(self): self.assertIsNone(spearman((0,0,0),(0,0,0)))
 def test_average_ranks_with_ties(self):
  from src.qc.replicates import ranks
  self.assertEqual([1.5,1.5,3.0,4.0],ranks((1,1,2,3)))
  self.assertAlmostEqual(0.8333333333333334,spearman((1,1,2,3),(1,2,2,3)))
 def test_condition_group_three_replicates(self):
  r=condition_results((sample("s1","R1"),sample("s2","R2"),sample("s3","R3")),((1,2,3),(2,3,4),(3,4,5)),"v1")[0]; self.assertEqual((3,3,3,"spearman"),(r.sample_count,r.biological_replicate_count,r.evaluable_pairs,r.correlation_method)); self.assertEqual((1,1,1),r.correlation_summary)
 def test_one_replicate(self): self.assertIn("fewer_than_two_biological_replicates",condition_results((sample("s","R1"),),((1,2),),"v1")[0].unevaluable_reasons)
 def test_technical_replicates_not_aggregated(self):
  r=condition_results((sample("s1","R1","T1"),sample("s2","R1","T2"),sample("s3","R2")),((1,2),(1,2),(2,3)),"v1")[0]; self.assertEqual(2,r.biological_replicate_count); self.assertEqual(0,r.evaluable_pairs); self.assertIn("technical_replicates_require_unapproved_aggregation",r.unevaluable_reasons); self.assertEqual("REVIEW",r.findings[0].severity)
 def test_shared_control_is_one_condition(self):
  results=condition_results((sample("c1","R1",c="control"),sample("c2","R2",c="control")),((1,2),(2,3)),"v1"); self.assertEqual(1,len(results))
