import unittest
from src.qc.rules import evaluate_thresholds
from src.qc.result import SampleResult,Finding,Severity,status_for,Status

S=SampleResult("s","c","r",None,10,2,1,.5,0,1,2,(0,1,2),1,1,())
class RuleTests(unittest.TestCase):
 def test_null_inactive(self): self.assertFalse(evaluate_thresholds(S,{"minimum_library_size":None}))
 def test_warning_configured(self): self.assertEqual(Severity.WARNING,evaluate_thresholds(S,{"minimum_library_size":11})[0].severity)
 def test_status_precedence(self):
  def f(s): return Finding(s,"r","d",None,"p","m")
  self.assertEqual(Status.PASS,status_for([])); self.assertEqual(Status.PASS_WITH_WARNINGS,status_for([f(Severity.WARNING)])); self.assertEqual(Status.NEEDS_REVIEW,status_for([f(Severity.WARNING),f(Severity.REVIEW)])); self.assertEqual(Status.FAIL,status_for([f(Severity.ERROR),f(Severity.REVIEW)])); self.assertEqual(Status.PASS,status_for([f(Severity.INFO)]))
