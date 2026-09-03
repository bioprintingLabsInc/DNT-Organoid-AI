import unittest
from src.qc.metrics import basic_metrics,ratio

class MetricsTests(unittest.TestCase):
 def test_exact_metrics(self):
  m=basic_metrics((0,1,3,6)); self.assertEqual((10,3,1,.25),(m["library_size"],m["detected_genes"],m["zero_count_genes"],m["zero_fraction"])); self.assertEqual((0,2.0,6),(m["minimum"],m["median"],m["maximum"])); self.assertEqual((.75,2.0,3.75),m["quantiles"])
 def test_ratios_and_zero(self): self.assertEqual(2,ratio(10,5)); self.assertIsNone(ratio(0,0))
