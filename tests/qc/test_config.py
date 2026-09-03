import tempfile,unittest
from pathlib import Path
from src.qc.config import load_rules,rules_checksum

class ConfigTests(unittest.TestCase):
 def test_defaults_and_nulls(self):
  r=load_rules(); self.assertEqual("spearman",r["replicate_correlation"]["method"]); self.assertTrue(all(v is None for v in r["warning_thresholds"].values()))
 def test_malformed_and_unknown(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"x.yaml"; p.write_text("- bad",encoding="utf-8")
   with self.assertRaises(ValueError): load_rules(p)
   p.write_text("unknown: true",encoding="utf-8")
   with self.assertRaises(ValueError): load_rules(p)
 def test_checksum_semantic(self):
  a={"b":2,"a":1}; b={"a":1,"b":2}; self.assertEqual(rules_checksum(a),rules_checksum(b))
