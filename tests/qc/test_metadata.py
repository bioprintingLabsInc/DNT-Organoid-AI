import tempfile,unittest
from copy import deepcopy
from pathlib import Path
from src.qc.metadata import reconcile
from .fixtures import metadata,write_matrix

class MetadataTests(unittest.TestCase):
 def test_complete_mapping_and_nonmutation(self):
  with tempfile.TemporaryDirectory() as d: p=write_matrix(Path(d)/"m.tsv"); data=metadata(p); before=deepcopy(data); c,f=reconcile(data,tuple(x["sample_id"] for x in data["samples"]),"counts_asset","bulk_assay")
  self.assertEqual(6,len(c)); self.assertEqual(before,data); self.assertFalse(f)
 def test_unexpected_and_missing(self):
  with tempfile.TemporaryDirectory() as d: data=metadata(write_matrix(Path(d)/"m.tsv")); _,f=reconcile(data,("sample_c1","alien"),"counts_asset","bulk_assay")
  rules={x.rule_id for x in f}; self.assertIn("unexpected_matrix_samples",rules); self.assertIn("missing_expected_samples",rules)
 def test_ambiguous_and_assay_asset_mismatch(self):
  with tempfile.TemporaryDirectory() as d: data=metadata(write_matrix(Path(d)/"m.tsv")); data["samples"].append(dict(data["samples"][0])); data["sequencing_assays"][0]["sample_ids"].remove("sample_c2"); _,f=reconcile(data,("sample_c1","sample_c2"),"counts_asset","wrong")
  rules={x.rule_id for x in f}; self.assertIn("ambiguous_sample_mapping",rules); self.assertIn("assay_missing",rules); self.assertIn("asset_assay_mismatch",rules)
 def test_no_filename_biology(self):
  with tempfile.TemporaryDirectory() as d: data=metadata(write_matrix(Path(d)/"DMSO_treatment.tsv")); c,_=reconcile(data,("sample_t1",),"counts_asset","bulk_assay")
  self.assertEqual("treatment_10um_24h_d45",c["sample_t1"]["condition_id"])
