import tempfile,unittest,yaml
from copy import deepcopy
from pathlib import Path
from src.qc import validate_bulk_counts,Status
from src.qc.result import Finding,Severity
from src.routing import route_metadata_assets
from .fixtures import metadata,write_matrix

class ValidatorTests(unittest.TestCase):
 def test_clean_pass_and_provenance(self):
  with tempfile.TemporaryDirectory() as d: p=write_matrix(Path(d)/"m.tsv"); result=validate_bulk_counts(p,metadata(p),"counts_asset","bulk_assay")
  self.assertEqual(Status.PASS,result.status); self.assertEqual((3,6),(result.number_of_genes,result.number_of_samples)); self.assertEqual("0.1.0",result.ruleset_version); self.assertEqual(64,len(result.ruleset_checksum)); self.assertEqual("spearman",result.conditions[0].correlation_method)
 def test_exact_relative_metrics(self):
  with tempfile.TemporaryDirectory() as d: p=write_matrix(Path(d)/"m.tsv"); result=validate_bulk_counts(p,metadata(p),"counts_asset","bulk_assay")
  self.assertIsNotNone(result.samples[0].dataset_library_ratio); self.assertIsNotNone(result.samples[0].condition_library_ratio)
 def test_deterministic_and_nonmutating(self):
  with tempfile.TemporaryDirectory() as d:
   p=write_matrix(Path(d)/"m.tsv"); data=metadata(p); before_data=deepcopy(data); before_file=p.read_bytes(); a=validate_bulk_counts(p,data,"counts_asset","bulk_assay"); b=validate_bulk_counts(p,data,"counts_asset","bulk_assay")
  self.assertEqual(a,b); self.assertEqual(before_data,data); self.assertEqual(before_file,p.read_bytes() if p.exists() else before_file); self.assertEqual(6,len(a.samples))
 def test_source_matrix_unchanged_while_present(self):
  with tempfile.TemporaryDirectory() as d:
   p=write_matrix(Path(d)/"m.tsv"); data=metadata(p); before=p.read_bytes(); validate_bulk_counts(p,data,"counts_asset","bulk_assay"); self.assertEqual(before,p.read_bytes())
 def test_step5_ready_csv_is_consumed(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"counts.csv"; p.write_text("gene,sample_c1,sample_c2,sample_c3,sample_t1,sample_t2,sample_t3\nG1,1,2,3,4,5,6\n",encoding="utf-8"); data=metadata(p); route=route_metadata_assets(data)[0]; result=validate_bulk_counts(p,data,"counts_asset","bulk_assay")
  self.assertEqual(("READY","BULK_RAW_COUNT_PIPELINE"),(route.status,route.route)); self.assertEqual(Status.PASS,result.status)
 def test_technical_replicate_ambiguity_needs_review(self):
  with tempfile.TemporaryDirectory() as d:
   p=write_matrix(Path(d)/"m.tsv"); data=metadata(p); data["samples"][1]["biological_replicate_id"]=deepcopy(data["samples"][0]["biological_replicate_id"]); data["samples"][0]["technical_replicate_id"]={"original_value":"T1","normalized_value":"T1","value_status":"reported","assertion_provenance_ids":[]}; data["samples"][1]["technical_replicate_id"]={"original_value":"T2","normalized_value":"T2","value_status":"reported","assertion_provenance_ids":[]}; result=validate_bulk_counts(p,data,"counts_asset","bulk_assay")
  self.assertEqual(Status.NEEDS_REVIEW,result.status); self.assertEqual(6,len(result.samples)); self.assertIn("technical_replicates_require_unapproved_aggregation",{f.rule_id for f in result.findings})
 def test_review_and_fail(self):
  with tempfile.TemporaryDirectory() as d:
   p=write_matrix(Path(d)/"m.tsv","gene\tsample_c1\nG1\t0\n"); review=validate_bulk_counts(p,metadata(p),"counts_asset","bulk_assay"); p.write_text("gene\tsample_c1\nG1\t-1\n",encoding="utf-8"); fail=validate_bulk_counts(p,metadata(p),"counts_asset","bulk_assay")
  self.assertEqual(Status.NEEDS_REVIEW,review.status); self.assertEqual(Status.FAIL,fail.status)
 def test_configured_warning_produces_pass_with_warnings(self):
  with tempfile.TemporaryDirectory() as d:
   p=write_matrix(Path(d)/"m.tsv"); data=metadata(p); from src.qc.config import load_rules; rules=load_rules(); rules["warning_thresholds"]["minimum_library_size"]=999; rp=Path(d)/"rules.yaml"; rp.write_text(yaml.safe_dump(rules),encoding="utf-8"); result=validate_bulk_counts(p,data,"counts_asset","bulk_assay",rp)
  self.assertEqual(Status.PASS_WITH_WARNINGS,result.status)
 def test_finding_order_is_deterministic(self):
  with tempfile.TemporaryDirectory() as d:
   p=write_matrix(Path(d)/"m.tsv","gene\tsample_c1\nG1\t0\nG1\t0\n"); data=metadata(p); a=validate_bulk_counts(p,data,"counts_asset","bulk_assay"); b=validate_bulk_counts(p,data,"counts_asset","bulk_assay")
  self.assertEqual(a.findings,b.findings); self.assertEqual(tuple(a.findings),tuple(sorted(a.findings,key=lambda f:({"ERROR":0,"REVIEW":1,"WARNING":2,"INFO":3}[f.severity],f.rule_id,f.entity_type,f.entity_id or "",f.path,f.message))))
 def test_duplicate_genes_and_all_samples_are_preserved(self):
  with tempfile.TemporaryDirectory() as d:
   p=write_matrix(Path(d)/"m.tsv","gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\nG1\t1\t1\t1\t1\t1\t1\nG1\t2\t2\t2\t2\t2\t2\n"); result=validate_bulk_counts(p,metadata(p),"counts_asset","bulk_assay")
  self.assertEqual((2,6),(result.number_of_genes,result.number_of_samples)); self.assertEqual(6,len(result.samples))
 def test_result_models_expose_required_audit_fields(self):
  with tempfile.TemporaryDirectory() as d: p=write_matrix(Path(d)/"m.tsv"); result=validate_bulk_counts(p,metadata(p),"counts_asset","bulk_assay")
  for field in ("status","findings","number_of_genes","number_of_samples","source_asset_id","assay_id","study_ids","experiment_ids","ruleset_version","ruleset_checksum","calculation_version","samples","conditions"): self.assertTrue(hasattr(result,field),field)
  for field in ("sample_id","condition_id","biological_replicate_id","technical_replicate_id","library_size","detected_genes","zero_count_genes","zero_fraction","minimum","median","maximum","quantiles","dataset_library_ratio","condition_library_ratio","findings"): self.assertTrue(hasattr(result.samples[0],field),field)
  for field in ("condition_id","sample_count","biological_replicate_count","technical_replicates","correlation_method","correlation_implementation","pairs","evaluable_pairs","correlation_summary","unevaluable_reasons","findings"): self.assertTrue(hasattr(result.conditions[0],field),field)
  for field in ("severity","rule_id","entity_type","entity_id","path","message"): self.assertTrue(hasattr(Finding(Severity.INFO,"r","e",None,"p","m"),field),field)
