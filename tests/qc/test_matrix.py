import tempfile,unittest
from pathlib import Path
from src.qc.config import load_rules
from src.qc.matrix import read_matrix
from src.qc.result import Severity

class MatrixTests(unittest.TestCase):
 def check(self,text,rule,severity=Severity.ERROR):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"m.tsv"; p.write_text(text,encoding="utf-8"); _,f=read_matrix(p,load_rules())
  self.assertIn((rule,severity),{(x.rule_id,x.severity) for x in f})
 def test_valid(self):
  self.check("gene\ts1\nG1\t1\n","x",Severity.INFO) if False else None
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"m.tsv"; p.write_text("gene\ts1\nG1\t1\n",encoding="utf-8"); m,f=read_matrix(p,load_rules()); self.assertEqual((1,1),(len(m.gene_ids),len(m.sample_ids))); self.assertFalse(f)
 def test_valid_csv(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"m.csv"; p.write_text("gene,s1\nG1,1\n",encoding="utf-8"); m,f=read_matrix(p,load_rules()); self.assertEqual(("s1",),m.sample_ids); self.assertFalse(f)
 def test_missing_gene_column(self): self.check("\ts1\nG1\t1\n","gene_identifier_column_required")
 def test_no_samples(self): self.check("gene\nG1\n","sample_columns_required")
 def test_nonnumeric(self): self.check("gene\ts1\nG1\tx\n","numeric_counts_required")
 def test_negative(self): self.check("gene\ts1\nG1\t-1\n","non_negative_counts_required")
 def test_fractional(self): self.check("gene\ts1\nG1\t1.2\n","integer_like_counts_required")
 def test_duplicate_samples(self): self.check("gene\ts1\ts1\nG1\t1\t2\n","unique_sample_identifiers_required")
 def test_duplicate_genes(self): self.check("gene\ts1\nG1\t1\nG1\t2\n","duplicate_gene_identifiers",Severity.REVIEW)
 def test_missing_count(self): self.check("gene\ts1\nG1\t\n","missing_counts_forbidden")
 def test_invalid_cells_are_not_imputed(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"m.tsv"; p.write_text("gene\ts1\nG1\t\n",encoding="utf-8"); matrix,_=read_matrix(p,load_rules())
  self.assertIsNone(matrix)
 def test_empty_gene(self): self.check("gene\ts1\n\t1\n","empty_gene_identifiers")
 def test_empty_and_malformed(self): self.check("","malformed_matrix"); self.check("gene\ts1\nG1\n","malformed_row")
 def test_zero_gene_and_sample(self):
  text="gene\ts1\ts2\nG1\t0\t0\nG2\t0\t1\n"; self.check(text,"zero_only_genes",Severity.INFO); self.check(text,"zero_only_samples",Severity.REVIEW)
