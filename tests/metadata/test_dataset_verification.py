"""Tests for DNT AI dataset verification and master metadata deliverables."""

import json
import pathlib
import unittest
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
INVENTORY_CSV = REPO_ROOT / "data" / "metadata" / "dnt_master_dataset_inventory.csv"
INVENTORY_XLSX = REPO_ROOT / "data" / "metadata" / "dnt_master_dataset_inventory.xlsx"
SAMPLE_METADATA_CSV = REPO_ROOT / "data" / "metadata" / "dnt_sample_metadata_v1.csv"
SAMPLE_METADATA_XLSX = REPO_ROOT / "data" / "metadata" / "dnt_sample_metadata_v1.xlsx"
VERIFICATION_LOG_CSV = REPO_ROOT / "outputs" / "verification" / "master_verification_log.csv"
VERIFICATION_MANIFEST = REPO_ROOT / "outputs" / "verification" / "master_verification_manifest.json"

EXPECTED_22_FIELDS = [
    "Dataset ID", "Sample ID", "Species", "Cell Type / Organoid Type", "Cell Line / Donor",
    "Organoid Age / Differentiation Day", "Age at Exposure", "Age at RNA Collection",
    "Compound Name", "Compound ID", "Compound Concentration", "Exposure Duration",
    "Exposure Pattern", "Treatment / Control Group", "Biological Replicates",
    "RNA-seq Type", "Sequencing Platform", "Study / Batch ID", "Gene Expression Data",
    "DNT Outcome Label", "DNT Evidence Source", "DNT Label Confidence"
]


class DatasetVerificationTests(unittest.TestCase):
    def test_master_inventory_files_exist_and_match(self):
        self.assertTrue(INVENTORY_CSV.exists(), "Master inventory CSV must exist")
        self.assertTrue(INVENTORY_XLSX.exists(), "Master inventory XLSX must exist")
        
        df_csv = pd.read_csv(INVENTORY_CSV)
        self.assertEqual(len(df_csv), 1666, "Master inventory must contain exactly 1,666 records")
        self.assertIn("gse", df_csv.columns)
        self.assertIn("decision", df_csv.columns)
        self.assertIn("biological_model", df_csv.columns)
        self.assertIn("compounds", df_csv.columns)
        
        # Verify decision vocabulary
        valid_decisions = {"INCLUDE", "NEEDS_REVIEW", "EXCLUDE"}
        actual_decisions = set(df_csv["decision"].unique())
        self.assertTrue(actual_decisions.issubset(valid_decisions), f"Invalid decisions: {actual_decisions}")
        
    def test_sample_metadata_conforms_to_22_fields(self):
        self.assertTrue(SAMPLE_METADATA_CSV.exists(), "Sample metadata CSV must exist")
        self.assertTrue(SAMPLE_METADATA_XLSX.exists(), "Sample metadata XLSX must exist")
        
        df = pd.read_csv(SAMPLE_METADATA_CSV)
        self.assertGreater(len(df), 1000, "Should have extracted sample metadata for included datasets")
        self.assertEqual(list(df.columns), EXPECTED_22_FIELDS, "Must match the exact 22 metadata column names in order")
        
        # Check required columns have no null values
        self.assertEqual(df["Dataset ID"].isna().sum(), 0)
        self.assertEqual(df["Sample ID"].isna().sum(), 0)
        self.assertEqual(df["Species"].isna().sum(), 0)
        self.assertTrue((df["Species"] == "Homo sapiens").all())
        
        # Group validity
        valid_groups = {"Control", "Treatment"}
        self.assertTrue(set(df["Treatment / Control Group"].unique()).issubset(valid_groups))
        
        # Biological replicates
        self.assertTrue((df["Biological Replicates"] == "YES").all())
        
        # Outcome labels
        valid_outcomes = {"CONTROL_BASELINE", "DNT_RELEVANT_NEUROTOXICANT", "POTENTIAL_NEUROTOXICANT"}
        self.assertTrue(set(df["DNT Outcome Label"].unique()).issubset(valid_outcomes))
        
    def test_manifest_integrity_and_decision_accounting(self):
        self.assertTrue(VERIFICATION_MANIFEST.exists(), "Manifest must exist")
        with open(VERIFICATION_MANIFEST, "r", encoding="utf-8") as f:
            manifest = json.load(f)
            
        self.assertEqual(manifest["total_datasets_reviewed"], 1666)
        decisions = manifest["decisions"]
        total_decisions = sum(decisions.values())
        self.assertEqual(total_decisions, 1666, "Sum of decision counts must equal 1,666")
        self.assertGreater(decisions["INCLUDE"], 50, "Should include at least 50 high-confidence DNT datasets")
        self.assertGreater(manifest["unique_compounds_identified"], 50, "Should identify diverse compounds")


if __name__ == "__main__":
    unittest.main()
