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

TASK1_CSV = REPO_ROOT / "outputs" / "verification" / "task1_reverified_70_datasets.csv"
TASK2_CSV = REPO_ROOT / "outputs" / "verification" / "task2_investigated_346_datasets.csv"
TASK3_CSV = REPO_ROOT / "outputs" / "verification" / "task3_audited_1250_datasets.csv"
PHASE2_DECISION_LOG_CSV = REPO_ROOT / "outputs" / "verification" / "phase2_decision_log.csv"
PHASE2_DECISION_LOG_XLSX = REPO_ROOT / "outputs" / "verification" / "phase2_decision_log.xlsx"
PHASE2_COMPLETENESS_JSON = REPO_ROOT / "outputs" / "verification" / "phase2_metadata_completeness.json"
PHASE2_COMPLETENESS_CSV = REPO_ROOT / "outputs" / "verification" / "phase2_metadata_completeness.csv"

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
        
        # Verify decision vocabulary (Phase 3 adds REFERENCE_DATASET)
        valid_decisions = {"INCLUDE", "NEEDS_REVIEW", "EXCLUDE", "REFERENCE_DATASET"}
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
        self.assertGreaterEqual(decisions["INCLUDE"], 45, "Should include at least 45 high-confidence DNT datasets in Phase 2")
        self.assertGreater(manifest["unique_compounds_identified"], 50, "Should identify diverse compounds")

    def test_phase2_task_breakdowns_and_accounting(self):
        self.assertTrue(TASK1_CSV.exists(), "Task 1 reverified 70 CSV must exist")
        self.assertTrue(TASK2_CSV.exists(), "Task 2 investigated 346 CSV must exist")
        self.assertTrue(TASK3_CSV.exists(), "Task 3 audited 1250 CSV must exist")
        
        df1 = pd.read_csv(TASK1_CSV)
        df2 = pd.read_csv(TASK2_CSV)
        df3 = pd.read_csv(TASK3_CSV)
        
        self.assertEqual(len(df1), 70, "Task 1 must re-verify exactly 70 datasets")
        self.assertEqual(len(df2), 346, "Task 2 must investigate exactly 346 datasets")
        self.assertEqual(len(df3), 1250, "Task 3 must audit exactly 1,250 datasets")
        self.assertEqual(len(df1) + len(df2) + len(df3), 1666, "Total tasks must sum to exactly 1,666 datasets")

    def test_phase2_metadata_completeness_and_zero_fabrication(self):
        self.assertTrue(PHASE2_COMPLETENESS_JSON.exists(), "Completeness JSON must exist")
        self.assertTrue(PHASE2_COMPLETENESS_CSV.exists(), "Completeness CSV must exist")
        
        with open(PHASE2_COMPLETENESS_JSON, "r", encoding="utf-8") as f:
            comp = json.load(f)
            
        # Verify all 22 canonical fields are evaluated in completeness report
        for field in EXPECTED_22_FIELDS:
            self.assertIn(field, comp, f"Field {field} must be evaluated in completeness report")
        self.assertGreaterEqual(len(comp), 22, "Must evaluate completeness for at least all 22 metadata fields")
        
        # Zero fabrication test: Compound ID must be 0% since PubChem IDs were not invented
        self.assertEqual(comp["Compound ID"]["actual_completeness_pct"], 0.0)
        self.assertEqual(comp["Compound ID"]["valid_count"], 0)
        
        # High-completeness fields
        self.assertEqual(comp["Dataset ID"]["actual_completeness_pct"], 100.0)
        self.assertEqual(comp["Sample ID"]["actual_completeness_pct"], 100.0)
        self.assertEqual(comp["Species"]["actual_completeness_pct"], 100.0)
        self.assertEqual(comp["Compound Name"]["actual_completeness_pct"], 100.0)

    def test_phase2_decision_log(self):
        self.assertTrue(PHASE2_DECISION_LOG_CSV.exists(), "Phase 2 decision log CSV must exist")
        self.assertTrue(PHASE2_DECISION_LOG_XLSX.exists(), "Phase 2 decision log XLSX must exist")
        
        df_log = pd.read_csv(PHASE2_DECISION_LOG_CSV)
        self.assertEqual(len(df_log), 1666, "Decision log must cover all 1,666 records")
        self.assertIn("phase2_decision", df_log.columns)
        self.assertIn("phase2_category", df_log.columns)
        self.assertIn("previous_decision", df_log.columns)


if __name__ == "__main__":
    unittest.main()
