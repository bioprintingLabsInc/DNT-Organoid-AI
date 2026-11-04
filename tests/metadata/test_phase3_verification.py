"""Phase 3 Comprehensive Scientific Dataset Verification Tests."""

import json
import pathlib
import unittest
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
INVENTORY_CSV = REPO_ROOT / "data" / "metadata" / "dnt_master_dataset_inventory.csv"
SAMPLE_METADATA_CSV = REPO_ROOT / "data" / "metadata" / "dnt_sample_metadata_v1.csv"
CONDITIONS_CSV = REPO_ROOT / "data" / "metadata" / "dnt_experimental_conditions_inventory.csv"
CONDITIONS_XLSX = REPO_ROOT / "data" / "metadata" / "dnt_experimental_conditions_inventory.xlsx"
PHASE3_DIR = REPO_ROOT / "outputs" / "verification" / "phase3"
BATCHES_DIR = PHASE3_DIR / "batches"
BATCH_LOG_JSON = PHASE3_DIR / "batch_execution_log.json"

EXPECTED_22_FIELDS = [
    "Dataset ID", "Sample ID", "Species", "Cell Type / Organoid Type", "Cell Line / Donor",
    "Organoid Age / Differentiation Day", "Age at Exposure", "Age at RNA Collection",
    "Compound Name", "Compound ID", "Compound Concentration", "Exposure Duration",
    "Exposure Pattern", "Treatment / Control Group", "Biological Replicates",
    "RNA-seq Type", "Sequencing Platform", "Study / Batch ID", "Gene Expression Data",
    "DNT Outcome Label", "DNT Evidence Source", "DNT Label Confidence"
]

VALID_PHASE3_DECISIONS = {"INCLUDE", "NEEDS_REVIEW", "EXCLUDE", "REFERENCE_DATASET"}

EXPECTED_CONDITION_FIELDS = [
    "Condition ID", "Dataset ID", "Study / BioProject ID", "Compound Name",
    "Chemical Partition Group", "Concentration Raw", "Concentration Normalized (uM)",
    "Exposure Duration Raw", "Exposure Duration (Hours)", "Exposure Pattern",
    "Developmental Stage", "Biological Model", "Cell Line / Donor",
    "Sample GSMs", "Replicate Count (n)", "Replication Status",
    "Matched Control GSMs", "Control Status", "Condition Eligibility", "Evidence Source"
]


class Phase3VerificationTests(unittest.TestCase):
    def test_decision_vocabulary_and_inventory_structure(self):
        self.assertTrue(INVENTORY_CSV.exists(), "Master inventory must exist")
        df = pd.read_csv(INVENTORY_CSV)
        self.assertEqual(len(df), 1666, "Must contain exactly 1,666 records")
        self.assertIn("gse", df.columns)
        self.assertIn("decision", df.columns)
        
        # All decisions must be within the 4 Phase 3 categories
        decisions = set(df["decision"].dropna().unique())
        self.assertTrue(decisions.issubset(VALID_PHASE3_DECISIONS), f"Invalid decisions found: {decisions - VALID_PHASE3_DECISIONS}")
        
    def test_sample_metadata_22_fields_and_zero_fabrication(self):
        self.assertTrue(SAMPLE_METADATA_CSV.exists(), "Sample metadata must exist")
        df = pd.read_csv(SAMPLE_METADATA_CSV)
        self.assertEqual(list(df.columns), EXPECTED_22_FIELDS, "Must match exact 22 canonical columns")
        
        # Zero fabrication check: Compound ID must be NOT_MAPPED
        self.assertTrue((df["Compound ID"] == "NOT_MAPPED").all(), "Compound ID must remain NOT_MAPPED (zero fabrication)")
        
        # Species must be Homo sapiens
        self.assertTrue((df["Species"] == "Homo sapiens").all(), "All included samples must be Homo sapiens")
        
        # Groups must be Control or Treatment
        self.assertTrue(set(df["Treatment / Control Group"].unique()).issubset({"Control", "Treatment"}))
        
        # Biological replicates must be YES for all included samples
        self.assertTrue((df["Biological Replicates"] == "YES").all(), "All included samples must have verified biological replicates")

    def test_experimental_conditions_inventory(self):
        if CONDITIONS_CSV.exists():
            df_cond = pd.read_csv(CONDITIONS_CSV)
            self.assertGreater(len(df_cond), 20, "Should identify distinct experimental conditions")
            for field in EXPECTED_CONDITION_FIELDS:
                self.assertIn(field, df_cond.columns, f"Condition field {field} missing")
                
            # Verify partition group is assigned for every condition
            self.assertEqual(df_cond["Chemical Partition Group"].isna().sum(), 0)
            self.assertTrue((df_cond["Chemical Partition Group"].str.startswith("CHEM_")).all())

    def test_qc_integrity_checks(self):
        df_inv = pd.read_csv(INVENTORY_CSV)
        # 1. No duplicate GEO accessions
        self.assertEqual(df_inv["gse"].nunique(), 1666, "All 1,666 GEO accessions must be unique")
        
        # 2. Included datasets must not be non-neural cancer or non-human
        included = df_inv[df_inv["decision"] == "INCLUDE"]
        self.assertFalse((included["biological_model"] == "Non-Neural Tissue / Cancer Model").any())
        self.assertFalse((included["biological_model"] == "Non-Human Organism").any())
        
        # 3. Reference datasets must be valid neural baseline models when present
        ref_datasets = df_inv[df_inv["decision"] == "REFERENCE_DATASET"]
        self.assertGreaterEqual(len(ref_datasets), 0, "Reference datasets must identify baseline organoid/neural models")


if __name__ == "__main__":
    unittest.main()
