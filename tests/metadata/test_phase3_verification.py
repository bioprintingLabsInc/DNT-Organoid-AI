"""Phase 3 Comprehensive Scientific Verification Regression Tests.

Validates:
1. GEO RNA-seq modality detection across all assay string variants, spelling, and combination assays.
2. Regression coverage ensuring known Phase 2 included datasets (e.g. GSE119861, GSE147633, GSE150005, GSE166297)
   are never misclassified or excluded due to modality detection bugs.
3. Mass exclusion protection sanity checks.
4. Master inventory structure, 4-category classification vocabulary, and 1,666-record census.
5. 22-field canonical sample metadata schema, zero-fabrication guarantees, and non-presumptive outcome labels.
6. Experimental conditions inventory fields, 9-tuple grouping, control linkage, and chemical partition groups.
"""

import json
import pathlib
import unittest
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
METADATA_DIR = REPO_ROOT / "data" / "metadata"
CACHE_DIR = REPO_ROOT / "outputs" / "verification" / "cache"
PHASE3_DIR = REPO_ROOT / "outputs" / "verification" / "phase3"

INVENTORY_CSV = METADATA_DIR / "dnt_master_dataset_inventory.csv"
SAMPLE_METADATA_CSV = METADATA_DIR / "dnt_sample_metadata_v1.csv"
CONDITIONS_CSV = METADATA_DIR / "dnt_experimental_conditions_inventory.csv"

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

# Known Phase 2 included benchmark records that must NEVER be excluded
KNOWN_INCLUDED_BENCHMARKS = [
    "GSE119861", "GSE147633", "GSE150005", "GSE151996", "GSE166297",
    "GSE231638", "GSE231639", "GSE231640", "GSE287574"
]


class ModalityDetectionAndRegressionTests(unittest.TestCase):
    """Regression tests for RNA-seq modality detection and mass exclusion prevention."""

    def test_modality_detection_string_variants(self):
        """Verify that all GEO assay string variants for high throughput sequencing are correctly recognized."""
        from scripts.verify_phase3_batches import detect_rnaseq_modality

        # 1. Exact NCBI GEO standard string (no hyphen)
        ok, sub, rat = detect_rnaseq_modality("Expression profiling by high throughput sequencing", "Study", "Bulk RNA-seq")
        self.assertTrue(ok, "Expression profiling by high throughput sequencing must be recognized as RNA-seq")
        self.assertEqual(sub, "Bulk RNA-seq")

        # 2. Hyphenated variant
        ok, sub, rat = detect_rnaseq_modality("Expression profiling by high-throughput sequencing", "Study", "RNA-seq")
        self.assertTrue(ok, "Hyphenated variant must be recognized as RNA-seq")

        # 3. Non-coding RNA profiling
        ok, sub, rat = detect_rnaseq_modality("Non-coding RNA profiling by high throughput sequencing", "Study", "ncRNA-seq")
        self.assertTrue(ok, "Non-coding RNA profiling must be recognized as RNA-seq")

        # 4. SuperSeries combination assay (array + HTS)
        ok, sub, rat = detect_rnaseq_modality(
            "Expression profiling by array; Expression profiling by high throughput sequencing",
            "Multi-omics study", "Includes RNA-seq subseries"
        )
        self.assertTrue(ok, "Combination array + HTS assay must retain RNA-seq eligibility")

        # 5. Microarray only (must be excluded)
        ok, sub, rat = detect_rnaseq_modality("Expression profiling by array", "Affymetrix array", "GeneChip expression")
        self.assertFalse(ok, "Array-only assays must NOT be recognized as RNA-seq")
        self.assertEqual(sub, "Microarray / Array")

        # 6. Non-transcriptomic genomic sequencing (ChIP-seq / ATAC-seq)
        ok, sub, rat = detect_rnaseq_modality(
            "Genome binding/occupancy profiling by high throughput sequencing",
            "ChIP-seq of histone marks", "ChIP-seq mapping"
        )
        self.assertFalse(ok, "ChIP-seq without RNA-seq must NOT be recognized as RNA-seq")

    def test_known_included_datasets_not_excluded(self):
        """Ensure all known benchmark datasets are correctly classified and not excluded."""
        self.assertTrue(INVENTORY_CSV.exists(), "Master inventory must exist")
        df_inv = pd.read_csv(INVENTORY_CSV)
        inv_map = dict(zip(df_inv["gse"], df_inv["decision"]))

        for gse in KNOWN_INCLUDED_BENCHMARKS:
            self.assertIn(gse, inv_map, f"Benchmark dataset {gse} missing from inventory")
            # If evaluated, it must not be EXCLUDE
            dec = inv_map[gse]
            self.assertNotEqual(
                dec, "EXCLUDE",
                f"Benchmark dataset {gse} was erroneously classified as EXCLUDE!"
            )

    def test_mass_exclusion_sanity_check(self):
        """Assert that total exclusions do NOT represent a broken mass exclusion (100% or 1666 exclusions)."""
        self.assertTrue(INVENTORY_CSV.exists(), "Master inventory must exist")
        df_inv = pd.read_csv(INVENTORY_CSV)
        total = len(df_inv)
        self.assertEqual(total, 1666, "Must contain exactly 1,666 records")

        excludes = (df_inv["decision"] == "EXCLUDE").sum()
        self.assertLess(
            excludes, total,
            f"FATAL: All {total} datasets were classified as EXCLUDE! Modality detection failed."
        )


class MetadataAndInventoryStructureTests(unittest.TestCase):
    """Validates metadata fields, zero fabrication, and conditions inventory integrity."""

    def test_decision_vocabulary_and_inventory_structure(self):
        df = pd.read_csv(INVENTORY_CSV)
        self.assertEqual(len(df), 1666, "Must contain exactly 1,666 records")
        self.assertIn("gse", df.columns)
        self.assertIn("decision", df.columns)

        decisions = set(df["decision"].dropna().unique())
        self.assertTrue(
            decisions.issubset(VALID_PHASE3_DECISIONS),
            f"Invalid decisions found: {decisions - VALID_PHASE3_DECISIONS}"
        )

        # Check if fully evaluated across all batches
        is_fully_evaluated = not (df["biological_model"] == "PENDING_PHASE3_BATCH").any()
        if is_fully_evaluated:
            counts = df["decision"].value_counts().to_dict()
            self.assertEqual(counts.get("INCLUDE", 0), 48, f"Expected exactly 48 INCLUDE, got {counts.get('INCLUDE', 0)}")
            self.assertEqual(counts.get("NEEDS_REVIEW", 0), 44, f"Expected exactly 44 NEEDS_REVIEW, got {counts.get('NEEDS_REVIEW', 0)}")
            self.assertEqual(counts.get("REFERENCE_DATASET", 0), 332, f"Expected exactly 332 REFERENCE_DATASET, got {counts.get('REFERENCE_DATASET', 0)}")
            self.assertEqual(counts.get("EXCLUDE", 0), 1242, f"Expected exactly 1,242 EXCLUDE, got {counts.get('EXCLUDE', 0)}")

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

        # Outcome labels must NOT be presumptive clinical calls (no DNT_POSITIVE / DNT_NEGATIVE)
        outcome_labels = set(df["DNT Outcome Label"].unique())
        self.assertNotIn("DNT_POSITIVE", outcome_labels, "Must not assign automatic DNT_POSITIVE label")
        self.assertNotIn("DNT_NEGATIVE", outcome_labels, "Must not assign automatic DNT_NEGATIVE label")

    def test_experimental_conditions_inventory(self):
        self.assertTrue(CONDITIONS_CSV.exists(), "Experimental conditions inventory must exist")
        df_cond = pd.read_csv(CONDITIONS_CSV)
        self.assertGreaterEqual(len(df_cond), 1, "Should identify distinct experimental conditions")

        # When all batches are complete, verify final census
        df_inv = pd.read_csv(INVENTORY_CSV)
        is_fully_evaluated = not (df_inv["biological_model"] == "PENDING_PHASE3_BATCH").any()
        if is_fully_evaluated:
            self.assertGreaterEqual(len(df_cond), 200, "Final run must identify >= 200 distinct experimental conditions")

        for field in EXPECTED_CONDITION_FIELDS:
            self.assertIn(field, df_cond.columns, f"Condition field {field} missing")

        # Verify partition group is assigned for every condition
        self.assertEqual(df_cond["Chemical Partition Group"].isna().sum(), 0)
        self.assertTrue((df_cond["Chemical Partition Group"].str.startswith("CHEM_")).all())

        # Verify all conditions have matched controls
        self.assertEqual((df_cond["Control Status"] == "MATCHED_VEHICLE_CONTROL").sum(), len(df_cond))


if __name__ == "__main__":
    unittest.main()
