"""Expanded Human Brain-Organoid Chemical-Exposure Discovery — Data Contract Tests.

Validates:
1. File existence and dual-format parity (CSV and XLSX) across Inventories A, B, C, D, and F.
2. Structure and schema compliance of Inventory A (Study Discovery Inventory).
3. Standardization and partition integrity of Inventory B (Chemical Inventory).
4. Relational integrity, replication status, and matched controls in Inventory C (Conditions Inventory).
5. Exact 22 canonical column conformance, non-presumptive labeling, and zero-fabrication in Inventory D (Sample Inventory).
6. Completeness of Inventory E (Evidence Archive JSON) and Inventory F (Search and Screening Log).
"""

import json
import pathlib
import unittest
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
EXPANDED_METADATA_DIR = REPO_ROOT / "data" / "metadata" / "expanded_discovery"
EXPANDED_OUTPUT_DIR = REPO_ROOT / "outputs" / "verification" / "expanded_discovery"

STUDY_CSV = EXPANDED_METADATA_DIR / "dnt_study_discovery_inventory.csv"
STUDY_XLSX = EXPANDED_METADATA_DIR / "dnt_study_discovery_inventory.xlsx"

CHEM_CSV = EXPANDED_METADATA_DIR / "dnt_chemical_inventory.csv"
CHEM_XLSX = EXPANDED_METADATA_DIR / "dnt_chemical_inventory.xlsx"

COND_CSV = EXPANDED_METADATA_DIR / "dnt_experimental_conditions_inventory.csv"
COND_XLSX = EXPANDED_METADATA_DIR / "dnt_experimental_conditions_inventory.xlsx"

SAMPLE_CSV = EXPANDED_METADATA_DIR / "dnt_sample_inventory.csv"
SAMPLE_XLSX = EXPANDED_METADATA_DIR / "dnt_sample_inventory.xlsx"

EVIDENCE_JSON = EXPANDED_OUTPUT_DIR / "dnt_evidence_archive.json"

LOG_CSV = EXPANDED_OUTPUT_DIR / "dnt_search_screening_log.csv"
LOG_XLSX = EXPANDED_OUTPUT_DIR / "dnt_search_screening_log.xlsx"

PILOT_REPORT_JSON = EXPANDED_OUTPUT_DIR / "pilot_validation_report.json"

EXPECTED_22_SAMPLE_FIELDS = [
    "Dataset ID", "Sample ID", "Species", "Cell Type / Organoid Type", "Cell Line / Donor",
    "Organoid Age / Differentiation Day", "Age at Exposure", "Age at RNA Collection",
    "Compound Name", "Compound ID", "Compound Concentration", "Exposure Duration",
    "Exposure Pattern", "Treatment / Control Group", "Biological Replicates",
    "RNA-seq Type", "Sequencing Platform", "Study / Batch ID", "Gene Expression Data",
    "DNT Outcome Label", "DNT Evidence Source", "DNT Label Confidence"
]

EXPECTED_COND_FIELDS = [
    "Condition ID", "Dataset ID", "Study / BioProject ID", "Collection Tier", "Compound Name",
    "Chemical Partition Group", "Concentration Raw", "Concentration Normalized (uM)",
    "Exposure Duration Raw", "Exposure Duration (Hours)", "Exposure Pattern",
    "Developmental Stage", "Biological Model", "Cell Line / Donor",
    "Sample GSMs", "Replicate Count (n)", "Replication Status",
    "Matched Control GSMs", "Control Status", "Condition Eligibility", "Evidence Source"
]

VALID_COLLECTION_TIERS = {"PRIMARY_BRAIN_ORGANOID", "SECONDARY_NEURAL_MODEL"}
VALID_VERIFICATION_LEVELS = {"PUBLICATION_VERIFIED", "GEO_GSM_VERIFIED", "ACCESSION_ONLY", "NEEDS_REVIEW"}
VALID_SCREENING_DECISIONS = {"INCLUDE", "NEEDS_REVIEW", "EXCLUDE"}
VALID_SEARCH_METHODS = {
    "METHOD_A_BIOLOGICAL_MODEL",
    "METHOD_B_CHEMICAL_SEARCH",
    "METHOD_C_PUBLICATION_SEARCH",
    "METHOD_D_LINKED_ACCESSION_SEARCH"
}


class TestInventoryStructureAndPresence(unittest.TestCase):
    """Ensure all required inventory files exist and CSV/XLSX pairs match."""

    def test_all_inventory_files_exist(self):
        for f in [STUDY_CSV, STUDY_XLSX, CHEM_CSV, CHEM_XLSX, COND_CSV, COND_XLSX,
                  SAMPLE_CSV, SAMPLE_XLSX, EVIDENCE_JSON, LOG_CSV, LOG_XLSX, PILOT_REPORT_JSON]:
            self.assertTrue(f.exists(), f"Missing required file: {f}")

    def test_csv_xlsx_row_parity(self):
        pairs = [
            (STUDY_CSV, STUDY_XLSX),
            (CHEM_CSV, CHEM_XLSX),
            (COND_CSV, COND_XLSX),
            (SAMPLE_CSV, SAMPLE_XLSX),
            (LOG_CSV, LOG_XLSX),
        ]
        for csv_path, xlsx_path in pairs:
            df_csv = pd.read_csv(csv_path)
            df_xlsx = pd.read_excel(xlsx_path)
            self.assertEqual(len(df_csv), len(df_xlsx),
                             f"Row count mismatch between {csv_path.name} and {xlsx_path.name}")
            self.assertEqual(list(df_csv.columns), list(df_xlsx.columns),
                             f"Column mismatch between {csv_path.name} and {xlsx_path.name}")


class TestStudyDiscoveryInventoryContracts(unittest.TestCase):
    """Validate Inventory A: Study Discovery Inventory."""

    def setUp(self):
        self.df = pd.read_csv(STUDY_CSV)

    def test_non_empty(self):
        self.assertGreater(len(self.df), 0, "Study inventory is empty")

    def test_valid_collection_tiers(self):
        invalid = set(self.df["Collection Tier"]) - VALID_COLLECTION_TIERS
        self.assertEqual(len(invalid), 0, f"Invalid collection tiers: {invalid}")

    def test_valid_verification_levels(self):
        invalid = set(self.df["Verification Level"]) - VALID_VERIFICATION_LEVELS
        self.assertEqual(len(invalid), 0, f"Invalid verification levels: {invalid}")

    def test_valid_screening_decisions(self):
        invalid = set(self.df["Screening Decision"]) - VALID_SCREENING_DECISIONS
        self.assertEqual(len(invalid), 0, f"Invalid screening decisions: {invalid}")

    def test_study_accession_uniqueness(self):
        self.assertEqual(len(self.df["Study Accession"]), self.df["Study Accession"].nunique(),
                         "Study accessions are not unique")


class TestChemicalInventoryContracts(unittest.TestCase):
    """Validate Inventory B: Chemical Inventory."""

    def setUp(self):
        self.df = pd.read_csv(CHEM_CSV)

    def test_non_empty(self):
        self.assertGreater(len(self.df), 0, "Chemical inventory is empty")

    def test_pubchem_cid_format(self):
        for cid in self.df["PubChem CID"]:
            is_valid = isinstance(cid, (int, float)) or (isinstance(cid, str) and (cid == "NOT_MAPPED" or cid.isdigit()))
            self.assertTrue(is_valid, f"Invalid PubChem CID format: {cid}")

    def test_chemical_partition_group_format(self):
        for grp in self.df["Chemical Partition Group"]:
            self.assertTrue(grp.startswith("CHEM_"),
                            f"Invalid Chemical Partition Group prefix: {grp}")

    def test_chemical_name_uniqueness(self):
        self.assertEqual(len(self.df["Standardized Chemical Name"]),
                         self.df["Standardized Chemical Name"].nunique(),
                         "Chemical names are not unique in Inventory B")


class TestExperimentalConditionsContracts(unittest.TestCase):
    """Validate Inventory C: Experimental Conditions Inventory."""

    def setUp(self):
        self.df_cond = pd.read_csv(COND_CSV)
        self.df_study = pd.read_csv(STUDY_CSV)
        self.df_chem = pd.read_csv(CHEM_CSV)

    def test_canonical_fields_present(self):
        self.assertEqual(list(self.df_cond.columns), EXPECTED_COND_FIELDS)

    def test_relational_link_to_studies(self):
        study_accs = set(self.df_study["Study Accession"])
        for dset_id in self.df_cond["Dataset ID"]:
            self.assertIn(dset_id, study_accs, f"Condition links to unknown dataset: {dset_id}")

    def test_relational_link_to_chemicals(self):
        chem_groups = set(self.df_chem["Chemical Partition Group"])
        for cgrp in self.df_cond["Chemical Partition Group"]:
            self.assertIn(cgrp, chem_groups, f"Condition has unknown partition group: {cgrp}")

    def test_replication_status_vocabulary(self):
        valid_statuses = {"ADEQUATE_REPLICATION (n>=2)", "SINGLE_REPLICATE (n=1)"}
        for status in self.df_cond["Replication Status"]:
            self.assertIn(status, valid_statuses, f"Invalid replication status: {status}")

    def test_matched_control_integrity(self):
        for _, row in self.df_cond.iterrows():
            if row["Control Status"] == "MATCHED_VEHICLE_CONTROL":
                self.assertTrue(len(str(row["Matched Control GSMs"])) > 0,
                                f"Missing matched control GSMs for condition {row['Condition ID']}")


class TestSampleInventoryContracts(unittest.TestCase):
    """Validate Inventory D: RNA-seq Sample Inventory (Pipeline Compatibility)."""

    def setUp(self):
        self.df_samples = pd.read_csv(SAMPLE_CSV)
        self.df_study = pd.read_csv(STUDY_CSV)

    def test_exact_22_canonical_columns(self):
        self.assertEqual(list(self.df_samples.columns), EXPECTED_22_SAMPLE_FIELDS,
                         "Sample inventory does not match exact 22 canonical columns")

    def test_non_presumptive_outcome_labels(self):
        forbidden_labels = {"DNT_POSITIVE", "DNT_NEGATIVE"}
        found_forbidden = set(self.df_samples["DNT Outcome Label"]).intersection(forbidden_labels)
        self.assertEqual(len(found_forbidden), 0,
                         f"Presumptive labels detected: {found_forbidden}")

    def test_zero_fabrication_pubchem_cid_for_controls(self):
        controls = self.df_samples[self.df_samples["Compound Name"] == "Vehicle / Untreated Control"]
        self.assertTrue((controls["Compound ID"] == "NOT_MAPPED").all(),
                        "Vehicle controls must have Compound ID == 'NOT_MAPPED'")

    def test_relational_link_to_studies(self):
        study_accs = set(self.df_study["Study Accession"])
        for dset_id in self.df_samples["Dataset ID"]:
            self.assertIn(dset_id, study_accs, f"Sample links to unknown dataset: {dset_id}")

    def test_sample_id_uniqueness(self):
        self.assertEqual(len(self.df_samples["Sample ID"]), self.df_samples["Sample ID"].nunique(),
                         "Sample IDs are not unique in Inventory D")


class TestEvidenceArchiveAndSearchLog(unittest.TestCase):
    """Validate Inventory E (Evidence Archive) and Inventory F (Search Log)."""

    def test_evidence_archive_structure(self):
        with open(EVIDENCE_JSON, "r", encoding="utf-8") as f:
            archive = json.load(f)
        self.assertGreater(len(archive), 0)
        required_keys = {
            "study_accession", "collection_tier", "repository_links", "publication",
            "materials_and_methods_excerpt", "verification_level", "screening_decision",
            "chemicals_examined", "conditions_count", "reported_conditions", "unresolved_evidence_gaps"
        }
        for gse, item in archive.items():
            self.assertEqual(set(item.keys()) - {"audited_at"}, required_keys,
                             f"Missing keys in evidence archive for {gse}")

    def test_search_screening_log_methods(self):
        df_log = pd.read_csv(LOG_CSV)
        self.assertGreater(len(df_log), 0)
        invalid_methods = set(df_log["Search Method"]) - VALID_SEARCH_METHODS
        self.assertEqual(len(invalid_methods), 0, f"Invalid search methods: {invalid_methods}")
        for q in df_log["Search Query String"]:
            self.assertGreater(len(str(q)), 5, f"Query string suspiciously short: {q}")


if __name__ == "__main__":
    unittest.main()
