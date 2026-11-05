"""Integration and unit tests for Phase U3.2 Processed Matrix Quality Control."""

from pathlib import Path
import tempfile
import unittest

from tests.metadata.fixtures import canonical_bulk_study, reported, with_provenance
from src.ingestion.processed_matrix_reader import ProcessedExpressionMatrix, read_processed_matrix
from src.qc.processed_matrix_qc import (
    ProcessedConditionMetrics,
    ProcessedMatrixQCResult,
    ProcessedSampleMetrics,
    reconcile_processed_metadata,
    validate_processed_matrix,
)
from src.qc.result import Severity, Status


def _make_metadata(
    samples: list[dict],
    assays: list[dict] | None = None,
    study_id: str = "study_u3_test",
) -> dict:
    doc = canonical_bulk_study()
    doc["studies"][0]["study_id"] = study_id
    doc["samples"] = samples
    if assays is not None:
        doc["sequencing_assays"] = assays
    else:
        sample_ids = [s["sample_id"] for s in samples]
        doc["sequencing_assays"][0]["sample_ids"] = sample_ids
        doc["input_data_assets"][0]["assay_ids"] = [doc["sequencing_assays"][0]["assay_id"]]
    return with_provenance(doc)


class ProcessedMatrixQCTests(unittest.TestCase):
    """Test Processed Matrix QC: continuous values, metrics, metadata reconciliation, and replicate correlations."""

    def test_continuous_floating_point_matrix_passes_without_integer_errors(self):
        """CRITICAL: Continuous floating-point matrix must NOT trigger raw-count integer errors."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("ENSG00000141510", "ENSG00000012048", "ENSG00000139618"),
            sample_ids=("sample_ctrl_1", "sample_ctrl_2"),
            columns=(
                (12.456, 0.0, 55.789),
                (14.123, 0.512, 53.901),
            ),
        )
        meta = _make_metadata(
            samples=[
                {
                    "sample_id": "sample_ctrl_1",
                    "condition_id": "cond_control",
                    "biological_replicate_id": reported("rep1"),
                    "technical_replicate_id": None,
                },
                {
                    "sample_id": "sample_ctrl_2",
                    "condition_id": "cond_control",
                    "biological_replicate_id": reported("rep2"),
                    "technical_replicate_id": None,
                },
            ]
        )

        qc_res = validate_processed_matrix(matrix, meta)

        self.assertEqual(Status.PASS, qc_res.status)
        self.assertTrue(qc_res.is_valid)
        self.assertEqual(3, qc_res.number_of_genes)
        self.assertEqual(2, qc_res.number_of_samples)

        # Zero integer-count errors
        self.assertFalse(any("integer" in f.rule_id for f in qc_res.findings))
        self.assertFalse(any(f.severity == Severity.ERROR for f in qc_res.findings))

    def test_sample_level_continuous_metrics_calculation(self):
        """Verify accurate computation of continuous metrics (library expression, detected genes, quantiles)."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2", "G3", "G4"),
            sample_ids=("s1",),
            columns=((10.0, 0.0, 20.0, 30.0),),
        )
        meta = _make_metadata(
            samples=[
                {
                    "sample_id": "s1",
                    "condition_id": "c1",
                    "biological_replicate_id": reported("rep1"),
                }
            ]
        )

        qc_res = validate_processed_matrix(matrix, meta)
        s_res = qc_res.samples[0]

        self.assertEqual("s1", s_res.sample_id)
        self.assertEqual(60.0, s_res.total_expression)
        self.assertEqual(3, s_res.detected_genes)
        self.assertEqual(1, s_res.zero_expression_genes)
        self.assertEqual(0.25, s_res.zero_fraction)
        self.assertEqual(0.0, s_res.minimum)
        self.assertEqual(15.0, s_res.median)
        self.assertEqual(30.0, s_res.maximum)
        self.assertEqual(15.0, s_res.mean)

    def test_biological_replicate_spearman_rank_correlation(self):
        """Verify Spearman correlation calculation between biological replicates within condition."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2", "G3", "G4", "G5"),
            sample_ids=("s1", "s2"),
            columns=(
                (1.0, 2.0, 3.0, 4.0, 5.0),
                (2.0, 4.0, 6.0, 8.0, 10.0),  # Perfect monotonic ranking
            ),
        )
        meta = _make_metadata(
            samples=[
                {
                    "sample_id": "s1",
                    "condition_id": "c1",
                    "biological_replicate_id": reported("rep1"),
                },
                {
                    "sample_id": "s2",
                    "condition_id": "c1",
                    "biological_replicate_id": reported("rep2"),
                },
            ]
        )

        qc_res = validate_processed_matrix(matrix, meta)
        self.assertEqual(1, len(qc_res.conditions))
        cond = qc_res.conditions[0]
        self.assertEqual("c1", cond.condition_id)
        self.assertEqual(1, cond.evaluable_pairs)
        self.assertIsNotNone(cond.correlation_summary)
        # Perfect rank correlation
        self.assertAlmostEqual(1.0, cond.correlation_summary[0], places=5)
        self.assertAlmostEqual(1.0, cond.correlation_summary[1], places=5)

    def test_metadata_reconciliation_unexpected_samples(self):
        """Matrix samples absent from canonical metadata emit unexpected_matrix_samples REVIEW finding."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("known_sample", "phantom_sample"),
            columns=((1.0,), (2.0,)),
        )
        meta = _make_metadata(
            samples=[
                {
                    "sample_id": "known_sample",
                    "condition_id": "c1",
                    "biological_replicate_id": reported("rep1"),
                }
            ]
        )

        qc_res = validate_processed_matrix(matrix, meta)
        self.assertEqual(Status.NEEDS_REVIEW, qc_res.status)
        self.assertTrue(
            any(
                f.rule_id == "unexpected_matrix_samples"
                and f.entity_id == "phantom_sample"
                and f.severity == Severity.REVIEW
                for f in qc_res.findings
            )
        )

    def test_metadata_reconciliation_missing_expected_samples(self):
        """Assay sample absent from matrix emits missing_expected_samples REVIEW finding."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("sample1",),
            columns=((1.0,),),
        )
        meta = _make_metadata(
            samples=[
                {"sample_id": "sample1", "condition_id": "c1", "biological_replicate_id": reported("rep1")},
                {"sample_id": "sample2", "condition_id": "c1", "biological_replicate_id": reported("rep2")},
            ]
        )
        # Both sample1 and sample2 are expected in assay
        meta["sequencing_assays"][0]["sample_ids"] = ["sample1", "sample2"]

        qc_res = validate_processed_matrix(matrix, meta)
        self.assertEqual(Status.NEEDS_REVIEW, qc_res.status)
        self.assertTrue(
            any(
                f.rule_id == "missing_expected_samples"
                and f.entity_id == "sample2"
                and f.severity == Severity.REVIEW
                for f in qc_res.findings
            )
        )

    def test_metadata_reconciliation_ambiguous_sample_mapping(self):
        """Duplicate sample records in canonical metadata emit ambiguous_sample_mapping ERROR finding."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("dup_sample",),
            columns=((1.0,),),
        )
        meta = _make_metadata(
            samples=[
                {"sample_id": "dup_sample", "condition_id": "c1", "biological_replicate_id": reported("rep1")},
                {"sample_id": "dup_sample", "condition_id": "c2", "biological_replicate_id": reported("rep2")},
            ]
        )

        qc_res = validate_processed_matrix(matrix, meta)
        self.assertEqual(Status.FAIL, qc_res.status)
        self.assertTrue(
            any(
                f.rule_id == "ambiguous_sample_mapping"
                and f.entity_id == "dup_sample"
                and f.severity == Severity.ERROR
                for f in qc_res.findings
            )
        )

    def test_metadata_reconciliation_assay_mismatch(self):
        """Sample mapped in metadata but absent from assay sample_ids emits assay_sample_mismatch ERROR."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("sample_assigned", "sample_unassigned"),
            columns=((1.0,), (2.0,)),
        )
        meta = _make_metadata(
            samples=[
                {"sample_id": "sample_assigned", "condition_id": "c1", "biological_replicate_id": reported("rep1")},
                {"sample_id": "sample_unassigned", "condition_id": "c2", "biological_replicate_id": reported("rep2")},
            ]
        )
        # Only sample_assigned is in assay
        meta["sequencing_assays"][0]["sample_ids"] = ["sample_assigned"]

        qc_res = validate_processed_matrix(matrix, meta)
        self.assertEqual(Status.FAIL, qc_res.status)
        self.assertTrue(
            any(
                f.rule_id == "assay_sample_mismatch"
                and f.entity_id == "sample_unassigned"
                and f.severity == Severity.ERROR
                for f in qc_res.findings
            )
        )

    def test_zero_only_genes_and_samples(self):
        """Zero-only genes emit INFO and zero-only samples emit REVIEW findings."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("ACTIVE_GENE", "SILENT_GENE"),
            sample_ids=("normal_sample", "blank_sample"),
            columns=(
                (15.2, 0.0),
                (0.0, 0.0),  # blank_sample has 0 for all genes
            ),
        )
        meta = _make_metadata(
            samples=[
                {"sample_id": "normal_sample", "condition_id": "c1", "biological_replicate_id": reported("rep1")},
                {"sample_id": "blank_sample", "condition_id": "c2", "biological_replicate_id": reported("rep2")},
            ]
        )

        qc_res = validate_processed_matrix(matrix, meta)
        self.assertEqual(Status.NEEDS_REVIEW, qc_res.status)
        self.assertTrue(any(f.rule_id == "zero_only_genes" and f.entity_id == "SILENT_GENE" for f in qc_res.findings))
        self.assertTrue(any(f.rule_id == "zero_only_samples" and f.entity_id == "blank_sample" for f in qc_res.findings))

    def test_real_gse232218_human_brain_organoid_tpm_qc(self):
        """Validate real human GSE232218 brain organoid TPM matrix against canonical metadata."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        self.assertTrue(fixture_path.exists(), "GSE232218 TPM fixture missing")

        meta = _make_metadata(
            samples=[
                {"sample_id": "GSM7321024", "condition_id": "cond_vehicle", "biological_replicate_id": reported("rep1")},
                {"sample_id": "GSM7321026", "condition_id": "cond_vpa200", "biological_replicate_id": reported("rep1")},
                {"sample_id": "GSM7321027", "condition_id": "cond_vpa400", "biological_replicate_id": reported("rep1")},
            ],
            study_id="GSE232218",
        )
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026", "GSM7321027"]

        qc_res = validate_processed_matrix(fixture_path, meta)

        self.assertEqual(Status.PASS, qc_res.status)
        self.assertTrue(qc_res.is_valid)
        self.assertEqual(15, qc_res.number_of_genes)
        self.assertEqual(3, qc_res.number_of_samples)
        self.assertEqual(3, len(qc_res.samples))

        # Check sample GSM7321024
        s1 = next(s for s in qc_res.samples if s.sample_id == "GSM7321024")
        self.assertAlmostEqual(933.5086, s1.total_expression, places=3)
        self.assertEqual(15, s1.detected_genes)
        self.assertEqual(0, s1.zero_expression_genes)


if __name__ == "__main__":
    unittest.main()
