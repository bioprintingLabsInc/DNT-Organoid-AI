"""Acceptance and integration tests for Phase U3.2: Processed Expression Matrix Intake & Provenance Gate."""

from pathlib import Path
import tempfile
import unittest

from src.ingestion import inspect_asset
from src.ingestion.models import Confidence, RegisteredAsset, Route
from src.normalization.provenance_gate import (
    ExpressionUnit,
    ProcessedExpressionDataset,
    TransformationScale,
    evaluate_provenance_gate,
    ingest_processed_expression_matrix,
)
from src.qc.matrix import read_matrix
from src.qc.processed_matrix_qc import validate_processed_matrix
from src.qc.result import Status as QCStatus
from src.qc.validator import validate_bulk_counts
from src.routing import dispatch_route_execution, plan_route_execution
from src.routing.detector import detect_format
from src.routing.router import _select_route
from tests.routing.fixtures import routing_metadata, write_asset


class PhaseU32ProcessedExpressionGateAcceptanceTests(unittest.TestCase):
    """Verify Phase U3.2 Processed Expression Matrix Intake & Provenance Gate acceptance requirements."""

    def test_detection_and_routing_of_processed_expression_matrix(self):
        """Tabular file with floating-point expression values is detected as processed_expression_matrix."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        self.assertTrue(fixture_path.exists())

        insp = inspect_asset(fixture_path)
        asset = RegisteredAsset(
            "asset_gse232218",
            str(fixture_path),
            "processed_expression_matrix",
            ("assay_gse232218",),
            ("bulk_rna_seq",),
        )
        det = detect_format(asset, insp)

        self.assertEqual("processed_expression_matrix", det.detected_format)
        self.assertEqual(Confidence.HIGH, det.confidence)

        route, findings = _select_route(asset, det)
        self.assertEqual(Route.PROCESSED_EXPRESSION_PIPELINE, route)
        self.assertEqual((), findings)

    def test_u3_2_boundary_plan_route_execution_stops_at_execution_not_implemented(self):
        """CRITICAL: U3.2 scope boundary leaves dispatcher execution at EXECUTION_NOT_IMPLEMENTED (U3.4 boundary)."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = routing_metadata(fixture_path, "processed_expression_matrix")

        plan = plan_route_execution(meta)

        self.assertEqual(Route.PROCESSED_EXPRESSION_PIPELINE, plan.route)
        self.assertEqual("EXECUTION_NOT_IMPLEMENTED", plan.status)
        self.assertIsNone(plan.executor_id)
        self.assertIsNone(plan.downstream_route)

    def test_prevention_processed_matrix_cannot_enter_raw_count_pipeline(self):
        """CRITICAL: Floating-point processed matrices are rejected by raw-count QC (Step 6A) with integer error."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = routing_metadata(fixture_path, "raw_counts")

        raw_qc_res = validate_bulk_counts(
            fixture_path,
            meta,
            "asset_gse232218",
            meta["sequencing_assays"][0]["assay_id"],
        )

        self.assertEqual(QCStatus.FAIL, raw_qc_res.status)
        self.assertTrue(any(f.rule_id == "integer_like_counts_required" for f in raw_qc_res.findings))

    def test_u3_2_full_intake_qc_and_provenance_gating_success(self):
        """U3.2 successfully accepts real human GEO brain organoid TPM matrix and establishes full provenance."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = routing_metadata(fixture_path, "processed_expression_matrix")

        # Configure canonical metadata samples matching fixture
        meta["samples"] = [
            {"sample_id": "GSM7321024", "condition_id": "cond_vehicle", "biological_replicate_id": "rep1"},
            {"sample_id": "GSM7321026", "condition_id": "cond_vpa200", "biological_replicate_id": "rep1"},
            {"sample_id": "GSM7321027", "condition_id": "cond_vpa400", "biological_replicate_id": "rep1"},
        ]
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026", "GSM7321027"]
        meta["processing_provenance"] = {
            "expression_unit": "TPM",
            "transformation_scale": "linear",
            "normalization_method": "tpm",
            "software_pipeline": "Kallisto v0.48.0",
            "reference_annotation": "Ensembl Release 112",
        }

        dataset = ingest_processed_expression_matrix(fixture_path, meta)

        self.assertEqual(QCStatus.PASS, dataset.status)
        self.assertTrue(dataset.is_inferentially_eligible)
        self.assertEqual("limma_trend", dataset.target_analysis_engine)
        self.assertEqual(("DESeq2",), dataset.blocked_analysis_engines)

        # Numerical preservation
        self.assertEqual(42.1582, dataset.matrix.get_value("ENSG00000141510", "GSM7321024"))
        self.assertEqual(15, dataset.matrix.number_of_genes)
        self.assertEqual(3, dataset.matrix.number_of_samples)

        # Gene conservation law
        self.assertTrue(dataset.is_gene_universe_conserved)

    def test_u3_2_blocks_inferential_processing_when_provenance_is_undocumented(self):
        """If unit or transformation scale is undocumented, inferential processing is blocked."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = routing_metadata(fixture_path, "processed_expression_matrix")
        meta["samples"] = [
            {"sample_id": "GSM7321024", "condition_id": "cond_vehicle"},
            {"sample_id": "GSM7321026", "condition_id": "cond_vpa200"},
            {"sample_id": "GSM7321027", "condition_id": "cond_vpa400"},
        ]
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026", "GSM7321027"]
        # Explicitly omit processing_provenance
        meta.pop("processing_provenance", None)
        meta["input_data_assets"][0].pop("processing_provenance", None)

        dataset = ingest_processed_expression_matrix(fixture_path, meta)

        # Must fail / block with INSUFFICIENT_PROCESSING_PROVENANCE
        self.assertEqual(QCStatus.FAIL, dataset.status)
        self.assertFalse(dataset.is_inferentially_eligible)
        self.assertIsNone(dataset.target_analysis_engine)
        self.assertTrue(any(f.rule_id == "INSUFFICIENT_PROCESSING_PROVENANCE" for f in dataset.gate_result.findings))


if __name__ == "__main__":
    unittest.main()
