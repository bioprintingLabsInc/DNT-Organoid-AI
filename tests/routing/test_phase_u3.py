"""Integration and acceptance tests for Universal RNA-seq Ingestion v2 — Phase U3.1."""

from pathlib import Path
import tempfile
import unittest

from src.ingestion import inspect_asset
from src.ingestion.models import (
    Confidence,
    RegisteredAsset,
    Route,
)
from src.qc.validator import validate_bulk_counts
from src.routing import (
    RouteExecutionPlan,
    RouteExecutionResult,
    dispatch_route_execution,
    plan_route_execution,
)
from src.routing.detector import detect_format
from src.routing.router import _select_route
from tests.routing.fixtures import routing_metadata, write_asset


class UniversalIngestionPhaseU3Tests(unittest.TestCase):
    """Verify Phase U3.1 format detection, routing, execution planning, and dispatching."""

    def test_xlsx_deg_table_detection(self):
        """Excel workbook containing a DEG sheet is detected as differential_expression_table."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.xlsx")
        insp = inspect_asset(fixture_path)
        asset = RegisteredAsset(
            "asset_vpa_xlsx",
            str(fixture_path),
            "differential_expression_table",
            ("assay_vpa",),
            ("bulk_rna_seq",),
        )
        det = detect_format(asset, insp)

        self.assertEqual("differential_expression_table", det.detected_format)
        self.assertEqual(Confidence.HIGH, det.confidence)

        route, findings = _select_route(asset, det)
        self.assertEqual(Route.DIFFERENTIAL_EXPRESSION_PIPELINE, route)
        self.assertEqual((), findings)

    def test_tsv_deg_table_detection(self):
        """TSV file containing real published DEG data is detected as differential_expression_table."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        insp = inspect_asset(fixture_path)
        asset = RegisteredAsset(
            "asset_vpa_tsv",
            str(fixture_path),
            "differential_expression_table",
            ("assay_vpa",),
            ("bulk_rna_seq",),
        )
        det = detect_format(asset, insp)

        self.assertEqual("differential_expression_table", det.detected_format)
        self.assertEqual(Confidence.HIGH, det.confidence)

        route, findings = _select_route(asset, det)
        self.assertEqual(Route.DIFFERENTIAL_EXPRESSION_PIPELINE, route)
        self.assertEqual((), findings)

    def test_deg_table_execution_planning(self):
        """When U3 execution is enabled, plan_route_execution assigns deg_table_pipeline_runner."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        meta = routing_metadata(fixture_path, "differential_expression_table")

        # When enable_u3_deg_execution=True
        plan = plan_route_execution(meta, enable_u3_deg_execution=True)

        self.assertEqual(Route.DIFFERENTIAL_EXPRESSION_PIPELINE, plan.route)
        self.assertEqual("READY_FOR_EXECUTION", plan.status)
        self.assertEqual("deg_table_pipeline_runner", plan.executor_id)
        self.assertEqual("differential_expression_contrast_result", plan.target_format)
        self.assertIsNone(plan.downstream_route)

        prov = plan.provenance_chain
        self.assertEqual("differential_expression_table", prov["input_format"])
        self.assertEqual("DIFFERENTIAL_EXPRESSION_PIPELINE", prov["selected_route"])
        self.assertEqual("deg_table_pipeline_runner", prov["executor"])

    def test_u1_backward_compatibility_preserved(self):
        """Without U3 execution flag, plan_route_execution preserves U1 EXECUTION_NOT_IMPLEMENTED status."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        meta = routing_metadata(fixture_path, "differential_expression_table")

        plan = plan_route_execution(meta)  # Default: U1 behavior
        self.assertEqual(Route.DIFFERENTIAL_EXPRESSION_PIPELINE, plan.route)
        self.assertEqual("EXECUTION_NOT_IMPLEMENTED", plan.status)
        self.assertIsNone(plan.executor_id)

    def test_deg_table_dispatch_execution_end_to_end(self):
        """End-to-end dispatch of a published DEG table runs adapter and emits contrast result artifact."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        meta = routing_metadata(fixture_path, "differential_expression_table")

        plan = plan_route_execution(meta, enable_u3_deg_execution=True)
        res = dispatch_route_execution(plan)

        self.assertTrue(res.is_success)
        self.assertEqual("COMPLETED", res.execution_state)
        self.assertEqual("differential_expression_contrast_result", res.emitted_format)

        contrast = res.artifacts.get("contrast_result")
        self.assertIsNotNone(contrast)
        self.assertEqual(500, contrast.total_input_genes_count)
        self.assertTrue(contrast.eligible_canonical_genes_count >= 490)
        self.assertTrue(contrast.is_gene_universe_conserved)

    def test_deg_table_custom_executor_dispatch(self):
        """dispatch_route_execution respects a caller-supplied custom executor callable."""
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        meta = routing_metadata(fixture_path, "differential_expression_table")
        plan = plan_route_execution(meta, enable_u3_deg_execution=True)

        executed_assets = []

        def mock_executor(assets, plan, **kwargs):
            executed_assets.extend(assets)
            return RouteExecutionResult(
                plan=plan,
                is_success=True,
                execution_state="CUSTOM_EXECUTOR_COMPLETED",
                emitted_asset_path=assets[0].local_path,
                emitted_format="custom_result",
                provenance_record={"custom": True},
            )

        res = dispatch_route_execution(plan, executor=mock_executor)
        self.assertTrue(res.is_success)
        self.assertEqual("CUSTOM_EXECUTOR_COMPLETED", res.execution_state)
        self.assertEqual(1, len(executed_assets))

    def test_prevention_deg_table_cannot_enter_step6a_count_validator(self):
        """CRITICAL: Published DEG tables are non-integer summaries and cannot pass Step 6A count QC."""
        from src.qc.result import Status as QCStatus
        fixture_path = Path("tests/fixtures/vpa_organoid_deg_published.tsv")
        meta = routing_metadata(fixture_path, "raw_counts")
        qc_result = validate_bulk_counts(
            fixture_path,
            meta,
            "asset_deg",
            "assay_deg",
        )
        # Must fail bulk count QC with FAIL status because it is not an integer count matrix
        self.assertEqual(QCStatus.FAIL, qc_result.status)
        self.assertTrue(any(f.rule_id == "integer_like_counts_required" for f in qc_result.findings))


if __name__ == "__main__":
    unittest.main()
