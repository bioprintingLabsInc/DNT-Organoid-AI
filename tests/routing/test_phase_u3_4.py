"""Acceptance and integration tests for Phase U3.4: Dispatcher Activation & Full Integration."""

from pathlib import Path
import tempfile
import unittest

from tests.metadata.fixtures import reported, with_provenance
from tests.routing.fixtures import routing_metadata
from src.differential_expression.models import DifferentialExpressionDataset
from src.ingestion.models import Route
from src.normalization.provenance_gate import ProcessedExpressionDataset
from src.routing import dispatch_route_execution, plan_route_execution
from src.routing.dispatcher import RouteExecutionPlan, RouteExecutionResult


class PhaseU34DispatcherActivationTests(unittest.TestCase):
    """Verify Phase U3.4 Processed Expression Pipeline Dispatcher Activation requirements."""

    def test_plan_route_execution_activates_processed_expression_pipeline(self):
        """1. When enabled, processed expression matrices plan READY_FOR_EXECUTION to processed_expression_pipeline_runner."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        self.assertTrue(fixture_path.exists())
        meta = routing_metadata(fixture_path, "processed_expression_matrix")

        plan = plan_route_execution(meta, enable_u3_processed_execution=True)

        self.assertEqual(Route.PROCESSED_EXPRESSION_PIPELINE, plan.route)
        self.assertEqual("READY_FOR_EXECUTION", plan.status)
        self.assertEqual("processed_expression_pipeline_runner", plan.executor_id)
        self.assertEqual("differential_expression_dataset", plan.target_format)

    def test_dispatch_blocks_under_replicated_data_without_differential_inference(self):
        """2. Single-replicate datasets (n=1 per condition, e.g. GSE232218) halt at EXECUTION_BLOCKED_UNDER_REPLICATED."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = routing_metadata(fixture_path, "processed_expression_matrix")
        meta["samples"] = [
            {
                "sample_id": "GSM7321024",
                "experiment_id": "experiment_1",
                "condition_id": "control",
                "biological_source_id": "line_1",
                "organoid_context_id": "cortical_1",
                "biological_replicate_id": reported("rep1"),
                "developmental_age_or_stage_at_collection": reported("day 46"),
                "collection_time_relative_to_exposure": reported("24 hours"),
                "metadata_validation_status": "not_validated",
            },
            {
                "sample_id": "GSM7321026",
                "experiment_id": "experiment_1",
                "condition_id": "treatment_10um_24h_d45",
                "biological_source_id": "line_1",
                "organoid_context_id": "cortical_1",
                "biological_replicate_id": reported("rep1"),
                "developmental_age_or_stage_at_collection": reported("day 46"),
                "collection_time_relative_to_exposure": reported("24 hours"),
                "metadata_validation_status": "not_validated",
            },
        ]
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026"]
        meta = with_provenance(meta)
        meta["processing_provenance"] = {
            "expression_unit": "TPM",
            "transformation_scale": "linear",
            "normalization_method": "tpm",
            "software_pipeline": "Kallisto v0.48.0",
            "reference_annotation": "Ensembl Release 112",
        }

        plan = plan_route_execution(meta, enable_u3_processed_execution=True)
        self.assertEqual("READY_FOR_EXECUTION", plan.status)

        res = dispatch_route_execution(plan)

        self.assertFalse(res.is_success)
        self.assertEqual("EXECUTION_BLOCKED_UNDER_REPLICATED", res.execution_state)
        self.assertIn("dataset", res.artifacts)
        self.assertNotIn("differential_expression_dataset", res.artifacts)
        self.assertIn("fewer than 2 biological replicates", res.error_message.lower())

    def test_dispatch_blocks_insufficient_provenance(self):
        """3. Processed expression matrices with undocumented unit or scale halt at EXECUTION_BLOCKED_INSUFFICIENT_PROVENANCE."""
        fixture_path = Path("tests/fixtures/gse232218_human_brain_organoid_tpm.tsv")
        meta = routing_metadata(fixture_path, "processed_expression_matrix")
        meta["samples"] = [
            {
                "sample_id": "GSM7321024",
                "experiment_id": "experiment_1",
                "condition_id": "control",
                "biological_source_id": "line_1",
                "organoid_context_id": "cortical_1",
                "biological_replicate_id": reported("rep1"),
                "developmental_age_or_stage_at_collection": reported("day 46"),
                "collection_time_relative_to_exposure": reported("24 hours"),
                "metadata_validation_status": "not_validated",
            },
            {
                "sample_id": "GSM7321026",
                "experiment_id": "experiment_1",
                "condition_id": "treatment_10um_24h_d45",
                "biological_source_id": "line_1",
                "organoid_context_id": "cortical_1",
                "biological_replicate_id": reported("rep1"),
                "developmental_age_or_stage_at_collection": reported("day 46"),
                "collection_time_relative_to_exposure": reported("24 hours"),
                "metadata_validation_status": "not_validated",
            },
        ]
        meta["sequencing_assays"][0]["sample_ids"] = ["GSM7321024", "GSM7321026"]
        meta = with_provenance(meta)
        meta.pop("processing_provenance", None)
        meta["input_data_assets"][0].pop("processing_provenance", None)

        plan = plan_route_execution(meta, enable_u3_processed_execution=True)
        self.assertEqual("READY_FOR_EXECUTION", plan.status)

        res = dispatch_route_execution(plan)

        self.assertFalse(res.is_success)
        self.assertEqual("EXECUTION_BLOCKED_INSUFFICIENT_PROVENANCE", res.execution_state)
        self.assertIn("dataset", res.artifacts)

    def test_dispatch_full_success_with_replicated_continuous_dataset(self):
        """4. Replicated (>= 2 reps/arm) continuous processed dataset dispatches full limma-trend pipeline."""
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write("gene\tsample_c1\tsample_c2\tsample_t1\tsample_t2\n")
            f.write("ENSG00000141510\t2.0\t2.2\t32.0\t31.8\n")
            f.write("ENSG00000012048\t16.0\t15.8\t1.0\t1.2\n")
            tmp_tsv = Path(f.name)

        try:
            meta = routing_metadata(tmp_tsv, "processed_expression_matrix")
            meta["samples"] = [s for s in meta["samples"] if s["sample_id"] in {"sample_c1", "sample_c2", "sample_t1", "sample_t2"}]
            meta["sequencing_assays"][0]["sample_ids"] = ["sample_c1", "sample_c2", "sample_t1", "sample_t2"]
            meta = with_provenance(meta)
            meta["processing_provenance"] = {
                "expression_unit": "TPM",
                "transformation_scale": "linear",
                "normalization_method": "tpm",
                "software_pipeline": "RSEM v1.3.3",
                "reference_annotation": "Ensembl Release 112",
            }

            plan = plan_route_execution(meta, enable_u3_processed_execution=True)
            self.assertEqual("READY_FOR_EXECUTION", plan.status)

            r_bin = r"C:\Users\chari\.r\R-4.4.3\bin\x64\Rscript.exe"
            res = dispatch_route_execution(plan, r_binary_path=r_bin)

            self.assertTrue(res.is_success)
            self.assertEqual("COMPLETED", res.execution_state)
            self.assertEqual("differential_expression_dataset", res.emitted_format)
            self.assertIn("dataset", res.artifacts)
            self.assertIn("differential_expression_dataset", res.artifacts)

            de_ds = res.artifacts["differential_expression_dataset"]
            self.assertIsInstance(de_ds, DifferentialExpressionDataset)
            self.assertEqual(1, len(de_ds.eligible_results))

            contrast_res = de_ds.eligible_results[0]
            self.assertEqual(2, len(contrast_res.gene_results))
            tp53 = next(g for g in contrast_res.gene_results if g.canonical_gene_id == "ENSG00000141510")
            self.assertGreater(tp53.log2_fold_change, 3.0)
            self.assertEqual("OK", tp53.result_status)
            self.assertIn("moderated_t", tp53.result_notes)
        finally:
            tmp_tsv.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
