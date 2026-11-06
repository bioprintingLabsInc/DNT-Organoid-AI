"""Unit and integration tests for Phase U3.3 limma-trend differential inference engine."""

from pathlib import Path
import tempfile
import unittest

from tests.metadata.fixtures import canonical_bulk_study, reported, with_provenance
from src.differential_expression.limma_runner import (
    execute_limma_trend_r,
    prepare_continuous_expression_matrix,
    run_limma_trend_contrast,
    run_limma_trend_pipeline,
)
from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
    GeneDifferentialExpressionResult,
)
from src.gene_harmonization.models import GeneIdentifierType, GeneMappingStatus, HarmonizedGene
from src.ingestion.processed_matrix_reader import ProcessedExpressionMatrix
from src.normalization.provenance_gate import (
    ExpressionUnit,
    ProcessedExpressionDataset,
    ProcessedNormalizationMethod,
    ProcessingProvenance,
    ProvenanceGateResult,
    TransformationScale,
    ingest_processed_expression_matrix,
)
from src.qc.processed_matrix_qc import ProcessedMatrixQCResult
from src.qc.result import Severity, Status
from src.response_builder.models import ContrastStatus


def _create_mock_dataset(
    gene_ids: tuple[str, ...],
    sample_ids: tuple[str, ...],
    columns: tuple[tuple[float, ...], ...],
    provenance: ProcessingProvenance,
    harmonized_genes: tuple[HarmonizedGene, ...] | None = None,
    dataset_status: Status = Status.PASS,
) -> ProcessedExpressionDataset:
    matrix = ProcessedExpressionMatrix(gene_ids=gene_ids, sample_ids=sample_ids, columns=columns)
    
    if harmonized_genes is None:
        h_list = []
        for idx, g_id in enumerate(gene_ids, start=1):
            h_list.append(
                HarmonizedGene(
                    source_index=idx,
                    original_gene_id=g_id,
                    identifier_type="ensembl_gene_id",
                    canonical_gene_id=g_id,
                    approved_symbol=f"SYM_{g_id}",
                    mapping_status=GeneMappingStatus.UNIQUELY_MAPPED,
                    mapping_reason=None,
                    candidate_canonical_ids=(),
                    reference_id="ensembl_human_genes_v112",
                    reference_version="112",
                    source_asset_id="asset_mock",
                )
            )
        harmonized_genes = tuple(h_list)

    gate_result = ProvenanceGateResult(
        status=dataset_status,
        findings=(),
        provenance=provenance,
        is_valid_intake=True,
        is_provenance_valid=dataset_status != Status.FAIL,
        is_u3_3_evaluable=dataset_status != Status.FAIL,
        is_statistically_eligible=True,
        statistical_eligibility_reason=None,
        is_inferentially_eligible=dataset_status != Status.FAIL,
        deseq2_eligible=False,
        deseq2_ineligibility_reason="Continuous data requires limma-trend",
        limma_trend_eligible=True,
        target_analysis_engine="limma_trend",
    )

    qc_res = ProcessedMatrixQCResult(
        status=dataset_status,
        findings=(),
        number_of_genes=len(gene_ids),
        number_of_samples=len(sample_ids),
        source_asset_id="asset_mock",
        assay_id="assay_mock",
        study_ids=("study_mock",),
        experiment_ids=("exp_mock",),
        matrix=matrix,
        samples=(),
        conditions=(),
        ruleset_version="processed_matrix_qc_v1",
        calculation_version="processed_matrix_qc_engine_v1",
    )

    return ProcessedExpressionDataset(
        dataset_id="dataset_mock_u3_3",
        matrix=matrix,
        provenance=provenance,
        qc_result=qc_res,
        gate_result=gate_result,
        harmonized_genes=harmonized_genes,
        eligible_canonical_genes_count=len(harmonized_genes),
        excluded_genes_count=0,
        total_input_genes_count=len(gene_ids),
        status=dataset_status,
        is_valid_intake=True,
        is_provenance_valid=dataset_status != Status.FAIL,
        is_u3_3_evaluable=dataset_status != Status.FAIL,
        is_statistically_eligible=True,
        statistical_eligibility_reason=None,
        is_inferentially_eligible=dataset_status != Status.FAIL,
        target_analysis_engine="limma_trend",
        blocked_analysis_engines=("DESeq2",),
    )


class LimmaTrendRunnerTests(unittest.TestCase):
    """Test Phase U3.3 limma-trend differential inference engine."""

    def test_continuous_scale_preparation_linear_to_log2(self):
        """1. Linear scale continuous values are transformed to log2(x + 1) with documented offset."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1", "G2"),
            sample_ids=("s1", "s2"),
            columns=(
                (0.0, 3.0),   # log2(0+1) = 0.0, log2(3+1) = 2.0
                (7.0, 15.0),  # log2(7+1) = 3.0, log2(15+1) = 4.0
            ),
        )
        prov = ProcessingProvenance(
            expression_unit=ExpressionUnit.TPM,
            transformation_scale=TransformationScale.LINEAR,
            normalization_method=ProcessedNormalizationMethod.TPM,
        )

        trans_cols, updated_prov, findings = prepare_continuous_expression_matrix(matrix, prov)
        self.assertEqual([], findings)
        self.assertEqual(TransformationScale.LOG2, updated_prov.transformation_scale)
        self.assertTrue(updated_prov.is_log_transformed)
        self.assertEqual(1.0, updated_prov.pseudocount)

        # Check values
        self.assertAlmostEqual(0.0, trans_cols[0][0])
        self.assertAlmostEqual(2.0, trans_cols[0][1])
        self.assertAlmostEqual(3.0, trans_cols[1][0])
        self.assertAlmostEqual(4.0, trans_cols[1][1])

    def test_continuous_scale_preparation_preserves_already_log2_values(self):
        """2. Already log2-transformed values are preserved without double logging."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("s1",),
            columns=((5.5,),),
        )
        prov = ProcessingProvenance(
            expression_unit=ExpressionUnit.LOG_CPM,
            transformation_scale=TransformationScale.LOG2,
            normalization_method=ProcessedNormalizationMethod.CPM,
            is_log_transformed=True,
        )

        trans_cols, updated_prov, findings = prepare_continuous_expression_matrix(matrix, prov)
        self.assertEqual([], findings)
        self.assertEqual(5.5, trans_cols[0][0])

    def test_continuous_scale_preparation_rejects_negative_linear_values(self):
        """3. Negative values in linear scale emit NEGATIVE_VALUES_IN_LINEAR_SCALE ERROR."""
        matrix = ProcessedExpressionMatrix(
            gene_ids=("G1",),
            sample_ids=("s1",),
            columns=((-1.5,),),
        )
        prov = ProcessingProvenance(
            expression_unit=ExpressionUnit.TPM,
            transformation_scale=TransformationScale.LINEAR,
            normalization_method=ProcessedNormalizationMethod.TPM,
        )

        trans_cols, updated_prov, findings = prepare_continuous_expression_matrix(matrix, prov)
        self.assertTrue(any(f.rule_id == "NEGATIVE_VALUES_IN_LINEAR_SCALE" and f.severity == Severity.ERROR for f in findings))

    def test_insufficient_biological_replicates_blocked_from_inference(self):
        """4. Contrasts with fewer than 2 biological replicates per arm are BLOCKED."""
        prov = ProcessingProvenance(
            expression_unit=ExpressionUnit.TPM,
            transformation_scale=TransformationScale.LINEAR,
            normalization_method=ProcessedNormalizationMethod.TPM,
        )
        # 1 treatment sample, 2 control samples
        dataset = _create_mock_dataset(
            gene_ids=("G1", "G2"),
            sample_ids=("s_trt1", "s_ctrl1", "s_ctrl2"),
            columns=((10.0, 5.0), (2.0, 1.0), (3.0, 2.0)),
            provenance=prov,
        )
        meta = {
            "study_id": "test_under_replicated",
            "samples": [
                {"sample_id": "s_trt1", "condition_id": "cond_trt", "biological_replicate_id": reported("rep1")},
                {"sample_id": "s_ctrl1", "condition_id": "cond_ctrl", "biological_replicate_id": reported("rep1")},
                {"sample_id": "s_ctrl2", "condition_id": "cond_ctrl", "biological_replicate_id": reported("rep2")},
            ],
            "treatment_control_relationships": [
                {
                    "relationship_id": "rel_1",
                    "relationship_type": "treatment_control",
                    "treatment_condition_id": "cond_trt",
                    "control_condition_ids": ["cond_ctrl"],
                }
            ],
        }

        res = run_limma_trend_contrast(dataset, meta)
        self.assertEqual(ContrastStatus.BLOCKED, res.status)
        self.assertFalse(res.is_inferentially_eligible)
        self.assertTrue(any(f.rule_id == "insufficient_biological_replicates" for f in res.findings))

    def test_full_limma_trend_differential_inference_execution(self):
        """5. Replicated 3-vs-3 continuous dataset executes official limma-trend inference with real R."""
        prov = ProcessingProvenance(
            expression_unit=ExpressionUnit.TPM,
            transformation_scale=TransformationScale.LINEAR,
            normalization_method=ProcessedNormalizationMethod.TPM,
            software_pipeline="Salmon v1.10.0",
        )
        # 3 control, 3 treated samples
        # Gene 1 is strongly upregulated in treatment: ctrl around 2.0 (linear), trt around 64.0 (linear)
        # Gene 2 is strongly downregulated: ctrl around 32.0, trt around 1.0
        # Gene 3 is unchanged: ctrl around 10.0, trt around 10.0
        samples = ("c1", "c2", "c3", "t1", "t2", "t3")
        gene_ids = ("ENSG00000141510", "ENSG00000012048", "ENSG00000139618")
        columns = (
            (2.1, 31.5, 9.8),   # c1
            (1.9, 32.8, 10.2),  # c2
            (2.0, 33.1, 10.1),  # c3
            (63.8, 1.1, 9.9),   # t1
            (64.2, 0.9, 10.3),  # t2
            (64.0, 1.0, 10.0),  # t3
        )
        dataset = _create_mock_dataset(gene_ids, samples, columns, prov)

        meta = {
            "study_id": "test_limma_3v3",
            "samples": [
                {"sample_id": "c1", "condition_id": "cond_ctrl", "biological_replicate_id": reported("rep1")},
                {"sample_id": "c2", "condition_id": "cond_ctrl", "biological_replicate_id": reported("rep2")},
                {"sample_id": "c3", "condition_id": "cond_ctrl", "biological_replicate_id": reported("rep3")},
                {"sample_id": "t1", "condition_id": "cond_trt", "biological_replicate_id": reported("rep1")},
                {"sample_id": "t2", "condition_id": "cond_trt", "biological_replicate_id": reported("rep2")},
                {"sample_id": "t3", "condition_id": "cond_trt", "biological_replicate_id": reported("rep3")},
            ],
            "treatment_control_relationships": [
                {
                    "relationship_id": "contrast_trt_vs_ctrl",
                    "relationship_type": "treatment_control",
                    "treatment_condition_id": "cond_trt",
                    "control_condition_ids": ["cond_ctrl"],
                }
            ],
            "agent_name": "Valproic Acid",
            "concentration_or_dose": 200.0,
            "concentration_or_dose_unit": "uM",
            "exposure_duration": 48.0,
            "exposure_duration_unit": "hours",
        }

        res = run_limma_trend_contrast(dataset, meta)

        self.assertEqual(ContrastStatus.ELIGIBLE, res.status)
        self.assertTrue(res.is_inferentially_eligible)
        self.assertEqual(3, len(res.gene_results))
        self.assertEqual("limma_trend_v1", res.de_version)
        self.assertEqual("~ 0 + condition", res.design_formula)
        self.assertIn("limma_trend", res.design_rationale)

        # Gene 1 upregulation: log2(65) - log2(3) ~ 4.43
        g1 = next(g for g in res.gene_results if g.canonical_gene_id == "ENSG00000141510")
        self.assertIsNotNone(g1.log2_fold_change)
        self.assertGreater(g1.log2_fold_change, 3.5)
        self.assertIsNotNone(g1.wald_statistic)  # Moderated t-statistic
        self.assertGreater(g1.wald_statistic, 0.0)
        self.assertIsNotNone(g1.p_value)
        self.assertLess(g1.p_value, 0.01)
        self.assertIsNotNone(g1.adjusted_p_value_bh)
        # SE uninvented
        self.assertIsNone(g1.lfc_standard_error)
        self.assertIn("moderated_t", g1.result_notes)

        # Gene 2 downregulation: log2(2) - log2(33) < -3.5
        g2 = next(g for g in res.gene_results if g.canonical_gene_id == "ENSG00000012048")
        self.assertLess(g2.log2_fold_change, -3.5)
        self.assertLess(g2.wald_statistic, 0.0)

        # Law of Conservation of the Gene Universe
        self.assertTrue(res.is_gene_universe_conserved)
        self.assertEqual(3, res.total_input_genes_count)
        self.assertEqual(3, res.eligible_canonical_genes_count)
        self.assertEqual(0, res.excluded_genes_count)

    def test_limma_trend_pipeline_multi_contrast_aggregation(self):
        """6. Full pipeline wrapper aggregates contrast results into DifferentialExpressionDataset."""
        prov = ProcessingProvenance(
            expression_unit=ExpressionUnit.TPM,
            transformation_scale=TransformationScale.LOG2,
            normalization_method=ProcessedNormalizationMethod.TPM,
            is_log_transformed=True,
        )
        samples = ("c1", "c2", "t1", "t2")
        gene_ids = ("G1", "G2")
        columns = ((1.0, 2.0), (1.1, 2.1), (3.0, 0.5), (3.2, 0.6))
        dataset = _create_mock_dataset(gene_ids, samples, columns, prov)

        meta = {
            "study_id": "test_pipeline",
            "samples": [
                {"sample_id": "c1", "condition_id": "cond_ctrl", "biological_replicate_id": reported("rep1")},
                {"sample_id": "c2", "condition_id": "cond_ctrl", "biological_replicate_id": reported("rep2")},
                {"sample_id": "t1", "condition_id": "cond_trt", "biological_replicate_id": reported("rep1")},
                {"sample_id": "t2", "condition_id": "cond_trt", "biological_replicate_id": reported("rep2")},
            ],
            "treatment_control_relationships": [
                {
                    "relationship_id": "contrast_1",
                    "relationship_type": "treatment_control",
                    "treatment_condition_id": "cond_trt",
                    "control_condition_ids": ["cond_ctrl"],
                }
            ],
        }

        pipeline_res = run_limma_trend_pipeline(dataset, meta)
        self.assertIsInstance(pipeline_res, DifferentialExpressionDataset)
        self.assertEqual(1, len(pipeline_res.contrast_results))
        self.assertEqual(1, len(pipeline_res.eligible_results))
        self.assertEqual(ContrastStatus.ELIGIBLE, pipeline_res.contrast_results[0].status)


if __name__ == "__main__":
    unittest.main()
