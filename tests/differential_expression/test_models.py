"""Unit tests for Step 7B Differential Expression v1 data models."""

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
    ExcludedGeneAudit,
    GeneDifferentialExpressionResult,
)
from src.response_builder.models import (
    ContrastExposure,
    ContrastStatus,
    Finding,
    Severity,
    Status,
    ordered_findings,
    status_for,
)


def test_gene_differential_expression_result_model():
    """Test GeneDifferentialExpressionResult instantiation and serialization."""
    gene_res = GeneDifferentialExpressionResult(
        canonical_gene_id="ENSG00000141510",
        approved_symbol="TP53",
        base_mean=254.5,
        log2_fold_change=1.45,
        lfc_standard_error=0.22,
        wald_statistic=6.59,
        p_value=4.3e-11,
        adjusted_p_value_bh=1.2e-9,
        result_status="OK",
        result_notes=None,
    )
    d = gene_res.to_dict()
    assert d["canonical_gene_id"] == "ENSG00000141510"
    assert d["approved_symbol"] == "TP53"
    assert d["base_mean"] == 254.5
    assert d["log2_fold_change"] == 1.45
    assert d["result_status"] == "OK"
    assert d["result_notes"] is None


def test_gene_differential_expression_result_outlier():
    """Test GeneDifferentialExpressionResult outlier preservation with NA statistics."""
    outlier_gene = GeneDifferentialExpressionResult(
        canonical_gene_id="ENSG00000139618",
        approved_symbol="BRCA2",
        base_mean=14500.0,
        log2_fold_change=8.45,
        lfc_standard_error=2.1,
        wald_statistic=4.02,
        p_value=None,
        adjusted_p_value_bh=None,
        result_status="OUTLIER",
        result_notes="p-value is NA due to Cook's distance outlier detection.",
    )
    d = outlier_gene.to_dict()
    assert d["p_value"] is None
    assert d["adjusted_p_value_bh"] is None
    assert d["result_status"] == "OUTLIER"
    assert "Cook's distance" in d["result_notes"]


def test_differential_expression_contrast_result_model():
    """Test DifferentialExpressionContrastResult fields and properties."""
    exposure = ContrastExposure(
        exposure_id="exp_01",
        agent_name="Valproic acid",
        concentration_or_dose=200.0,
        concentration_or_dose_unit="uM",
    )
    gene_res = GeneDifferentialExpressionResult(
        canonical_gene_id="ENSG00000141510",
        approved_symbol="TP53",
        base_mean=100.0,
        log2_fold_change=1.2,
        lfc_standard_error=0.3,
        wald_statistic=4.0,
        p_value=0.001,
        adjusted_p_value_bh=0.01,
    )
    contrast_res = DifferentialExpressionContrastResult(
        contrast_id="contrast_001",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="treat_01",
        matched_control_condition_ids=("ctrl_01",),
        treatment_sample_ids=("s1", "s2", "s3"),
        control_sample_ids=("s4", "s5", "s6"),
        treatment_biological_replicate_ids=("rep1", "rep2", "rep3"),
        control_biological_replicate_ids=("rep4", "rep5", "rep6"),
        treatment_replicate_count=3,
        control_replicate_count=3,
        treatment_exposures=(exposure,),
        design_formula="~ condition",
        design_rationale="Single biological source.",
        design_matrix_rank=2,
        residual_degrees_of_freedom=4,
        size_factors_used={"s1": 1.0, "s2": 1.0, "s3": 1.0, "s4": 1.0, "s5": 1.0, "s6": 1.0},
        gene_results=(gene_res,),
        status=ContrastStatus.ELIGIBLE,
    )
    assert contrast_res.is_inferentially_eligible is True
    d = contrast_res.to_dict()
    assert d["contrast_id"] == "contrast_001"
    assert d["status"] == "ELIGIBLE"
    assert len(d["gene_results"]) == 1


def test_differential_expression_dataset_model():
    """Test DifferentialExpressionDataset container aggregation."""
    exposure = ContrastExposure(exposure_id="exp_01", agent_name="VPA")
    c1 = DifferentialExpressionContrastResult(
        contrast_id="c1",
        study_id="s1",
        experiment_id="e1",
        cohort_id="co1",
        relationship_id="r1",
        treatment_condition_id="t1",
        matched_control_condition_ids=("c0",),
        treatment_sample_ids=("s1", "s2"),
        control_sample_ids=("s3", "s4"),
        treatment_biological_replicate_ids=("r1", "r2"),
        control_biological_replicate_ids=("r3", "r4"),
        treatment_replicate_count=2,
        control_replicate_count=2,
        status=ContrastStatus.ELIGIBLE,
        gene_results=(
            GeneDifferentialExpressionResult(
                canonical_gene_id="ENSG00000000001",
                approved_symbol="GENE1",
                base_mean=50.0,
                log2_fold_change=0.5,
                lfc_standard_error=0.1,
                wald_statistic=5.0,
                p_value=1e-5,
                adjusted_p_value_bh=1e-4,
            ),
        ),
    )
    dataset = DifferentialExpressionDataset(
        status=Status.PASS,
        findings=(),
        contrast_results=(c1,),
        source_asset_id="asset_01",
        reference_identity="Ensembl_112",
        builder_version="1.0.0",
        de_version="1.0.0",
    )
    assert len(dataset.eligible_results) == 1
    assert len(dataset.blocked_results) == 0
    d = dataset.to_dict()
    assert d["number_of_contrasts"] == 1
    assert d["number_of_eligible"] == 1


def test_excluded_gene_audit_model():
    """Test ExcludedGeneAudit model and gene universe conservation property."""
    audit = ExcludedGeneAudit(
        original_gene_id="OLD_GENE_1",
        source_index=3,
        mapping_status="UNMAPPED",
        canonical_gene_id=None,
        exclusion_reason="Not found in Ensembl reference.",
    )
    d = audit.to_dict()
    assert d["original_gene_id"] == "OLD_GENE_1"
    assert d["source_index"] == 3
    assert d["mapping_status"] == "UNMAPPED"
    assert d["canonical_gene_id"] is None
    assert "Not found" in d["exclusion_reason"]

    # Test conservation logic on contrast result
    gene_res = GeneDifferentialExpressionResult(
        canonical_gene_id="ENSG00000141510",
        approved_symbol="TP53",
        base_mean=100.0,
        log2_fold_change=1.2,
        lfc_standard_error=0.3,
        wald_statistic=4.0,
        p_value=0.001,
        adjusted_p_value_bh=0.01,
    )
    cr = DifferentialExpressionContrastResult(
        contrast_id="c_test",
        study_id="s1",
        experiment_id="e1",
        cohort_id="co1",
        relationship_id="r1",
        treatment_condition_id="t1",
        matched_control_condition_ids=("c0",),
        treatment_sample_ids=("s1", "s2"),
        control_sample_ids=("s3", "s4"),
        treatment_biological_replicate_ids=("r1", "r2"),
        control_biological_replicate_ids=("r3", "r4"),
        treatment_replicate_count=2,
        control_replicate_count=2,
        gene_results=(gene_res,),
        excluded_genes=(audit,),
        total_input_genes_count=2,
    )
    assert cr.eligible_canonical_genes_count == 1
    assert cr.excluded_genes_count == 1
    assert cr.is_gene_universe_conserved is True
