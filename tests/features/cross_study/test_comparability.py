"""Unit and integration tests for Stage 8A.2 Cross-Study Comparability Evaluation.

Governed Scientific Tests:
1. Direction concordance logic (CONCORDANT_UP, CONCORDANT_DOWN, DISCORDANT, NO_DIRECTION, UNDETERMINED).
2. Comparability gating across compounds, doses, stages, durations, models, and replicate verification.
3. Cross-study independence: blocks declaration of intra-study conditions/cell lines as independent studies.
4. Explicit replicate verification: missing verification does NOT default to approved.
5. Organoid model compatibility: differing organoid platforms block synthesis.
6. Missing metadata rejection: missing/unknown metadata is NOT interpreted as verified equivalence.
7. Strict enforcement of synthesis boundary: zero pooled estimates when conditions differ.
8. Comparison QC auditor integrity.
9. Integration verification with real GSE232218 vs GSE187006 artifacts.
"""

from pathlib import Path
import pytest
import pandas as pd

from src.features.cross_study.comparability import (
    compare_cross_study_gene_responses,
    derive_direction_concordance,
    evaluate_study_comparability,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyComparisonDataset,
    CrossStudyStatus,
    DirectionConcordance,
    EstimateStatus,
    GeneCrossStudyComparisonRecord,
    ResponseDirection,
    StudyComparabilityAssessment,
)
from src.features.cross_study.qc import (
    ComparisonQCSummary,
    audit_comparison_qc,
)
from src.features.cross_study.serialization import (
    CANONICAL_COMPARISON_COLUMNS,
    export_comparison_to_csv,
    export_comparison_to_dataframe,
)


def _make_dummy_record(
    gene_id: str = "ENSG00000141510",
    study_id: str = "study_1",
    study_family_id: str | None = None,
    independent_study_id: str | None = None,
    contrast_id: str = "contrast_1",
    compound_id: str = "CID:3121",
    compound_name: str = "Valproic Acid",
    dose_value: float | None = 200.0,
    dose_unit: str | None = "uM",
    developmental_stage: str | None = "Day 6",
    exposure_duration: float | None = 48.0,
    exposure_duration_unit: str | None = "h",
    model_type: str | None = "Organoid A",
    treatment_replicate_count: int = 3,
    control_replicate_count: int = 3,
    log2_fc: float = 1.5,
    p_val: float = 0.01,
    direction: ResponseDirection = ResponseDirection.UP,
    estimate_status: EstimateStatus = EstimateStatus.OK,
) -> ConditionGeneResponseRecord:
    fam_id = study_family_id or study_id
    ind_id = independent_study_id or study_id
    return ConditionGeneResponseRecord(
        canonical_ensembl_gene_id=gene_id,
        gene_symbol="DUMMY",
        study_id=study_id,
        study_family_id=fam_id,
        independent_study_id=ind_id,
        experiment_id="exp_1",
        contrast_id=contrast_id,
        treatment_condition_id="trt",
        matched_control_condition_ids=("ctrl",),
        compound_id=compound_id,
        compound_name=compound_name,
        dose_value=dose_value,
        dose_unit=dose_unit,
        dose_original=f"{dose_value} {dose_unit}" if dose_value is not None else None,
        exposure_duration=exposure_duration,
        exposure_duration_unit=exposure_duration_unit,
        developmental_stage=developmental_stage,
        model_type=model_type,
        biological_source_ids=("donor_1",),
        organoid_context_ids=("ctx_1",),
        treatment_replicate_count=treatment_replicate_count,
        control_replicate_count=control_replicate_count,
        treatment_sample_ids=("s1", "s2", "s3"),
        control_sample_ids=("s4", "s5", "s6"),
        treatment_biological_replicate_ids=("r1", "r2", "r3"),
        control_biological_replicate_ids=("r4", "r5", "r6"),
        base_mean=100.0,
        log2_fold_change=log2_fc,
        log2_fold_change_se=0.2,
        wald_statistic=log2_fc / 0.2,
        p_value=p_val,
        adjusted_p_value=p_val,
        response_direction=direction,
        estimate_status=estimate_status,
        missing_reason=None,
        cross_study_status=CrossStudyStatus.SINGLE_STUDY_ONLY,
    )


def test_derive_direction_concordance() -> None:
    """Test direction concordance mapping."""
    assert derive_direction_concordance(ResponseDirection.UP, ResponseDirection.UP) == DirectionConcordance.CONCORDANT_UP
    assert derive_direction_concordance(ResponseDirection.DOWN, ResponseDirection.DOWN) == DirectionConcordance.CONCORDANT_DOWN
    assert derive_direction_concordance(ResponseDirection.UP, ResponseDirection.DOWN) == DirectionConcordance.DISCORDANT
    assert derive_direction_concordance(ResponseDirection.DOWN, ResponseDirection.UP) == DirectionConcordance.DISCORDANT
    assert derive_direction_concordance(ResponseDirection.NO_DIRECTION, ResponseDirection.UP) == DirectionConcordance.NO_DIRECTION
    assert derive_direction_concordance(ResponseDirection.UP, ResponseDirection.UNDETERMINED) == DirectionConcordance.UNDETERMINED


def test_study_comparability_detection_of_mismatches() -> None:
    """Test comparability gating catches dose, stage, duration, compound, and replicate mismatches."""
    r1 = _make_dummy_record(study_id="study_1")

    # 1. Dose mismatch
    r2_dose = _make_dummy_record(study_id="study_2", dose_value=1000.0, dose_unit="uM")
    eval_dose = evaluate_study_comparability(r1, r2_dose, replicates_verified1=True, replicates_verified2=True)
    assert eval_dose.is_synthesis_justified is False
    assert eval_dose.cross_study_status == CrossStudyStatus.NONCOMPARABLE_CONTEXT
    assert any("Dose mismatch" in d for d in eval_dose.differences)

    # 2. Stage mismatch
    r2_stage = _make_dummy_record(study_id="study_2", developmental_stage="Day 42")
    eval_stage = evaluate_study_comparability(r1, r2_stage, replicates_verified1=True, replicates_verified2=True)
    assert eval_stage.is_synthesis_justified is False
    assert any("Developmental stage mismatch" in d for d in eval_stage.differences)

    # 3. Duration mismatch
    r2_dur = _make_dummy_record(study_id="study_2", exposure_duration=72.0)
    eval_dur = evaluate_study_comparability(r1, r2_dur, replicates_verified1=True, replicates_verified2=True)
    assert eval_dur.is_synthesis_justified is False
    assert any("Exposure duration mismatch" in d for d in eval_dur.differences)

    # 4. Compound mismatch
    r2_chem = _make_dummy_record(study_id="study_2", compound_id="CID:9999", compound_name="Other Chemical")
    eval_chem = evaluate_study_comparability(r1, r2_chem, replicates_verified1=True, replicates_verified2=True)
    assert eval_chem.is_synthesis_justified is False
    assert any("Compound mismatch" in d for d in eval_chem.differences)

    # 5. Under-replicated study (< 2 replicates)
    r2_rep = _make_dummy_record(study_id="study_2", treatment_replicate_count=1)
    eval_rep = evaluate_study_comparability(r1, r2_rep, replicates_verified1=True, replicates_verified2=True)
    assert eval_rep.is_synthesis_justified is False
    assert any("Replicate deficiency" in d for d in eval_rep.differences)

    # 6. Identical conditions with explicit replicate verification -> Synthesis justified
    r2_identical = _make_dummy_record(study_id="study_2")
    eval_identical = evaluate_study_comparability(r1, r2_identical, replicates_verified1=True, replicates_verified2=True)
    assert eval_identical.is_synthesis_justified is True
    assert eval_identical.cross_study_status == CrossStudyStatus.COMPARABLE_INDEPENDENT_STUDIES


def test_cross_study_independence_enforcement() -> None:
    """Test explicit requirement for independent studies; intra-study cell lines cannot count as independent."""
    # Two cell lines from same study family (e.g. U1M vs U2F in GSE187006)
    r_u1m = _make_dummy_record(study_id="GSE187006", study_family_id="GSE187006", contrast_id="U1M")
    r_u2f = _make_dummy_record(study_id="GSE187006", study_family_id="GSE187006", contrast_id="U2F")

    eval_same_study = evaluate_study_comparability(r_u1m, r_u2f, replicates_verified1=True, replicates_verified2=True)
    assert eval_same_study.is_independent_studies is False
    assert eval_same_study.is_synthesis_justified is False
    assert eval_same_study.cross_study_status == CrossStudyStatus.SINGLE_STUDY_ONLY
    assert any("same study or study family" in d for d in eval_same_study.differences)

    # Missing study identifier
    r_missing_study = _make_dummy_record(study_id="", independent_study_id="", contrast_id="c2")
    eval_missing = evaluate_study_comparability(r_u1m, r_missing_study, replicates_verified1=True, replicates_verified2=True)
    assert eval_missing.is_independent_studies is False
    assert eval_missing.is_synthesis_justified is False
    assert any("Study identifiers missing or unverified" in d for d in eval_missing.differences)


def test_biological_replicates_explicit_verification_required() -> None:
    """Test replicate verification does NOT default to approved; explicit status is required."""
    r1 = _make_dummy_record(study_id="study_1")
    r2 = _make_dummy_record(study_id="study_2")

    # Omitted verification parameters (default None)
    eval_default = evaluate_study_comparability(r1, r2)
    assert eval_default.are_replicates_verified is False
    assert eval_default.study1_replicates_verified is False
    assert eval_default.study2_replicates_verified is False
    assert eval_default.is_synthesis_justified is False
    assert any("replicates not explicitly verified" in d for d in eval_default.differences)

    # Only one study verified
    eval_partial = evaluate_study_comparability(r1, r2, replicates_verified1=True, replicates_verified2=False)
    assert eval_partial.are_replicates_verified is False
    assert eval_partial.is_synthesis_justified is False


def test_organoid_model_compatibility_gating() -> None:
    """Test organoid model platform compatibility is strictly enforced."""
    r1 = _make_dummy_record(study_id="study_1", model_type="Cerebral Cortex Organoid (SOSR-CO)")
    r2 = _make_dummy_record(study_id="study_2", model_type="Human Forebrain Organoid")

    eval_model = evaluate_study_comparability(r1, r2, replicates_verified1=True, replicates_verified2=True)
    assert eval_model.is_organoid_model_matched is False
    assert eval_model.is_synthesis_justified is False
    assert eval_model.cross_study_status == CrossStudyStatus.NONCOMPARABLE_CONTEXT
    assert any("Organoid model mismatch" in d for d in eval_model.differences)


def test_missing_metadata_blocks_synthesis() -> None:
    """Test missing or placeholder metadata is strictly not interpreted as verified equivalence."""
    r1 = _make_dummy_record(study_id="study_1")

    # Missing dose
    r2_no_dose = _make_dummy_record(study_id="study_2", dose_value=None, dose_unit=None)
    eval_no_dose = evaluate_study_comparability(r1, r2_no_dose, replicates_verified1=True, replicates_verified2=True)
    assert eval_no_dose.is_dose_matched is False
    assert eval_no_dose.is_synthesis_justified is False

    # Placeholder / unknown stage
    r2_unknown_stage = _make_dummy_record(study_id="study_2", developmental_stage="UNKNOWN")
    eval_stage = evaluate_study_comparability(r1, r2_unknown_stage, replicates_verified1=True, replicates_verified2=True)
    assert eval_stage.is_developmental_stage_matched is False
    assert eval_stage.is_synthesis_justified is False

    # Missing model
    r2_no_model = _make_dummy_record(study_id="study_2", model_type="unspecified")
    eval_model = evaluate_study_comparability(r1, r2_no_model, replicates_verified1=True, replicates_verified2=True)
    assert eval_model.is_organoid_model_matched is False
    assert eval_model.is_synthesis_justified is False


def test_cross_study_comparison_rejection_of_pooling() -> None:
    """Test compare_cross_study_gene_responses blocks pooling when synthesis is not justified."""
    r1 = _make_dummy_record(gene_id="ENSG00000001", dose_value=200.0, developmental_stage="Day 6")
    r2 = _make_dummy_record(gene_id="ENSG00000001", study_id="study_2", dose_value=1000.0, developmental_stage="Day 42")

    reg1 = ConditionGeneResponseRegistry(records=(r1,), registry_id="reg1")
    reg2 = ConditionGeneResponseRegistry(records=(r2,), registry_id="reg2")

    cmp = compare_cross_study_gene_responses(reg1, reg2, replicates_verified1=True, replicates_verified2=True)
    assert cmp.assessment.is_synthesis_justified is False
    assert cmp.assessment.cross_study_status == CrossStudyStatus.NONCOMPARABLE_CONTEXT

    rec = cmp.gene_comparisons[0]
    assert rec.pooled_log2_fold_change is None
    assert rec.pooled_standard_error is None
    assert rec.direction_concordance == DirectionConcordance.CONCORDANT_UP
    assert rec.is_concordant is True


def test_comparison_qc_auditor() -> None:
    """Test audit_comparison_qc validates clean comparison and flags boundary violations."""
    r1 = _make_dummy_record(gene_id="ENSG00000001", dose_value=200.0)
    r2 = _make_dummy_record(gene_id="ENSG00000001", study_id="study_2", dose_value=1000.0)

    reg1 = ConditionGeneResponseRegistry(records=(r1,), registry_id="reg1")
    reg2 = ConditionGeneResponseRegistry(records=(r2,), registry_id="reg2")

    cmp = compare_cross_study_gene_responses(reg1, reg2, replicates_verified1=True, replicates_verified2=True)
    qc = audit_comparison_qc(cmp, reg1, reg2)
    assert qc.is_valid is True
    assert qc.numerical_mismatches_count == 0
    assert qc.synthesis_boundary_violations_count == 0

    # Inject artificial boundary violation: manufactured pooled estimate when not justified
    corrupted_rec = GeneCrossStudyComparisonRecord(
        canonical_ensembl_gene_id=r1.canonical_ensembl_gene_id,
        gene_symbol=r1.gene_symbol,
        study1_id=r1.study_id,
        study1_contrast_id=r1.contrast_id,
        study1_log2_fold_change=r1.log2_fold_change,
        study1_log2_fold_change_se=r1.log2_fold_change_se,
        study1_wald_statistic=r1.wald_statistic,
        study1_p_value=r1.p_value,
        study1_adjusted_p_value=r1.adjusted_p_value,
        study1_response_direction=r1.response_direction,
        study1_estimate_status=r1.estimate_status,
        study2_id=r2.study_id,
        study2_contrast_id=r2.contrast_id,
        study2_log2_fold_change=r2.log2_fold_change,
        study2_log2_fold_change_se=r2.log2_fold_change_se,
        study2_wald_statistic=r2.wald_statistic,
        study2_p_value=r2.p_value,
        study2_adjusted_p_value=r2.adjusted_p_value,
        study2_response_direction=r2.response_direction,
        study2_estimate_status=r2.estimate_status,
        direction_concordance=DirectionConcordance.CONCORDANT_UP,
        is_concordant=True,
        cross_study_status=CrossStudyStatus.NONCOMPARABLE_CONTEXT,
        is_synthesis_justified=False,
        pooled_log2_fold_change=1.234,  # VIOLATION
        pooled_standard_error=0.123,     # VIOLATION
    )
    corrupted_cmp = CrossStudyComparisonDataset(
        comparison_id=cmp.comparison_id,
        assessment=cmp.assessment,
        gene_comparisons=(corrupted_rec,),
    )
    corrupted_qc = audit_comparison_qc(corrupted_cmp, reg1, reg2)
    assert corrupted_qc.is_valid is False
    assert corrupted_qc.synthesis_boundary_violations_count > 0


def test_real_gse232218_vs_gse187006_comparison_artifacts() -> None:
    """Integration test verifying real comparison artifacts between GSE232218 and GSE187006."""
    u1m_csv = Path("data/processed/GSE232218_vs_GSE187006_U1M_cross_study_comparison.csv")
    u2f_csv = Path("data/processed/GSE232218_vs_GSE187006_U2F_cross_study_comparison.csv")

    if not u1m_csv.exists() or not u2f_csv.exists():
        pytest.skip("Comparison CSV artifacts not found. Run scripts/run_gse187006_and_comparability.py first.")

    df_u1m = pd.read_csv(u1m_csv)
    df_u2f = pd.read_csv(u2f_csv)

    # 1. Total rows match 60,221 union gene universe
    assert len(df_u1m) == 60221
    assert len(df_u2f) == 60221

    # 2. Check canonical governed comparison columns
    assert list(df_u1m.columns) == CANONICAL_COMPARISON_COLUMNS
    assert list(df_u2f.columns) == CANONICAL_COMPARISON_COLUMNS

    # 3. Synthesis boundary: pooled estimates strictly absent
    assert df_u1m["pooled_log2_fold_change"].isna().all()
    assert df_u1m["pooled_standard_error"].isna().all()
    assert (df_u1m["cross_study_status"] == "NONCOMPARABLE_CONTEXT").all()
    assert (df_u1m["is_synthesis_justified"] == False).all()

    # 4. Direction concordance breakdown
    conc_up = (df_u1m["direction_concordance"] == "CONCORDANT_UP").sum()
    conc_down = (df_u1m["direction_concordance"] == "CONCORDANT_DOWN").sum()
    disc = (df_u1m["direction_concordance"] == "DISCORDANT").sum()

    assert conc_up == 4877
    assert conc_down == 4739
    assert disc == 5736
    concordance_rate = (conc_up + conc_down) / (conc_up + conc_down + disc)
    assert abs(concordance_rate - 0.6264) < 0.001
