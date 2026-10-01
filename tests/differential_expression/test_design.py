"""Unit tests for Step 7B DesignResolver covering Cases A-E and estimability."""

from src.differential_expression.design import DesignResolver, matrix_rank
from src.response_builder.models import (
    ContrastStatus,
    MolecularResponseContrast,
    Severity,
)


def make_test_contrast(
    treatment_samples: tuple[str, ...],
    control_samples: tuple[str, ...],
    treatment_bio_reps: tuple[str, ...],
    control_bio_reps: tuple[str, ...],
    biological_sources: tuple[str, ...] = ("source_1",),
) -> MolecularResponseContrast:
    return MolecularResponseContrast(
        contrast_id="test_contrast",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="treat_01",
        matched_control_condition_ids=("ctrl_01",),
        treatment_sample_ids=treatment_samples,
        control_sample_ids=control_samples,
        treatment_biological_replicate_ids=treatment_bio_reps,
        control_biological_replicate_ids=control_bio_reps,
        biological_source_ids=biological_sources,
        status=ContrastStatus.ELIGIBLE,
    )


def test_matrix_rank_calculation():
    """Verify pure python Gaussian elimination matrix rank."""
    # Identity 3x3
    id3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    assert matrix_rank(id3) == 3

    # Rank deficient 3x3
    rank2 = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [0.0, 1.0, 1.0]]
    assert matrix_rank(rank2) == 2


def test_case_a_single_source_condition_design():
    """Case A: All samples from one biological source -> ~ condition."""
    contrast = make_test_contrast(
        treatment_samples=("t1", "t2", "t3"),
        control_samples=("c1", "c2", "c3"),
        treatment_bio_reps=("rep_t1", "rep_t2", "rep_t3"),
        control_bio_reps=("rep_c1", "rep_c2", "rep_c3"),
        biological_sources=("cell_line_A",),
    )
    res = DesignResolver.resolve(contrast)
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.formula == "~ condition"
    assert "Single biological source" in res.rationale
    assert res.rank == 2
    assert res.residual_degrees_of_freedom == 4
    assert res.is_estimable is True


def test_case_b_multiple_paired_biological_sources():
    """Case B: Multiple paired biological sources represented in both arms -> ~ biological_source + condition."""
    contrast = make_test_contrast(
        treatment_samples=("t1", "t2"),
        control_samples=("c1", "c2"),
        treatment_bio_reps=("rep_t1", "rep_t2"),
        control_bio_reps=("rep_c1", "rep_c2"),
        biological_sources=("donor_A", "donor_B"),
    )
    meta = {
        "t1": {"sample_id": "t1", "biological_source_id": "donor_A", "biological_replicate_id": "rep_t1"},
        "c1": {"sample_id": "c1", "biological_source_id": "donor_A", "biological_replicate_id": "rep_c1"},
        "t2": {"sample_id": "t2", "biological_source_id": "donor_B", "biological_replicate_id": "rep_t2"},
        "c2": {"sample_id": "c2", "biological_source_id": "donor_B", "biological_replicate_id": "rep_c2"},
    }
    res = DesignResolver.resolve(contrast, meta)
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.formula == "~ biological_source + condition"
    assert "Multiple paired biological sources" in res.rationale
    assert res.rank == 3
    assert res.residual_degrees_of_freedom == 1
    assert res.is_estimable is True


def test_case_c_independent_unpaired_unique_sources():
    """Case C: Independent unpaired unique sources -> ~ condition without source fixed effects."""
    contrast = make_test_contrast(
        treatment_samples=("t1", "t2", "t3"),
        control_samples=("c1", "c2", "c3"),
        treatment_bio_reps=("rep_t1", "rep_t2", "rep_t3"),
        control_bio_reps=("rep_c1", "rep_c2", "rep_c3"),
        biological_sources=("donor_1", "donor_2", "donor_3", "donor_4", "donor_5", "donor_6"),
    )
    meta = {
        "t1": {"sample_id": "t1", "biological_source_id": "donor_1", "biological_replicate_id": "rep_t1"},
        "t2": {"sample_id": "t2", "biological_source_id": "donor_2", "biological_replicate_id": "rep_t2"},
        "t3": {"sample_id": "t3", "biological_source_id": "donor_3", "biological_replicate_id": "rep_t3"},
        "c1": {"sample_id": "c1", "biological_source_id": "donor_4", "biological_replicate_id": "rep_c1"},
        "c2": {"sample_id": "c2", "biological_source_id": "donor_5", "biological_replicate_id": "rep_c2"},
        "c3": {"sample_id": "c3", "biological_source_id": "donor_6", "biological_replicate_id": "rep_c3"},
    }
    res = DesignResolver.resolve(contrast, meta)
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.formula == "~ condition"
    assert "Independent unpaired biological sources" in res.rationale
    assert res.rank == 2
    assert res.residual_degrees_of_freedom == 4
    assert res.is_estimable is True


def test_case_d_treatment_source_confounding_blocked():
    """Case D: Perfect treatment/source confounding -> BLOCKED with treatment_source_confounding."""
    contrast = make_test_contrast(
        treatment_samples=("t1", "t2", "t3"),
        control_samples=("c1", "c2", "c3"),
        treatment_bio_reps=("rep_t1", "rep_t2", "rep_t3"),
        control_bio_reps=("rep_c1", "rep_c2", "rep_c3"),
        biological_sources=("source_A", "source_B"),
    )
    # All treatments from source_A, all controls from source_B
    meta = {
        "t1": {"sample_id": "t1", "biological_source_id": "source_A", "biological_replicate_id": "rep_t1"},
        "t2": {"sample_id": "t2", "biological_source_id": "source_A", "biological_replicate_id": "rep_t2"},
        "t3": {"sample_id": "t3", "biological_source_id": "source_A", "biological_replicate_id": "rep_t3"},
        "c1": {"sample_id": "c1", "biological_source_id": "source_B", "biological_replicate_id": "rep_c1"},
        "c2": {"sample_id": "c2", "biological_source_id": "source_B", "biological_replicate_id": "rep_c2"},
        "c3": {"sample_id": "c3", "biological_source_id": "source_B", "biological_replicate_id": "rep_c3"},
    }
    res = DesignResolver.resolve(contrast, meta)
    assert res.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "treatment_source_confounding" for f in res.findings)


def test_case_e_partially_crossed_sources_needs_review():
    """Case E: Partially crossed biological sources -> NEEDS_REVIEW."""
    contrast = make_test_contrast(
        treatment_samples=("t1", "t2", "t3"),
        control_samples=("c1", "c2"),
        treatment_bio_reps=("rep_t1", "rep_t2", "rep_t3"),
        control_bio_reps=("rep_c1", "rep_c2"),
        biological_sources=("source_A", "source_B"),
    )
    # Source A in treatment and control, but source B only in treatment
    meta = {
        "t1": {"sample_id": "t1", "biological_source_id": "source_A", "biological_replicate_id": "rep_t1"},
        "t2": {"sample_id": "t2", "biological_source_id": "source_A", "biological_replicate_id": "rep_t2"},
        "t3": {"sample_id": "t3", "biological_source_id": "source_B", "biological_replicate_id": "rep_t3"},
        "c1": {"sample_id": "c1", "biological_source_id": "source_A", "biological_replicate_id": "rep_c1"},
        "c2": {"sample_id": "c2", "biological_source_id": "source_A", "biological_replicate_id": "rep_c2"},
    }
    res = DesignResolver.resolve(contrast, meta)
    assert res.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "partially_crossed_biological_sources" for f in res.findings)


def test_technical_replicates_trigger_needs_review():
    """Duplicate sample columns with same biological_replicate_id trigger NEEDS_REVIEW."""
    contrast = make_test_contrast(
        treatment_samples=("t1", "t2", "t3"),
        control_samples=("c1", "c2"),
        treatment_bio_reps=("rep_1", "rep_1", "rep_2"),  # t1 and t2 share rep_1
        control_bio_reps=("rep_c1", "rep_c2"),
        biological_sources=("cell_line_A",),
    )
    meta = {
        "t1": {"sample_id": "t1", "biological_replicate_id": "rep_1"},
        "t2": {"sample_id": "t2", "biological_replicate_id": "rep_1"},
        "t3": {"sample_id": "t3", "biological_replicate_id": "rep_2"},
        "c1": {"sample_id": "c1", "biological_replicate_id": "rep_c1"},
        "c2": {"sample_id": "c2", "biological_replicate_id": "rep_c2"},
    }
    res = DesignResolver.resolve(contrast, meta)
    assert res.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "unresolved_technical_replicates" for f in res.findings)


def test_insufficient_biological_replicates_blocked():
    """Single biological replicate in treatment or control -> BLOCKED."""
    contrast = make_test_contrast(
        treatment_samples=("t1",),
        control_samples=("c1", "c2"),
        treatment_bio_reps=("rep_t1",),
        control_bio_reps=("rep_c1", "rep_c2"),
        biological_sources=("cell_line_A",),
    )
    res = DesignResolver.resolve(contrast)
    assert res.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "insufficient_biological_replicates" for f in res.findings)


def test_zero_residual_degrees_of_freedom_blocked():
    """Design with 0 residual degrees of freedom -> BLOCKED."""
    # 2 paired sources with only 1 sample each (total 2 samples, but needs 3 parameters)
    contrast = make_test_contrast(
        treatment_samples=("t1",),
        control_samples=("c1",),
        treatment_bio_reps=("rep_t1",),
        control_bio_reps=("rep_c1",),
        biological_sources=("source_A",),
    )
    res = DesignResolver.resolve(contrast)
    assert res.status == ContrastStatus.BLOCKED
