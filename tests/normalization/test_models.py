"""Unit tests for Normalization v1 data models, enums, and dispositions."""

import pytest

from src.normalization.models import (
    CohortNormalizationMetrics,
    CollisionDisposition,
    Finding,
    NormalizationCohort,
    NormalizationMethod,
    NormalizationReviewDisposition,
    NormalizedCountMatrix,
    ResponseEligibility,
    SampleSizeFactor,
    Severity,
    Status,
    ordered,
    status_for,
)


def test_response_eligibility_states() -> None:
    assert ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL == "ELIGIBLE_WITH_MATCHED_CONTROL"
    assert ResponseEligibility.INELIGIBLE_MISSING_MATCHED_CONTROL == "INELIGIBLE_MISSING_MATCHED_CONTROL"
    assert ResponseEligibility.CONTROL_BASELINE == "CONTROL_BASELINE"
    assert ResponseEligibility.UNASSIGNED_OR_AMBIGUOUS == "UNASSIGNED_OR_AMBIGUOUS"


def test_normalization_method_enum() -> None:
    assert NormalizationMethod.RATIO == "ratio"
    assert NormalizationMethod.POSCOUNTS == "poscounts"


def test_sample_size_factor_model() -> None:
    sf = SampleSizeFactor(
        sample_id="sample_01",
        size_factor=1.234,
        raw_library_size=1000000,
        normalized_library_size=810372.77,
        cohort_id="cohort_01",
        response_eligibility=ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL,
        matched_control_sample_ids=("ctrl_01", "ctrl_02"),
        treatment_control_status="treatment",
    )
    d = sf.to_dict()
    assert d["sample_id"] == "sample_01"
    assert d["size_factor"] == 1.234
    assert d["response_eligibility"] == "ELIGIBLE_WITH_MATCHED_CONTROL"
    assert d["matched_control_sample_ids"] == ["ctrl_01", "ctrl_02"]


def test_normalized_count_matrix() -> None:
    mat = NormalizedCountMatrix(
        gene_ids=("ENSG001", "ENSG002"),
        sample_ids=("s1", "s2"),
        columns=((10.5, 20.0), (15.2, 30.1)),
    )
    assert mat.get_count("ENSG001", "s1") == 10.5
    assert mat.get_count("ENSG002", "s2") == 30.1
    assert mat.get_count("NONEXISTENT", "s1") is None
    d = mat.to_dict()
    assert d["number_of_genes"] == 2
    assert d["number_of_samples"] == 2


def test_normalization_review_disposition() -> None:
    disp = NormalizationReviewDisposition(
        decision="AUTHORIZE_POSCOUNTS",
        rationale="Sparse dataset where all genes contain zero; poscounts scientifically justified.",
        reviewer="curator_01",
        review_date="2026-09-24",
        affected_cohort_id="cohort_exp1_1",
    )
    assert disp.decision == "AUTHORIZE_POSCOUNTS"
    d = disp.to_dict()
    assert d["reviewer"] == "curator_01"

    coerced = NormalizationReviewDisposition.from_any({
        "decision": "AUTHORIZE_POSCOUNTS",
        "rationale": "Justified",
        "reviewer": "curator_02",
        "review_date": "2026-09-24",
    })
    assert coerced is not None
    assert coerced.reviewer == "curator_02"

    with pytest.raises(ValueError):
        NormalizationReviewDisposition(decision="", rationale="r", reviewer="rev", review_date="d")


def test_collision_disposition() -> None:
    cdisp = CollisionDisposition(
        decision="DROP_COLLISIONS_FROM_ESTIMATION",
        rationale="Collisions represent duplicated annotations; dropping from size factor estimation.",
        reviewer="curator_01",
        review_date="2026-09-24",
        resolved_collision_ids=("ENSG00000123456",),
    )
    assert cdisp.decision == "DROP_COLLISIONS_FROM_ESTIMATION"
    assert cdisp.resolved_collision_ids == ("ENSG00000123456",)

    coerced = CollisionDisposition.from_any({
        "decision": "DROP",
        "rationale": "Valid reason",
        "reviewer": "curator_01",
        "review_date": "2026-09-24",
    })
    assert coerced is not None

    with pytest.raises(ValueError):
        CollisionDisposition(decision="", rationale="r", reviewer="rev", review_date="d")


def test_finding_ordering_and_status() -> None:
    f_info = Finding(Severity.INFO, "info_rule", "sample", "s1", "path", "info")
    f_warn = Finding(Severity.WARNING, "warn_rule", "sample", "s1", "path", "warn")
    f_rev = Finding(Severity.REVIEW, "rev_rule", "sample", "s1", "path", "review")
    f_err = Finding(Severity.ERROR, "err_rule", "sample", "s1", "path", "error")

    ordered_findings = ordered([f_info, f_err, f_warn, f_rev])
    assert ordered_findings[0].severity == Severity.ERROR
    assert ordered_findings[1].severity == Severity.REVIEW
    assert ordered_findings[2].severity == Severity.WARNING
    assert ordered_findings[3].severity == Severity.INFO

    assert status_for([f_info]) == Status.PASS
    assert status_for([f_warn, f_info]) == Status.PASS_WITH_WARNINGS
    assert status_for([f_rev, f_warn]) == Status.NEEDS_REVIEW
    assert status_for([f_err, f_rev]) == Status.FAIL
