"""Unit tests for response_builder models and dataclasses."""

import pytest

from src.response_builder.models import (
    ContrastStatus,
    Finding,
    MolecularResponseContrast,
    ResponseContrastDataset,
    Severity,
    Status,
    ordered_findings,
    status_for,
)


def test_contrast_status_enum() -> None:
    assert ContrastStatus.ELIGIBLE == "ELIGIBLE"
    assert ContrastStatus.NEEDS_REVIEW == "NEEDS_REVIEW"
    assert ContrastStatus.BLOCKED == "BLOCKED"


def test_ordered_findings() -> None:
    f1 = Finding(Severity.INFO, "r1", "entity", "e1", "p1", "info msg")
    f2 = Finding(Severity.ERROR, "r2", "entity", "e2", "p2", "error msg")
    f3 = Finding(Severity.REVIEW, "r3", "entity", "e3", "p3", "review msg")
    f4 = Finding(Severity.WARNING, "r4", "entity", "e4", "p4", "warn msg")

    sorted_f = ordered_findings([f1, f2, f3, f4])
    assert [f.severity for f in sorted_f] == [
        Severity.ERROR,
        Severity.REVIEW,
        Severity.WARNING,
        Severity.INFO,
    ]


def test_status_for() -> None:
    assert status_for([Finding(Severity.INFO, "r", "e", "id", "p", "m")]) == Status.PASS
    assert status_for([Finding(Severity.WARNING, "r", "e", "id", "p", "m")]) == Status.PASS_WITH_WARNINGS
    assert status_for([Finding(Severity.REVIEW, "r", "e", "id", "p", "m")]) == Status.NEEDS_REVIEW
    assert status_for([
        Finding(Severity.REVIEW, "r", "e", "id", "p", "m"),
        Finding(Severity.ERROR, "r2", "e", "id", "p", "m2"),
    ]) == Status.FAIL


def test_molecular_response_contrast_replication_and_eligibility() -> None:
    contrast = MolecularResponseContrast(
        contrast_id="contrast_01",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="cond_treat",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("s_t1", "s_t2", "s_t3"),
        control_sample_ids=("s_c1", "s_c2"),
        treatment_biological_replicate_ids=("rep1", "rep2", "rep3"),
        control_biological_replicate_ids=("c_rep1", "c_rep2"),
        agent_name="Valproic Acid",
        agent_identifier="CID_107",
        vehicle="DMSO",
        concentration_or_dose=200,
        concentration_or_dose_unit="uM",
        exposure_start_time_or_stage="Day 25",
        developmental_age_or_stage_at_exposure="Day 25",
        exposure_duration=5,
        exposure_duration_unit="days",
        washout_or_recovery_duration=0,
        washout_or_recovery_duration_unit="days",
        treatment_collection_age_or_stage="Day 30",
        control_collection_age_or_stage="Day 30",
        biological_source_ids=("bio_01",),
        organoid_context_ids=("org_01",),
        sample_comparison_exception_ids=(),
        status=ContrastStatus.ELIGIBLE,
        findings=(),
    )

    assert contrast.treatment_replicate_count == 3
    assert contrast.control_replicate_count == 2
    assert contrast.is_inferentially_eligible is True

    d = contrast.to_dict()
    assert d["contrast_id"] == "contrast_01"
    assert d["treatment_replicate_count"] == 3
    assert d["control_replicate_count"] == 2
    assert d["is_inferentially_eligible"] is True
    assert d["status"] == "ELIGIBLE"


def test_molecular_response_contrast_under_replicated_is_not_eligible() -> None:
    contrast = MolecularResponseContrast(
        contrast_id="contrast_under_replicated",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_01",
        treatment_condition_id="cond_treat",
        matched_control_condition_ids=("cond_ctrl",),
        treatment_sample_ids=("s_t1",),
        control_sample_ids=("s_c1", "s_c2"),
        treatment_biological_replicate_ids=("rep1",),
        control_biological_replicate_ids=("c_rep1", "c_rep2"),
        agent_name="Valproic Acid",
        agent_identifier="CID_107",
        vehicle="DMSO",
        concentration_or_dose=200,
        concentration_or_dose_unit="uM",
        exposure_start_time_or_stage="Day 25",
        developmental_age_or_stage_at_exposure="Day 25",
        exposure_duration=5,
        exposure_duration_unit="days",
        washout_or_recovery_duration=0,
        washout_or_recovery_duration_unit="days",
        treatment_collection_age_or_stage="Day 30",
        control_collection_age_or_stage="Day 30",
        biological_source_ids=("bio_01",),
        organoid_context_ids=("org_01",),
        sample_comparison_exception_ids=(),
        status=ContrastStatus.NEEDS_REVIEW,
        findings=(),
    )

    assert contrast.treatment_replicate_count == 1
    assert contrast.control_replicate_count == 2
    assert contrast.is_inferentially_eligible is False
