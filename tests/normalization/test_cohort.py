"""Unit tests for metadata-derived Normalization Cohort construction."""

from typing import Any
import pytest

from src.normalization.cohort import NormalizationCohortBuilder
from src.normalization.models import ResponseEligibility, Severity


def create_mock_metadata(
    study_id: str = "study_001",
    exp_id: str = "exp_001",
    assay_id_t: str = "assay_lane_1",
    assay_id_c: str = "assay_lane_2",
    modality: str = "bulk_rna_seq",
    strat_t: str = "polyA",
    strat_c: str = "polyA",
) -> dict[str, Any]:
    return {
        "studies": [{"study_id": study_id}],
        "experiments": [{"experiment_id": exp_id, "study_id": study_id}],
        "conditions": [
            {
                "condition_id": "cond_treat",
                "experiment_id": exp_id,
                "treatment_control_status": {"original_value": "treatment", "value_status": "reported"},
            },
            {
                "condition_id": "cond_ctrl",
                "experiment_id": exp_id,
                "treatment_control_status": {"original_value": "control", "value_status": "reported"},
            },
        ],
        "samples": [
            {"sample_id": "s_t1", "experiment_id": exp_id, "condition_id": "cond_treat"},
            {"sample_id": "s_t2", "experiment_id": exp_id, "condition_id": "cond_treat"},
            {"sample_id": "s_c1", "experiment_id": exp_id, "condition_id": "cond_ctrl"},
            {"sample_id": "s_c2", "experiment_id": exp_id, "condition_id": "cond_ctrl"},
        ],
        "sequencing_assays": [
            {
                "assay_id": assay_id_t,
                "sample_ids": ["s_t1", "s_t2"],
                "rna_seq_modality": {"original_value": modality, "value_status": "reported"},
                "library_strategy": {"original_value": strat_t, "value_status": "reported"},
                "sequencing_platform": {"original_value": "Illumina NovaSeq", "value_status": "reported"},
                "strandedness": {"original_value": "reverse", "value_status": "reported"},
                "read_layout": {"original_value": "paired_end", "value_status": "reported"},
            },
            {
                "assay_id": assay_id_c,
                "sample_ids": ["s_c1", "s_c2"],
                "rna_seq_modality": {"original_value": modality, "value_status": "reported"},
                "library_strategy": {"original_value": strat_c, "value_status": "reported"},
                "sequencing_platform": {"original_value": "Illumina NovaSeq", "value_status": "reported"},
                "strandedness": {"original_value": "reverse", "value_status": "reported"},
                "read_layout": {"original_value": "paired_end", "value_status": "reported"},
            },
        ],
        "treatment_control_relationships": [
            {
                "relationship_id": "rel_01",
                "treatment_condition_ids": ["cond_treat"],
                "matched_control_condition_ids": ["cond_ctrl"],
            }
        ],
    }


def test_cohort_preserves_treatment_and_matched_controls_across_assays() -> None:
    # Treatment is on assay_lane_1, control is on assay_lane_2
    meta = create_mock_metadata(assay_id_t="assay_lane_1", assay_id_c="assay_lane_2")
    builder = NormalizationCohortBuilder(meta)
    sample_ids = ("s_t1", "s_t2", "s_c1", "s_c2")

    cohorts, findings = builder.build_cohorts(sample_ids)
    assert not any(f.severity == Severity.ERROR for f in findings)
    assert len(cohorts) == 1

    cohort = cohorts[0]
    assert set(cohort.sample_ids) == set(sample_ids)
    assert set(cohort.assay_ids) == {"assay_lane_1", "assay_lane_2"}
    assert set(cohort.treatment_sample_ids) == {"s_t1", "s_t2"}
    assert set(cohort.control_sample_ids) == {"s_c1", "s_c2"}

    # Resolve controls
    controls = builder.resolve_sample_controls(cohort)
    assert controls["s_t1"][0] == ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL
    assert set(controls["s_t1"][1]) == {"s_c1", "s_c2"}
    assert controls["s_c1"][0] == ResponseEligibility.CONTROL_BASELINE


def test_cohort_strictly_rejects_cross_study_pooling() -> None:
    meta = create_mock_metadata()
    # Add a second study and sample
    meta["studies"].append({"study_id": "study_002"})
    meta["experiments"].append({"experiment_id": "exp_002", "study_id": "study_002"})
    meta["conditions"].append({"condition_id": "cond_002", "experiment_id": "exp_002", "treatment_control_status": "treatment"})
    meta["samples"].append({"sample_id": "s_foreign", "experiment_id": "exp_002", "condition_id": "cond_002"})
    meta["sequencing_assays"].append({
        "assay_id": "assay_foreign",
        "sample_ids": ["s_foreign"],
        "rna_seq_modality": "bulk_rna_seq",
    })

    builder = NormalizationCohortBuilder(meta)
    cohorts, findings = builder.build_cohorts(("s_t1", "s_foreign"))
    assert any(f.rule_id == "cross_study_pooling_forbidden" for f in findings)
    assert len(cohorts) == 0


def test_cohort_rejects_non_bulk_modality() -> None:
    meta = create_mock_metadata(modality="scRNA_seq")
    builder = NormalizationCohortBuilder(meta)
    cohorts, findings = builder.build_cohorts(("s_t1", "s_c1"))
    assert any(f.rule_id == "incompatible_modality" for f in findings)


def test_cohort_detects_conflicting_library_strategies() -> None:
    meta = create_mock_metadata(strat_t="polyA", strat_c="ribo-depletion")
    builder = NormalizationCohortBuilder(meta)
    cohorts, findings = builder.build_cohorts(("s_t1", "s_t2", "s_c1", "s_c2"))
    assert any(f.rule_id == "cohort_technical_context_ambiguous" for f in findings)
    assert len(cohorts) == 1  # Formed with review finding


def test_cohort_tags_missing_matched_controls() -> None:
    meta = create_mock_metadata()
    # Remove relationship so treatment condition has no control
    meta["treatment_control_relationships"] = []
    builder = NormalizationCohortBuilder(meta)
    cohorts, findings = builder.build_cohorts(("s_t1", "s_t2"))
    assert len(cohorts) == 1
    cohort = cohorts[0]
    assert len(cohort.control_sample_ids) == 0

    controls = builder.resolve_sample_controls(cohort)
    assert controls["s_t1"][0] == ResponseEligibility.INELIGIBLE_MISSING_MATCHED_CONTROL
    assert controls["s_t1"][1] == ()
