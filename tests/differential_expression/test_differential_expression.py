"""Comprehensive acceptance test suite for Step 7B Differential Expression v1.

Exercises all 28 required acceptance criteria:
1. Valid 3-treatment vs 3-control raw-count comparison -> successful DE result.
2. Positive-control synthetic fixture where a known gene is higher in treatment -> positive log2FoldChange.
3. Known lower-expression fixture -> negative log2FoldChange.
4. Reversing sample/input order -> identical contrast direction and equivalent results.
5. Day30 treatment vs Day30 control -> allowed.
6. Day30 treatment vs Day60 control cannot enter Step 7B because Step 7A blocks it.
7. Two dose conditions -> two separate DE results.
8. Shared explicit control -> separate treatment-specific DE results.
9. Single biological replicate -> DE execution blocked.
10. Technical replicates cannot inflate biological N.
11. Duplicate biological replicate columns without an approved aggregation rule -> NEEDS_REVIEW.
12. One cell line/source with replicated treatment and control -> design ~ condition.
13. Multiple paired biological sources represented in both arms -> design ~ biological_source + condition.
14. Perfect treatment/source confounding -> BLOCKED (treatment_source_confounding).
15. Independent unpaired unique biological sources -> valid ~ condition design.
16. Rank-deficient design -> blocked before DESeq2.
17. Locked Step 6C size factors are used exactly.
18. Missing or invalid size factor -> blocked.
19. Diagnostic normalized matrix accidentally supplied as DE input -> rejected.
20. Raw integer counts remain unchanged after DE.
21. DNT labels/evidence added to metadata -> identical DE result.
22. Multi-agent treatment retains all exposure metadata.
23. No log2FC/FDR threshold removes genes.
24. Gene with NA p/FDR remains represented.
25. Ambiguous/unmapped/colliding genes never become falsely mapped canonical gene-response features.
26. Exact R/DESeq2 runtime versions verified (strict checking).
27. Python wrapper output numerically agrees with direct R DESeq2 reference execution on the same fixture.
28. Full existing repository regression remains green.
"""

import copy
import math
import subprocess
import tempfile
from pathlib import Path
import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
    GeneDifferentialExpressionResult,
)
from src.differential_expression.pipeline import (
    run_contrast_differential_expression,
    run_differential_expression,
)
from src.differential_expression.r_bridge import (
    LOCKED_BIOC_VERSION,
    LOCKED_DESEQ2_VERSION,
    LOCKED_R_VERSION,
    probe_r_environment,
    resolve_r_binary,
)
from src.normalization.models import (
    NormalizedCohortResult,
    NormalizedCountMatrix,
    NormalizedDataset,
    SampleSizeFactor,
    Status as NormalizationStatus,
)
from src.response_builder.models import (
    ContrastExposure,
    ContrastStatus,
    Finding,
    MolecularResponseContrast,
    ResponseContrastDataset,
    Severity,
    Status,
)

from tests.differential_expression.conftest import create_synthetic_de_fixtures


def test_acceptance_1_valid_3trt_vs_3ctrl_raw_count_comparison():
    """1. Valid 3-treatment vs 3-control raw-count comparison -> successful DE result."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.is_inferentially_eligible is True
    assert len(res.gene_results) == 25
    assert res.design_formula == "~ condition"
    assert res.design_matrix_rank == 2
    assert res.residual_degrees_of_freedom == 4
    for gr in res.gene_results:
        assert gr.canonical_gene_id.startswith("ENSG")
        assert gr.base_mean > 0
        assert gr.log2_fold_change is not None
        assert gr.wald_statistic is not None
        assert gr.p_value is not None


def test_acceptance_2_positive_log2fc_orientation():
    """2. Positive-control synthetic fixture where a known gene is higher in treatment -> positive log2FoldChange."""
    custom = {
        "GENE_1": ([50, 50, 50], [200, 200, 200])  # control=50, treatment=200
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        custom_gene_counts=custom,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    gene1 = next(g for g in res.gene_results if g.canonical_gene_id == "ENSG00000000001")
    assert gene1.log2_fold_change is not None
    assert gene1.log2_fold_change > 0.0
    # Expected log2(200 / 50) = 2.0
    assert pytest.approx(gene1.log2_fold_change, abs=0.2) == 2.0


def test_acceptance_3_negative_log2fc_orientation():
    """3. Known lower-expression fixture -> negative log2FoldChange."""
    custom = {
        "GENE_1": ([200, 200, 200], [50, 50, 50])  # control=200, treatment=50
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        custom_gene_counts=custom,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    gene1 = next(g for g in res.gene_results if g.canonical_gene_id == "ENSG00000000001")
    assert gene1.log2_fold_change is not None
    assert gene1.log2_fold_change < 0.0
    # Expected log2(50 / 200) = -2.0
    assert pytest.approx(gene1.log2_fold_change, abs=0.2) == -2.0


def test_acceptance_4_input_and_sample_ordering_invariance():
    """4. Reversing sample/input order -> identical contrast direction and equivalent results."""
    custom = {
        "GENE_1": ([50, 50, 50], [200, 200, 200]),
        "GENE_2": ([150, 150, 150], [75, 75, 75]),
    }
    contrast1, norm_ds1, harm_ds1, meta1 = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        custom_gene_counts=custom,
    )
    res1 = run_contrast_differential_expression(
        contrast=contrast1,
        normalized_dataset=norm_ds1,
        harmonized_dataset=harm_ds1,
        sample_metadata=meta1,
    )

    # Reversed sample ordering in contrast
    contrast2 = MolecularResponseContrast(
        contrast_id=contrast1.contrast_id,
        study_id=contrast1.study_id,
        experiment_id=contrast1.experiment_id,
        cohort_id=contrast1.cohort_id,
        relationship_id=contrast1.relationship_id,
        treatment_condition_id=contrast1.treatment_condition_id,
        matched_control_condition_ids=contrast1.matched_control_condition_ids,
        treatment_sample_ids=tuple(reversed(contrast1.treatment_sample_ids)),
        control_sample_ids=tuple(reversed(contrast1.control_sample_ids)),
        treatment_biological_replicate_ids=tuple(reversed(contrast1.treatment_biological_replicate_ids)),
        control_biological_replicate_ids=tuple(reversed(contrast1.control_biological_replicate_ids)),
        biological_source_ids=contrast1.biological_source_ids,
        status=contrast1.status,
    )
    res2 = run_contrast_differential_expression(
        contrast=contrast2,
        normalized_dataset=norm_ds1,
        harmonized_dataset=harm_ds1,
        sample_metadata=meta1,
    )

    assert len(res1.gene_results) == len(res2.gene_results)
    for g1, g2 in zip(res1.gene_results, res2.gene_results):
        assert g1.canonical_gene_id == g2.canonical_gene_id
        if g1.log2_fold_change is not None and g2.log2_fold_change is not None:
            assert pytest.approx(g1.log2_fold_change, abs=1e-5) == g2.log2_fold_change
        if g1.p_value is not None and g2.p_value is not None:
            assert pytest.approx(g1.p_value, abs=1e-5) == g2.p_value


def test_acceptance_5_day30_treatment_vs_day30_control_allowed():
    """5. Day30 treatment vs Day30 control -> allowed."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        developmental_age="Day 30",
        control_developmental_age="Day 30",
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.developmental_age_or_stage_at_exposure == "Day 30"
    assert res.treatment_collection_age_normalized == "Day 30"
    assert res.control_collection_age_normalized == "Day 30"


def test_acceptance_6_day30_treatment_vs_day60_control_blocked_upstream():
    """6. Day30 treatment vs Day60 control cannot enter Step 7B because Step 7A blocks it."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        developmental_age="Day 30",
        control_developmental_age="Day 60",
        contrast_status=ContrastStatus.BLOCKED,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    # DESeq2 must NOT have run; result must be BLOCKED with 0 gene results
    assert res.status == ContrastStatus.BLOCKED
    assert len(res.gene_results) == 0
    assert any("skipped" in f.message.lower() for f in res.findings)


def test_acceptance_7_two_dose_conditions_produce_two_separate_de_results():
    """7. Two dose conditions -> two separate DE results."""
    contrast_low, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        exposures=[ContrastExposure(exposure_id="e1", agent_name="VPA", concentration_or_dose=100.0, concentration_or_dose_unit="uM")],
    )
    # Create second contrast with dose=200
    contrast_high = MolecularResponseContrast(
        contrast_id="contrast_high_dose",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_02",
        treatment_condition_id="treat_high",
        matched_control_condition_ids=("ctrl_01",),
        treatment_sample_ids=contrast_low.treatment_sample_ids,
        control_sample_ids=contrast_low.control_sample_ids,
        treatment_biological_replicate_ids=contrast_low.treatment_biological_replicate_ids,
        control_biological_replicate_ids=contrast_low.control_biological_replicate_ids,
        treatment_exposures=(
            ContrastExposure(exposure_id="e2", agent_name="VPA", concentration_or_dose=200.0, concentration_or_dose_unit="uM"),
        ),
        biological_source_ids=contrast_low.biological_source_ids,
        status=ContrastStatus.ELIGIBLE,
    )
    contrast_dataset = ResponseContrastDataset(
        status=Status.PASS,
        findings=(),
        contrasts=(contrast_low, contrast_high),
        treatment_conditions_encountered=("treat_01", "treat_high"),
        source_asset_id="asset_01",
        reference_identity="Ensembl_112",
        builder_version="1.0.0",
        metadata_hash="dummy_hash",
    )
    de_dataset = run_differential_expression(
        contrast_dataset=contrast_dataset,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert len(de_dataset.contrast_results) == 2
    assert de_dataset.contrast_results[0].contrast_id == "contrast_001"
    assert de_dataset.contrast_results[1].contrast_id == "contrast_high_dose"
    assert de_dataset.contrast_results[0].treatment_exposures[0].concentration_or_dose == 100.0
    assert de_dataset.contrast_results[1].treatment_exposures[0].concentration_or_dose == 200.0


def test_acceptance_8_shared_explicit_control_separate_treatment_results():
    """8. Shared explicit control -> separate treatment-specific DE results."""
    contrast1, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
    )
    contrast2 = MolecularResponseContrast(
        contrast_id="contrast_cmp2",
        study_id="study_01",
        experiment_id="exp_01",
        cohort_id="cohort_01",
        relationship_id="rel_02",
        treatment_condition_id="treat_cmp2",
        matched_control_condition_ids=("ctrl_01",),
        treatment_sample_ids=contrast1.treatment_sample_ids,
        control_sample_ids=contrast1.control_sample_ids,  # Shared control!
        treatment_biological_replicate_ids=contrast1.treatment_biological_replicate_ids,
        control_biological_replicate_ids=contrast1.control_biological_replicate_ids,
        biological_source_ids=contrast1.biological_source_ids,
        status=ContrastStatus.ELIGIBLE,
    )
    contrast_dataset = ResponseContrastDataset(
        status=Status.PASS,
        findings=(),
        contrasts=(contrast1, contrast2),
        treatment_conditions_encountered=("treat_01", "treat_cmp2"),
        source_asset_id="asset_01",
        reference_identity="Ensembl_112",
        builder_version="1.0.0",
        metadata_hash="dummy_hash",
    )
    de_dataset = run_differential_expression(
        contrast_dataset=contrast_dataset,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert len(de_dataset.contrast_results) == 2
    assert de_dataset.contrast_results[0].control_sample_ids == de_dataset.contrast_results[1].control_sample_ids


def test_acceptance_9_single_biological_replicate_blocked():
    """9. Single biological replicate -> DE execution blocked."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=1,  # Single replicate!
        n_control=3,
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.BLOCKED
    assert len(res.gene_results) == 0
    assert any(f.rule_id == "insufficient_biological_replicates" or "skipped" in f.message.lower() for f in res.findings)


def test_acceptance_10_technical_replicates_cannot_inflate_biological_n():
    """10. Technical replicates cannot inflate biological N."""
    # 2 treatment samples, but both share rep_1
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=2,
        n_control=2,
        treatment_bio_reps=["rep_1", "rep_1"],  # technical reps!
        control_bio_reps=["rep_c1", "rep_c2"],
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    # Cannot be ELIGIBLE
    assert res.status != ContrastStatus.ELIGIBLE
    assert len(res.gene_results) == 0


def test_acceptance_11_duplicate_biological_replicate_columns_needs_review():
    """11. Duplicate biological replicate columns without an approved aggregation rule -> NEEDS_REVIEW."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        treatment_bio_reps=["rep_1", "rep_1", "rep_2"],  # duplicate rep_1
        control_bio_reps=["rep_c1", "rep_c2", "rep_c3"],
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.NEEDS_REVIEW
    assert any(f.rule_id == "unresolved_technical_replicates" for f in res.findings)
    assert len(res.gene_results) == 0


def test_acceptance_12_single_source_replicated_design_condition():
    """12. One cell line/source with replicated treatment and control -> design ~ condition."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        treatment_source="cell_line_1",
        control_source="cell_line_1",
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.design_formula == "~ condition"


def test_acceptance_13_paired_biological_sources_design_source_plus_condition():
    """13. Multiple paired biological sources represented in both arms -> design ~ biological_source + condition."""
    sample_sources = {
        "ctrl_s1": "donor_A", "ctrl_s2": "donor_B",
        "trt_s1": "donor_A", "trt_s2": "donor_B",
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=2,
        n_control=2,
        sample_sources=sample_sources,
        biological_sources=["donor_A", "donor_B"],
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.design_formula == "~ biological_source + condition"
    assert res.design_matrix_rank == 3
    assert res.residual_degrees_of_freedom == 1


def test_acceptance_14_perfect_treatment_source_confounding_blocked():
    """14. Perfect treatment/source confounding -> BLOCKED (treatment_source_confounding); DESeq2 not executed."""
    sample_sources = {
        "ctrl_s1": "source_A", "ctrl_s2": "source_A", "ctrl_s3": "source_A",
        "trt_s1": "source_B", "trt_s2": "source_B", "trt_s3": "source_B",
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        sample_sources=sample_sources,
        biological_sources=["source_A", "source_B"],
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "treatment_source_confounding" for f in res.findings)
    assert len(res.gene_results) == 0


def test_acceptance_15_independent_unpaired_unique_sources():
    """15. Independent unpaired unique biological sources -> valid ~ condition design."""
    sample_sources = {
        "ctrl_s1": "donor_1", "ctrl_s2": "donor_2", "ctrl_s3": "donor_3",
        "trt_s1": "donor_4", "trt_s2": "donor_5", "trt_s3": "donor_6",
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        sample_sources=sample_sources,
        biological_sources=[f"donor_{i}" for i in range(1, 7)],
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert res.design_formula == "~ condition"


def test_acceptance_16_rank_deficient_design_blocked_before_deseq2():
    """16. Rank-deficient design -> blocked before DESeq2."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=1,
        n_control=1,
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.BLOCKED
    assert len(res.gene_results) == 0


def test_acceptance_17_locked_step6c_size_factors_used_exactly():
    """17. Locked Step 6C size factors are used exactly."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        size_factor_val=1.234,
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    for sid, sf in res.size_factors_used.items():
        assert sf == 1.234


def test_acceptance_18_missing_or_invalid_size_factor_blocked():
    """18. Missing or invalid size factor -> blocked."""
    # Negative size factor
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        invalid_size_factor=-1.0,
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "invalid_size_factors" for f in res.findings)


def test_acceptance_19_diagnostic_normalized_matrix_accidentally_supplied_rejected():
    """19. Diagnostic normalized matrix accidentally supplied as DE input -> rejected."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        use_diagnostic_as_raw=True,
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.BLOCKED
    assert any(f.rule_id == "invalid_count_matrix_type" for f in res.findings)


def test_acceptance_20_raw_integer_counts_immutable():
    """20. Raw integer counts remain unchanged after DE."""
    import hashlib

    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
    )
    cohort = norm_ds.cohort_results[0]
    raw_mat = cohort.raw_count_matrix

    orig_gene_ids = copy.deepcopy(raw_mat.gene_ids)
    orig_sample_ids = copy.deepcopy(raw_mat.sample_ids)
    orig_cols = copy.deepcopy(raw_mat.columns)
    orig_hash = hashlib.sha256(str((orig_gene_ids, orig_sample_ids, orig_cols)).encode()).hexdigest()

    run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )

    post_hash = hashlib.sha256(str((raw_mat.gene_ids, raw_mat.sample_ids, raw_mat.columns)).encode()).hexdigest()
    assert raw_mat.gene_ids == orig_gene_ids
    assert raw_mat.sample_ids == orig_sample_ids
    assert raw_mat.columns == orig_cols
    assert post_hash == orig_hash


def test_acceptance_21_dnt_labels_have_zero_effect_on_de_output():
    """21. DNT labels/evidence added to metadata -> identical DE result."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
    )
    res_unlabelled = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )

    # Add DNT labels / OECD evidence to metadata
    meta_with_dnt = copy.deepcopy(meta)
    for sid in meta_with_dnt:
        meta_with_dnt[sid]["dnt_label"] = "DNT_POSITIVE"
        meta_with_dnt[sid]["oecd_guideline"] = "TG_426"
        meta_with_dnt[sid]["dnt_target_list"] = ["GENE_1", "GENE_5"]

    res_labelled = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta_with_dnt,
    )

    assert len(res_unlabelled.gene_results) == len(res_labelled.gene_results)
    for g1, g2 in zip(res_unlabelled.gene_results, res_labelled.gene_results):
        assert g1.canonical_gene_id == g2.canonical_gene_id
        if g1.log2_fold_change is not None:
            assert pytest.approx(g1.log2_fold_change, abs=1e-6) == g2.log2_fold_change
        if g1.p_value is not None:
            assert pytest.approx(g1.p_value, abs=1e-6) == g2.p_value


def test_acceptance_22_multi_agent_treatment_retains_all_exposure_metadata():
    """22. Multi-agent treatment retains all exposure metadata."""
    exp1 = ContrastExposure(exposure_id="e1", agent_name="CompoundA", concentration_or_dose=10.0, concentration_or_dose_unit="uM")
    exp2 = ContrastExposure(exposure_id="e2", agent_name="CompoundB", concentration_or_dose=50.0, concentration_or_dose_unit="nM")
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        exposures=[exp1, exp2],
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert len(res.treatment_exposures) == 2
    assert res.treatment_exposures[0].agent_name == "CompoundA"
    assert res.treatment_exposures[1].agent_name == "CompoundB"


def test_acceptance_23_no_log2fc_or_fdr_threshold_removes_genes():
    """23. No log2FC/FDR threshold removes genes."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    # Exactly all 25 genes must be returned even if FDR > 0.05
    assert len(res.gene_results) == 25


def test_acceptance_24_gene_with_na_statistics_preserved():
    """24. Gene with NA p/FDR remains represented."""
    custom = {
        "GENE_1": ([0, 0, 0], [0, 0, 0]),  # all zeros
        "GENE_2": ([10, 10, 10], [10, 10, 100000]),  # extreme outlier
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        custom_gene_counts=custom,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE

    # Zero count gene
    gene_zero = next(g for g in res.gene_results if g.canonical_gene_id == "ENSG00000000001")
    assert gene_zero.base_mean == 0.0
    assert gene_zero.log2_fold_change is None
    assert gene_zero.p_value is None
    assert gene_zero.result_status == "ZERO_COUNTS"

    # Outlier gene
    gene_outlier = next(g for g in res.gene_results if g.canonical_gene_id == "ENSG00000000002")
    assert gene_outlier.p_value is None
    assert gene_outlier.result_status == "OUTLIER"


def test_acceptance_25_ambiguous_and_colliding_genes_excluded_from_canonical_features():
    """25. Ambiguous/unmapped/colliding genes never become falsely mapped canonical gene-response features."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        include_unmapped_gene=True,
        include_colliding_gene=True,
        include_ambiguous_gene=True,
        include_invalid_gene=True,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    # Exactly the 25 uniquely mapped genes should be in output
    assert len(res.gene_results) == 25
    can_ids = {g.canonical_gene_id for g in res.gene_results}
    assert "UNMAPPED_GENE" not in can_ids
    assert "ENSG_COLLISION_01" not in can_ids
    assert "AMBIGUOUS_GENE" not in can_ids
    assert "INVALID_GENE" not in can_ids

    # Verify excluded_genes audit
    assert len(res.excluded_genes) == 5
    assert res.total_input_genes_count == 30
    assert res.eligible_canonical_genes_count == 25
    assert res.excluded_genes_count == 5
    assert res.is_gene_universe_conserved is True


def test_gene_universe_case_a_unmapped_gene_in_audit():
    """Case A: unmapped gene is excluded from canonical DE but appears in audit."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        include_unmapped_gene=True,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert "UNMAPPED_GENE" not in {g.canonical_gene_id for g in res.gene_results}
    audit = next(eg for eg in res.excluded_genes if eg.original_gene_id == "UNMAPPED_GENE")
    assert audit.mapping_status == "UNMAPPED"
    assert audit.canonical_gene_id is None
    assert audit.source_index == 25
    assert audit.exclusion_reason != ""


def test_gene_universe_case_b_ambiguous_gene_in_audit():
    """Case B: ambiguous gene is excluded but appears in audit."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        include_ambiguous_gene=True,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert "AMBIGUOUS_GENE" not in {g.canonical_gene_id for g in res.gene_results}
    audit = next(eg for eg in res.excluded_genes if eg.original_gene_id == "AMBIGUOUS_GENE")
    assert audit.mapping_status == "AMBIGUOUS"
    assert audit.source_index == 25
    assert "multiple candidate" in audit.exclusion_reason


def test_gene_universe_case_c_collision_genes_in_audit():
    """Case C: collision gene is excluded but appears in audit."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        include_colliding_gene=True,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    assert "ENSG_COLLISION_01" not in {g.canonical_gene_id for g in res.gene_results}
    coll1 = next(eg for eg in res.excluded_genes if eg.original_gene_id == "COLLIDING_GENE_1")
    coll2 = next(eg for eg in res.excluded_genes if eg.original_gene_id == "COLLIDING_GENE_2")
    assert coll1.mapping_status == "COLLISION"
    assert coll1.canonical_gene_id == "ENSG_COLLISION_01"
    assert coll2.mapping_status == "COLLISION"
    assert coll2.canonical_gene_id == "ENSG_COLLISION_01"


def test_gene_universe_case_d_conservation_equation():
    """Case D: conservation: total input rows = DE genes + all harmonization-excluded rows."""
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        include_unmapped_gene=True,
        include_colliding_gene=True,
        include_ambiguous_gene=True,
        include_invalid_gene=True,
    )
    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    raw_mat = norm_ds.cohort_results[0].raw_count_matrix
    assert len(raw_mat.gene_ids) == 30
    assert res.total_input_genes_count == 30
    assert res.eligible_canonical_genes_count == 25
    assert res.excluded_genes_count == 5
    assert res.total_input_genes_count == (res.eligible_canonical_genes_count + res.excluded_genes_count)
    assert res.is_gene_universe_conserved is True


def test_gene_universe_exclusion_independent_of_outcome_and_dnt_labels():
    """Exclusion decision must be independent of treatment outcome, log2FC, p-value, FDR, and DNT label."""
    custom = {
        "COLLIDING_GENE_1": ([50, 50, 50], [5000, 5000, 5000]),
        "UNMAPPED_GENE": ([100, 100, 100], [0, 0, 0]),
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        include_unmapped_gene=True,
        include_colliding_gene=True,
        custom_gene_counts=custom,
    )
    for sid in meta:
        meta[sid]["dnt_label"] = "DNT_POSITIVE"

    res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert res.status == ContrastStatus.ELIGIBLE
    can_ids = {g.canonical_gene_id for g in res.gene_results}
    assert "UNMAPPED_GENE" not in can_ids
    assert "ENSG_COLLISION_01" not in can_ids
    assert len(res.excluded_genes) == 3
    assert res.is_gene_universe_conserved is True


def test_acceptance_26_exact_r_deseq2_runtime_versions_verified():
    """26. Exact R/DESeq2 runtime versions verified (strict checking)."""
    env_info = probe_r_environment()
    assert env_info.is_available is True
    assert env_info.r_version == LOCKED_R_VERSION
    assert env_info.bioc_version == LOCKED_BIOC_VERSION
    assert env_info.deseq2_version == LOCKED_DESEQ2_VERSION


def test_acceptance_27_python_wrapper_agrees_with_direct_r_deseq2():
    """27. Python wrapper output numerically agrees with direct R DESeq2 reference execution on the same fixture."""
    custom = {
        "GENE_1": ([50, 55, 45], [200, 220, 190]),
    }
    contrast, norm_ds, harm_ds, meta = create_synthetic_de_fixtures(
        n_treatment=3,
        n_control=3,
        n_genes=25,
        custom_gene_counts=custom,
    )
    py_res = run_contrast_differential_expression(
        contrast=contrast,
        normalized_dataset=norm_ds,
        harmonized_dataset=harm_ds,
        sample_metadata=meta,
    )
    assert py_res.status == ContrastStatus.ELIGIBLE
    py_gene1 = next(g for g in py_res.gene_results if g.canonical_gene_id == "ENSG00000000001")

    # Direct R execution on the exact same matrix
    raw_mat = norm_ds.cohort_results[0].raw_count_matrix
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        counts_tsv = tmp / "counts.tsv"
        with counts_tsv.open("w") as f:
            f.write("gene_id\t" + "\t".join(raw_mat.sample_ids) + "\n")
            for i, gid in enumerate(raw_mat.gene_ids):
                can_id = f"ENSG{i+1:011d}"
                f.write(can_id + "\t" + "\t".join(str(raw_mat.columns[j][i]) for j in range(len(raw_mat.sample_ids))) + "\n")

        r_bin = resolve_r_binary()
        r_code = f"""
        suppressPackageStartupMessages(library(DESeq2))
        counts <- as.matrix(read.delim('{counts_tsv.as_posix()}', row.names=1, sep='\\t'))
        storage.mode(counts) <- 'integer'
        coldata <- data.frame(
            condition=factor(c('control', 'control', 'control', 'treatment', 'treatment', 'treatment'), levels=c('control', 'treatment')),
            row.names=colnames(counts)
        )
        dds <- DESeqDataSetFromMatrix(counts, coldata, ~ condition)
        sizeFactors(dds) <- rep(1, 6)
        dds <- DESeq(dds, test='Wald', betaPrior=FALSE, minReplicatesForReplace=Inf, quiet=TRUE)
        res <- results(dds, contrast=c('condition', 'treatment', 'control'), independentFiltering=FALSE)
        cat(sprintf('%.6f,%.6f,%.8e', res['ENSG00000000001', 'log2FoldChange'], res['ENSG00000000001', 'stat'], res['ENSG00000000001', 'pvalue']))
        """
        proc = subprocess.run([r_bin, "--vanilla", "-e", r_code], capture_output=True, text=True)
        assert proc.returncode == 0
        parts = proc.stdout.strip().split(",")
        r_lfc = float(parts[0])
        r_stat = float(parts[1])
        r_pval = float(parts[2])

        assert pytest.approx(py_gene1.log2_fold_change, abs=1e-4) == r_lfc
        assert pytest.approx(py_gene1.wald_statistic, abs=1e-4) == r_stat
        assert pytest.approx(py_gene1.p_value, rel=1e-3) == r_pval
