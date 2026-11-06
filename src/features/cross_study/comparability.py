"""Stage 8A.2: Cross-Study Comparability Evaluation and Alignment Framework.

Governing Scientific Invariants:
1. Comparability Gating: Independent studies testing the same chemical are evaluated for
   regimen equivalence (dose, developmental stage, duration, organoid model, biological replication).
2. Non-Comparable Context: When experimental regimens differ (e.g. Day 6 vs Day 28 organoids,
   200 uM vs 250 uM VPA), formal statistical pooling or meta-analysis is strictly prohibited
   (CrossStudyStatus.NONCOMPARABLE_CONTEXT).
3. No Unweighted Averaging: Unweighted fold change averaging or synthetic variance manufacturing
   is strictly prohibited across independent studies.
4. Lossless Record Preservation: Condition-level records from each study are preserved with 100%
   numerical and identity fidelity.
5. Descriptive Direction Alignment: Matched canonical Ensembl genes are aligned descriptively
   (CONCORDANT_UP, CONCORDANT_DOWN, DISCORDANT, NO_DIRECTION, UNDETERMINED) without claiming
   toxicity outcomes or consuming DNT labels.
6. Law of Conservation of the Gene Universe: All genes evaluated in Study 1 and Study 2 are conserved.
"""

from datetime import datetime, timezone
import math
from typing import Any

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


def derive_direction_concordance(
    dir1: ResponseDirection,
    dir2: ResponseDirection,
) -> DirectionConcordance:
    """Derive descriptive concordance of response directions between two independent studies.

    This is non-evaluative: does NOT imply biological significance, mechanism, or DNT outcome.
    """
    if dir1 == ResponseDirection.UNDETERMINED or dir2 == ResponseDirection.UNDETERMINED:
        return DirectionConcordance.UNDETERMINED
    if dir1 == ResponseDirection.NO_DIRECTION or dir2 == ResponseDirection.NO_DIRECTION:
        return DirectionConcordance.NO_DIRECTION
    if dir1 == ResponseDirection.UP and dir2 == ResponseDirection.UP:
        return DirectionConcordance.CONCORDANT_UP
    if dir1 == ResponseDirection.DOWN and dir2 == ResponseDirection.DOWN:
        return DirectionConcordance.CONCORDANT_DOWN
    if (dir1 == ResponseDirection.UP and dir2 == ResponseDirection.DOWN) or (
        dir1 == ResponseDirection.DOWN and dir2 == ResponseDirection.UP
    ):
        return DirectionConcordance.DISCORDANT
    return DirectionConcordance.UNDETERMINED


def _is_valid_metadata_value(val: Any) -> bool:
    """Check if metadata value is present, non-empty, and not an uninformative placeholder."""
    if val is None:
        return False
    if isinstance(val, str):
        s = val.strip().lower()
        if not s or s in ("unknown", "unspecified", "none", "n/a", "na", "null", "missing"):
            return False
    return True


def evaluate_study_comparability(
    record1: ConditionGeneResponseRecord,
    record2: ConditionGeneResponseRecord,
    replicates_verified1: bool | None = None,
    replicates_verified2: bool | None = None,
) -> StudyComparabilityAssessment:
    """Evaluate biological, chemical, and statistical comparability between two study conditions.

    Governed Scientific Invariants:
    1. Independent study verification: Compared records must have explicit evidence of coming
       from independent studies or study families. Different cell lines, donors, or conditions
       within the same study cannot count as independent studies.
    2. Biological replication: Missing replicate verification must NOT default to approved;
       both studies require explicit validated status (replicates_verified is True) and >= 2 reps.
    3. Regimen & model matching: Dose, developmental stage, exposure duration, AND organoid model
       system must match.
    4. Missing metadata: Missing, unknown, or unspecified metadata is strictly NOT interpreted as
       verified equivalence.
    5. Final synthesis eligibility: True ONLY if independent studies, compound, dose, stage,
       duration, organoid model match, and replicates are verified.
    """
    differences: list[str] = []

    # 1. Study Independence Verification
    s1_id = (record1.study_id or "").strip()
    s2_id = (record2.study_id or "").strip()
    s1_fam = (record1.study_family_id or "").strip()
    s2_fam = (record2.study_family_id or "").strip()
    s1_ind = (record1.independent_study_id or "").strip()
    s2_ind = (record2.independent_study_id or "").strip()

    has_study_identifiers = (
        _is_valid_metadata_value(s1_id)
        and _is_valid_metadata_value(s2_id)
        and _is_valid_metadata_value(s1_fam)
        and _is_valid_metadata_value(s2_fam)
        and _is_valid_metadata_value(s1_ind)
        and _is_valid_metadata_value(s2_ind)
    )

    is_same_study = bool(
        s1_id.lower() == s2_id.lower()
        or s1_fam.lower() == s2_fam.lower()
        or s1_ind.lower() == s2_ind.lower()
    )

    is_independent_studies = has_study_identifiers and not is_same_study

    if not has_study_identifiers:
        differences.append(
            f"Study identifiers missing or unverified: '{record1.study_id}' vs '{record2.study_id}'"
        )
    elif is_same_study:
        differences.append(
            f"Records belong to the same study or study family: '{s1_fam}' "
            "(different cell lines, donors, or conditions within the same study are not independent studies)"
        )

    # 2. Compound Identity Match
    c1_id = (record1.compound_id or "").strip()
    c2_id = (record2.compound_id or "").strip()
    c1_name = (record1.compound_name or "").strip().lower()
    c2_name = (record2.compound_name or "").strip().lower()

    has_compound_meta = (
        (_is_valid_metadata_value(c1_id) or _is_valid_metadata_value(c1_name))
        and (_is_valid_metadata_value(c2_id) or _is_valid_metadata_value(c2_name))
    )
    is_compound_matched = bool(
        has_compound_meta
        and (
            (c1_id and c2_id and c1_id.lower() == c2_id.lower())
            or (c1_name and c2_name and c1_name == c2_name)
        )
    )
    if not has_compound_meta:
        differences.append("Compound metadata missing, unknown, or unspecified")
    elif not is_compound_matched:
        differences.append(
            f"Compound mismatch: '{record1.compound_name}' ({record1.compound_id}) vs '{record2.compound_name}' ({record2.compound_id})"
        )

    # 3. Dose Match
    d1 = record1.dose_value
    d2 = record2.dose_value
    u1 = (record1.dose_unit or "").strip().lower()
    u2 = (record2.dose_unit or "").strip().lower()
    has_dose_meta = (
        d1 is not None
        and d2 is not None
        and _is_valid_metadata_value(u1)
        and _is_valid_metadata_value(u2)
    )
    is_dose_matched = bool(
        has_dose_meta
        and abs(d1 - d2) < 1e-6
        and u1 == u2
    )
    if not has_dose_meta:
        differences.append("Dose or concentration metadata missing, unknown, or unspecified")
    elif not is_dose_matched:
        differences.append(f"Dose mismatch: {d1} {record1.dose_unit} vs {d2} {record2.dose_unit}")

    # 4. Developmental Stage Match
    st1 = (record1.developmental_stage or "").strip().lower()
    st2 = (record2.developmental_stage or "").strip().lower()
    has_stage_meta = _is_valid_metadata_value(st1) and _is_valid_metadata_value(st2)
    is_developmental_stage_matched = bool(has_stage_meta and st1 == st2)
    if not has_stage_meta:
        differences.append("Developmental stage metadata missing, unknown, or unspecified")
    elif not is_developmental_stage_matched:
        differences.append(
            f"Developmental stage mismatch: '{record1.developmental_stage}' vs '{record2.developmental_stage}'"
        )

    # 5. Exposure Duration Match
    dur1 = record1.exposure_duration
    dur2 = record2.exposure_duration
    du1 = (record1.exposure_duration_unit or "").strip().lower()
    du2 = (record2.exposure_duration_unit or "").strip().lower()
    has_dur_meta = (
        dur1 is not None
        and dur2 is not None
        and _is_valid_metadata_value(du1)
        and _is_valid_metadata_value(du2)
    )
    is_duration_matched = bool(
        has_dur_meta
        and abs(dur1 - dur2) < 1e-6
        and du1 == du2
    )
    if not has_dur_meta:
        differences.append("Exposure duration metadata missing, unknown, or unspecified")
    elif not is_duration_matched:
        differences.append(
            f"Exposure duration mismatch: {dur1} {record1.exposure_duration_unit} vs {dur2} {record2.exposure_duration_unit}"
        )

    # 6. Organoid Model Match
    m1 = (record1.model_type or "").strip().lower()
    m2 = (record2.model_type or "").strip().lower()
    has_model_meta = _is_valid_metadata_value(m1) and _is_valid_metadata_value(m2)
    is_organoid_model_matched = bool(has_model_meta and m1 == m2)
    if not has_model_meta:
        differences.append("Organoid model metadata missing, unknown, or unspecified")
    elif not is_organoid_model_matched:
        differences.append(f"Organoid model mismatch: '{record1.model_type}' vs '{record2.model_type}'")

    # 7. Biological Replicates Verification
    rep_adequacy1 = record1.treatment_replicate_count >= 2 and record1.control_replicate_count >= 2
    rep_adequacy2 = record2.treatment_replicate_count >= 2 and record2.control_replicate_count >= 2
    study1_reps_verified = bool(replicates_verified1 is True)
    study2_reps_verified = bool(replicates_verified2 is True)
    are_replicates_verified = study1_reps_verified and study2_reps_verified and rep_adequacy1 and rep_adequacy2
    if not are_replicates_verified:
        diff_str: list[str] = []
        if not study1_reps_verified:
            diff_str.append(f"Study 1 ({record1.study_id}) replicates not explicitly verified")
        if not study2_reps_verified:
            diff_str.append(f"Study 2 ({record2.study_id}) replicates not explicitly verified")
        if not rep_adequacy1:
            diff_str.append(f"Study 1 under-replicated (trt={record1.treatment_replicate_count}, ctrl={record1.control_replicate_count})")
        if not rep_adequacy2:
            diff_str.append(f"Study 2 under-replicated (trt={record2.treatment_replicate_count}, ctrl={record2.control_replicate_count})")
        differences.append("Replicate deficiency: " + "; ".join(diff_str))

    # Synthesis Decision
    is_synthesis_justified = (
        is_independent_studies
        and is_compound_matched
        and is_dose_matched
        and is_developmental_stage_matched
        and is_duration_matched
        and is_organoid_model_matched
        and are_replicates_verified
    )

    if not is_compound_matched:
        status = CrossStudyStatus.NONCOMPARABLE_CONTEXT
        rationale = "Formal statistical synthesis rejected: distinct chemical compounds cannot be pooled."
    elif not is_independent_studies:
        status = CrossStudyStatus.SINGLE_STUDY_ONLY
        rationale = (
            "Formal cross-study synthesis rejected: compared records do not come from independent studies or study families. "
            f"Both records belong to study family '{record1.study_family_id}' (different cell lines, donors, or conditions within the same study are not independent studies)."
        )
    elif not is_synthesis_justified:
        status = CrossStudyStatus.NONCOMPARABLE_CONTEXT
        rationale = (
            f"Formal statistical synthesis rejected: experimental conditions differ ({'; '.join(differences)}). "
            "Descriptive cross-study alignment is permitted, but combining effect sizes into a pooled estimate "
            "or unweighted averaging is strictly prohibited."
        )
    else:
        status = CrossStudyStatus.COMPARABLE_INDEPENDENT_STUDIES
        rationale = (
            "Independent studies share verified comparable experimental regimens, developmental stages, "
            "organoid model systems, and verified biological replication. Formally eligible for prespecified meta-analysis."
        )

    return StudyComparabilityAssessment(
        study1_id=record1.study_id,
        study2_id=record2.study_id,
        contrast1_id=record1.contrast_id,
        contrast2_id=record2.contrast_id,
        compound_id=record1.compound_id,
        compound_name=record1.compound_name,
        is_independent_studies=is_independent_studies,
        is_compound_matched=is_compound_matched,
        is_dose_matched=is_dose_matched,
        study1_dose=record1.dose_value,
        study2_dose=record2.dose_value,
        dose_unit=record1.dose_unit,
        is_developmental_stage_matched=is_developmental_stage_matched,
        study1_stage=record1.developmental_stage,
        study2_stage=record2.developmental_stage,
        is_duration_matched=is_duration_matched,
        study1_duration=record1.exposure_duration,
        study2_duration=record2.exposure_duration,
        is_organoid_model_matched=is_organoid_model_matched,
        study1_model=record1.model_type,
        study2_model=record2.model_type,
        are_replicates_verified=are_replicates_verified,
        study1_replicates_verified=study1_reps_verified,
        study2_replicates_verified=study2_reps_verified,
        is_synthesis_justified=is_synthesis_justified,
        cross_study_status=status,
        rationale=rationale,
        differences=tuple(differences),
    )


def compare_cross_study_gene_responses(
    registry1: ConditionGeneResponseRegistry,
    registry2: ConditionGeneResponseRegistry,
    contrast1_id: str | None = None,
    contrast2_id: str | None = None,
    replicates_verified1: bool | None = None,
    replicates_verified2: bool | None = None,
    comparison_id: str | None = None,
) -> CrossStudyComparisonDataset:
    """Compare and align gene responses across two independent studies.

    Preserves 100% of underlying condition records.
    When experimental conditions conflict, descriptive alignment is performed but
    statistical pooling (pooled log2FC, SE) is strictly prohibited.
    """
    recs1 = registry1.get_records_for_contrast(contrast1_id) if contrast1_id else registry1.records
    recs2 = registry2.get_records_for_contrast(contrast2_id) if contrast2_id else registry2.records

    if not recs1 or not recs2:
        raise ValueError("Both registries must contain records for cross-study comparison.")

    rep1 = recs1[0]
    rep2 = recs2[0]
    assessment = evaluate_study_comparability(
        record1=rep1,
        record2=rep2,
        replicates_verified1=replicates_verified1,
        replicates_verified2=replicates_verified2,
    )

    map1 = {r.canonical_ensembl_gene_id: r for r in recs1}
    map2 = {r.canonical_ensembl_gene_id: r for r in recs2}

    all_gene_ids = sorted(list(set(map1.keys()) | set(map2.keys())))
    gene_comparisons: list[GeneCrossStudyComparisonRecord] = []

    for gid in all_gene_ids:
        r1 = map1.get(gid)
        r2 = map2.get(gid)

        symbol = r1.gene_symbol if (r1 and r1.gene_symbol) else (r2.gene_symbol if r2 else None)

        if r1 is not None and r2 is not None:
            concordance = derive_direction_concordance(r1.response_direction, r2.response_direction)
            is_conc: bool | None = None
            if concordance in (DirectionConcordance.CONCORDANT_UP, DirectionConcordance.CONCORDANT_DOWN):
                is_conc = True
            elif concordance == DirectionConcordance.DISCORDANT:
                is_conc = False

            rec = GeneCrossStudyComparisonRecord(
                canonical_ensembl_gene_id=gid,
                gene_symbol=symbol,
                study1_id=rep1.study_id,
                study1_contrast_id=rep1.contrast_id,
                study1_log2_fold_change=r1.log2_fold_change,
                study1_log2_fold_change_se=r1.log2_fold_change_se,
                study1_wald_statistic=r1.wald_statistic,
                study1_p_value=r1.p_value,
                study1_adjusted_p_value=r1.adjusted_p_value,
                study1_response_direction=r1.response_direction,
                study1_estimate_status=r1.estimate_status,
                study2_id=rep2.study_id,
                study2_contrast_id=rep2.contrast_id,
                study2_log2_fold_change=r2.log2_fold_change,
                study2_log2_fold_change_se=r2.log2_fold_change_se,
                study2_wald_statistic=r2.wald_statistic,
                study2_p_value=r2.p_value,
                study2_adjusted_p_value=r2.adjusted_p_value,
                study2_response_direction=r2.response_direction,
                study2_estimate_status=r2.estimate_status,
                direction_concordance=concordance,
                is_concordant=is_conc,
                cross_study_status=assessment.cross_study_status,
                is_synthesis_justified=assessment.is_synthesis_justified,
                pooled_log2_fold_change=None,
                pooled_standard_error=None,
            )
        elif r1 is not None:
            rec = GeneCrossStudyComparisonRecord(
                canonical_ensembl_gene_id=gid,
                gene_symbol=symbol,
                study1_id=rep1.study_id,
                study1_contrast_id=rep1.contrast_id,
                study1_log2_fold_change=r1.log2_fold_change,
                study1_log2_fold_change_se=r1.log2_fold_change_se,
                study1_wald_statistic=r1.wald_statistic,
                study1_p_value=r1.p_value,
                study1_adjusted_p_value=r1.adjusted_p_value,
                study1_response_direction=r1.response_direction,
                study1_estimate_status=r1.estimate_status,
                study2_id=rep2.study_id,
                study2_contrast_id=rep2.contrast_id,
                study2_log2_fold_change=None,
                study2_log2_fold_change_se=None,
                study2_wald_statistic=None,
                study2_p_value=None,
                study2_adjusted_p_value=None,
                study2_response_direction=ResponseDirection.UNDETERMINED,
                study2_estimate_status=EstimateStatus.INVALID_STATISTIC,
                direction_concordance=DirectionConcordance.UNDETERMINED,
                is_concordant=None,
                cross_study_status=assessment.cross_study_status,
                is_synthesis_justified=assessment.is_synthesis_justified,
                pooled_log2_fold_change=None,
                pooled_standard_error=None,
            )
        else:
            assert r2 is not None
            rec = GeneCrossStudyComparisonRecord(
                canonical_ensembl_gene_id=gid,
                gene_symbol=symbol,
                study1_id=rep1.study_id,
                study1_contrast_id=rep1.contrast_id,
                study1_log2_fold_change=None,
                study1_log2_fold_change_se=None,
                study1_wald_statistic=None,
                study1_p_value=None,
                study1_adjusted_p_value=None,
                study1_response_direction=ResponseDirection.UNDETERMINED,
                study1_estimate_status=EstimateStatus.INVALID_STATISTIC,
                study2_id=rep2.study_id,
                study2_contrast_id=rep2.contrast_id,
                study2_log2_fold_change=r2.log2_fold_change,
                study2_log2_fold_change_se=r2.log2_fold_change_se,
                study2_wald_statistic=r2.wald_statistic,
                study2_p_value=r2.p_value,
                study2_adjusted_p_value=r2.adjusted_p_value,
                study2_response_direction=r2.response_direction,
                study2_estimate_status=r2.estimate_status,
                direction_concordance=DirectionConcordance.UNDETERMINED,
                is_concordant=None,
                cross_study_status=assessment.cross_study_status,
                is_synthesis_justified=assessment.is_synthesis_justified,
                pooled_log2_fold_change=None,
                pooled_standard_error=None,
            )

        gene_comparisons.append(rec)

    comp_id = comparison_id or f"cmp_{rep1.study_id}_{rep2.study_id}_{rep1.compound_id}"
    return CrossStudyComparisonDataset(
        comparison_id=comp_id,
        assessment=assessment,
        gene_comparisons=tuple(gene_comparisons),
        created_at=datetime.now(timezone.utc).isoformat(),
        provenance={
            "study1_id": rep1.study_id,
            "study2_id": rep2.study_id,
            "contrast1_id": rep1.contrast_id,
            "contrast2_id": rep2.contrast_id,
            "total_genes": len(gene_comparisons),
        },
    )
