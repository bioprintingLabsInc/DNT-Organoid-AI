"""Scientific contrast construction and validation engine for Treatment-versus-Matched-Control Contrast Builder v1.

Enforces:
1. Contrasts strictly originate from explicit treatment_control_relationships.
2. Separate contrast per treatment condition (no automated pooling of treatment conditions).
3. Preservation of study and experiment boundaries (strictly no cross-study or cross-experiment pooling).
4. Developmental age strict matching (blocking on collection age mismatches, review on ambiguity).
5. Organoid context and biological source compatibility.
6. True biological replication requirements for inferential downstream DE (technical replicates do not inflate counts).
7. Preservation and non-silent handling of sample comparison exceptions.
8. Complete independence from DNT labels/evidence.
9. Zero gene filtering.
10. Fully deterministic contrast IDs and output sorting.
"""

from collections import defaultdict
from collections.abc import Iterable
import hashlib
import json
import re
from typing import Any

from src.ingestion.registry import reported_value
from src.normalization.models import (
    NormalizedCohortResult,
    NormalizedDataset,
    Status as NormalizationStatus,
)
from .models import (
    ContrastStatus,
    Finding,
    MolecularResponseContrast,
    ResponseContrastDataset,
    Severity,
    Status,
    ordered_findings,
    status_for,
)

BUILDER_VERSION = "1.0.0"


def extract_metadata_value(val_rec: Any) -> Any:
    """Extract reported or normalized value from scalar or value_record dict."""
    if val_rec is None:
        return None
    if isinstance(val_rec, (str, int, float, bool)):
        return val_rec
    if isinstance(val_rec, dict):
        rep = reported_value(val_rec)
        if rep is not None:
            return rep
        norm = val_rec.get("normalized_value")
        if norm is not None:
            return norm
        return val_rec.get("original_value")
    return None


def normalize_age_string(val: Any) -> str | None:
    """Normalize a developmental age or stage string to a canonical token for strict comparison.

    Maps e.g. 'Day 30', 'day 30', 'd30', '30 days', 30 -> 'day_30'.
    Returns None if missing, unknown, or empty.
    """
    if val is None:
        return None
    s = str(val).strip().lower()
    if not s or s in ("unknown", "missing", "missing_not_reported", "none", "null", "not_reported"):
        return None

    clean = re.sub(r"[\s\-_]+", " ", s)

    # e.g., 'day 30' or 'd 30'
    m_day = re.match(r"^(?:day|d)\s*([0-9]+(?:\.[0-9]+)?)$", clean)
    if m_day:
        return f"day_{m_day.group(1)}"

    # e.g., '30 days' or '30 day' or '30 d'
    m_days = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*(?:days|day|d)$", clean)
    if m_days:
        return f"day_{m_days.group(1)}"

    # e.g., plain number '30'
    m_num = re.match(r"^([0-9]+(?:\.[0-9]+)?)$", clean)
    if m_num:
        return f"day_{m_num.group(1)}"

    # e.g., 'gw 10' or 'gw_10'
    m_gw = re.match(r"^(?:gw|gestational week)\s*([0-9]+(?:\.[0-9]+)?)$", clean)
    if m_gw:
        return f"gw_{m_gw.group(1)}"

    return clean.replace(" ", "_")


def sanitize_id_token(token: str) -> str:
    """Sanitize identifier token to safe alphanumeric and underscore characters."""
    return re.sub(r"[^A-Za-z0-9_]+", "_", token).strip("_")


class ResponseContrastBuilder:
    """Orchestrates validation, control matching, and molecular-response contrast generation."""

    def __init__(self, metadata: dict[str, Any], normalized: NormalizedDataset) -> None:
        self.metadata = metadata
        self.normalized = normalized
        self._index_data()

    def _index_data(self) -> None:
        # 1. Index metadata collections
        self.samples_by_id = {
            s["sample_id"]: s
            for s in self.metadata.get("samples", [])
            if isinstance(s, dict) and "sample_id" in s
        }
        self.conditions_by_id = {
            c["condition_id"]: c
            for c in self.metadata.get("conditions", [])
            if isinstance(c, dict) and "condition_id" in c
        }
        self.experiments_by_id = {
            e["experiment_id"]: e
            for e in self.metadata.get("experiments", [])
            if isinstance(e, dict) and "experiment_id" in e
        }
        self.studies_by_id = {
            s["study_id"]: s
            for s in self.metadata.get("studies", [])
            if isinstance(s, dict) and "study_id" in s
        }
        self.organoid_contexts_by_id = {
            o["organoid_context_id"]: o
            for o in self.metadata.get("organoid_contexts", [])
            if isinstance(o, dict) and "organoid_context_id" in o
        }
        self.biological_sources_by_id = {
            b["biological_source_id"]: b
            for b in self.metadata.get("biological_sources", [])
            if isinstance(b, dict) and "biological_source_id" in b
        }

        # Map condition_id -> list of sample_ids
        self.samples_by_condition: dict[str, list[str]] = defaultdict(list)
        for sid, s in self.samples_by_id.items():
            cid = s.get("condition_id")
            if cid:
                self.samples_by_condition[cid].append(sid)

        # Index exposures by condition_id (prefer condition-level exposure)
        self.exposures_by_condition: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for exp in self.metadata.get("exposures", []):
            if isinstance(exp, dict) and "condition_id" in exp:
                self.exposures_by_condition[exp["condition_id"]].append(exp)

        # Index sample comparison exceptions by relationship_id
        self.exceptions_by_relationship: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for exc in self.metadata.get("sample_comparison_exceptions", []):
            if isinstance(exc, dict) and "relationship_id" in exc:
                self.exceptions_by_relationship[exc["relationship_id"]].append(exc)

        # 2. Index NormalizedDataset
        # Map sample_id -> NormalizedCohortResult
        self.cohort_result_by_sample: dict[str, NormalizedCohortResult] = {}
        for cr in self.normalized.cohort_results:
            for sid in cr.cohort.sample_ids:
                self.cohort_result_by_sample[sid] = cr

    def build_dataset(self) -> ResponseContrastDataset:
        """Evaluate all declared treatment-control relationships and construct contrasts."""
        dataset_findings: list[Finding] = []
        contrasts: list[MolecularResponseContrast] = []

        declared_relationships = [
            r
            for r in self.metadata.get("treatment_control_relationships", [])
            if isinstance(r, dict) and "relationship_id" in r
        ]

        if not declared_relationships:
            dataset_findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="no_treatment_control_relationships",
                    entity_type="metadata",
                    entity_id=None,
                    path="treatment_control_relationships",
                    message="Metadata contains no explicit treatment_control_relationships. No contrasts can be constructed.",
                )
            )

        for rel in declared_relationships:
            rel_contrasts, rel_findings = self._evaluate_relationship(rel)
            contrasts.extend(rel_contrasts)
            dataset_findings.extend(rel_findings)

        # Sort contrasts deterministically
        sorted_contrasts = tuple(
            sorted(
                contrasts,
                key=lambda c: (c.study_id, c.experiment_id, c.cohort_id, c.treatment_condition_id, c.contrast_id),
            )
        )

        all_findings = list(dataset_findings)
        for c in sorted_contrasts:
            all_findings.extend(c.findings)

        sorted_dataset_findings = ordered_findings(dataset_findings)
        overall_status = status_for(all_findings)

        metadata_serialized = json.dumps(self.metadata, sort_keys=True, default=str)
        metadata_hash = hashlib.sha256(metadata_serialized.encode("utf-8")).hexdigest()

        return ResponseContrastDataset(
            status=overall_status,
            findings=sorted_dataset_findings,
            contrasts=sorted_contrasts,
            source_asset_id=self.normalized.source_asset_id,
            reference_identity=self.normalized.reference_identity,
            builder_version=BUILDER_VERSION,
            metadata_hash=metadata_hash,
        )

    def _evaluate_relationship(
        self,
        rel: dict[str, Any],
    ) -> tuple[list[MolecularResponseContrast], list[Finding]]:
        """Evaluate an explicit treatment_control_relationship record into condition-level contrasts."""
        rel_id = str(rel["relationship_id"])
        rel_findings: list[Finding] = []
        contrasts: list[MolecularResponseContrast] = []

        treatment_cond_ids = tuple(str(c) for c in rel.get("treatment_condition_ids", ()) if str(c).strip())
        control_cond_ids = tuple(str(c) for c in rel.get("matched_control_condition_ids", ()) if str(c).strip())

        if not treatment_cond_ids:
            rel_findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="relationship_missing_treatment_conditions",
                    entity_type="treatment_control_relationship",
                    entity_id=rel_id,
                    path=f"treatment_control_relationships.{rel_id}.treatment_condition_ids",
                    message=f"Relationship '{rel_id}' declares zero treatment conditions.",
                )
            )
            return [], rel_findings

        # Scientific Invariant: Each treatment condition becomes its OWN molecular-response contrast
        for treat_cid in sorted(set(treatment_cond_ids)):
            contrast = self._build_single_contrast(rel_id, treat_cid, control_cond_ids)
            contrasts.append(contrast)

        return contrasts, rel_findings

    def _build_single_contrast(
        self,
        rel_id: str,
        treat_cid: str,
        control_cond_ids: tuple[str, ...],
    ) -> MolecularResponseContrast:
        """Construct and validate a single molecular-response contrast for one treatment condition."""
        findings: list[Finding] = []
        status = ContrastStatus.ELIGIBLE

        treat_cond = self.conditions_by_id.get(treat_cid)
        if not treat_cond:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="treatment_condition_not_found",
                    entity_type="condition",
                    entity_id=treat_cid,
                    path=f"conditions.{treat_cid}",
                    message=f"Treatment condition '{treat_cid}' referenced in relationship '{rel_id}' does not exist in metadata.",
                )
            )
            return self._build_empty_contrast(
                rel_id=rel_id,
                treat_cid=treat_cid,
                control_cond_ids=control_cond_ids,
                status=ContrastStatus.BLOCKED,
                findings=findings,
            )

        treat_status = extract_metadata_value(treat_cond.get("treatment_control_status"))
        if treat_status != "treatment":
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="treatment_control_role_inconsistency",
                    entity_type="condition",
                    entity_id=treat_cid,
                    path=f"conditions.{treat_cid}.treatment_control_status",
                    message=(
                        f"Condition '{treat_cid}' is declared as a treatment arm in relationship '{rel_id}', "
                        f"but its treatment_control_status is '{treat_status}' (expected 'treatment')."
                    ),
                )
            )
            status = ContrastStatus.BLOCKED

        treat_exp_id = treat_cond.get("experiment_id") or "unknown_exp"
        treat_exp = self.experiments_by_id.get(treat_exp_id, {})
        treat_study_id = treat_exp.get("study_id") or "unknown_study"

        # 1. Matched Control Condition Verification
        if not control_cond_ids:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="missing_matched_control",
                    entity_type="treatment_control_relationship",
                    entity_id=rel_id,
                    path=f"treatment_control_relationships.{rel_id}.matched_control_condition_ids",
                    message=f"Relationship '{rel_id}' has no matched control conditions declared for treatment '{treat_cid}'.",
                )
            )
            status = ContrastStatus.BLOCKED

        valid_control_conds: list[dict[str, Any]] = []
        for ctrl_cid in control_cond_ids:
            ctrl_cond = self.conditions_by_id.get(ctrl_cid)
            if not ctrl_cond:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="control_condition_not_found",
                        entity_type="condition",
                        entity_id=ctrl_cid,
                        path=f"conditions.{ctrl_cid}",
                        message=f"Control condition '{ctrl_cid}' referenced in relationship '{rel_id}' does not exist in metadata.",
                    )
                )
                status = ContrastStatus.BLOCKED
                continue

            ctrl_status = extract_metadata_value(ctrl_cond.get("treatment_control_status"))
            if ctrl_status != "control":
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="treatment_control_role_inconsistency",
                        entity_type="condition",
                        entity_id=ctrl_cid,
                        path=f"conditions.{ctrl_cid}.treatment_control_status",
                        message=(
                            f"Condition '{ctrl_cid}' is declared as matched control in relationship '{rel_id}', "
                            f"but its treatment_control_status is '{ctrl_status}' (expected 'control')."
                        ),
                    )
                )
                status = ContrastStatus.BLOCKED

            ctrl_exp_id = ctrl_cond.get("experiment_id") or "unknown_exp"
            ctrl_exp = self.experiments_by_id.get(ctrl_exp_id, {})
            ctrl_study_id = ctrl_exp.get("study_id") or "unknown_study"

            # 2. Study and Experiment Boundary Enforcement
            if ctrl_exp_id != treat_exp_id:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="cross_experiment_comparison_blocked",
                        entity_type="contrast",
                        entity_id=f"{treat_cid}_vs_{ctrl_cid}",
                        path=f"conditions.{ctrl_cid}.experiment_id",
                        message=(
                            f"Cross-experiment comparison strictly blocked: treatment '{treat_cid}' belongs to "
                            f"experiment '{treat_exp_id}', while control '{ctrl_cid}' belongs to experiment '{ctrl_exp_id}'."
                        ),
                    )
                )
                status = ContrastStatus.BLOCKED

            if ctrl_study_id != treat_study_id:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="cross_study_comparison_blocked",
                        entity_type="contrast",
                        entity_id=f"{treat_cid}_vs_{ctrl_cid}",
                        path=f"conditions.{ctrl_cid}.study_id",
                        message=(
                            f"Cross-study comparison strictly blocked: treatment '{treat_cid}' belongs to "
                            f"study '{treat_study_id}', while control '{ctrl_cid}' belongs to study '{ctrl_study_id}'."
                        ),
                    )
                )
                status = ContrastStatus.BLOCKED

            valid_control_conds.append(ctrl_cond)

        # 3. Sample Resolution
        treat_sample_ids = sorted(self.samples_by_condition.get(treat_cid, []))
        if not treat_sample_ids:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="missing_treatment_samples",
                    entity_type="condition",
                    entity_id=treat_cid,
                    path=f"samples",
                    message=f"Treatment condition '{treat_cid}' has zero samples in metadata.",
                )
            )
            status = ContrastStatus.BLOCKED

        raw_control_sample_ids: list[str] = []
        for ctrl_cid in control_cond_ids:
            raw_control_sample_ids.extend(self.samples_by_condition.get(ctrl_cid, []))
        raw_control_sample_ids = sorted(set(raw_control_sample_ids))

        if not raw_control_sample_ids and control_cond_ids:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="missing_matched_control_samples",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path=f"samples",
                    message=f"Matched control condition(s) {control_cond_ids} have zero samples in metadata.",
                )
            )
            status = ContrastStatus.BLOCKED

        # 4. Sample Comparison Exceptions Evaluation
        applied_exception_ids: list[str] = []
        active_control_sample_ids = list(raw_control_sample_ids)
        rel_exceptions = self.exceptions_by_relationship.get(rel_id, [])

        for exc in rel_exceptions:
            exc_id = str(exc.get("exception_id", "unnamed_exception"))
            affected_sids = set(exc.get("affected_sample_ids", []))
            all_contrast_sids = set(treat_sample_ids) | set(active_control_sample_ids)

            if not affected_sids.isdisjoint(all_contrast_sids):
                applied_exception_ids.append(exc_id)
                replacements = exc.get("replacement_control_sample_ids")

                if replacements and isinstance(replacements, list):
                    # Validate replacement control samples
                    valid_replacements: list[str] = []
                    for rep_sid in replacements:
                        rep_s = self.samples_by_id.get(rep_sid)
                        if rep_s and rep_s.get("condition_id") in control_cond_ids:
                            valid_replacements.append(rep_sid)

                    if valid_replacements:
                        # Replace affected control samples
                        active_control_sample_ids = [
                            sid for sid in active_control_sample_ids if sid not in affected_sids
                        ]
                        active_control_sample_ids.extend(valid_replacements)
                        active_control_sample_ids = sorted(set(active_control_sample_ids))
                        findings.append(
                            Finding(
                                severity=Severity.INFO,
                                rule_id="sample_comparison_exception_applied",
                                entity_type="sample_comparison_exception",
                                entity_id=exc_id,
                                path=f"sample_comparison_exceptions.{exc_id}",
                                message=(
                                    f"Applied sample comparison exception '{exc_id}': replaced affected control sample(s) "
                                    f"{sorted(affected_sids)} with replacement control sample(s) {sorted(valid_replacements)}."
                                ),
                            )
                        )
                    else:
                        findings.append(
                            Finding(
                                severity=Severity.REVIEW,
                                rule_id="unresolved_sample_comparison_exception",
                                entity_type="sample_comparison_exception",
                                entity_id=exc_id,
                                path=f"sample_comparison_exceptions.{exc_id}",
                                message=(
                                    f"Sample comparison exception '{exc_id}' affects samples {sorted(affected_sids)}, "
                                    "but supplied replacement controls are invalid or absent."
                                ),
                            )
                        )
                        if status != ContrastStatus.BLOCKED:
                            status = ContrastStatus.NEEDS_REVIEW
                else:
                    findings.append(
                        Finding(
                            severity=Severity.REVIEW,
                            rule_id="unresolved_sample_comparison_exception",
                            entity_type="sample_comparison_exception",
                            entity_id=exc_id,
                            path=f"sample_comparison_exceptions.{exc_id}",
                            message=(
                                f"Sample comparison exception '{exc_id}' affects contrast samples {sorted(affected_sids)} "
                                "without declared replacement controls. Comparison requires expert review."
                            ),
                        )
                    )
                    if status != ContrastStatus.BLOCKED:
                        status = ContrastStatus.NEEDS_REVIEW

        treat_samples = [self.samples_by_id[sid] for sid in treat_sample_ids if sid in self.samples_by_id]
        control_samples = [self.samples_by_id[sid] for sid in active_control_sample_ids if sid in self.samples_by_id]

        # 5. Normalization Cohort Linkage and Integrity
        cohort_ids: set[str] = set()
        for sid in treat_sample_ids + active_control_sample_ids:
            cr = self.cohort_result_by_sample.get(sid)
            if cr:
                cohort_ids.add(cr.cohort.cohort_id)
                if cr.status == NormalizationStatus.FAIL:
                    findings.append(
                        Finding(
                            severity=Severity.ERROR,
                            rule_id="upstream_normalization_failure",
                            entity_type="cohort",
                            entity_id=cr.cohort.cohort_id,
                            path=f"cohort_results.{cr.cohort.cohort_id}",
                            message=f"Sample '{sid}' belongs to cohort '{cr.cohort.cohort_id}' which failed normalization.",
                        )
                    )
                    status = ContrastStatus.BLOCKED
                elif cr.status == NormalizationStatus.NEEDS_REVIEW:
                    findings.append(
                        Finding(
                            severity=Severity.REVIEW,
                            rule_id="upstream_normalization_needs_review",
                            entity_type="cohort",
                            entity_id=cr.cohort.cohort_id,
                            path=f"cohort_results.{cr.cohort.cohort_id}",
                            message=f"Sample '{sid}' belongs to cohort '{cr.cohort.cohort_id}' which is at NEEDS_REVIEW in normalization.",
                        )
                    )
                    if status != ContrastStatus.BLOCKED:
                        status = ContrastStatus.NEEDS_REVIEW
            else:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="sample_missing_from_normalization",
                        entity_type="sample",
                        entity_id=sid,
                        path=f"samples.{sid}",
                        message=f"Sample '{sid}' is not present in any cohort of the normalized dataset.",
                    )
                )
                status = ContrastStatus.BLOCKED

        if len(cohort_ids) == 1:
            cohort_id = next(iter(cohort_ids))
        elif len(cohort_ids) > 1:
            cohort_id = "_vs_".join(sorted(cohort_ids))
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="cross_cohort_comparison_blocked",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="cohort_id",
                    message=(
                        f"Contrast samples span multiple normalization cohorts: {sorted(cohort_ids)}. "
                        "Contrasts must originate within a single coherent normalization cohort."
                    ),
                )
            )
            status = ContrastStatus.BLOCKED
        else:
            cohort_id = "unknown_cohort"

        # 6. Developmental Age / Stage Matching Invariants
        treat_raw_ages = [extract_metadata_value(s.get("developmental_age_or_stage_at_collection")) for s in treat_samples]
        control_raw_ages = [extract_metadata_value(s.get("developmental_age_or_stage_at_collection")) for s in control_samples]

        treat_norm_ages = {normalize_age_string(a) for a in treat_raw_ages}
        control_norm_ages = {normalize_age_string(a) for a in control_raw_ages}

        treat_age_rep = treat_raw_ages[0] if len(set(treat_raw_ages)) == 1 and treat_raw_ages[0] is not None else (
            ", ".join(sorted({str(a) for a in treat_raw_ages if a is not None})) or None
        )
        control_age_rep = control_raw_ages[0] if len(set(control_raw_ages)) == 1 and control_raw_ages[0] is not None else (
            ", ".join(sorted({str(a) for a in control_raw_ages if a is not None})) or None
        )

        treat_age_ambiguous = None in treat_norm_ages or len(treat_norm_ages) != 1
        control_age_ambiguous = None in control_norm_ages or len(control_norm_ages) != 1

        if treat_samples and control_samples:
            if treat_age_ambiguous or control_age_ambiguous:
                findings.append(
                    Finding(
                        severity=Severity.REVIEW,
                        rule_id="developmental_age_ambiguous",
                        entity_type="contrast",
                        entity_id=treat_cid,
                        path="samples.developmental_age_or_stage_at_collection",
                        message=(
                            f"Developmental age or stage at collection is ambiguous, missing, or divergent across samples "
                            f"(treatment: {treat_norm_ages}, control: {control_norm_ages}). "
                            "Do not guess or silently reconcile; requires expert review."
                        ),
                    )
                )
                if status != ContrastStatus.BLOCKED:
                    status = ContrastStatus.NEEDS_REVIEW
            else:
                t_single_age = next(iter(treat_norm_ages))
                c_single_age = next(iter(control_norm_ages))
                if t_single_age != c_single_age:
                    findings.append(
                        Finding(
                            severity=Severity.ERROR,
                            rule_id="developmental_age_mismatch",
                            entity_type="contrast",
                            entity_id=treat_cid,
                            path="samples.developmental_age_or_stage_at_collection",
                            message=(
                                f"Developmental age mismatch: treatment collected at '{treat_age_rep}' (normalized '{t_single_age}'), "
                                f"while control collected at '{control_age_rep}' (normalized '{c_single_age}'). "
                                "Comparing across different developmental ages is strictly invalid."
                            ),
                        )
                    )
                    status = ContrastStatus.BLOCKED

        # 7. Organoid Context and Biological Source Compatibility
        treat_org_contexts = tuple(sorted({
            str(extract_metadata_value(s.get("organoid_context_id")))
            for s in treat_samples
            if extract_metadata_value(s.get("organoid_context_id"))
        }))
        control_org_contexts = tuple(sorted({
            str(extract_metadata_value(s.get("organoid_context_id")))
            for s in control_samples
            if extract_metadata_value(s.get("organoid_context_id"))
        }))
        all_org_contexts = tuple(sorted(set(treat_org_contexts) | set(control_org_contexts)))

        if treat_org_contexts and control_org_contexts:
            if set(treat_org_contexts).isdisjoint(set(control_org_contexts)):
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="organoid_context_mismatch",
                        entity_type="contrast",
                        entity_id=treat_cid,
                        path="samples.organoid_context_id",
                        message=(
                            f"Incompatible organoid contexts: treatment has {treat_org_contexts} "
                            f"while control has {control_org_contexts}."
                        ),
                    )
                )
                status = ContrastStatus.BLOCKED

        treat_bio_sources = tuple(sorted({
            str(extract_metadata_value(s.get("biological_source_id")))
            for s in treat_samples
            if extract_metadata_value(s.get("biological_source_id"))
        }))
        control_bio_sources = tuple(sorted({
            str(extract_metadata_value(s.get("biological_source_id")))
            for s in control_samples
            if extract_metadata_value(s.get("biological_source_id"))
        }))
        all_bio_sources = tuple(sorted(set(treat_bio_sources) | set(control_bio_sources)))

        # 8. Biological Replicate Counting (Technical replicates do NOT inflate counts)
        treat_bio_reps = tuple(sorted({
            str(extract_metadata_value(s.get("biological_replicate_id")))
            for s in treat_samples
            if extract_metadata_value(s.get("biological_replicate_id")) is not None
        }))
        control_bio_reps = tuple(sorted({
            str(extract_metadata_value(s.get("biological_replicate_id")))
            for s in control_samples
            if extract_metadata_value(s.get("biological_replicate_id")) is not None
        }))

        treat_rep_count = len(set(treat_bio_reps))
        control_rep_count = len(set(control_bio_reps))

        if treat_samples and control_samples:
            if treat_rep_count < 2 or control_rep_count < 2:
                findings.append(
                    Finding(
                        severity=Severity.REVIEW,
                        rule_id="insufficient_biological_replication",
                        entity_type="contrast",
                        entity_id=treat_cid,
                        path="samples.biological_replicate_id",
                        message=(
                            f"Insufficient biological replication for inferential DE: treatment has {treat_rep_count} "
                            f"biological replicate(s), control has {control_rep_count} biological replicate(s). "
                            "Minimum 2 required; single biological replicates cannot be presented as model-ready inferential comparisons."
                        ),
                    )
                )
                if status != ContrastStatus.BLOCKED:
                    status = ContrastStatus.NEEDS_REVIEW

        # 9. Exposure Metadata Extraction (Never influenced by DNT labels)
        exp_records = self.exposures_by_condition.get(treat_cid, [])
        # Prefer condition-level exposure
        chosen_exp = next((e for e in exp_records if not e.get("sample_id")), exp_records[0] if exp_records else {})

        agent_name = extract_metadata_value(chosen_exp.get("agent_name"))
        agent_id = extract_metadata_value(chosen_exp.get("agent_identifier"))
        vehicle = extract_metadata_value(chosen_exp.get("vehicle"))
        dose = extract_metadata_value(chosen_exp.get("concentration_or_dose"))
        dose_unit = extract_metadata_value(chosen_exp.get("concentration_or_dose_unit"))
        exp_start = extract_metadata_value(chosen_exp.get("exposure_start_time_or_stage"))
        exp_age = extract_metadata_value(chosen_exp.get("developmental_age_or_stage_at_exposure"))
        exp_dur = extract_metadata_value(chosen_exp.get("exposure_duration"))
        exp_dur_unit = extract_metadata_value(chosen_exp.get("exposure_duration_unit"))
        washout_dur = extract_metadata_value(chosen_exp.get("washout_or_recovery_duration"))
        washout_unit = extract_metadata_value(chosen_exp.get("washout_or_recovery_duration_unit"))

        # 10. Construct Deterministic Contrast ID
        ctrl_part = "_and_".join(sorted(control_cond_ids)) if control_cond_ids else "no_control"
        contrast_id_raw = f"contrast_{treat_study_id}_{treat_exp_id}_{treat_cid}_vs_{ctrl_part}"
        contrast_id = sanitize_id_token(contrast_id_raw)

        return MolecularResponseContrast(
            contrast_id=contrast_id,
            study_id=treat_study_id,
            experiment_id=treat_exp_id,
            cohort_id=cohort_id,
            relationship_id=rel_id,
            treatment_condition_id=treat_cid,
            matched_control_condition_ids=tuple(sorted(control_cond_ids)),
            treatment_sample_ids=tuple(treat_sample_ids),
            control_sample_ids=tuple(active_control_sample_ids),
            treatment_biological_replicate_ids=treat_bio_reps,
            control_biological_replicate_ids=control_bio_reps,
            agent_name=str(agent_name) if agent_name is not None else None,
            agent_identifier=str(agent_id) if agent_id is not None else None,
            vehicle=str(vehicle) if vehicle is not None else None,
            concentration_or_dose=dose,
            concentration_or_dose_unit=str(dose_unit) if dose_unit is not None else None,
            exposure_start_time_or_stage=str(exp_start) if exp_start is not None else None,
            developmental_age_or_stage_at_exposure=str(exp_age) if exp_age is not None else None,
            exposure_duration=exp_dur,
            exposure_duration_unit=str(exp_dur_unit) if exp_dur_unit is not None else None,
            washout_or_recovery_duration=washout_dur,
            washout_or_recovery_duration_unit=str(washout_unit) if washout_unit is not None else None,
            treatment_collection_age_or_stage=str(treat_age_rep) if treat_age_rep is not None else None,
            control_collection_age_or_stage=str(control_age_rep) if control_age_rep is not None else None,
            biological_source_ids=all_bio_sources,
            organoid_context_ids=all_org_contexts,
            sample_comparison_exception_ids=tuple(sorted(applied_exception_ids)),
            status=status,
            findings=ordered_findings(findings),
        )

    def _build_empty_contrast(
        self,
        rel_id: str,
        treat_cid: str,
        control_cond_ids: tuple[str, ...],
        status: ContrastStatus,
        findings: list[Finding],
    ) -> MolecularResponseContrast:
        ctrl_part = "_and_".join(sorted(control_cond_ids)) if control_cond_ids else "no_control"
        contrast_id = sanitize_id_token(f"contrast_unresolved_{treat_cid}_vs_{ctrl_part}")
        return MolecularResponseContrast(
            contrast_id=contrast_id,
            study_id="unknown_study",
            experiment_id="unknown_exp",
            cohort_id="unknown_cohort",
            relationship_id=rel_id,
            treatment_condition_id=treat_cid,
            matched_control_condition_ids=tuple(sorted(control_cond_ids)),
            treatment_sample_ids=(),
            control_sample_ids=(),
            treatment_biological_replicate_ids=(),
            control_biological_replicate_ids=(),
            agent_name=None,
            agent_identifier=None,
            vehicle=None,
            concentration_or_dose=None,
            concentration_or_dose_unit=None,
            exposure_start_time_or_stage=None,
            developmental_age_or_stage_at_exposure=None,
            exposure_duration=None,
            exposure_duration_unit=None,
            washout_or_recovery_duration=None,
            washout_or_recovery_duration_unit=None,
            treatment_collection_age_or_stage=None,
            control_collection_age_or_stage=None,
            biological_source_ids=(),
            organoid_context_ids=(),
            sample_comparison_exception_ids=(),
            status=status,
            findings=ordered_findings(findings),
        )


def build_response_contrasts(
    normalized: NormalizedDataset,
    metadata: dict[str, Any],
) -> ResponseContrastDataset:
    """Public functional entrypoint to build molecular-response contrasts from normalized data and metadata."""
    builder = ResponseContrastBuilder(metadata=metadata, normalized=normalized)
    return builder.build_dataset()
