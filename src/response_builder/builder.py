"""Scientific contrast construction and validation engine for Treatment-versus-Matched-Control Contrast Builder v1.

Enforces:
1. Every treatment condition represented in an eligible cohort has an auditable disposition (never silently dropped).
2. Contrasts strictly originate from explicit treatment_control_relationships.
3. Separate contrast per treatment condition (no automated pooling of treatment conditions).
4. Preservation of all condition-planned exposure records without flattening or silent truncation.
5. Preservation of study and experiment boundaries (strictly no cross-study or cross-experiment comparisons).
6. Developmental age strict matching (preferring governed normalized values, blocking on true mismatches, review on ambiguity).
7. Organoid context and biological source compatibility evaluated on scientific attributes rather than identifier-only.
8. True biological replication requirements for inferential downstream DE (technical replicates do not inflate counts).
9. Preservation and non-silent handling of sample comparison exceptions.
10. Complete independence from DNT labels/evidence.
11. Zero gene filtering.
12. Fully deterministic contrast IDs, exposure ordering, and output sorting.
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
    ContrastExposure,
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
    """Extract governed normalized or reported value from scalar or value_record dict."""
    if val_rec is None:
        return None
    if isinstance(val_rec, (str, int, float, bool)):
        return val_rec
    if isinstance(val_rec, dict):
        norm = val_rec.get("normalized_value")
        if norm is not None:
            return norm
        rep = reported_value(val_rec)
        if rep is not None:
            return rep
        return val_rec.get("original_value")
    return None


def normalize_age_string(val: Any) -> str | None:
    """Normalize a developmental age or stage string to a canonical token for strict comparison.

    Recognizes standard age expressions (e.g. 'Day 30', 'day 30', 'd30', '30 days', 30 -> 'day_30';
    'GW 10', '10 gw' -> 'gw_10'; 'PCW 8' -> 'pcw_8'; 'week 4' -> 'week_4').
    Returns None if missing, unknown, or if no standard pattern matches (never invents a normalized age).
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

    # e.g., 'gw 10' or 'gw_10' or '10 gw'
    m_gw = re.match(r"^(?:gw|gestational week)\s*([0-9]+(?:\.[0-9]+)?)$", clean)
    if m_gw:
        return f"gw_{m_gw.group(1)}"
    m_gw_rev = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*(?:gw|gestational weeks?)$", clean)
    if m_gw_rev:
        return f"gw_{m_gw_rev.group(1)}"

    # e.g., 'pcw 8' or 'post conception week 8'
    m_pcw = re.match(r"^(?:pcw|post conception week)\s*([0-9]+(?:\.[0-9]+)?)$", clean)
    if m_pcw:
        return f"pcw_{m_pcw.group(1)}"

    # e.g., 'week 4' or '4 weeks'
    m_wk = re.match(r"^(?:week|wk)\s*([0-9]+(?:\.[0-9]+)?)$", clean)
    if m_wk:
        return f"week_{m_wk.group(1)}"
    m_wk_rev = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*(?:weeks|wks?)$", clean)
    if m_wk_rev:
        return f"week_{m_wk_rev.group(1)}"

    # If it is already a normalized canonical token like 'day_30' or 'gw_10':
    if re.match(r"^(?:day|gw|pcw|week)_[0-9]+(?:\.[0-9]+)?$", s):
        return s

    return None


def resolve_developmental_age(val_rec: Any) -> tuple[str | None, str | None, bool]:
    """Resolve developmental age/stage preferring governed normalized values.

    Returns:
        (canonical_normalized_age, reported_age_str, is_resolved)
    """
    if val_rec is None:
        return None, None, False

    if isinstance(val_rec, dict):
        gov_norm = val_rec.get("normalized_value")
        orig_val = val_rec.get("original_value")
        rep_val = reported_value(val_rec)
        reported_str = str(orig_val if orig_val is not None else (rep_val if rep_val is not None else gov_norm))

        if gov_norm is not None and str(gov_norm).strip() != "":
            norm_str = normalize_age_string(gov_norm)
            if norm_str is not None:
                return norm_str, reported_str, True
            return str(gov_norm).strip().lower().replace(" ", "_"), reported_str, True

        norm_from_orig = normalize_age_string(orig_val)
        if norm_from_orig is not None:
            return norm_from_orig, reported_str, True

        return None, reported_str, False

    norm_scalar = normalize_age_string(val_rec)
    rep_scalar = str(val_rec)
    if norm_scalar is not None:
        return norm_scalar, rep_scalar, True
    return None, rep_scalar, False


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

        # Index exposures: separate condition-planned exposures from sample-specific deviations
        self.exposures_by_condition: dict[str, list[dict[str, Any]]] = defaultdict(list)
        self.exposure_deviations_by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for exp in self.metadata.get("exposures", []):
            if not isinstance(exp, dict):
                continue
            cid = exp.get("condition_id")
            sid = exp.get("sample_id")
            if sid:
                self.exposure_deviations_by_sample[sid].append(exp)
            elif cid:
                self.exposures_by_condition[cid].append(exp)

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

        # 3. Discover all treatment conditions encountered in eligible cohorts
        self.treatment_conditions_encountered_set: set[str] = set()
        for cr in self.normalized.cohort_results:
            for sid in cr.cohort.sample_ids:
                sample = self.samples_by_id.get(sid)
                if not sample:
                    continue
                cid = sample.get("condition_id")
                if not cid:
                    continue
                cond = self.conditions_by_id.get(cid)
                if not cond:
                    continue
                role = extract_metadata_value(cond.get("treatment_control_status"))
                if role == "treatment":
                    self.treatment_conditions_encountered_set.add(cid)

        # Also include any treatment condition declared in relationships if represented in normalized cohorts
        for rel in self.metadata.get("treatment_control_relationships", []):
            if isinstance(rel, dict):
                for tcid in rel.get("treatment_condition_ids", ()):
                    t_str = str(tcid).strip()
                    if t_str and any(sid in self.cohort_result_by_sample for sid in self.samples_by_condition.get(t_str, [])):
                        self.treatment_conditions_encountered_set.add(t_str)

    def build_dataset(self) -> ResponseContrastDataset:
        """Evaluate all declared relationships and unlinked treatment conditions into auditable contrasts."""
        dataset_findings: list[Finding] = []
        contrasts: list[MolecularResponseContrast] = []
        evaluated_treatment_conditions: set[str] = set()

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
                    message="Metadata contains no explicit treatment_control_relationships.",
                )
            )

        for rel in declared_relationships:
            rel_contrasts, rel_findings = self._evaluate_relationship(rel)
            contrasts.extend(rel_contrasts)
            dataset_findings.extend(rel_findings)
            for c in rel_contrasts:
                evaluated_treatment_conditions.add(c.treatment_condition_id)

        # Conservation Invariant: No treatment condition may be silently dropped.
        # Every treatment condition represented in an eligible cohort must have an auditable disposition.
        unlinked_conditions = sorted(self.treatment_conditions_encountered_set - evaluated_treatment_conditions)
        for unlinked_cid in unlinked_conditions:
            blocked_contrast = self._build_unlinked_treatment_contrast(unlinked_cid)
            contrasts.append(blocked_contrast)

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
            treatment_conditions_encountered=tuple(sorted(self.treatment_conditions_encountered_set)),
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

    def _extract_treatment_exposures(
        self,
        treat_cid: str,
        treat_sample_ids: list[str],
    ) -> tuple[tuple[ContrastExposure, ...], tuple[dict[str, Any], ...], list[Finding]]:
        """Preserve all condition-planned exposures and sample-specific deviations deterministically."""
        findings: list[Finding] = []
        planned_raw = self.exposures_by_condition.get(treat_cid, [])

        planned_exps: list[ContrastExposure] = []
        for idx, exp in enumerate(planned_raw):
            exp_id = str(exp.get("exposure_id") or f"exp_{treat_cid}_{idx+1}")
            agent_name = extract_metadata_value(exp.get("agent_name"))
            agent_id = extract_metadata_value(exp.get("agent_identifier"))
            vehicle = extract_metadata_value(exp.get("vehicle"))
            dose = extract_metadata_value(exp.get("concentration_or_dose"))
            dose_unit = extract_metadata_value(exp.get("concentration_or_dose_unit"))
            exp_start = extract_metadata_value(exp.get("exposure_start_time_or_stage"))
            exp_age = extract_metadata_value(exp.get("developmental_age_or_stage_at_exposure"))
            exp_dur = extract_metadata_value(exp.get("exposure_duration"))
            exp_dur_unit = extract_metadata_value(exp.get("exposure_duration_unit"))
            washout_dur = extract_metadata_value(exp.get("washout_or_recovery_duration"))
            washout_unit = extract_metadata_value(exp.get("washout_or_recovery_duration_unit"))

            planned_exps.append(
                ContrastExposure(
                    exposure_id=exp_id,
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
                )
            )

        # Sort exposures deterministically by exposure_id, agent_name, agent_identifier
        sorted_planned = tuple(
            sorted(
                planned_exps,
                key=lambda e: (e.exposure_id, e.agent_name or "", e.agent_identifier or ""),
            )
        )

        # Collect sample-level exposure deviations for samples in this treatment condition
        deviations: list[dict[str, Any]] = []
        for sid in sorted(treat_sample_ids):
            for dev in self.exposure_deviations_by_sample.get(sid, []):
                deviations.append(dev)

        if deviations:
            findings.append(
                Finding(
                    severity=Severity.INFO,
                    rule_id="sample_level_exposure_deviation_present",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="exposures",
                    message=(
                        f"Treatment condition '{treat_cid}' has {len(deviations)} sample-specific exposure deviation record(s). "
                        "Condition-planned exposures remain authoritative and preserved."
                    ),
                )
            )

        return sorted_planned, tuple(deviations), findings

    def _build_unlinked_treatment_contrast(self, treat_cid: str) -> MolecularResponseContrast:
        """Construct an auditable BLOCKED contrast for a treatment condition lacking matched-control relationship."""
        treat_cond = self.conditions_by_id.get(treat_cid, {})
        treat_exp_id = treat_cond.get("experiment_id") or "unknown_exp"
        treat_exp = self.experiments_by_id.get(treat_exp_id, {})
        treat_study_id = treat_exp.get("study_id") or "unknown_study"

        treat_sample_ids = sorted(self.samples_by_condition.get(treat_cid, []))
        treat_samples = [self.samples_by_id[sid] for sid in treat_sample_ids if sid in self.samples_by_id]

        cohort_ids = {
            self.cohort_result_by_sample[sid].cohort.cohort_id
            for sid in treat_sample_ids
            if sid in self.cohort_result_by_sample
        }
        cohort_id = next(iter(cohort_ids)) if len(cohort_ids) == 1 else (
            "_vs_".join(sorted(cohort_ids)) if cohort_ids else "unknown_cohort"
        )

        treat_bio_reps = tuple(sorted({
            str(extract_metadata_value(s.get("biological_replicate_id")))
            for s in treat_samples
            if extract_metadata_value(s.get("biological_replicate_id")) is not None
        }))

        treat_exposures, deviations, exp_findings = self._extract_treatment_exposures(treat_cid, treat_sample_ids)

        treat_norm_ages: set[str] = set()
        treat_rep_ages: list[str] = []
        for s in treat_samples:
            norm_a, rep_a, _ = resolve_developmental_age(s.get("developmental_age_or_stage_at_collection"))
            if norm_a:
                treat_norm_ages.add(norm_a)
            if rep_a:
                treat_rep_ages.append(rep_a)

        treat_age_rep = treat_rep_ages[0] if len(set(treat_rep_ages)) == 1 and treat_rep_ages else (
            ", ".join(sorted(set(treat_rep_ages))) or None
        )
        treat_age_norm = next(iter(treat_norm_ages)) if len(treat_norm_ages) == 1 else None

        treat_org_contexts = tuple(sorted({
            str(extract_metadata_value(s.get("organoid_context_id")))
            for s in treat_samples
            if extract_metadata_value(s.get("organoid_context_id"))
        }))
        treat_bio_sources = tuple(sorted({
            str(extract_metadata_value(s.get("biological_source_id")))
            for s in treat_samples
            if extract_metadata_value(s.get("biological_source_id"))
        }))

        findings = [
            Finding(
                severity=Severity.ERROR,
                rule_id="missing_explicit_treatment_control_relationship",
                entity_type="condition",
                entity_id=treat_cid,
                path=f"conditions.{treat_cid}",
                message=(
                    f"Treatment condition '{treat_cid}' exists in an eligible normalized cohort but has no "
                    "explicit treatment_control_relationship declared in metadata."
                ),
            )
        ] + exp_findings

        contrast_id = sanitize_id_token(f"contrast_{treat_study_id}_{treat_exp_id}_{treat_cid}_vs_no_control")

        single_exp = treat_exposures[0] if len(treat_exposures) == 1 else None

        return MolecularResponseContrast(
            contrast_id=contrast_id,
            study_id=treat_study_id,
            experiment_id=treat_exp_id,
            cohort_id=cohort_id,
            relationship_id="missing_relationship",
            treatment_condition_id=treat_cid,
            matched_control_condition_ids=(),
            treatment_sample_ids=tuple(treat_sample_ids),
            control_sample_ids=(),
            treatment_biological_replicate_ids=treat_bio_reps,
            control_biological_replicate_ids=(),
            treatment_exposures=treat_exposures,
            sample_level_exposure_deviations=deviations,
            agent_name=single_exp.agent_name if single_exp else None,
            agent_identifier=single_exp.agent_identifier if single_exp else None,
            vehicle=single_exp.vehicle if single_exp else None,
            concentration_or_dose=single_exp.concentration_or_dose if single_exp else None,
            concentration_or_dose_unit=single_exp.concentration_or_dose_unit if single_exp else None,
            exposure_start_time_or_stage=single_exp.exposure_start_time_or_stage if single_exp else None,
            developmental_age_or_stage_at_exposure=single_exp.developmental_age_or_stage_at_exposure if single_exp else None,
            exposure_duration=single_exp.exposure_duration if single_exp else None,
            exposure_duration_unit=single_exp.exposure_duration_unit if single_exp else None,
            washout_or_recovery_duration=single_exp.washout_or_recovery_duration if single_exp else None,
            washout_or_recovery_duration_unit=single_exp.washout_or_recovery_duration_unit if single_exp else None,
            treatment_collection_age_or_stage=treat_age_rep,
            control_collection_age_or_stage=None,
            treatment_collection_age_normalized=treat_age_norm,
            control_collection_age_normalized=None,
            biological_source_ids=treat_bio_sources,
            organoid_context_ids=treat_org_contexts,
            sample_comparison_exception_ids=(),
            status=ContrastStatus.BLOCKED,
            findings=ordered_findings(findings),
        )

    def _evaluate_organoid_context_compatibility(
        self,
        treat_cid: str,
        treat_samples: list[dict[str, Any]],
        control_samples: list[dict[str, Any]],
    ) -> tuple[tuple[str, ...], list[Finding], ContrastStatus]:
        """Evaluate biological compatibility based on governed scientific attributes, not ID identity alone."""
        findings: list[Finding] = []
        status = ContrastStatus.ELIGIBLE

        treat_ctx_ids = tuple(sorted({
            str(extract_metadata_value(s.get("organoid_context_id")))
            for s in treat_samples
            if extract_metadata_value(s.get("organoid_context_id"))
        }))
        control_ctx_ids = tuple(sorted({
            str(extract_metadata_value(s.get("organoid_context_id")))
            for s in control_samples
            if extract_metadata_value(s.get("organoid_context_id"))
        }))
        all_ctx_ids = tuple(sorted(set(treat_ctx_ids) | set(control_ctx_ids)))

        if not treat_ctx_ids or not control_ctx_ids:
            findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="unresolved_organoid_context_compatibility",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="samples.organoid_context_id",
                    message=(
                        f"Organoid context identifier missing or unspecified in treatment ({treat_ctx_ids}) "
                        f"or control ({control_ctx_ids}). Cannot establish biological compatibility."
                    ),
                )
            )
            return all_ctx_ids, findings, ContrastStatus.NEEDS_REVIEW

        treat_records = [self.organoid_contexts_by_id.get(cid) for cid in treat_ctx_ids]
        control_records = [self.organoid_contexts_by_id.get(cid) for cid in control_ctx_ids]

        if any(r is None for r in treat_records) or any(r is None for r in control_records):
            unresolved = [
                cid
                for cid, rec in zip(treat_ctx_ids + control_ctx_ids, treat_records + control_records, strict=False)
                if rec is None
            ]
            findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="unresolved_organoid_context_compatibility",
                    entity_type="organoid_context",
                    entity_id=treat_cid,
                    path="metadata.organoid_contexts",
                    message=f"Organoid context record(s) {sorted(set(unresolved))} referenced by samples do not exist in metadata.",
                )
            )
            return all_ctx_ids, findings, ContrastStatus.NEEDS_REVIEW

        # Extract governed scientific attributes
        treat_types = {extract_metadata_value(r.get("organoid_type")) for r in treat_records if r and extract_metadata_value(r.get("organoid_type"))}
        control_types = {extract_metadata_value(r.get("organoid_type")) for r in control_records if r and extract_metadata_value(r.get("organoid_type"))}

        treat_regions = {extract_metadata_value(r.get("brain_region_or_model_identity")) for r in treat_records if r and extract_metadata_value(r.get("brain_region_or_model_identity"))}
        control_regions = {extract_metadata_value(r.get("brain_region_or_model_identity")) for r in control_records if r and extract_metadata_value(r.get("brain_region_or_model_identity"))}

        treat_protocols = {extract_metadata_value(r.get("differentiation_protocol_identifier")) for r in treat_records if r and extract_metadata_value(r.get("differentiation_protocol_identifier")) is not None}
        control_protocols = {extract_metadata_value(r.get("differentiation_protocol_identifier")) for r in control_records if r and extract_metadata_value(r.get("differentiation_protocol_identifier")) is not None}

        treat_versions = {extract_metadata_value(r.get("differentiation_protocol_version")) for r in treat_records if r and extract_metadata_value(r.get("differentiation_protocol_version")) is not None}
        control_versions = {extract_metadata_value(r.get("differentiation_protocol_version")) for r in control_records if r and extract_metadata_value(r.get("differentiation_protocol_version")) is not None}

        # Check for biological model / organoid type conflict
        if treat_types and control_types and treat_types.isdisjoint(control_types):
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="organoid_context_mismatch",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="organoid_contexts.organoid_type",
                    message=(
                        f"Conflicting organoid types between treatment ({sorted(treat_types)}) "
                        f"and control ({sorted(control_types)}). Biologically incompatible comparison."
                    ),
                )
            )
            status = ContrastStatus.BLOCKED

        # Check for brain region / model identity conflict
        if treat_regions and control_regions and treat_regions.isdisjoint(control_regions):
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="organoid_context_mismatch",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="organoid_contexts.brain_region_or_model_identity",
                    message=(
                        f"Conflicting brain region/model identity between treatment ({sorted(treat_regions)}) "
                        f"and control ({sorted(control_regions)}). Biologically incompatible comparison."
                    ),
                )
            )
            status = ContrastStatus.BLOCKED

        # Check differentiation protocol conflict when known
        if treat_protocols and control_protocols and treat_protocols.isdisjoint(control_protocols):
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="differentiation_protocol_mismatch",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="organoid_contexts.differentiation_protocol_identifier",
                    message=(
                        f"Conflicting differentiation protocols between treatment ({sorted(treat_protocols)}) "
                        f"and control ({sorted(control_protocols)})."
                    ),
                )
            )
            status = ContrastStatus.BLOCKED

        # Check protocol version conflict when known
        if treat_versions and control_versions and treat_versions.isdisjoint(control_versions):
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="differentiation_protocol_version_mismatch",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="organoid_contexts.differentiation_protocol_version",
                    message=(
                        f"Conflicting differentiation protocol versions between treatment ({sorted(treat_versions)}) "
                        f"and control ({sorted(control_versions)})."
                    ),
                )
            )
            status = ContrastStatus.BLOCKED

        # Check for missing critical attributes when records exist but fields are blank
        if not treat_types and not control_types and not treat_regions and not control_regions:
            findings.append(
                Finding(
                    severity=Severity.REVIEW,
                    rule_id="unresolved_organoid_context_compatibility",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="organoid_contexts",
                    message="Organoid context attributes (organoid_type, brain_region_or_model_identity) are missing or blank.",
                )
            )
            if status != ContrastStatus.BLOCKED:
                status = ContrastStatus.NEEDS_REVIEW

        return all_ctx_ids, findings, status

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
                    path="samples",
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
                    path="samples",
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
                    valid_replacements: list[str] = []
                    for rep_sid in replacements:
                        rep_s = self.samples_by_id.get(rep_sid)
                        if rep_s and rep_s.get("condition_id") in control_cond_ids:
                            valid_replacements.append(rep_sid)

                    if valid_replacements:
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
        treat_norm_ages: set[str] = set()
        treat_rep_ages: list[str] = []
        treat_unresolved_flags: list[bool] = []
        for s in treat_samples:
            norm_a, rep_a, is_res = resolve_developmental_age(s.get("developmental_age_or_stage_at_collection"))
            if norm_a:
                treat_norm_ages.add(norm_a)
            if rep_a:
                treat_rep_ages.append(rep_a)
            treat_unresolved_flags.append(not is_res)

        control_norm_ages: set[str] = set()
        control_rep_ages: list[str] = []
        control_unresolved_flags: list[bool] = []
        for s in control_samples:
            norm_a, rep_a, is_res = resolve_developmental_age(s.get("developmental_age_or_stage_at_collection"))
            if norm_a:
                control_norm_ages.add(norm_a)
            if rep_a:
                control_rep_ages.append(rep_a)
            control_unresolved_flags.append(not is_res)

        treat_age_rep = treat_rep_ages[0] if len(set(treat_rep_ages)) == 1 and treat_rep_ages else (
            ", ".join(sorted(set(treat_rep_ages))) or None
        )
        control_age_rep = control_rep_ages[0] if len(set(control_rep_ages)) == 1 and control_rep_ages else (
            ", ".join(sorted(set(control_rep_ages))) or None
        )

        treat_norm_age_single = next(iter(treat_norm_ages)) if len(treat_norm_ages) == 1 else None
        control_norm_age_single = next(iter(control_norm_ages)) if len(control_norm_ages) == 1 else None

        if treat_samples and control_samples:
            # Check for unresolved or ambiguous age
            if any(treat_unresolved_flags) or any(control_unresolved_flags) or not treat_rep_ages or not control_rep_ages:
                # If both have identical reported text (e.g. both 'stage_A'), they match identically
                if set(treat_rep_ages) == set(control_rep_ages) and len(set(treat_rep_ages)) == 1:
                    pass
                else:
                    findings.append(
                        Finding(
                            severity=Severity.REVIEW,
                            rule_id="unresolved_developmental_age_equivalence",
                            entity_type="contrast",
                            entity_id=treat_cid,
                            path="samples.developmental_age_or_stage_at_collection",
                            message=(
                                f"Developmental age equivalence cannot be established from metadata "
                                f"(treatment reported: {sorted(set(treat_rep_ages))}, control reported: {sorted(set(control_rep_ages))}). "
                                "Never invent a normalized age; requires expert review."
                            ),
                        )
                    )
                    if status != ContrastStatus.BLOCKED:
                        status = ContrastStatus.NEEDS_REVIEW
            elif len(treat_norm_ages) != 1 or len(control_norm_ages) != 1:
                findings.append(
                    Finding(
                        severity=Severity.REVIEW,
                        rule_id="developmental_age_ambiguous",
                        entity_type="contrast",
                        entity_id=treat_cid,
                        path="samples.developmental_age_or_stage_at_collection",
                        message=(
                            f"Developmental age or stage at collection is divergent within arm across samples "
                            f"(treatment normalized: {sorted(treat_norm_ages)}, control normalized: {sorted(control_norm_ages)}). "
                            "Requires expert review."
                        ),
                    )
                )
                if status != ContrastStatus.BLOCKED:
                    status = ContrastStatus.NEEDS_REVIEW
            else:
                t_age = next(iter(treat_norm_ages))
                c_age = next(iter(control_norm_ages))
                if t_age != c_age:
                    findings.append(
                        Finding(
                            severity=Severity.ERROR,
                            rule_id="developmental_age_mismatch",
                            entity_type="contrast",
                            entity_id=treat_cid,
                            path="samples.developmental_age_or_stage_at_collection",
                            message=(
                                f"Developmental age mismatch: treatment collected at '{treat_age_rep}' (normalized '{t_age}'), "
                                f"while control collected at '{control_age_rep}' (normalized '{c_age}'). "
                                "Comparing across different developmental ages is strictly invalid."
                            ),
                        )
                    )
                    status = ContrastStatus.BLOCKED

        # 7. Organoid Context Compatibility (Attribute-based scientific compatibility)
        all_org_contexts, org_findings, org_status = self._evaluate_organoid_context_compatibility(
            treat_cid, treat_samples, control_samples
        )
        findings.extend(org_findings)
        if org_status == ContrastStatus.BLOCKED:
            status = ContrastStatus.BLOCKED
        elif org_status == ContrastStatus.NEEDS_REVIEW and status != ContrastStatus.BLOCKED:
            status = ContrastStatus.NEEDS_REVIEW

        # 8. Biological Source Compatibility
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

        # Check biological source records for species consistency
        treat_sources = [self.biological_sources_by_id.get(sid) for sid in treat_bio_sources]
        ctrl_sources = [self.biological_sources_by_id.get(sid) for sid in control_bio_sources]
        treat_species = {extract_metadata_value(r.get("organism_or_species")) for r in treat_sources if r and extract_metadata_value(r.get("organism_or_species"))}
        ctrl_species = {extract_metadata_value(r.get("organism_or_species")) for r in ctrl_sources if r and extract_metadata_value(r.get("organism_or_species"))}

        if treat_species and ctrl_species and treat_species.isdisjoint(ctrl_species):
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="cross_species_comparison_blocked",
                    entity_type="contrast",
                    entity_id=treat_cid,
                    path="biological_sources.organism_or_species",
                    message=f"Cross-species comparison blocked: treatment ({sorted(treat_species)}) vs control ({sorted(ctrl_species)}).",
                )
            )
            status = ContrastStatus.BLOCKED

        # 9. Biological Replicate Counting (Technical replicates do NOT inflate counts)
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

        # 10. Exposure Metadata Extraction (Never influenced by DNT labels; all planned records preserved)
        treat_exposures, deviations, exp_findings = self._extract_treatment_exposures(treat_cid, treat_sample_ids)
        findings.extend(exp_findings)

        single_exp = treat_exposures[0] if len(treat_exposures) == 1 else None

        # 11. Construct Deterministic Contrast ID
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
            treatment_exposures=treat_exposures,
            sample_level_exposure_deviations=deviations,
            agent_name=single_exp.agent_name if single_exp else None,
            agent_identifier=single_exp.agent_identifier if single_exp else None,
            vehicle=single_exp.vehicle if single_exp else None,
            concentration_or_dose=single_exp.concentration_or_dose if single_exp else None,
            concentration_or_dose_unit=single_exp.concentration_or_dose_unit if single_exp else None,
            exposure_start_time_or_stage=single_exp.exposure_start_time_or_stage if single_exp else None,
            developmental_age_or_stage_at_exposure=single_exp.developmental_age_or_stage_at_exposure if single_exp else None,
            exposure_duration=single_exp.exposure_duration if single_exp else None,
            exposure_duration_unit=single_exp.exposure_duration_unit if single_exp else None,
            washout_or_recovery_duration=single_exp.washout_or_recovery_duration if single_exp else None,
            washout_or_recovery_duration_unit=single_exp.washout_or_recovery_duration_unit if single_exp else None,
            treatment_collection_age_or_stage=str(treat_age_rep) if treat_age_rep is not None else None,
            control_collection_age_or_stage=str(control_age_rep) if control_age_rep is not None else None,
            treatment_collection_age_normalized=treat_norm_age_single,
            control_collection_age_normalized=control_norm_age_single,
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
            treatment_exposures=(),
            sample_level_exposure_deviations=(),
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
            treatment_collection_age_normalized=None,
            control_collection_age_normalized=None,
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
