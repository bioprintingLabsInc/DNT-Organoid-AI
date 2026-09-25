"""Metadata-derived Normalization Cohort construction.

Preserves study and experiment boundaries, bulk RNA-seq modality, compatible technical context,
and ensures biologically matched treatment and control samples remain in the same cohort.
"""

from collections import defaultdict
import hashlib
import json
from typing import Any

from src.ingestion.registry import reported_value
from .models import (
    Finding,
    NormalizationCohort,
    ResponseEligibility,
    Severity,
)


def _extract_val(val_rec: Any) -> Any:
    """Extract reported or normalized value from scalar or value_record dict."""
    if val_rec is None:
        return None
    if isinstance(val_rec, (str, int, float, bool)):
        return val_rec
    if isinstance(val_rec, dict):
        rep = reported_value(val_rec)
        if rep is not None:
            return rep
        return val_rec.get("normalized_value") or val_rec.get("original_value")
    return None


class NormalizationCohortBuilder:
    """Constructs Normalization Cohorts from canonical metadata and matrix samples."""

    def __init__(self, metadata: dict[str, Any]) -> None:
        self.metadata = metadata
        self._index_metadata()

    def _index_metadata(self) -> None:
        self.samples_by_id = {
            s["sample_id"]: s
            for s in self.metadata.get("samples", [])
            if isinstance(s, dict) and "sample_id" in s
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
        self.conditions_by_id = {
            c["condition_id"]: c
            for c in self.metadata.get("conditions", [])
            if isinstance(c, dict) and "condition_id" in c
        }
        self.assays_by_id = {
            a["assay_id"]: a
            for a in self.metadata.get("sequencing_assays", [])
            if isinstance(a, dict) and "assay_id" in a
        }
        self.relationships = [
            r
            for r in self.metadata.get("treatment_control_relationships", [])
            if isinstance(r, dict) and "relationship_id" in r
        ]

        # Map sample_id to assay_ids
        self.assays_by_sample: dict[str, list[str]] = defaultdict(list)
        for assay_id, assay in self.assays_by_id.items():
            for sid in assay.get("sample_ids", []):
                self.assays_by_sample[sid].append(assay_id)

    def build_cohorts(
        self,
        sample_ids: tuple[str, ...],
        source_asset_id: str | None = None,
    ) -> tuple[tuple[NormalizationCohort, ...], list[Finding]]:
        """Group samples into Normalization Cohorts obeying all scientific invariants."""
        findings: list[Finding] = []

        # 1. Validate sample presence and basic metadata integrity
        sample_records: list[dict[str, Any]] = []
        for sid in sample_ids:
            s_rec = self.samples_by_id.get(sid)
            if not s_rec:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="sample_missing_from_metadata",
                        entity_type="sample",
                        entity_id=sid,
                        path=f"samples.{sid}",
                        message=f"Sample '{sid}' in expression matrix is absent from canonical metadata.",
                    )
                )
                continue
            sample_records.append(s_rec)

        if any(f.severity == Severity.ERROR for f in findings):
            return (), findings

        # 2. Reconcile study and experiment per sample
        exp_study_groups: dict[tuple[str, str], list[str]] = defaultdict(list)
        sample_to_exp: dict[str, str] = {}
        sample_to_study: dict[str, str] = {}
        sample_to_cond: dict[str, str] = {}

        for s_rec in sample_records:
            sid = s_rec["sample_id"]
            exp_id = s_rec.get("experiment_id")
            cond_id = s_rec.get("condition_id")
            sample_to_cond[sid] = cond_id

            if not exp_id or exp_id not in self.experiments_by_id:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="sample_experiment_missing",
                        entity_type="sample",
                        entity_id=sid,
                        path=f"samples.{sid}.experiment_id",
                        message=f"Sample '{sid}' references missing experiment '{exp_id}'.",
                    )
                )
                continue

            exp_rec = self.experiments_by_id[exp_id]
            study_id = exp_rec.get("study_id")
            if not study_id:
                findings.append(
                    Finding(
                        severity=Severity.ERROR,
                        rule_id="experiment_study_missing",
                        entity_type="experiment",
                        entity_id=exp_id,
                        path=f"experiments.{exp_id}.study_id",
                        message=f"Experiment '{exp_id}' has no study_id.",
                    )
                )
                continue

            sample_to_exp[sid] = exp_id
            sample_to_study[sid] = study_id
            exp_study_groups[(study_id, exp_id)].append(sid)

        if any(f.severity == Severity.ERROR for f in findings):
            return (), findings

        # Check for cross-study pooling attempt across the input asset
        unique_studies = {study for study, _ in exp_study_groups.keys()}
        if len(unique_studies) > 1:
            findings.append(
                Finding(
                    severity=Severity.ERROR,
                    rule_id="cross_study_pooling_forbidden",
                    entity_type="dataset",
                    entity_id=source_asset_id,
                    path="metadata.studies",
                    message=(
                        f"Matrix samples span multiple distinct studies: {sorted(unique_studies)}. "
                        "Pooling samples across unrelated public studies into a shared normalization cohort is strictly prohibited."
                    ),
                )
            )
            return (), findings

        cohorts: list[NormalizationCohort] = []

        # 3. For each (study_id, experiment_id), construct cohorts
        for (study_id, exp_id), exp_samples in exp_study_groups.items():
            # Check modality for all samples
            for sid in exp_samples:
                assigned_assays = self.assays_by_sample.get(sid, [])
                if not assigned_assays:
                    findings.append(
                        Finding(
                            severity=Severity.REVIEW,
                            rule_id="sample_lacks_assay",
                            entity_type="sample",
                            entity_id=sid,
                            path=f"samples.{sid}",
                            message=f"Sample '{sid}' has no assigned sequencing assay in canonical metadata.",
                        )
                    )
                    continue

                for aid in assigned_assays:
                    a_rec = self.assays_by_id.get(aid, {})
                    mod = _extract_val(a_rec.get("rna_seq_modality"))
                    if mod != "bulk_rna_seq":
                        findings.append(
                            Finding(
                                severity=Severity.ERROR,
                                rule_id="incompatible_modality",
                                entity_type="sequencing_assay",
                                entity_id=aid,
                                path=f"sequencing_assays.{aid}.rna_seq_modality",
                                message=(
                                    f"Sequencing assay '{aid}' has modality '{mod}'. "
                                    "Bulk RNA-seq Normalization v1 strictly requires 'bulk_rna_seq'."
                                ),
                            )
                        )

            if any(f.severity == Severity.ERROR for f in findings):
                continue

            # Resolve technical context across samples in this experiment
            sample_tech_context: dict[str, tuple[Any, Any, Any, Any]] = {}
            for sid in exp_samples:
                assigned_assays = self.assays_by_sample.get(sid, [])
                if assigned_assays:
                    a_rec = self.assays_by_id[assigned_assays[0]]
                    strat = _extract_val(a_rec.get("library_strategy"))
                    plat = _extract_val(a_rec.get("sequencing_platform"))
                    strand = _extract_val(a_rec.get("strandedness"))
                    layout = _extract_val(a_rec.get("read_layout"))
                    sample_tech_context[sid] = (strat, plat, strand, layout)
                else:
                    sample_tech_context[sid] = (None, None, None, None)

            # Map conditions to relationships within this experiment
            exp_relationships = [
                r for r in self.relationships
                if any(
                    self.conditions_by_id.get(cid, {}).get("experiment_id") == exp_id
                    for cid in r.get("treatment_condition_ids", []) + r.get("matched_control_condition_ids", [])
                )
            ]

            # Build connected components of conditions linked by treatment-control relationships
            # to guarantee that treatment conditions and matched control conditions stay together
            condition_adjacency: dict[str, set[str]] = defaultdict(set)
            condition_to_rel_ids: dict[str, set[str]] = defaultdict(set)

            for r in exp_relationships:
                rid = r["relationship_id"]
                t_cids = set(r.get("treatment_condition_ids", []))
                c_cids = set(r.get("matched_control_condition_ids", []))
                all_rel_cids = t_cids | c_cids
                for c1 in all_rel_cids:
                    condition_to_rel_ids[c1].add(rid)
                    for c2 in all_rel_cids:
                        if c1 != c2:
                            condition_adjacency[c1].add(c2)

            visited_conditions: set[str] = set()
            condition_components: list[set[str]] = []

            exp_conditions = {sample_to_cond[sid] for sid in exp_samples if sample_to_cond[sid]}
            for cid in exp_conditions:
                if cid not in visited_conditions:
                    comp = set()
                    queue = [cid]
                    visited_conditions.add(cid)
                    while queue:
                        curr = queue.pop(0)
                        comp.add(curr)
                        for neighbor in condition_adjacency.get(curr, []):
                            if neighbor not in visited_conditions:
                                visited_conditions.add(neighbor)
                                queue.append(neighbor)
                    condition_components.append(comp)

            # Assign samples to components
            for comp_idx, comp_cids in enumerate(condition_components, 1):
                cohort_samples = [sid for sid in exp_samples if sample_to_cond.get(sid) in comp_cids]
                if not cohort_samples:
                    continue

                # Check technical context compatibility within the planned comparison cohort
                tech_contexts = {sample_tech_context[sid] for sid in cohort_samples}
                # If library_strategy differs across matched treatment and controls: conflict!
                library_strategies = {tc[0] for tc in tech_contexts if tc[0] is not None}
                if len(library_strategies) > 1:
                    findings.append(
                        Finding(
                            severity=Severity.REVIEW,
                            rule_id="cohort_technical_context_ambiguous",
                            entity_type="cohort",
                            entity_id=f"cohort_{exp_id}_{comp_idx}",
                            path="library_strategy",
                            message=(
                                f"Samples in comparison cohort have conflicting library strategies: {library_strategies}. "
                                "Biologically matched samples cannot be reliably co-normalized across differing library strategies."
                            ),
                        )
                    )

                first_tc = next(iter(tech_contexts)) if tech_contexts else (None, None, None, None)
                strat, plat, strand, layout = first_tc

                # Collect all assay_ids represented in this cohort
                cohort_assays = sorted({aid for sid in cohort_samples for aid in self.assays_by_sample.get(sid, [])})

                # Determine treatment vs control vs unassigned
                treatment_sids: list[str] = []
                control_sids: list[str] = []
                unassigned_sids: list[str] = []

                for sid in cohort_samples:
                    cid = sample_to_cond.get(sid)
                    c_rec = self.conditions_by_id.get(cid, {})
                    tc_status = _extract_val(c_rec.get("treatment_control_status"))
                    if tc_status == "treatment":
                        treatment_sids.append(sid)
                    elif tc_status == "control":
                        control_sids.append(sid)
                    else:
                        unassigned_sids.append(sid)

                # Collect relationship IDs
                comp_rel_ids = sorted({rid for cid in comp_cids for rid in condition_to_rel_ids.get(cid, set())})

                # Create deterministic metadata hash for the cohort
                cohort_hash_payload = {
                    "study_id": study_id,
                    "experiment_id": exp_id,
                    "sample_ids": sorted(cohort_samples),
                    "assay_ids": cohort_assays,
                    "library_strategy": strat,
                    "sequencing_platform": plat,
                    "conditions": sorted(comp_cids),
                    "relationship_ids": comp_rel_ids,
                }
                cohort_hash = hashlib.sha256(
                    json.dumps(cohort_hash_payload, sort_keys=True).encode("utf-8")
                ).hexdigest()

                cohort_id = f"cohort_{study_id}_{exp_id}_{comp_idx}"

                cohort = NormalizationCohort(
                    cohort_id=cohort_id,
                    study_id=study_id,
                    experiment_id=exp_id,
                    sample_ids=tuple(cohort_samples),
                    assay_ids=tuple(cohort_assays),
                    rna_seq_modality="bulk_rna_seq",
                    library_strategy=strat,
                    sequencing_platform=plat,
                    strandedness=strand,
                    read_layout=layout,
                    treatment_sample_ids=tuple(treatment_sids),
                    control_sample_ids=tuple(control_sids),
                    unassigned_sample_ids=tuple(unassigned_sids),
                    treatment_control_relationship_ids=tuple(comp_rel_ids),
                    metadata_hash=cohort_hash,
                )
                cohorts.append(cohort)

        return tuple(cohorts), findings

    def resolve_sample_controls(
        self,
        cohort: NormalizationCohort,
    ) -> dict[str, tuple[ResponseEligibility, tuple[str, ...], str]]:
        """Resolve matched control sample IDs and downstream ResponseEligibility for each sample in a cohort.

        Returns:
            dict mapping sample_id -> (ResponseEligibility, matched_control_sample_ids, treatment_control_status)
        """
        results: dict[str, tuple[ResponseEligibility, tuple[str, ...], str]] = {}

        # Build condition -> matched control conditions mapping for this cohort's relationships
        condition_to_controls: dict[str, set[str]] = defaultdict(set)
        for rid in cohort.treatment_control_relationship_ids:
            rel = next((r for r in self.relationships if r.get("relationship_id") == rid), None)
            if rel:
                t_cids = rel.get("treatment_condition_ids", [])
                c_cids = rel.get("matched_control_condition_ids", [])
                for t_cid in t_cids:
                    condition_to_controls[t_cid].update(c_cids)

        # Index control samples in the cohort by condition_id
        control_samples_by_condition: dict[str, list[str]] = defaultdict(list)
        for sid in cohort.control_sample_ids:
            s_rec = self.samples_by_id.get(sid, {})
            cid = s_rec.get("condition_id")
            if cid:
                control_samples_by_condition[cid].append(sid)

        for sid in cohort.sample_ids:
            s_rec = self.samples_by_id.get(sid, {})
            cid = s_rec.get("condition_id")
            c_rec = self.conditions_by_id.get(cid, {})
            tc_status = _extract_val(c_rec.get("treatment_control_status")) or "unassigned"

            if tc_status == "control":
                results[sid] = (ResponseEligibility.CONTROL_BASELINE, (), "control")
            elif tc_status == "treatment":
                matched_cids = condition_to_controls.get(cid, set())
                matched_sids: list[str] = []
                for m_cid in sorted(matched_cids):
                    matched_sids.extend(control_samples_by_condition.get(m_cid, []))
                matched_sids = sorted(set(matched_sids))

                if matched_sids:
                    elig = ResponseEligibility.ELIGIBLE_WITH_MATCHED_CONTROL
                else:
                    elig = ResponseEligibility.INELIGIBLE_MISSING_MATCHED_CONTROL

                results[sid] = (elig, tuple(matched_sids), "treatment")
            else:
                results[sid] = (ResponseEligibility.UNASSIGNED_OR_AMBIGUOUS, (), "unassigned")

        return results
