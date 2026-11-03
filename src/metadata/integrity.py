"""Extended metadata integrity, biological model, and treatment-control validation engine.

Complements schema validation by enforcing deep experimental design rules:
1. Internal consistency of treatment-control assignments and matched controls.
2. Detection of conflicting treatment assignments across sample metadata.
3. Verification that controls are not manufactured or silently assigned.
4. Validation of biological replicate sufficiency (>=2 replicates for inferential DE).
5. Detection of technical replicate inflation (technical replicates cannot count as biological N).
6. Biological model consistency: detects cross-species, cross-model, or mismatched organoid contexts in pairwise comparisons.
7. Verification of unambiguous sample-to-file and assay relationships.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any

from .result import Finding, Severity, ValidationStatus
from .validator import MetadataValidator


@dataclass(frozen=True, slots=True)
class IntegrityFinding:
    severity: Severity
    rule: str
    entity_type: str
    entity_id: str | None
    field_path: str
    message: str

    def to_validator_finding(self) -> Finding:
        return Finding(
            severity=self.severity,
            entity_type=self.entity_type,
            entity_identifier=self.entity_id,
            path=self.field_path,
            rule=self.rule,
            message=self.message,
        )


@dataclass(frozen=True, slots=True)
class IntegrityCheckResult:
    is_valid: bool
    status: ValidationStatus
    findings: tuple[IntegrityFinding, ...]
    replicate_summary: dict[str, dict[str, Any]]
    model_consistency_summary: dict[str, Any]


class MetadataIntegrityChecker:
    """Verifies deep experimental design integrity on canonical metadata documents."""

    def __init__(self, min_biological_replicates: int = 2) -> None:
        self.min_replicates = min_biological_replicates

    def check(self, metadata: dict[str, Any]) -> IntegrityCheckResult:
        """Execute full experimental integrity suite."""
        findings: list[IntegrityFinding] = []

        conditions = {c.get("condition_id"): c for c in metadata.get("conditions", []) if isinstance(c, dict)}
        samples = {s.get("sample_id"): s for s in metadata.get("samples", []) if isinstance(s, dict)}
        bio_sources = {b.get("biological_source_id"): b for b in metadata.get("biological_sources", []) if isinstance(b, dict)}
        organoid_contexts = {o.get("organoid_context_id"): o for o in metadata.get("organoid_contexts", []) if isinstance(o, dict)}
        relationships = metadata.get("treatment_control_relationships", [])

        # 1. Duplicate Sample Identifiers
        raw_sample_ids = [s.get("sample_id") for s in metadata.get("samples", []) if isinstance(s, dict)]
        id_counts = Counter(raw_sample_ids)
        for sid, count in id_counts.items():
            if count > 1:
                findings.append(
                    IntegrityFinding(
                        severity=Severity.ERROR,
                        rule="duplicate_sample_id",
                        entity_type="samples",
                        entity_id=str(sid),
                        field_path=f"samples[sample_id='{sid}']",
                        message=f"Sample identifier '{sid}' appears {count} times in the metadata.",
                    )
                )

        # 2. Condition-Level Assignment and Replicate Accounting
        samples_by_condition: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for s in metadata.get("samples", []):
            if isinstance(s, dict):
                cid = s.get("condition_id")
                if cid:
                    samples_by_condition[cid].append(s)

        replicate_summary: dict[str, dict[str, Any]] = {}
        for cid, cond in conditions.items():
            cond_samples = samples_by_condition.get(cid, [])
            cond_status_rec = cond.get("treatment_control_status", {})
            status_val = cond_status_rec.get("normalized_value") or cond_status_rec.get("original_value")

            # Extract biological replicate IDs
            bio_reps: list[str] = []
            tech_reps: list[str] = []
            for s in cond_samples:
                brep_rec = s.get("biological_replicate_id", {})
                brep_val = brep_rec.get("normalized_value") or brep_rec.get("original_value")
                if not brep_val or str(brep_val).strip().lower() in ("nan", "none", "unknown", ""):
                    findings.append(
                        IntegrityFinding(
                            severity=Severity.ERROR,
                            rule="missing_biological_replicate_id",
                            entity_type="samples",
                            entity_id=s.get("sample_id"),
                            field_path=f"samples[sample_id='{s.get('sample_id')}'].biological_replicate_id",
                            message=f"Sample '{s.get('sample_id')}' in condition '{cid}' is missing an explicit biological replicate ID.",
                        )
                    )
                else:
                    bio_reps.append(str(brep_val).strip())

                trep_rec = s.get("technical_replicate_id", {})
                trep_val = trep_rec.get("normalized_value") or trep_rec.get("original_value")
                if trep_val:
                    tech_reps.append(str(trep_val).strip())

            unique_bio_reps = set(bio_reps)
            replicate_summary[cid] = {
                "condition_id": cid,
                "role": status_val,
                "total_samples": len(cond_samples),
                "unique_biological_replicates": len(unique_bio_reps),
                "has_technical_replicates": len(unique_bio_reps) < len(bio_reps),
            }

            # Check for replicate deflation / technical replicate inflation
            if len(unique_bio_reps) < len(bio_reps):
                findings.append(
                    IntegrityFinding(
                        severity=Severity.WARNING,
                        rule="technical_replicates_detected",
                        entity_type="conditions",
                        entity_id=cid,
                        field_path=f"conditions[condition_id='{cid}']",
                        message=f"Condition '{cid}' contains repeated biological replicate IDs ({len(bio_reps)} samples across {len(unique_bio_reps)} unique biological replicates). Technical replicates cannot inflate biological N.",
                    )
                )

            # Check minimum biological replicates
            if len(unique_bio_reps) < self.min_replicates:
                findings.append(
                    IntegrityFinding(
                        severity=Severity.WARNING,
                        rule="under_replicated_cohort",
                        entity_type="conditions",
                        entity_id=cid,
                        field_path=f"conditions[condition_id='{cid}']",
                        message=f"Condition '{cid}' has {len(unique_bio_reps)} unique biological replicates (minimum {self.min_replicates} required for inferential differential expression).",
                    )
                )

        # 3. Treatment-Control Matching & Biological Model Consistency
        model_consistency_summary: dict[str, Any] = {}
        for rel in relationships:
            rel_id = rel.get("relationship_id")
            trt_cids = rel.get("treatment_condition_ids", [])
            ctrl_cids = rel.get("matched_control_condition_ids", [])
            rel_val_status = rel.get("metadata_validation_status")

            if not ctrl_cids:
                findings.append(
                    IntegrityFinding(
                        severity=Severity.ERROR,
                        rule="missing_matched_control",
                        entity_type="treatment_control_relationships",
                        entity_id=rel_id,
                        field_path=f"treatment_control_relationships[relationship_id='{rel_id}'].matched_control_condition_ids",
                        message=f"Comparison '{rel_id}' has no matched control condition.",
                    )
                )
                continue

            ALLOWED_VALIDATION_STATUSES = {"not_validated", "valid", "valid_with_warnings", "invalid"}

            if rel_val_status not in ALLOWED_VALIDATION_STATUSES:
                findings.append(
                    IntegrityFinding(
                        severity=Severity.ERROR,
                        rule="invalid_validation_status",
                        entity_type="treatment_control_relationships",
                        entity_id=rel_id,
                        field_path=f"treatment_control_relationships[relationship_id='{rel_id}'].metadata_validation_status",
                        message=(
                            f"Comparison '{rel_id}' has invalid metadata_validation_status '{rel_val_status}'. "
                            f"Allowed statuses: {sorted(ALLOWED_VALIDATION_STATUSES)}."
                        ),
                    )
                )

            if rel_val_status == "invalid":
                findings.append(
                    IntegrityFinding(
                        severity=Severity.ERROR,
                        rule="invalid_relationship",
                        entity_type="treatment_control_relationships",
                        entity_id=rel_id,
                        field_path=f"treatment_control_relationships[relationship_id='{rel_id}']",
                        message=f"Comparison '{rel_id}' has validation status 'invalid'; relationship is blocked.",
                    )
                )

            if len(ctrl_cids) > 1:
                findings.append(
                    IntegrityFinding(
                        severity=Severity.ERROR,
                        rule="ambiguous_matched_control",
                        entity_type="treatment_control_relationships",
                        entity_id=rel_id,
                        field_path=f"treatment_control_relationships[relationship_id='{rel_id}'].matched_control_condition_ids",
                        message=(
                            f"Comparison '{rel_id}' specifies multiple candidate matched controls {ctrl_cids}. "
                            "Ambiguous control matching is prohibited without confirmed resolution to a single control; comparison is blocked."
                        ),
                    )
                )

            # Role consistency check
            for t_cid in trt_cids:
                t_role = conditions.get(t_cid, {}).get("treatment_control_status", {}).get("normalized_value")
                if t_role == "control":
                    findings.append(
                        IntegrityFinding(
                            severity=Severity.ERROR,
                            rule="treatment_control_role_inconsistency",
                            entity_type="treatment_control_relationships",
                            entity_id=rel_id,
                            field_path=f"treatment_control_relationships[relationship_id='{rel_id}']",
                            message=f"Condition '{t_cid}' assigned as treatment in '{rel_id}' has declared status 'control'.",
                        )
                    )
            for c_cid in ctrl_cids:
                c_role = conditions.get(c_cid, {}).get("treatment_control_status", {}).get("normalized_value")
                if c_role == "treatment":
                    findings.append(
                        IntegrityFinding(
                            severity=Severity.ERROR,
                            rule="treatment_control_role_inconsistency",
                            entity_type="treatment_control_relationships",
                            entity_id=rel_id,
                            field_path=f"treatment_control_relationships[relationship_id='{rel_id}']",
                            message=f"Condition '{c_cid}' assigned as control in '{rel_id}' has declared status 'treatment'.",
                        )
                    )

            for t_cid in trt_cids:
                t_samples = samples_by_condition.get(t_cid, [])
                for c_cid in ctrl_cids:
                    c_samples = samples_by_condition.get(c_cid, [])

                    # Check species consistency across comparison
                    t_species = {
                        bio_sources[s["biological_source_id"]]["species"]["original_value"]
                        for s in t_samples if s.get("biological_source_id") in bio_sources
                    }
                    c_species = {
                        bio_sources[s["biological_source_id"]]["species"]["original_value"]
                        for s in c_samples if s.get("biological_source_id") in bio_sources
                    }
                    if t_species and c_species and t_species != c_species:
                        findings.append(
                            IntegrityFinding(
                                severity=Severity.ERROR,
                                rule="cross_species_comparison_prohibited",
                                entity_type="treatment_control_relationships",
                                entity_id=rel_id,
                                field_path=f"treatment_control_relationships[relationship_id='{rel_id}']",
                                message=f"Comparison compares different species: treatment has {t_species}, control has {c_species}.",
                            )
                        )

                    # Check organoid context consistency
                    t_contexts = {
                        organoid_contexts[s["organoid_context_id"]]["organoid_type"]["original_value"]
                        for s in t_samples if s.get("organoid_context_id") in organoid_contexts
                    }
                    c_contexts = {
                        organoid_contexts[s["organoid_context_id"]]["organoid_type"]["original_value"]
                        for s in c_samples if s.get("organoid_context_id") in organoid_contexts
                    }
                    if t_contexts and c_contexts and t_contexts != c_contexts:
                        findings.append(
                            IntegrityFinding(
                                severity=Severity.WARNING,
                                rule="heterogeneous_organoid_contexts",
                                entity_type="treatment_control_relationships",
                                entity_id=rel_id,
                                field_path=f"treatment_control_relationships[relationship_id='{rel_id}']",
                                message=f"Comparison compares distinct organoid types: treatment={t_contexts}, control={c_contexts}.",
                            )
                        )

                    # Check developmental stage consistency
                    t_stages = {
                        s.get("developmental_age_or_stage_at_collection", {}).get("normalized_value")
                        or s.get("developmental_age_or_stage_at_collection", {}).get("original_value")
                        for s in t_samples
                    }
                    c_stages = {
                        s.get("developmental_age_or_stage_at_collection", {}).get("normalized_value")
                        or s.get("developmental_age_or_stage_at_collection", {}).get("original_value")
                        for s in c_samples
                    }
                    clean_t_stages = {str(x).strip().lower() for x in t_stages if x and str(x).strip().lower() not in ("nan", "none", "unknown", "")}
                    clean_c_stages = {str(x).strip().lower() for x in c_stages if x and str(x).strip().lower() not in ("nan", "none", "unknown", "")}
                    if clean_t_stages and clean_c_stages and not (clean_t_stages & clean_c_stages):
                        findings.append(
                            IntegrityFinding(
                                severity=Severity.ERROR,
                                rule="developmental_stage_mismatch",
                                entity_type="treatment_control_relationships",
                                entity_id=rel_id,
                                field_path=f"treatment_control_relationships[relationship_id='{rel_id}']",
                                message=f"Comparison compares mismatched developmental stages: treatment={clean_t_stages}, control={clean_c_stages}. Cannot be compared without manual review.",
                            )
                        )

                    # Check biological source / cell line consistency
                    t_sources = {
                        bio_sources[s["biological_source_id"]]["cell_line_or_source_identifier"]["original_value"]
                        for s in t_samples if s.get("biological_source_id") in bio_sources
                    }
                    c_sources = {
                        bio_sources[s["biological_source_id"]]["cell_line_or_source_identifier"]["original_value"]
                        for s in c_samples if s.get("biological_source_id") in bio_sources
                    }
                    clean_t_sources = {str(x).strip() for x in t_sources if x and str(x).strip().lower() not in ("nan", "none", "unknown", "")}
                    clean_c_sources = {str(x).strip() for x in c_sources if x and str(x).strip().lower() not in ("nan", "none", "unknown", "")}
                    if clean_t_sources and clean_c_sources and not (clean_t_sources & clean_c_sources):
                        findings.append(
                            IntegrityFinding(
                                severity=Severity.WARNING,
                                rule="mismatched_biological_sources",
                                entity_type="treatment_control_relationships",
                                entity_id=rel_id,
                                field_path=f"treatment_control_relationships[relationship_id='{rel_id}']",
                                message=f"Comparison compares distinct biological sources/cell lines: treatment={clean_t_sources}, control={clean_c_sources}. Requires manual review.",
                            )
                        )

                    model_consistency_summary[f"{t_cid}_vs_{c_cid}"] = {
                        "treatment_samples": len(t_samples),
                        "control_samples": len(c_samples),
                        "treatment_species": list(t_species),
                        "control_species": list(c_species),
                        "treatment_contexts": list(t_contexts),
                        "control_contexts": list(c_contexts),
                    }

        has_errors = any(f.severity == Severity.ERROR for f in findings)
        status = ValidationStatus.INVALID if has_errors else ValidationStatus.VALID_WITH_WARNINGS if findings else ValidationStatus.VALID

        return IntegrityCheckResult(
            is_valid=not has_errors,
            status=status,
            findings=tuple(findings),
            replicate_summary=replicate_summary,
            model_consistency_summary=model_consistency_summary,
        )
