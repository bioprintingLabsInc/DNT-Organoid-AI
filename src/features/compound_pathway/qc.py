"""Stage 8D: Strengthened Quality Control and Accounting Auditor for Evidence Relationships.

Verifies:
1. Zero DNT Outcome Label Leakage (circularity prevention).
2. Lossless Gene Universe Conservation per Study, Contrast, and Gene.
3. Exact Numerical Parity and Metadata Fidelity per Study, Contrast, and Gene.
4. Duplicate Evidence Identifier Detection.
5. Inconsistent Pathway Status Detection.
6. Unexpected Cross-Contrast Relationship Detection.
7. Scientific Causality & Direct Binding Safeguards.
"""

from dataclasses import dataclass
import math
from typing import Any

from src.features.compound_pathway.models import (
    ChemicalGenePathwayEvidenceDataset,
    EvidenceType,
    GenePathwayMembershipType,
    PathwayEvidenceStatus,
)


FORBIDDEN_LEAKAGE_TERMS = {
    "dnt_positive",
    "dnt_negative",
    "regulatory_positive",
    "regulatory_negative",
    "gold_standard_positive",
    "gold_standard_negative",
    "neurotoxic_label",
    "non_neurotoxic_label",
    "training_class",
    "test_class",
}


@dataclass(frozen=True, slots=True)
class EvidenceQCSummary:
    """Comprehensive QC outcome for a Stage 8D evidence dataset."""

    is_valid: bool
    total_evidence_records: int
    unique_genes_count: int
    input_genes_count: int
    unique_compounds_count: int
    unique_contrasts_count: int
    pathway_annotated_records_count: int
    unannotated_records_count: int
    leakage_violations_count: int
    numerical_mismatches_count: int
    causality_violations_count: int
    duplicate_ids_count: int
    cross_contrast_violations_count: int
    pathway_status_violations_count: int
    findings: tuple[str, ...]

    @property
    def notes(self) -> tuple[str, ...]:
        return self.findings

    def to_dict(self) -> dict[str, Any]:
        return {
            "is_valid": self.is_valid,
            "total_evidence_records": self.total_evidence_records,
            "unique_genes_count": self.unique_genes_count,
            "input_genes_count": self.input_genes_count,
            "unique_compounds_count": self.unique_compounds_count,
            "unique_contrasts_count": self.unique_contrasts_count,
            "pathway_annotated_records_count": self.pathway_annotated_records_count,
            "unannotated_records_count": self.unannotated_records_count,
            "leakage_violations_count": self.leakage_violations_count,
            "numerical_mismatches_count": self.numerical_mismatches_count,
            "causality_violations_count": self.causality_violations_count,
            "duplicate_ids_count": self.duplicate_ids_count,
            "cross_contrast_violations_count": self.cross_contrast_violations_count,
            "pathway_status_violations_count": self.pathway_status_violations_count,
            "findings": list(self.findings),
        }


def _approx_eq(v1: float | None, v2: float | None, tol: float = 1e-7) -> bool:
    if v1 is None and v2 is None:
        return True
    if v1 is None or v2 is None:
        return False
    if math.isnan(v1) and math.isnan(v2):
        return True
    return abs(v1 - v2) <= tol


def audit_evidence_dataset_qc(dataset: ChemicalGenePathwayEvidenceDataset) -> EvidenceQCSummary:
    """Execute exhaustive quality control and invariant auditing per study, contrast, and gene."""
    findings: list[str] = []
    leakage_violations = 0
    numerical_mismatches = 0
    causality_violations = 0
    duplicate_ids = 0
    cross_contrast_violations = 0
    pathway_status_violations = 0
    pathway_annotated = 0
    unannotated = 0

    # 1. Index input chemical_gene_records by (study_id, contrast_id, canonical_ensembl_gene_id)
    src_map = {
        (cg.exposure.study_id, cg.exposure.contrast_id, cg.canonical_ensembl_gene_id): cg
        for cg in dataset.chemical_gene_records
    }
    input_genes_by_contrast: dict[tuple[str, str], set[str]] = {}
    for cg in dataset.chemical_gene_records:
        key = (cg.exposure.study_id, cg.exposure.contrast_id)
        input_genes_by_contrast.setdefault(key, set()).add(cg.canonical_ensembl_gene_id)

    # Index input chemical_pathway_records by (contrast_id, pathway_id)
    pathway_ref_map = {
        (pw.exposure.contrast_id, pw.pathway_id): pw
        for pw in dataset.chemical_pathway_records
    }

    seen_evidence_ids: set[str] = set()
    observed_genes_by_contrast: dict[tuple[str, str], set[str]] = {}

    for rec in dataset.records:
        key = (rec.study_id, rec.contrast_id)
        observed_genes_by_contrast.setdefault(key, set()).add(rec.canonical_ensembl_gene_id)

        # 2. Duplicate Evidence ID Check
        if rec.evidence_id in seen_evidence_ids:
            duplicate_ids += 1
            findings.append(f"Duplicate evidence_id '{rec.evidence_id}' detected for {rec.canonical_ensembl_gene_id}.")
        else:
            seen_evidence_ids.add(rec.evidence_id)

        # 3. Leakage Audit
        text_corpus = f"{rec.scientific_limitations} {rec.provenance_notes or ''} {rec.compound_name}".lower()
        for term in FORBIDDEN_LEAKAGE_TERMS:
            if term in text_corpus:
                leakage_violations += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: Detected forbidden DNT supervisory term '{term}'.")

        # 4. Numerical Fidelity & Metadata Reconstruction Audit against Source
        src_key = (rec.study_id, rec.contrast_id, rec.canonical_ensembl_gene_id)
        src = src_map.get(src_key)
        if src is None:
            numerical_mismatches += 1
            findings.append(f"Gene {rec.canonical_ensembl_gene_id} in contrast {rec.contrast_id} not found in input chemical_gene_records.")
        else:
            if not _approx_eq(rec.log2_fold_change, src.log2_fold_change):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: log2FC mismatch ({rec.log2_fold_change} vs {src.log2_fold_change}).")
            if not _approx_eq(rec.log2_fold_change_se, src.log2_fold_change_se):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: SE mismatch ({rec.log2_fold_change_se} vs {src.log2_fold_change_se}).")
            if not _approx_eq(rec.wald_statistic, src.wald_statistic):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: statistic mismatch ({rec.wald_statistic} vs {src.wald_statistic}).")
            if not _approx_eq(rec.p_value, src.p_value):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: p_value mismatch ({rec.p_value} vs {src.p_value}).")
            if not _approx_eq(rec.adjusted_p_value, src.adjusted_p_value):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: FDR mismatch ({rec.adjusted_p_value} vs {src.adjusted_p_value}).")
            if not _approx_eq(rec.base_mean, src.base_mean):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: base_mean mismatch ({rec.base_mean} vs {src.base_mean}).")
            if rec.gene_response_direction != src.response_direction:
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: direction mismatch ({rec.gene_response_direction} vs {src.response_direction}).")
            if rec.gene_estimate_status != src.estimate_status:
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}, Gene {rec.canonical_ensembl_gene_id}: estimate_status mismatch ({rec.gene_estimate_status} vs {src.estimate_status}).")

            # Metadata reconstruction check
            if rec.treatment_condition_id != src.exposure.treatment_condition_id:
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}: treatment_condition_id mismatch ({rec.treatment_condition_id} vs {src.exposure.treatment_condition_id}).")
            if rec.compound_id != src.exposure.compound_id:
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}: compound_id mismatch ({rec.compound_id} vs {src.exposure.compound_id}).")
            if not _approx_eq(rec.dose_value, src.exposure.dose_value):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}: dose_value mismatch ({rec.dose_value} vs {src.exposure.dose_value}).")
            if not _approx_eq(rec.exposure_duration, src.exposure.exposure_duration):
                numerical_mismatches += 1
                findings.append(f"Contrast {rec.contrast_id}: exposure_duration mismatch ({rec.exposure_duration} vs {src.exposure.exposure_duration}).")

        # 5. Cross-Contrast Integrity Audit
        if rec.pathway_id is not None:
            pw_key = (rec.contrast_id, rec.pathway_id)
            if pw_key not in pathway_ref_map:
                cross_contrast_violations += 1
                findings.append(f"Contrast {rec.contrast_id}: Pathway '{rec.pathway_id}' does not belong to this contrast.")

        # 6. Causality & Direct Binding Audit
        if rec.is_direct_target_interaction and rec.evidence_type != EvidenceType.DIRECT_TARGET_INTERACTION:
            causality_violations += 1
            findings.append(f"Gene {rec.canonical_ensembl_gene_id}: Direct target interaction claimed without verified assay evidence.")

        # 7. Pathway Status Consistency Audit
        if rec.pathway_id is not None:
            pathway_annotated += 1
            if rec.pathway_evidence_status != PathwayEvidenceStatus.AVAILABLE:
                pathway_status_violations += 1
                findings.append(f"Gene {rec.canonical_ensembl_gene_id}: Pathway ID present ({rec.pathway_id}) but status is {rec.pathway_evidence_status}.")
        else:
            unannotated += 1
            if rec.pathway_evidence_status == PathwayEvidenceStatus.AVAILABLE:
                pathway_status_violations += 1
                findings.append(f"Gene {rec.canonical_ensembl_gene_id}: Unannotated record has status AVAILABLE.")

    # 8. Law of Conservation of the Gene Universe per Study & Contrast
    for contrast_key, in_genes in input_genes_by_contrast.items():
        obs_genes = observed_genes_by_contrast.get(contrast_key, set())
        missing = in_genes - obs_genes
        if missing:
            numerical_mismatches += len(missing)
            findings.append(f"Study/Contrast {contrast_key}: {len(missing)} input genes missing from evidence records.")
        extra = obs_genes - in_genes
        if extra:
            numerical_mismatches += len(extra)
            findings.append(f"Study/Contrast {contrast_key}: {len(extra)} unexpected genes present in evidence records.")

    is_valid = bool(
        leakage_violations == 0
        and numerical_mismatches == 0
        and causality_violations == 0
        and duplicate_ids == 0
        and cross_contrast_violations == 0
        and pathway_status_violations == 0
    )

    all_observed_genes = {rec.canonical_ensembl_gene_id for rec in dataset.records}
    all_input_genes = {cg.canonical_ensembl_gene_id for cg in dataset.chemical_gene_records}

    return EvidenceQCSummary(
        is_valid=is_valid,
        total_evidence_records=dataset.total_evidence_records,
        unique_genes_count=len(all_observed_genes),
        input_genes_count=len(all_input_genes),
        unique_compounds_count=len(dataset.unique_compounds),
        unique_contrasts_count=len(dataset.unique_contrasts),
        pathway_annotated_records_count=pathway_annotated,
        unannotated_records_count=unannotated,
        leakage_violations_count=leakage_violations,
        numerical_mismatches_count=numerical_mismatches,
        causality_violations_count=causality_violations,
        duplicate_ids_count=duplicate_ids,
        cross_contrast_violations_count=cross_contrast_violations,
        pathway_status_violations_count=pathway_status_violations,
        findings=tuple(findings),
    )
