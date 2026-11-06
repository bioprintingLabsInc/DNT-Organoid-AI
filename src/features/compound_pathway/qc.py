"""Stage 8D: Quality Control and Accounting Auditor for Evidence Relationships.

Verifies:
1. Zero DNT Outcome Label Leakage (circulatory prevention).
2. Lossless Gene Universe Conservation (100% of input genes preserved).
3. Exact Numerical Identity with Upstream Stage 8A.1 Statistics.
4. Scientific Causality & Direct Binding Safeguards.
5. Traceable Pathway Evidence Qualification.
"""

from dataclasses import dataclass
import math
from typing import Any

from src.features.compound_pathway.models import (
    ChemicalGenePathwayEvidenceDataset,
    EvidenceType,
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
    findings: tuple[str, ...]

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
            "findings": list(self.findings),
        }


def _approx_eq(v1: float | None, v2: float | None, tol: float = 1e-8) -> bool:
    if v1 is None and v2 is None:
        return True
    if v1 is None or v2 is None:
        return False
    if math.isnan(v1) and math.isnan(v2):
        return True
    return abs(v1 - v2) <= tol


def audit_evidence_dataset_qc(dataset: ChemicalGenePathwayEvidenceDataset) -> EvidenceQCSummary:
    """Execute exhaustive quality control and invariant auditing on an evidence dataset."""
    findings: list[str] = []
    leakage_violations = 0
    numerical_mismatches = 0
    causality_violations = 0
    pathway_annotated = 0
    unannotated = 0

    # Build reference lookup from input chemical_gene_records
    src_map = {cg.canonical_ensembl_gene_id: cg for cg in dataset.chemical_gene_records}
    input_gene_ids = set(src_map.keys())
    observed_gene_ids = set()

    for rec in dataset.records:
        gid = rec.canonical_ensembl_gene_id
        observed_gene_ids.add(gid)

        # 1. Leakage Audit
        text_corpus = f"{rec.scientific_limitations} {rec.provenance_notes or ''} {rec.compound_name}".lower()
        for term in FORBIDDEN_LEAKAGE_TERMS:
            if term in text_corpus:
                leakage_violations += 1
                findings.append(f"Gene {gid}: Detected forbidden DNT supervisory label term '{term}'.")

        # 2. Numerical Identity Audit against Stage 8A.1 Source
        src = src_map.get(gid)
        if src is None:
            numerical_mismatches += 1
            findings.append(f"Gene {gid}: Evidence record not found in input chemical_gene_records.")
        else:
            if not _approx_eq(rec.log2_fold_change, src.log2_fold_change):
                numerical_mismatches += 1
                findings.append(f"Gene {gid}: log2_fold_change mismatch ({rec.log2_fold_change} vs {src.log2_fold_change}).")
            if not _approx_eq(rec.wald_statistic, src.wald_statistic):
                numerical_mismatches += 1
                findings.append(f"Gene {gid}: wald_statistic mismatch ({rec.wald_statistic} vs {src.wald_statistic}).")
            if not _approx_eq(rec.p_value, src.p_value):
                numerical_mismatches += 1
                findings.append(f"Gene {gid}: p_value mismatch ({rec.p_value} vs {src.p_value}).")
            if rec.gene_response_direction != src.response_direction:
                numerical_mismatches += 1
                findings.append(f"Gene {gid}: direction mismatch ({rec.gene_response_direction} vs {src.response_direction}).")

        # 3. Causality & Direct Binding Audit
        if rec.is_direct_target_interaction and rec.evidence_type != EvidenceType.DIRECT_TARGET_INTERACTION:
            causality_violations += 1
            findings.append(f"Gene {gid}: Claimed direct target interaction without DIRECT_TARGET_INTERACTION evidence type.")

        # 4. Pathway Qualification Audit
        if rec.pathway_id is not None:
            pathway_annotated += 1
            if rec.pathway_evidence_status != PathwayEvidenceStatus.AVAILABLE:
                findings.append(f"Gene {gid}: Pathway ID present ({rec.pathway_id}) but status is {rec.pathway_evidence_status}.")
        else:
            unannotated += 1

    # 5. Law of Conservation of the Gene Universe
    missing_genes = input_gene_ids - observed_gene_ids
    if missing_genes:
        numerical_mismatches += len(missing_genes)
        findings.append(f"Law of Conservation of Gene Universe violated: {len(missing_genes)} input genes missing from evidence records.")

    is_valid = bool(
        leakage_violations == 0
        and numerical_mismatches == 0
        and causality_violations == 0
        and len(missing_genes) == 0
    )

    return EvidenceQCSummary(
        is_valid=is_valid,
        total_evidence_records=dataset.total_evidence_records,
        unique_genes_count=len(observed_gene_ids),
        input_genes_count=len(input_gene_ids),
        unique_compounds_count=len(dataset.unique_compounds),
        unique_contrasts_count=len(dataset.unique_contrasts),
        pathway_annotated_records_count=pathway_annotated,
        unannotated_records_count=unannotated,
        leakage_violations_count=leakage_violations,
        numerical_mismatches_count=numerical_mismatches,
        causality_violations_count=causality_violations,
        findings=tuple(findings),
    )
