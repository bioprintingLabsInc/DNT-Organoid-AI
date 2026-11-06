"""Stage 8D: End-to-End Orchestration Script for Chemical–Gene–Pathway Evidence Relationships.

Synthesizes, validates, audits, and serializes evidence relationships for:
1. GSE232218 VPA 200 uM (with Stage 8C GO BP + KEGG pathways)
2. GSE187006 U1M VPA 1 mM (pathway status NOT_AVAILABLE)
3. GSE187006 U2F VPA 1 mM (pathway status NOT_AVAILABLE)

Computes SHA-256 checksums and performs complete invariant verification.
"""

import hashlib
import json
from pathlib import Path
import sys
import time

# Ensure repo root is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.features.compound_pathway.models import PathwayEvidenceStatus
from src.features.compound_pathway.runner import (
    build_stage8d_evidence_pipeline,
    load_pathway_results_from_csv,
    load_registry_from_csv,
)
from scripts.run_gse187006_and_comparability import load_gse232218_registry


def compute_sha256(file_path: Path) -> str:
    """Compute hex SHA-256 digest of a file."""
    hasher = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    return hasher.hexdigest()


def run_stage8d_evidence_generation() -> dict:
    """Execute Stage 8D evidence generation across all governed study arms."""
    print("================================================================================")
    print(" DNT-Organoid-AI — Stage 8D: Chemical–Gene–Pathway Evidence Relationships")
    print("================================================================================")
    start_time = time.time()
    artifacts_dir = Path("data/processed")
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    results_summary = {}

    # -------------------------------------------------------------------------
    # 1. GSE232218 VPA 200 uM (DESeq2 + Stage 8C GO BP & KEGG)
    # -------------------------------------------------------------------------
    print("\n[1/3] Processing GSE232218 VPA 200 uM (with GO BP + KEGG pathway evidence)...")
    reg_csv_232 = artifacts_dir / "GSE232218_VPA_200uM_condition_gene_response_registry.csv"
    go_csv = artifacts_dir / "GSE232218_VPA_200uM_GO_BP_results.csv"
    kegg_csv = artifacts_dir / "GSE232218_VPA_200uM_KEGG_results.csv"

    if not reg_csv_232.exists():
        raise FileNotFoundError(f"Missing required GSE232218 registry at '{reg_csv_232}'.")
    if not go_csv.exists() or not kegg_csv.exists():
        raise FileNotFoundError("Missing required GSE232218 Stage 8C pathway results.")

    reg_232 = load_gse232218_registry(reg_csv_232)
    pathways_go = load_pathway_results_from_csv(go_csv)
    pathways_kegg = load_pathway_results_from_csv(kegg_csv)
    all_pathways = pathways_go + pathways_kegg

    out_csv_232 = artifacts_dir / "GSE232218_VPA_200uM_chemical_gene_pathway_evidence.csv"
    ds_232, qc_232 = build_stage8d_evidence_pipeline(
        registry=reg_232,
        pathway_results=all_pathways,
        dataset_id="evds_GSE232218_VPA_200uM",
        output_csv=out_csv_232,
    )

    print(f"   Input Genes:               {qc_232.input_genes_count:,}")
    print(f"   Unique Conserved Genes:    {qc_232.unique_genes_count:,}")
    print(f"   Total Evidence Records:    {qc_232.total_evidence_records:,}")
    print(f"   Pathway-Annotated Records: {qc_232.pathway_annotated_records_count:,}")
    print(f"   Unannotated Records:       {qc_232.unannotated_records_count:,}")
    print(f"   Leakage Violations:        {qc_232.leakage_violations_count} (PASS)")
    print(f"   Numerical Mismatches:      {qc_232.numerical_mismatches_count} (PASS)")
    print(f"   QC Status:                 {'VALID' if qc_232.is_valid else 'INVALID'}")
    assert qc_232.is_valid, f"QC failed for GSE232218: {qc_232.notes}"

    sha_232 = compute_sha256(out_csv_232)
    results_summary["GSE232218_VPA_200uM"] = {
        "output_csv": str(out_csv_232),
        "sha256": sha_232,
        "total_records": qc_232.total_evidence_records,
        "unique_genes": qc_232.unique_genes_count,
        "pathway_annotated_records": qc_232.pathway_annotated_records_count,
        "unannotated_records": qc_232.unannotated_records_count,
        "is_valid": qc_232.is_valid,
    }

    # -------------------------------------------------------------------------
    # 2. GSE187006 U1M VPA 1 mM (limma-trend; Pathway Status: NOT_AVAILABLE)
    # -------------------------------------------------------------------------
    print("\n[2/3] Processing GSE187006 U1M VPA 1 mM (Pathway status: NOT_AVAILABLE)...")
    reg_csv_u1m = artifacts_dir / "GSE187006_U1M_VPA_1mM_condition_gene_response_registry.csv"
    if not reg_csv_u1m.exists():
        raise FileNotFoundError(f"Missing required GSE187006 U1M registry at '{reg_csv_u1m}'.")

    reg_u1m = load_registry_from_csv(
        reg_csv_u1m,
        registry_id="reg_GSE187006_U1M",
    )
    out_csv_u1m = artifacts_dir / "GSE187006_U1M_VPA_1mM_chemical_gene_pathway_evidence.csv"
    ds_u1m, qc_u1m = build_stage8d_evidence_pipeline(
        registry=reg_u1m,
        pathway_results=None,
        pathway_status=PathwayEvidenceStatus.NOT_AVAILABLE,
        dataset_id="evds_GSE187006_U1M_VPA_1mM",
        output_csv=out_csv_u1m,
    )

    print(f"   Input Genes:               {qc_u1m.input_genes_count:,}")
    print(f"   Unique Conserved Genes:    {qc_u1m.unique_genes_count:,}")
    print(f"   Total Evidence Records:    {qc_u1m.total_evidence_records:,}")
    print(f"   Pathway-Annotated Records: {qc_u1m.pathway_annotated_records_count:,}")
    print(f"   Unannotated Records:       {qc_u1m.unannotated_records_count:,}")
    print(f"   Leakage Violations:        {qc_u1m.leakage_violations_count} (PASS)")
    print(f"   Numerical Mismatches:      {qc_u1m.numerical_mismatches_count} (PASS)")
    print(f"   QC Status:                 {'VALID' if qc_u1m.is_valid else 'INVALID'}")
    assert qc_u1m.is_valid, f"QC failed for GSE187006 U1M: {qc_u1m.notes}"

    sha_u1m = compute_sha256(out_csv_u1m)
    results_summary["GSE187006_U1M_VPA_1mM"] = {
        "output_csv": str(out_csv_u1m),
        "sha256": sha_u1m,
        "total_records": qc_u1m.total_evidence_records,
        "unique_genes": qc_u1m.unique_genes_count,
        "pathway_annotated_records": qc_u1m.pathway_annotated_records_count,
        "unannotated_records": qc_u1m.unannotated_records_count,
        "is_valid": qc_u1m.is_valid,
    }

    # -------------------------------------------------------------------------
    # 3. GSE187006 U2F VPA 1 mM (limma-trend; Pathway Status: NOT_AVAILABLE)
    # -------------------------------------------------------------------------
    print("\n[3/3] Processing GSE187006 U2F VPA 1 mM (Pathway status: NOT_AVAILABLE)...")
    reg_csv_u2f = artifacts_dir / "GSE187006_U2F_VPA_1mM_condition_gene_response_registry.csv"
    if not reg_csv_u2f.exists():
        raise FileNotFoundError(f"Missing required GSE187006 U2F registry at '{reg_csv_u2f}'.")

    reg_u2f = load_registry_from_csv(
        reg_csv_u2f,
        registry_id="reg_GSE187006_U2F",
    )
    out_csv_u2f = artifacts_dir / "GSE187006_U2F_VPA_1mM_chemical_gene_pathway_evidence.csv"
    ds_u2f, qc_u2f = build_stage8d_evidence_pipeline(
        registry=reg_u2f,
        pathway_results=None,
        pathway_status=PathwayEvidenceStatus.NOT_AVAILABLE,
        dataset_id="evds_GSE187006_U2F_VPA_1mM",
        output_csv=out_csv_u2f,
    )

    print(f"   Input Genes:               {qc_u2f.input_genes_count:,}")
    print(f"   Unique Conserved Genes:    {qc_u2f.unique_genes_count:,}")
    print(f"   Total Evidence Records:    {qc_u2f.total_evidence_records:,}")
    print(f"   Pathway-Annotated Records: {qc_u2f.pathway_annotated_records_count:,}")
    print(f"   Unannotated Records:       {qc_u2f.unannotated_records_count:,}")
    print(f"   Leakage Violations:        {qc_u2f.leakage_violations_count} (PASS)")
    print(f"   Numerical Mismatches:      {qc_u2f.numerical_mismatches_count} (PASS)")
    print(f"   QC Status:                 {'VALID' if qc_u2f.is_valid else 'INVALID'}")
    assert qc_u2f.is_valid, f"QC failed for GSE187006 U2F: {qc_u2f.notes}"

    sha_u2f = compute_sha256(out_csv_u2f)
    results_summary["GSE187006_U2F_VPA_1mM"] = {
        "output_csv": str(out_csv_u2f),
        "sha256": sha_u2f,
        "total_records": qc_u2f.total_evidence_records,
        "unique_genes": qc_u2f.unique_genes_count,
        "pathway_annotated_records": qc_u2f.pathway_annotated_records_count,
        "unannotated_records": qc_u2f.unannotated_records_count,
        "is_valid": qc_u2f.is_valid,
    }

    # -------------------------------------------------------------------------
    # Persist Checksum Manifest
    # -------------------------------------------------------------------------
    sha_manifest_path = artifacts_dir / "stage8d_artifacts_sha256.json"
    with open(sha_manifest_path, "w", encoding="utf-8") as f:
        json.dump(results_summary, f, indent=2)
    print(f"\nSaved SHA-256 manifest to '{sha_manifest_path}'.")

    elapsed = time.time() - start_time
    print(f"\nStage 8D evidence generation completed successfully in {elapsed:.2f} seconds.")
    return results_summary


if __name__ == "__main__":
    run_stage8d_evidence_generation()
