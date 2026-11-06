"""Stage 8A.1 & Stage 8A.2: GSE187006 Ingestion, limma-trend inference, and cross-study comparison against GSE232218.

Implements:
1. Universal RNA-seq U3.2 ingestion of real GSE187006 log2(TPM) continuous expression matrix.
2. U3.3 limma-trend differential inference via official R/Bioconductor limma for U1M and U2F contrasts.
3. Stage 8A.1 Condition-Level Gene Response Registry generation and QC audits.
4. Stage 8A.2 Cross-study comparability evaluation and descriptive direction concordance against GSE232218.
5. Strict non-comparable context enforcement (zero fold-change pooling, zero invented SE).
"""

import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from src.differential_expression.limma_runner import run_limma_trend_contrast
from src.features.cross_study.comparability import (
    compare_cross_study_gene_responses,
    evaluate_study_comparability,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
    DirectionConcordance,
    EstimateStatus,
    ResponseDirection,
)
from src.features.cross_study.qc import (
    audit_comparison_qc,
    audit_registry_qc,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)
from src.features.cross_study.serialization import (
    export_comparison_to_csv,
    export_registry_to_csv,
)
from src.normalization.provenance_gate import (
    ingest_processed_expression_matrix,
)

R_SCRIPT_PATH = "C:/Users/chari/.r/R-4.4.3/bin/x64/Rscript.exe"


def load_gse232218_registry(csv_path: Path) -> ConditionGeneResponseRegistry:
    """Load previously validated GSE232218 Stage 8A.1 registry from CSV."""
    df = pd.read_csv(csv_path)
    records = []
    for _, row in df.iterrows():
        rec = ConditionGeneResponseRecord(
            canonical_ensembl_gene_id=str(row["canonical_ensembl_gene_id"]),
            gene_symbol=str(row["gene_symbol"]) if pd.notna(row["gene_symbol"]) else None,
            study_id=str(row["study_id"]),
            study_family_id=str(row["study_family_id"]),
            independent_study_id=str(row["independent_study_id"]),
            experiment_id=str(row["experiment_id"]),
            contrast_id=str(row["contrast_id"]),
            treatment_condition_id="cond_vpa200",
            matched_control_condition_ids=("cond_vehicle",),
            compound_id=str(row["compound_id"]),
            compound_name=str(row["compound_name"]),
            dose_value=float(row["dose_value"]) if pd.notna(row["dose_value"]) else None,
            dose_unit=str(row["dose_unit"]) if pd.notna(row["dose_unit"]) else None,
            dose_original=str(row["dose_original"]) if pd.notna(row["dose_original"]) else None,
            exposure_duration=float(row["exposure_duration"]) if pd.notna(row["exposure_duration"]) else None,
            exposure_duration_unit="h",
            developmental_stage=str(row["developmental_stage"]) if pd.notna(row["developmental_stage"]) else None,
            model_type=str(row["model_type"]) if pd.notna(row["model_type"]) else None,
            biological_source_ids=("WA09",),
            organoid_context_ids=("organoid_ctx_01",),
            treatment_replicate_count=4,
            control_replicate_count=4,
            treatment_sample_ids=("Sample_3901-AT-3", "Sample_3901-AT-7", "Sample_3901-AT-12", "Sample_3901-AT-17"),
            control_sample_ids=("Sample_3901-AT-1", "Sample_3901-AT-5", "Sample_3901-AT-10", "Sample_3901-AT-15"),
            treatment_biological_replicate_ids=("rep_trt_3", "rep_trt_7", "rep_trt_12", "rep_trt_17"),
            control_biological_replicate_ids=("rep_ctrl_1", "rep_ctrl_5", "rep_ctrl_10", "rep_ctrl_15"),
            base_mean=float(row["base_mean"]) if pd.notna(row["base_mean"]) else 0.0,
            log2_fold_change=float(row["log2_fold_change"]) if pd.notna(row["log2_fold_change"]) else None,
            log2_fold_change_se=float(row["log2_fold_change_se"]) if pd.notna(row["log2_fold_change_se"]) else None,
            wald_statistic=float(row["wald_statistic"]) if pd.notna(row["wald_statistic"]) else None,
            p_value=float(row["p_value"]) if pd.notna(row["p_value"]) else None,
            adjusted_p_value=float(row["adjusted_p_value"]) if pd.notna(row["adjusted_p_value"]) else None,
            response_direction=ResponseDirection(str(row["response_direction"])),
            estimate_status=EstimateStatus(str(row["estimate_status"])),
            missing_reason=str(row["missing_reason"]) if pd.notna(row["missing_reason"]) else None,
            cross_study_status=CrossStudyStatus(str(row["cross_study_status"])),
            source_contrast_id=str(row["contrast_id"]),
            de_version="DESeq2_v1",
            registry_version="1.0.0",
        )
        records.append(rec)
    return ConditionGeneResponseRegistry(
        registry_id="reg_GSE232218_VPA_200uM_recovered",
        records=tuple(records),
        registry_version="1.0.0",
    )


def run_pipeline() -> None:
    print("=" * 80)
    print("STAGE 8A: GSE187006 INGESTION, LIMMA-TREND, & CROSS-STUDY COMPARISON")
    print("=" * 80)

    # 1. Prepare GSE187006 Metadata
    gse187006_raw = Path("data/raw/GSE187006/GSE187006_Processed_data.csv.gz")
    if not gse187006_raw.exists():
        raise FileNotFoundError(f"GSE187006 raw matrix not found at {gse187006_raw}")

    meta = {
        "dataset": {"dataset_id": "GSE187006", "study_id": "GSE187006"},
        "samples": [
            {"sample_id": "U1M.NC1", "condition_id": "U1M_NC", "biological_replicate_id": "U1M_NC1"},
            {"sample_id": "U1M.NC2", "condition_id": "U1M_NC", "biological_replicate_id": "U1M_NC2"},
            {"sample_id": "U1M.NC3", "condition_id": "U1M_NC", "biological_replicate_id": "U1M_NC3"},
            {"sample_id": "U1M.V1", "condition_id": "U1M_VPA", "biological_replicate_id": "U1M_V1"},
            {"sample_id": "U1M.V2", "condition_id": "U1M_VPA", "biological_replicate_id": "U1M_V2"},
            {"sample_id": "U1M.V3", "condition_id": "U1M_VPA", "biological_replicate_id": "U1M_V3"},
            {"sample_id": "U2F.NC2", "condition_id": "U2F_NC", "biological_replicate_id": "U2F_NC2"},
            {"sample_id": "U2F.NC3", "condition_id": "U2F_NC", "biological_replicate_id": "U2F_NC3"},
            {"sample_id": "U2F.NC4", "condition_id": "U2F_NC", "biological_replicate_id": "U2F_NC4"},
            {"sample_id": "U2F.V2", "condition_id": "U2F_VPA", "biological_replicate_id": "U2F_V2"},
            {"sample_id": "U2F.V3", "condition_id": "U2F_VPA", "biological_replicate_id": "U2F_V3"},
            {"sample_id": "U2F.V4", "condition_id": "U2F_VPA", "biological_replicate_id": "U2F_V4"},
        ],
        "conditions": [
            {
                "condition_id": "U1M_NC",
                "chemical_agent_name": None,
                "concentration_or_dose": 0.0,
                "concentration_or_dose_unit": "mM",
                "exposure_duration": 72.0,
                "exposure_duration_unit": "h",
                "developmental_age_or_stage_at_exposure": "Day 42 forebrain organoids",
                "model_type": "Human Forebrain Organoid (U1M)",
            },
            {
                "condition_id": "U1M_VPA",
                "chemical_agent_name": "Valproic Acid",
                "chemical_agent_id": "CID:3121",
                "concentration_or_dose": 1.0,
                "concentration_or_dose_unit": "mM",
                "exposure_duration": 72.0,
                "exposure_duration_unit": "h",
                "developmental_age_or_stage_at_exposure": "Day 42 forebrain organoids",
                "model_type": "Human Forebrain Organoid (U1M)",
            },
            {
                "condition_id": "U2F_NC",
                "chemical_agent_name": None,
                "concentration_or_dose": 0.0,
                "concentration_or_dose_unit": "mM",
                "exposure_duration": 72.0,
                "exposure_duration_unit": "h",
                "developmental_age_or_stage_at_exposure": "Day 42 forebrain organoids",
                "model_type": "Human Forebrain Organoid (U2F)",
            },
            {
                "condition_id": "U2F_VPA",
                "chemical_agent_name": "Valproic Acid",
                "chemical_agent_id": "CID:3121",
                "concentration_or_dose": 1.0,
                "concentration_or_dose_unit": "mM",
                "exposure_duration": 72.0,
                "exposure_duration_unit": "h",
                "developmental_age_or_stage_at_exposure": "Day 42 forebrain organoids",
                "model_type": "Human Forebrain Organoid (U2F)",
            },
        ],
        "relationships": [
            {
                "relationship_id": "GSE187006_U1M_VPA_1mM_72h_vs_NC",
                "relationship_type": "treatment_control",
                "treatment_condition_id": "U1M_VPA",
                "matched_control_condition_ids": ["U1M_NC"],
            },
            {
                "relationship_id": "GSE187006_U2F_VPA_1mM_72h_vs_NC",
                "relationship_type": "treatment_control",
                "treatment_condition_id": "U2F_VPA",
                "matched_control_condition_ids": ["U2F_NC"],
            },
        ],
        "processing_provenance": {
            "expression_unit": "TPM",
            "transformation_scale": "log2",
            "normalization_method": "tpm",
            "software_pipeline": "RSEM (1.3.3) / STAR (v2.7.6a)",
            "genome_build": "GRCh38.p13",
        },
    }

    # 2. Ingest via U3.2 Provenance Gate
    print("\n1. Ingesting GSE187006 via Universal RNA-seq U3.2 Provenance Gate...")
    t0 = time.time()
    ds = ingest_processed_expression_matrix(gse187006_raw, meta)
    print(f"   Dataset status:             {ds.status}")
    print(f"   Evaluability (U3.3):        {ds.is_u3_3_evaluable}")
    print(f"   Statistical eligibility:    {ds.is_statistically_eligible}")
    print(f"   Total input genes:          {ds.total_input_genes_count:,}")
    print(f"   Canonical eligible genes:   {ds.eligible_canonical_genes_count:,}")
    print(f"   Excluded genes:             {ds.excluded_genes_count:,}")
    print(f"   Universe conserved:         {ds.is_gene_universe_conserved}")

    if not ds.is_valid_intake or not ds.is_statistically_eligible:
        raise RuntimeError("GSE187006 intake failed or is not statistically eligible.")

    # 3. Differential Inference via U3.3 limma-trend for U1M
    print("\n2. Executing U3.3 limma-trend for Contrast 1: U1M VPA 1 mM vs NC...")
    t1 = time.time()
    res_u1m = run_limma_trend_contrast(
        dataset=ds,
        metadata=meta,
        contrast_id="GSE187006_U1M_VPA_1mM_72h_vs_NC",
        r_binary_path=R_SCRIPT_PATH,
    )
    print(f"   U1M Contrast status:        {res_u1m.status}")
    print(f"   Canonical genes evaluated:  {len(res_u1m.gene_results):,}")
    sig_u1m = [g for g in res_u1m.gene_results if g.adjusted_p_value_bh is not None and g.adjusted_p_value_bh < 0.05]
    print(f"   Genes with FDR < 0.05:      {len(sig_u1m):,}")

    # 4. Differential Inference via U3.3 limma-trend for U2F
    print("\n3. Executing U3.3 limma-trend for Contrast 2: U2F VPA 1 mM vs NC...")
    res_u2f = run_limma_trend_contrast(
        dataset=ds,
        metadata=meta,
        contrast_id="GSE187006_U2F_VPA_1mM_72h_vs_NC",
        r_binary_path=R_SCRIPT_PATH,
    )
    print(f"   U2F Contrast status:        {res_u2f.status}")
    print(f"   Canonical genes evaluated:  {len(res_u2f.gene_results):,}")
    sig_u2f = [g for g in res_u2f.gene_results if g.adjusted_p_value_bh is not None and g.adjusted_p_value_bh < 0.05]
    print(f"   Genes with FDR < 0.05:      {len(sig_u2f):,}")

    # 5. Build Condition Gene Response Registries (Stage 8A.1)
    print("\n4. Building Stage 8A.1 Registries for GSE187006...")
    model_map = {
        "GSE187006_U1M_VPA_1mM_72h_vs_NC": "Human Forebrain Organoid (U1M)",
        "GSE187006_U2F_VPA_1mM_72h_vs_NC": "Human Forebrain Organoid (U2F)",
    }

    reg_u1m = build_condition_gene_response_registry(
        contrasts=[res_u1m],
        registry_id="reg_GSE187006_U1M_VPA_1mM_v1",
        registry_version="1.0.0",
        model_type_map=model_map,
    )
    qc_u1m = audit_registry_qc(reg_u1m, [res_u1m])
    print(f"   U1M Registry QC valid:      {qc_u1m.is_valid} ({reg_u1m.total_records:,} records)")
    export_registry_to_csv(reg_u1m, "data/processed/GSE187006_U1M_VPA_1mM_condition_gene_response_registry.csv")

    reg_u2f = build_condition_gene_response_registry(
        contrasts=[res_u2f],
        registry_id="reg_GSE187006_U2F_VPA_1mM_v1",
        registry_version="1.0.0",
        model_type_map=model_map,
    )
    qc_u2f = audit_registry_qc(reg_u2f, [res_u2f])
    print(f"   U2F Registry QC valid:      {qc_u2f.is_valid} ({reg_u2f.total_records:,} records)")
    export_registry_to_csv(reg_u2f, "data/processed/GSE187006_U2F_VPA_1mM_condition_gene_response_registry.csv")

    reg_combined = build_condition_gene_response_registry(
        contrasts=[res_u1m, res_u2f],
        registry_id="reg_GSE187006_VPA_1mM_combined_v1",
        registry_version="1.0.0",
        model_type_map=model_map,
    )
    qc_combined = audit_registry_qc(reg_combined, [res_u1m, res_u2f])
    print(f"   Combined Registry QC valid: {qc_combined.is_valid} ({reg_combined.total_records:,} records)")
    export_registry_to_csv(reg_combined, "data/processed/GSE187006_VPA_1mM_condition_gene_response_registry.csv")

    # 6. Load GSE232218 Stage 8A.1 Baseline
    print("\n5. Loading GSE232218 Stage 8A.1 Baseline...")
    gse232218_csv = Path("data/processed/GSE232218_VPA_200uM_condition_gene_response_registry.csv")
    reg_gse232218 = load_gse232218_registry(gse232218_csv)
    print(f"   GSE232218 loaded:           {reg_gse232218.total_records:,} records")

    # 7. Stage 8A.2 Comparability Evaluation: GSE232218 vs GSE187006 U1M
    print("\n6. Executing Stage 8A.2 Cross-Study Comparability: GSE232218 vs GSE187006 U1M...")
    gse232218_cid = list(reg_gse232218.unique_contrasts)[0]
    cmp_u1m = compare_cross_study_gene_responses(
        registry1=reg_gse232218,
        registry2=reg_u1m,
        contrast1_id=gse232218_cid,
        contrast2_id="GSE187006_U1M_VPA_1mM_72h_vs_NC",
        comparison_id="cmp_GSE232218_vs_GSE187006_U1M",
    )
    print(f"   Synthesis Justified:        {cmp_u1m.assessment.is_synthesis_justified}")
    print(f"   Cross-Study Status:         {cmp_u1m.assessment.cross_study_status}")
    print(f"   Differences:                {cmp_u1m.assessment.differences}")
    print(f"   Total genes in comparison:  {cmp_u1m.total_genes_evaluated:,}")
    print(f"   Shared genes:               {cmp_u1m.shared_genes_count:,}")
    print(f"   Concordant UP:              {cmp_u1m.concordant_up_count:,}")
    print(f"   Concordant DOWN:            {cmp_u1m.concordant_down_count:,}")
    print(f"   Discordant:                 {cmp_u1m.discordant_count:,}")
    print(f"   Concordance rate:           {cmp_u1m.concordance_rate:.4f} ({(cmp_u1m.concordance_rate*100):.2f}%)")

    qc_cmp_u1m = audit_comparison_qc(cmp_u1m, reg_gse232218, reg_u1m)
    print(f"   Comparison QC valid:        {qc_cmp_u1m.is_valid}")
    export_comparison_to_csv(cmp_u1m, "data/processed/GSE232218_vs_GSE187006_U1M_cross_study_comparison.csv")

    # 8. Stage 8A.2 Comparability Evaluation: GSE232218 vs GSE187006 U2F
    print("\n7. Executing Stage 8A.2 Cross-Study Comparability: GSE232218 vs GSE187006 U2F...")
    cmp_u2f = compare_cross_study_gene_responses(
        registry1=reg_gse232218,
        registry2=reg_u2f,
        contrast1_id=gse232218_cid,
        contrast2_id="GSE187006_U2F_VPA_1mM_72h_vs_NC",
        comparison_id="cmp_GSE232218_vs_GSE187006_U2F",
    )
    print(f"   Synthesis Justified:        {cmp_u2f.assessment.is_synthesis_justified}")
    print(f"   Cross-Study Status:         {cmp_u2f.assessment.cross_study_status}")
    print(f"   Differences:                {cmp_u2f.assessment.differences}")
    print(f"   Total genes in comparison:  {cmp_u2f.total_genes_evaluated:,}")
    print(f"   Shared genes:               {cmp_u2f.shared_genes_count:,}")
    print(f"   Concordant UP:              {cmp_u2f.concordant_up_count:,}")
    print(f"   Concordant DOWN:            {cmp_u2f.concordant_down_count:,}")
    print(f"   Discordant:                 {cmp_u2f.discordant_count:,}")
    print(f"   Concordance rate:           {cmp_u2f.concordance_rate:.4f} ({(cmp_u2f.concordance_rate*100):.2f}%)")

    qc_cmp_u2f = audit_comparison_qc(cmp_u2f, reg_gse232218, reg_u2f)
    print(f"   Comparison QC valid:        {qc_cmp_u2f.is_valid}")
    export_comparison_to_csv(cmp_u2f, "data/processed/GSE232218_vs_GSE187006_U2F_cross_study_comparison.csv")

    # Summary of SHA256 checksums
    print("\n8. Generating SHA-256 Checksums of Processed Artifacts...")
    artifacts = [
        "data/processed/GSE187006_U1M_VPA_1mM_condition_gene_response_registry.csv",
        "data/processed/GSE187006_U2F_VPA_1mM_condition_gene_response_registry.csv",
        "data/processed/GSE187006_VPA_1mM_condition_gene_response_registry.csv",
        "data/processed/GSE232218_vs_GSE187006_U1M_cross_study_comparison.csv",
        "data/processed/GSE232218_vs_GSE187006_U2F_cross_study_comparison.csv",
    ]
    checksums = {}
    for a in artifacts:
        p = Path(a)
        if p.exists():
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            checksums[a] = h
            print(f"   {p.name}: {h}")

    with open("data/processed/stage8a2_artifacts_sha256.json", "w", encoding="utf-8") as f:
        json.dump(checksums, f, indent=2)

    print("\nStage 8A execution completed successfully.")


if __name__ == "__main__":
    run_pipeline()
