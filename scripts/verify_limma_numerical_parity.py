"""Script to independently verify limma-trend numerical parity on the harmonized gene universe.

Compares pipeline outputs with direct independent R execution on the exact same
canonical Ensembl gene universe and sample set for both U1M and U2F contrasts.
"""

from pathlib import Path
import subprocess
import json
import math
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import numpy as np

from src.normalization.provenance_gate import ingest_processed_expression_matrix
from src.differential_expression.limma_runner import run_limma_trend_contrast

R_BINARY = "C:/Users/chari/.r/R-4.4.3/bin/x64/Rscript.exe"

def main():
    print("=" * 80)
    print("INDEPENDENT LIMMA NUMERICAL PARITY VERIFICATION")
    print("=" * 80)

    gse187006_raw = Path("data/raw/GSE187006/GSE187006_Processed_data.csv.gz")
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
        "treatment_conditions": [
            {
                "condition_id": "U1M_VPA",
                "chemical_agent_id": "CID:3121",
                "agent_name": "Valproic Acid",
                "concentration_or_dose": 1000.0,
                "concentration_or_dose_unit": "uM",
                "exposure_duration": 72.0,
                "exposure_duration_unit": "h",
                "developmental_age_or_stage_at_exposure": "Day 42",
            },
            {
                "condition_id": "U2F_VPA",
                "chemical_agent_id": "CID:3121",
                "agent_name": "Valproic Acid",
                "concentration_or_dose": 1000.0,
                "concentration_or_dose_unit": "uM",
                "exposure_duration": 72.0,
                "exposure_duration_unit": "h",
                "developmental_age_or_stage_at_exposure": "Day 42",
            },
        ],
        "relationships": [
            {
                "relationship_id": "GSE187006_U1M_VPA_1mM_72h_vs_NC",
                "treatment_condition_id": "U1M_VPA",
                "control_condition_id": "U1M_NC",
                "comparison_type": "treatment_vs_vehicle_control",
                "is_active_for_analysis": True,
            },
            {
                "relationship_id": "GSE187006_U2F_VPA_1mM_72h_vs_NC",
                "treatment_condition_id": "U2F_VPA",
                "control_condition_id": "U2F_NC",
                "comparison_type": "treatment_vs_vehicle_control",
                "is_active_for_analysis": True,
            },
        ],
    }

    # 1. Ingest via U3.2
    ds = ingest_processed_expression_matrix(gse187006_raw, meta)
    print(f"Ingested {len(ds.matrix.gene_ids)} raw genes across {len(ds.matrix.sample_ids)} samples.")

    # 2. Extract eligible canonical genes exactly as in limma_runner.py
    eligible_genes = []
    seen_canonical_ids = set()
    for g in ds.harmonized_genes:
        if str(g.mapping_status) == "UNIQUELY_MAPPED" and g.canonical_gene_id:
            can_id = g.canonical_gene_id
            if can_id in seen_canonical_ids:
                continue
            seen_canonical_ids.add(can_id)
            eligible_genes.append(g)

    print(f"Eligible canonical genes: {len(eligible_genes)}")
    gene_row_indices = [g.source_index - 1 for g in eligible_genes]
    canonical_ids = [g.canonical_gene_id for g in eligible_genes]

    # Load existing registry outputs
    reg_u1m = pd.read_csv("data/processed/GSE187006_U1M_VPA_1mM_condition_gene_response_registry.csv")
    reg_u2f = pd.read_csv("data/processed/GSE187006_U2F_VPA_1mM_condition_gene_response_registry.csv")

    results = {}

    for arm_name, reg_df, trt_samples, ctrl_samples in [
        ("U1M", reg_u1m, ["U1M.V1", "U1M.V2", "U1M.V3"], ["U1M.NC1", "U1M.NC2", "U1M.NC3"]),
        ("U2F", reg_u2f, ["U2F.V2", "U2F.V3", "U2F.V4"], ["U2F.NC2", "U2F.NC3", "U2F.NC4"]),
    ]:
        print(f"\nEvaluating arm: {arm_name}")
        sample_subset = trt_samples + ctrl_samples
        sample_indices = [ds.matrix.sample_ids.index(s) for s in sample_subset]

        # Extract the exact sub-matrix rows for eligible canonical genes
        matrix_data = {}
        for s in sample_subset:
            s_idx = ds.matrix.sample_ids.index(s)
            full_col = ds.matrix.columns[s_idx]
            matrix_data[s] = [full_col[r_i] for r_i in gene_row_indices]

        sub_df = pd.DataFrame(matrix_data, index=canonical_ids)
        sub_df.index.name = "canonical_ensembl_gene_id"

        # Save to temp TSV
        tmp_expr = Path(f"data/interim/direct_r_{arm_name}_harmonized_expr.tsv")
        sub_df.to_csv(tmp_expr, sep="\t")

        # Independent R script
        r_code = f"""
        suppressPackageStartupMessages(library(limma))
        expr_df <- read.delim('{tmp_expr.as_posix()}', row.names=1, sep='\\t', check.names=FALSE)
        expr_mat <- as.matrix(expr_df)
        storage.mode(expr_mat) <- 'numeric'

        samples <- colnames(expr_mat)
        conditions <- factor(c(rep('treatment', {len(trt_samples)}), rep('control', {len(ctrl_samples)})))
        design <- model.matrix(~ 0 + conditions)
        colnames(design) <- c('control', 'treatment')

        fit <- lmFit(expr_mat, design)
        contrast <- makeContrasts(treatment - control, levels=design)
        fit2 <- contrasts.fit(fit, contrast)
        fit_eb <- eBayes(fit2, trend=TRUE)
        top <- topTable(fit_eb, number=Inf, adjust.method='BH', sort.by='none')

        top$canonical_ensembl_gene_id <- rownames(top)
        write.csv(top, 'data/interim/direct_r_{arm_name}_limma_harmonized_results.csv', row.names=FALSE)
        """
        r_script_file = Path(f"data/interim/run_direct_{arm_name}.R")
        r_script_file.write_text(r_code, encoding="utf-8")

        # Run direct R
        cmd = [R_BINARY, str(r_script_file)]
        subprocess.run(cmd, capture_output=True, text=True, check=True)

        direct_r_df = pd.read_csv(f"data/interim/direct_r_{arm_name}_limma_harmonized_results.csv")
        merged = pd.merge(
            reg_df,
            direct_r_df,
            on="canonical_ensembl_gene_id",
            suffixes=("_pipe", "_r"),
        )
        print(f"   Evaluated shared canonical loci: {len(merged):,}")

        # Compute differences
        diff_lfc = (merged["log2_fold_change"] - merged["logFC"]).abs()
        diff_t = (merged["wald_statistic"] - merged["t"]).abs()
        diff_p = (merged["p_value"] - merged["P.Value"]).abs()
        diff_fdr = (merged["adjusted_p_value"] - merged["adj.P.Val"]).abs()

        results[arm_name] = {
            "loci_count": len(merged),
            "max_abs_diff_log2fc": float(diff_lfc.max()),
            "mean_abs_diff_log2fc": float(diff_lfc.mean()),
            "max_abs_diff_moderated_t": float(diff_t.max()),
            "mean_abs_diff_moderated_t": float(diff_t.mean()),
            "max_abs_diff_p_value": float(diff_p.max()),
            "mean_abs_diff_p_value": float(diff_p.mean()),
            "max_abs_diff_fdr": float(diff_fdr.max()),
            "mean_abs_diff_fdr": float(diff_fdr.mean()),
        }

        print(f"   Max abs diff log2FC:     {diff_lfc.max():.2e}")
        print(f"   Max abs diff moderated t:{diff_t.max():.2e}")
        print(f"   Max abs diff p-value:    {diff_p.max():.2e}")
        print(f"   Max abs diff FDR:        {diff_fdr.max():.2e}")

    out_json = Path("data/processed/limma_numerical_parity_audit.json")
    out_json.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("\nParity audit saved to", out_json)

if __name__ == "__main__":
    main()
