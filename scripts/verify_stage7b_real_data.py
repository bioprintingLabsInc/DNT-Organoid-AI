"""Real-data validation script executing Stage 7B Differential Expression on GSE232218."""

import math
from pathlib import Path
import sys
import time

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.app.pipeline_runner import (
    ChemicalParams,
    parse_count_matrix,
    run_full_analysis_pipeline,
)


def main() -> None:
    # 1. Locate count matrix
    candidates = [
        Path("C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt/GSE232218_gene_expected_count.txt"),
        Path("C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt.gz"),
        Path("C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt"),
    ]
    raw_path = None
    for c in candidates:
        if c.exists():
            raw_path = c
            break

    if raw_path is None:
        raise FileNotFoundError(f"GSE232218 raw count file not found among candidates: {candidates}")

    print("=======================================================================")
    print("STAGE 7B REAL-DATA VALIDATION: GSE232218 VPA 200 µM in Human Organoids")
    print("=======================================================================")
    print(f"Loading raw count matrix from: {raw_path}")

    t0 = time.time()
    counts_df, sample_ids, excluded_annotations = parse_count_matrix(raw_path, return_excluded_columns=True)
    print(f"Parsed {len(counts_df):,} genes x {len(sample_ids)} samples in {time.time() - t0:.2f}s")
    print(f"Excluded annotations: {excluded_annotations}")

    # Set up sample assignments: 4 VPA 200 uM treatment vs 4 Vehicle control
    treatment_sids = {"Sample_3901-AT-3", "Sample_3901-AT-7", "Sample_3901-AT-12", "Sample_3901-AT-17"}
    control_sids = {"Sample_3901-AT-1", "Sample_3901-AT-5", "Sample_3901-AT-10", "Sample_3901-AT-15"}

    samples_rows = []
    for sid in sample_ids:
        if sid in treatment_sids:
            cond = "treatment"
            rep = sid.replace("Sample_3901-AT-", "rep_trt_")
        elif sid in control_sids:
            cond = "control"
            rep = sid.replace("Sample_3901-AT-", "rep_ctrl_")
        else:
            cond = "exclude"
            rep = ""
        samples_rows.append({"sample_id": sid, "condition_type": cond, "biological_replicate_id": rep})

    samples_df = pd.DataFrame(samples_rows)
    trt_count = sum(1 for r in samples_rows if r["condition_type"] == "treatment")
    ctrl_count = sum(1 for r in samples_rows if r["condition_type"] == "control")
    print(f"Configured cohort: {trt_count} treatment vs {ctrl_count} control samples.")

    chem_params = ChemicalParams(
        agent_name="Valproic Acid",
        agent_identifier="CID:3121",
        concentration_or_dose=200.0,
        concentration_or_dose_unit="uM",
        developmental_age="Day 6",
        exposure_duration=48.0,
        exposure_duration_unit="h",
        vehicle="water",
        organoid_context_id="organoid_ctx_01",
        organoid_type="Cerebral Cortex Organoid (SOSR-CO)",
        brain_region="cerebral_cortex",
        biological_source_id="WA09",
        species="Homo sapiens",
    )

    print("\nExecuting full Stage 1-7 pipeline (QC -> Harmonization -> Normalization -> Contrast -> DESeq2)...")
    t1 = time.time()
    pipe_res = run_full_analysis_pipeline(
        counts_input=counts_df,
        samples_input=samples_df,
        chemical_params=chem_params,
        source_accession="GSE232218",
    )
    duration = time.time() - t1
    print(f"Pipeline executed in {duration:.2f}s")

    assert pipe_res.is_success, f"Pipeline failed: {pipe_res.error_message}"
    contrast = pipe_res.contrast_result
    assert contrast is not None, "Contrast result is None"
    assert contrast.is_inferentially_eligible, "Contrast is not inferentially eligible"
    assert contrast.is_gene_universe_conserved, "Gene universe conservation failed"

    print("\nStage 7B Results:")
    print(f"  Contrast ID: {contrast.contrast_id}")
    print(f"  Design Formula: {contrast.design_formula}")
    print(f"  Total Input Rows: {contrast.total_input_genes_count:,}")
    print(f"  Evaluated Canonical Genes: {len(contrast.gene_results):,}")
    print(f"  Excluded Input Genes: {len(contrast.excluded_genes):,}")
    print(f"  Gene Universe Conserved: {contrast.is_gene_universe_conserved}")

    sig_up = sum(1 for g in contrast.gene_results if g.adjusted_p_value_bh and g.adjusted_p_value_bh < 0.05 and g.log2_fold_change and g.log2_fold_change > 0)
    sig_down = sum(1 for g in contrast.gene_results if g.adjusted_p_value_bh and g.adjusted_p_value_bh < 0.05 and g.log2_fold_change and g.log2_fold_change < 0)
    print(f"  Significant Upregulated (FDR < 0.05): {sig_up:,}")
    print(f"  Significant Downregulated (FDR < 0.05): {sig_down:,}")

    # Compare with verified registry
    baseline_path = Path("data/processed/GSE232218_VPA_200uM_condition_gene_response_registry.csv")
    if baseline_path.exists():
        print(f"\nComparing against verified baseline: {baseline_path}...")
        base_df = pd.read_csv(baseline_path)
        base_map = {row["canonical_ensembl_gene_id"]: row for _, row in base_df.iterrows()}

        matched = 0
        lfc_diffs = []
        stat_diffs = []
        for g in contrast.gene_results:
            if g.canonical_gene_id in base_map:
                matched += 1
                b_row = base_map[g.canonical_gene_id]
                b_lfc = b_row["log2_fold_change"]
                b_stat = b_row["wald_statistic"]
                if g.log2_fold_change is not None and not math.isnan(b_lfc):
                    lfc_diffs.append(abs(g.log2_fold_change - b_lfc))
                if g.wald_statistic is not None and not math.isnan(b_stat):
                    stat_diffs.append(abs(g.wald_statistic - b_stat))

        print(f"  Matched genes: {matched:,} / {len(contrast.gene_results):,}")
        max_lfc_diff = max(lfc_diffs) if lfc_diffs else 0.0
        max_stat_diff = max(stat_diffs) if stat_diffs else 0.0
        mean_lfc_diff = sum(lfc_diffs) / len(lfc_diffs) if lfc_diffs else 0.0
        mean_stat_diff = sum(stat_diffs) / len(stat_diffs) if stat_diffs else 0.0
        print(f"  Max log2FC difference: {max_lfc_diff:.6e} (Mean: {mean_lfc_diff:.6e})")
        print(f"  Max Wald stat difference: {max_stat_diff:.6e} (Mean: {mean_stat_diff:.6e})")

        assert max_lfc_diff < 1e-4, f"log2FC discrepancy exceeded tolerance: {max_lfc_diff}"
        assert max_stat_diff < 1e-4, f"Wald stat discrepancy exceeded tolerance: {max_stat_diff}"
        print("  -> NUMERICAL REPRODUCIBILITY VERIFIED (Max diff < 1e-4 across all genes)")

    print("\n=======================================================================")
    print("STAGE 7B REAL-DATA VALIDATION: ALL CHECKS PASSED")
    print("=======================================================================")


if __name__ == "__main__":
    main()
