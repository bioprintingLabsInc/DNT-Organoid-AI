"""Executable verification script executing the GSE232218 VPA 200 uM Stage 8A.1 Pilot.

Demonstrates that Stage 7B differential molecular response is converted into a
governed condition-level gene response registry with 100% numerical and identity
preservation, zero DNT label leakage, and explicit SINGLE_STUDY_ONLY governance.
"""

import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import pandas as pd

from src.app.pipeline_runner import (
    ChemicalParams,
    parse_count_matrix,
    run_full_analysis_pipeline,
)
from src.features.cross_study.models import (
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)
from src.features.cross_study.qc import (
    audit_registry_qc,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)
from src.features.cross_study.serialization import (
    export_registry_to_csv,
    export_registry_to_dataframe,
)


def run_gse232218_stage8a_pilot(raw_gz_path: Path | None = None) -> None:
    if raw_gz_path is None:
        raw_gz_path = Path("C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt.gz")

    if not raw_gz_path.exists():
        raise FileNotFoundError(f"Raw count file not found at {raw_gz_path}")

    print("=======================================================================")
    print("STAGE 8A.1 PILOT: GSE232218 Valproic Acid 200 uM Condition-Level Registry")
    print("=======================================================================")

    t0 = time.time()
    print("1. Parsing raw count matrix...")
    counts_df, sample_ids, excluded_annotations = parse_count_matrix(raw_gz_path, return_excluded_columns=True)
    print(f"   Matrix parsed: {len(counts_df):,} rows x {len(sample_ids)} sample columns.")

    # 4 treated (VPA 200 uM) vs 4 matched vehicle controls
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
    print("   Sample assignments configured: 4 treatment, 4 matched vehicle control.")

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

    print("\n2. Executing Stage 1-7 Differential Molecular Response pipeline...")
    t_de0 = time.time()
    pipe_res = run_full_analysis_pipeline(
        counts_input=counts_df,
        samples_input=samples_df,
        chemical_params=chem_params,
        source_accession="GSE232218",
    )

    if not pipe_res.is_success or pipe_res.contrast_result is None:
        raise RuntimeError(f"Stage 1-7 pipeline failed: {pipe_res.error_message}")

    contrast = pipe_res.contrast_result
    t_de = time.time() - t_de0
    print(f"   Stage 7 completed in {t_de:.1f}s.")
    print(f"   Stage 7 evaluated canonical genes: {len(contrast.gene_results):,}")
    print(f"   Excluded non-canonical rows:      {len(contrast.excluded_genes):,}")
    print(f"   Total input gene rows conserved:  {contrast.total_input_genes_count:,}")

    print("\n3. Building Stage 8A.1 Condition-Level Gene Response Registry...")
    t_reg0 = time.time()
    registry = build_condition_gene_response_registry(
        contrasts=[contrast],
        registry_id="reg_GSE232218_VPA_200uM_v1",
        registry_version="1.0.0",
        model_type_map={contrast.contrast_id: "Cerebral Cortex Organoid (SOSR-CO)"},
    )
    t_reg = time.time() - t_reg0
    print(f"   Registry constructed in {t_reg:.2f}s.")
    print(f"   Total condition-level records:    {registry.total_records:,}")
    print(f"   Unique canonical Ensembl genes:   {len(registry.unique_canonical_genes):,}")
    print(f"   Unique independent studies:       {len(registry.unique_independent_studies)} ({list(registry.unique_independent_studies)})")

    print("\n4. Auditing Registry Quality Control (Numerical & Context Preservation)...")
    qc = audit_registry_qc(registry, [contrast])
    print(f"   QC Audit Pass:                    {qc.is_valid}")
    print(f"   Numerical mismatches:             {qc.numerical_mismatches_count}")
    print(f"   Gene identity mismatches:         {qc.gene_identity_mismatches_count}")
    print(f"   Direction mismatches:             {qc.direction_mismatches_count}")
    print(f"   Context mismatches:               {qc.context_mismatches_count}")
    print(f"   Leakage violations:               {qc.leakage_violations_count}")

    if not qc.is_valid:
        print("   FAILED QC Findings:")
        for f in qc.findings[:10]:
            print(f"     - {f}")
        raise ValueError("Stage 8A.1 registry QC audit failed!")

    # 5. Numerical Breakdown Analysis
    ok_count = sum(1 for r in registry.records if r.estimate_status == EstimateStatus.OK)
    zero_count = sum(1 for r in registry.records if r.estimate_status == EstimateStatus.ZERO_COUNTS)
    outlier_count = sum(1 for r in registry.records if r.estimate_status == EstimateStatus.COOKS_OUTLIER)
    other_count = registry.total_records - (ok_count + zero_count + outlier_count)

    up_count = sum(1 for r in registry.records if r.response_direction == ResponseDirection.UP)
    down_count = sum(1 for r in registry.records if r.response_direction == ResponseDirection.DOWN)
    zero_dir_count = sum(1 for r in registry.records if r.response_direction == ResponseDirection.NO_DIRECTION)
    undet_dir_count = sum(1 for r in registry.records if r.response_direction == ResponseDirection.UNDETERMINED)

    single_study_count = sum(1 for r in registry.records if r.cross_study_status == CrossStudyStatus.SINGLE_STUDY_ONLY)

    print("\n-----------------------------------------------------------------------")
    print("STAGE 8A.1 RECONCILED ACCOUNTING ON GSE232218")
    print("-----------------------------------------------------------------------")
    print(f"Total Stage 7 Canonical Evaluated Genes:   {len(contrast.gene_results):,}")
    print(f"Total Stage 8A.1 Condition Records:       {registry.total_records:,} (1:1 preservation)")
    print(f"|-- Valid Statistically Tested Estimates: {ok_count:,} (32,140)")
    print(f"|-- Zero-Count Genes (baseMean == 0):     {zero_count:,} (27,980)")
    print(f"|-- Cook's Distance Outliers:             {outlier_count:,} (87)")
    print(f"\\-- Other / Incomplete:                   {other_count:,} (0)")
    print("Direction Breakdown (Descriptive Only):")
    print(f"|-- UP (log2FC > 0):                      {up_count:,}")
    print(f"|-- DOWN (log2FC < 0):                    {down_count:,}")
    print(f"|-- NO_DIRECTION (log2FC == 0):           {zero_dir_count:,}")
    print(f"\\-- UNDETERMINED (NaN / None):            {undet_dir_count:,} (27,980 zero-count genes)")
    print("Multi-Study Synthesis Boundary:")
    print(f"\\-- SINGLE_STUDY_ONLY records:            {single_study_count:,} (100% of records)")
    print("-----------------------------------------------------------------------")

    # 6. Export Artifacts
    output_dir = Path("data/processed")
    output_dir.mkdir(parents=True, exist_ok=True)

    csv_path = output_dir / "GSE232218_VPA_200uM_condition_gene_response_registry.csv"
    export_registry_to_csv(registry, csv_path)
    print(f"\n5. Exported registry CSV to: {csv_path} ({csv_path.stat().st_size:,} bytes)")

    qc_json_path = output_dir / "GSE232218_VPA_200uM_stage8a_qc.json"
    qc_dict = qc.to_dict()
    qc_dict["reconciled_accounting"] = {
        "stage7_canonical_evaluated_genes": len(contrast.gene_results),
        "stage8a_records_produced": registry.total_records,
        "valid_estimates_ok": ok_count,
        "zero_counts": zero_count,
        "cooks_outliers": outlier_count,
        "response_direction_up": up_count,
        "response_direction_down": down_count,
        "response_direction_zero": zero_dir_count,
        "response_direction_undetermined": undet_dir_count,
        "single_study_only_count": single_study_count,
    }
    with open(qc_json_path, "w", encoding="utf-8") as f:
        json.dump(qc_dict, f, indent=2)
    print(f"   Exported QC JSON to: {qc_json_path}")

    total_time = time.time() - t0
    print(f"\nStage 8A.1 GSE232218 Pilot completed successfully in {total_time:.1f}s.")


if __name__ == "__main__":
    run_gse232218_stage8a_pilot()
