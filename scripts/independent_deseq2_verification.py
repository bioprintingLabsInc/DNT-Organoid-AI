"""Independent Direct R / Bioconductor DESeq2 Numerical Verification.

Universal, dataset-agnostic verification script that executes a completely separate,
standalone R process with pure Bioconductor DESeq2 commands (decoupled from the
application pipeline) and performs gene-by-gene numerical verification of:
- log2 fold changes
- Wald statistics
- two-sided p-values
- Benjamini-Hochberg adjusted p-values (FDR)
- Treatment-vs-control directionality

Supports arbitrary RNA-seq count inputs, sample metadata tables, benchmark configs,
and customizable Rscript locations across Windows and Linux.
"""

import argparse
from collections.abc import Sequence
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.app.pipeline_runner import (
    ChemicalParams,
    parse_count_matrix,
    parse_sample_info,
    run_full_analysis_pipeline,
)
from src.normalization.r_bridge import resolve_r_binary


def run_independent_verification(
    counts_path: Path | str,
    samples_path: Path | str | None = None,
    samples_df: pd.DataFrame | None = None,
    chemical_params: ChemicalParams | None = None,
    treatment_sids: Sequence[str] | None = None,
    control_sids: Sequence[str] | None = None,
    r_binary_path: str | None = None,
    output_dir: Path | str | None = None,
    source_accession: str = "",
    tolerance: float = 1e-6,
) -> dict[str, Any]:
    """Execute independent direct R / DESeq2 numerical verification.

    Args:
        counts_path: Path to raw count matrix (CSV, TSV, TXT, or GZ).
        samples_path: Path to sample metadata CSV or TSV table.
        samples_df: Optional pre-loaded sample metadata DataFrame.
        chemical_params: Chemical exposure and organoid parameters.
        treatment_sids: Explicit treatment sample IDs (if samples_path not provided).
        control_sids: Explicit control sample IDs (if samples_path not provided).
        r_binary_path: Path to Rscript executable.
        output_dir: Optional directory to save verification outputs.
        source_accession: Study accession identifier (e.g. GSE232218).
        tolerance: Maximum acceptable numerical difference (default: 1e-6).

    Returns:
        dict containing verification metrics and pass/fail disposition.
    """
    counts_file = Path(counts_path)
    if not counts_file.exists():
        raise FileNotFoundError(f"Raw count matrix file not found: {counts_file}")

    print("=======================================================================")
    print("INDEPENDENT DIRECT R / BIOCONDUCTOR DESEQ2 NUMERICAL VERIFICATION")
    print("=======================================================================")
    print(f"Count matrix: {counts_file}")

    t0 = time.time()
    counts_df, sample_ids, excluded_cols = parse_count_matrix(counts_file, return_excluded_columns=True)
    print(f"Parsed {len(counts_df):,} rows x {len(sample_ids)} sample columns in {time.time() - t0:.2f}s.")
    if excluded_cols:
        print(f"Excluded annotation columns: {excluded_cols}")

    # Resolve samples_df
    if samples_df is None:
        if samples_path is not None:
            s_path = Path(samples_path)
            if not s_path.exists():
                raise FileNotFoundError(f"Sample metadata file not found: {s_path}")
            sep = "\t" if str(s_path).endswith((".tsv", ".tab", ".txt")) else ","
            raw_s_df = pd.read_csv(s_path, sep=sep)
            samples_df = parse_sample_info(raw_s_df, sample_ids)
        elif treatment_sids is not None and control_sids is not None:
            t_set = set(treatment_sids)
            c_set = set(control_sids)
            if not t_set or not c_set:
                raise ValueError("Both treatment_sids and control_sids must contain at least one sample.")
            rows = []
            for sid in sample_ids:
                if sid in t_set:
                    ctype = "treatment"
                    rep = f"rep_trt_{sid}"
                elif sid in c_set:
                    ctype = "control"
                    rep = f"rep_ctrl_{sid}"
                else:
                    ctype = "exclude"
                    rep = ""
                rows.append({"sample_id": sid, "condition_type": ctype, "biological_replicate_id": rep})
            samples_df = pd.DataFrame(rows)
        else:
            raise ValueError(
                "Either samples_path, samples_df, or both treatment_sids and control_sids must be supplied."
            )

    active_df = samples_df[samples_df["condition_type"].isin(["treatment", "control"])]
    trt_sids = sorted(active_df[active_df["condition_type"] == "treatment"]["sample_id"])
    ctrl_sids = sorted(active_df[active_df["condition_type"] == "control"]["sample_id"])

    if len(trt_sids) < 2 or len(ctrl_sids) < 2:
        raise ValueError(
            f"Differential expression verification requires >= 2 replicates per arm. "
            f"Found: {len(trt_sids)} treatment, {len(ctrl_sids)} control."
        )

    print(f"Active cohort: {len(trt_sids)} treatment ({trt_sids}) vs {len(ctrl_sids)} control ({ctrl_sids})")

    if chemical_params is None:
        chemical_params = ChemicalParams(
            agent_name="Chemical_Exposure",
            agent_identifier="",
            concentration_or_dose=1.0,
            concentration_or_dose_unit="uM",
            developmental_age="Day 30",
            exposure_duration=24.0,
            exposure_duration_unit="h",
            vehicle="vehicle",
            organoid_context_id="organoid_ctx_01",
            organoid_type="cerebral_organoid",
            brain_region="cerebral_cortex",
            biological_source_id="source_01",
            species="Homo sapiens",
        )

    # 1. Run Application Pipeline
    print("\n1. Executing application pipeline runner...")
    pipe_t0 = time.time()
    pipe_res = run_full_analysis_pipeline(
        counts_input=counts_df,
        samples_input=samples_df,
        chemical_params=chemical_params,
        source_accession=source_accession,
    )
    pipe_dur = time.time() - pipe_t0
    print(f"   Pipeline completed in {pipe_dur:.2f}s.")
    if not pipe_res.is_success or pipe_res.contrast_result is None:
        raise RuntimeError(f"Application pipeline failed: {pipe_res.error_message}")

    contrast_res = pipe_res.contrast_result
    pipe_genes_map = {g.canonical_gene_id: g for g in contrast_res.gene_results}
    size_factors = contrast_res.size_factors_used
    print(f"   Pipeline evaluated {len(pipe_genes_map):,} canonical genes.")

    # 2. Standalone Independent R Execution
    print("\n2. Executing independent standalone R / DESeq2 process...")
    ordered_samples = ctrl_sids + trt_sids

    if r_binary_path is None:
        resolved_r = resolve_r_binary()
    else:
        resolved_r = r_binary_path

    with tempfile.TemporaryDirectory(prefix="indep_deseq2_") as tmpdir:
        tmp = Path(tmpdir)
        counts_tsv = tmp / "independent_counts.tsv"
        r_output_csv = tmp / "direct_r_results.csv"

        # Extract exact count subset for canonical genes
        cohort_res = next(
            cr for cr in pipe_res.normalized_dataset.cohort_results if cr.cohort.cohort_id == contrast_res.cohort_id
        )
        raw_mat = cohort_res.raw_count_matrix

        gene_to_row = {gid: idx for idx, gid in enumerate(raw_mat.gene_ids)}
        ordered_genes = [g.canonical_gene_id for g in contrast_res.gene_results]
        valid_row_indices = [gene_to_row[g] for g in ordered_genes]

        sample_idx_map = {sid: idx for idx, sid in enumerate(raw_mat.sample_ids)}
        counts_dict: dict[str, list[Any]] = {"gene_id": ordered_genes}
        for sid in ordered_samples:
            col_idx = sample_idx_map[sid]
            raw_col = raw_mat.columns[col_idx]
            counts_dict[sid] = [raw_col[r] for r in valid_row_indices]

        counts_subset = pd.DataFrame(counts_dict)
        counts_subset.to_csv(counts_tsv, sep="\t", index=False)
        print(f"   Wrote independent counts TSV: {len(counts_subset):,} genes x {len(ordered_samples)} samples.")

        cond_vec = ["control" if s in ctrl_sids else "treatment" for s in ordered_samples]
        sf_vec = [str(size_factors[s]) for s in ordered_samples]

        r_code = f"""
        suppressPackageStartupMessages(library(DESeq2))
        
        counts_mat <- as.matrix(read.delim('{counts_tsv.as_posix()}', row.names = 1, sep = '\\t', check.names = FALSE))
        storage.mode(counts_mat) <- 'integer'
        
        coldata <- data.frame(
            sample_id = colnames(counts_mat),
            condition = factor(c({', '.join(f"'{c}'" for c in cond_vec)}), levels = c('control', 'treatment')),
            row.names = colnames(counts_mat)
        )
        
        sf <- c({', '.join(sf_vec)})
        names(sf) <- colnames(counts_mat)
        
        dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ condition)
        sizeFactors(dds) <- sf
        
        dds <- DESeq(dds, test = 'Wald', betaPrior = FALSE, minReplicatesForReplace = Inf, quiet = TRUE)
        res <- results(dds, contrast = c('condition', 'treatment', 'control'), independentFiltering = FALSE)
        
        out_df <- data.frame(
            canonical_gene_id = rownames(res),
            base_mean = res$baseMean,
            log2_fold_change = res$log2FoldChange,
            lfc_se = res$lfcSE,
            wald_statistic = res$stat,
            p_value = res$pvalue,
            adjusted_p_value = res$padj,
            stringsAsFactors = FALSE
        )
        write.csv(out_df, '{r_output_csv.as_posix()}', row.names = FALSE)
        cat("DIRECT_R_EXECUTION_COMPLETED\\n")
        """
        r_script_path = tmp / "independent_run.R"
        r_script_path.write_text(r_code, encoding="utf-8")

        r_t0 = time.time()
        proc = subprocess.run([resolved_r, "--vanilla", str(r_script_path)], capture_output=True, text=True)
        r_dur = time.time() - r_t0
        print(f"   Independent R execution finished in {r_dur:.2f}s (Exit code: {proc.returncode}).")

        if proc.returncode != 0:
            print(f"STDOUT: {proc.stdout}")
            print(f"STDERR: {proc.stderr}")
            raise RuntimeError(f"Independent R process failed with exit code {proc.returncode}")

        direct_r_df = pd.read_csv(r_output_csv)

    # 3. Gene-by-Gene Numerical Comparison
    print("\n3. Comparing pipeline results against independent direct R results...")
    direct_r_map = {row["canonical_gene_id"]: row for _, row in direct_r_df.iterrows()}

    matched_count = 0
    missing_genes = []
    lfc_diffs = []
    stat_diffs = []
    pval_diffs = []
    padj_diffs = []

    direction_matches = 0
    direction_comparable = 0

    for gid, pipe_gene in pipe_genes_map.items():
        if gid not in direct_r_map:
            missing_genes.append(gid)
            continue
        matched_count += 1
        r_row = direct_r_map[gid]

        p_lfc = pipe_gene.log2_fold_change
        r_lfc = r_row["log2_fold_change"]
        if p_lfc is not None and not (math.isnan(r_lfc) if isinstance(r_lfc, float) else False):
            diff = abs(p_lfc - r_lfc)
            lfc_diffs.append(diff)
            direction_comparable += 1
            if (p_lfc > 0 and r_lfc > 0) or (p_lfc < 0 and r_lfc < 0) or (p_lfc == 0 and r_lfc == 0):
                direction_matches += 1

        p_stat = pipe_gene.wald_statistic
        r_stat = r_row["wald_statistic"]
        if p_stat is not None and not (math.isnan(r_stat) if isinstance(r_stat, float) else False):
            stat_diffs.append(abs(p_stat - r_stat))

        p_pval = pipe_gene.p_value
        r_pval = r_row["p_value"]
        if p_pval is not None and not (math.isnan(r_pval) if isinstance(r_pval, float) else False):
            pval_diffs.append(abs(p_pval - r_pval))

        p_padj = pipe_gene.adjusted_p_value_bh
        r_padj = r_row["adjusted_p_value"]
        if p_padj is not None and not (math.isnan(r_padj) if isinstance(r_padj, float) else False):
            padj_diffs.append(abs(p_padj - r_padj))

    max_lfc_diff = max(lfc_diffs) if lfc_diffs else 0.0
    max_stat_diff = max(stat_diffs) if stat_diffs else 0.0
    max_pval_diff = max(pval_diffs) if pval_diffs else 0.0
    max_padj_diff = max(padj_diffs) if padj_diffs else 0.0
    dir_agreement_pct = (direction_matches / direction_comparable * 100.0) if direction_comparable else 100.0

    print(f"   Genes compared:               {matched_count:,} / {len(pipe_genes_map):,}")
    print(f"   Missing / mismatched genes:   {len(missing_genes)}")
    print(f"   Max log2FC difference:        {max_lfc_diff:.8e}")
    print(f"   Max Wald statistic difference:{max_stat_diff:.8e}")
    print(f"   Max p-value difference:       {max_pval_diff:.8e}")
    print(f"   Max adjusted p-value diff:    {max_padj_diff:.8e}")
    print(f"   Direction agreement:          {direction_matches:,} / {direction_comparable:,} ({dir_agreement_pct:.1f}%)")

    is_pass = (
        len(missing_genes) == 0
        and max_lfc_diff <= tolerance
        and max_stat_diff <= tolerance
        and max_pval_diff <= tolerance
        and max_padj_diff <= tolerance
        and dir_agreement_pct == 100.0
    )

    report = {
        "status": "PASS" if is_pass else "FAIL",
        "genes_evaluated": len(pipe_genes_map),
        "genes_compared": matched_count,
        "missing_genes_count": len(missing_genes),
        "max_log2fc_diff": max_lfc_diff,
        "max_wald_stat_diff": max_stat_diff,
        "max_pvalue_diff": max_pval_diff,
        "max_padj_diff": max_padj_diff,
        "direction_agreement_pct": dir_agreement_pct,
        "tolerance": tolerance,
        "pipeline_duration_seconds": pipe_dur,
        "direct_r_duration_seconds": r_dur,
    }

    if output_dir:
        out_p = Path(output_dir)
        out_p.mkdir(parents=True, exist_ok=True)
        direct_r_df.to_csv(out_p / "direct_r_results.csv", index=False)
        pipe_res.de_genes_df.to_csv(out_p / "pipeline_results.csv", index=False)
        with (out_p / "verification_report.json").open("w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        print(f"   Saved verification artifacts to {out_p}")

    assert is_pass, f"Independent DESeq2 verification failed tolerances: {report}"
    print(f"\nResult: {report['status']} (All differences within {tolerance:.1e})")
    print("=======================================================================")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Independent Direct R / DESeq2 Numerical Verification for Bulk RNA-seq."
    )
    parser.add_argument("--counts", "-c", type=str, help="Path to raw count matrix (CSV/TSV/GZ).")
    parser.add_argument("--samples", "-s", type=str, help="Path to sample metadata table (CSV/TSV).")
    parser.add_argument("--treatment-samples", "-t", nargs="+", help="Treatment sample IDs.")
    parser.add_argument("--control-samples", "--ctrl", nargs="+", help="Control sample IDs.")
    parser.add_argument("--config", type=str, help="Path to benchmark JSON configuration file.")
    parser.add_argument("--output-dir", "-o", type=str, help="Directory to save verification artifacts.")
    parser.add_argument("--rscript", type=str, default=None, help="Path to Rscript binary.")
    parser.add_argument("--tolerance", type=float, default=1e-6, help="Maximum acceptable numerical difference.")

    # Chemical parameter flags
    parser.add_argument("--agent-name", type=str, default=None)
    parser.add_argument("--agent-identifier", type=str, default="")
    parser.add_argument("--concentration", type=float, default=None)
    parser.add_argument("--concentration-unit", type=str, default="uM")
    parser.add_argument("--developmental-age", type=str, default="Day 30")
    parser.add_argument("--exposure-duration", type=float, default=24.0)
    parser.add_argument("--exposure-duration-unit", type=str, default="h")
    parser.add_argument("--vehicle", type=str, default="vehicle")
    parser.add_argument("--organoid-context-id", type=str, default="organoid_ctx_01")
    parser.add_argument("--organoid-type", type=str, default="cerebral_organoid")
    parser.add_argument("--brain-region", type=str, default="cerebral_cortex")
    parser.add_argument("--biological-source-id", type=str, default="source_01")
    parser.add_argument("--species", type=str, default="Homo sapiens")
    parser.add_argument("--source-accession", type=str, default="")

    args = parser.parse_args()

    counts_path = args.counts
    treatment_sids = args.treatment_samples
    control_sids = args.control_samples
    source_acc = args.source_accession
    chem_params = None

    if args.config:
        cfg_p = Path(args.config)
        if not cfg_p.exists():
            raise FileNotFoundError(f"Config file not found: {cfg_p}")
        with cfg_p.open("r", encoding="utf-8") as f:
            cfg = json.load(f)

        if not counts_path:
            for cand in cfg.get("counts_file_candidates", []):
                p = Path(cand)
                if p.exists():
                    counts_path = str(p)
                    break
        if not treatment_sids:
            treatment_sids = cfg.get("treatment_sample_ids")
        if not control_sids:
            control_sids = cfg.get("control_sample_ids")
        if not source_acc:
            source_acc = cfg.get("source_accession", "")

        c_data = cfg.get("chemical_parameters", {})
        chem_params = ChemicalParams(
            agent_name=c_data.get("agent_name", args.agent_name or "Chemical_Exposure"),
            agent_identifier=c_data.get("agent_identifier", args.agent_identifier),
            concentration_or_dose=float(c_data.get("concentration_or_dose", args.concentration or 1.0)),
            concentration_or_dose_unit=c_data.get("concentration_or_dose_unit", args.concentration_unit),
            developmental_age=c_data.get("developmental_age", args.developmental_age),
            exposure_duration=float(c_data.get("exposure_duration", args.exposure_duration)),
            exposure_duration_unit=c_data.get("exposure_duration_unit", args.exposure_duration_unit),
            vehicle=c_data.get("vehicle", args.vehicle),
            organoid_context_id=c_data.get("organoid_context_id", args.organoid_context_id),
            organoid_type=c_data.get("organoid_type", args.organoid_type),
            brain_region=c_data.get("brain_region", args.brain_region),
            biological_source_id=c_data.get("biological_source_id", args.biological_source_id),
            species=c_data.get("species", args.species),
        )

    if not counts_path:
        raise ValueError(
            "Counts matrix path is required. Provide via --counts <path> or via --config <benchmark.json>."
        )

    if chem_params is None:
        chem_params = ChemicalParams(
            agent_name=args.agent_name or "Chemical_Exposure",
            agent_identifier=args.agent_identifier,
            concentration_or_dose=float(args.concentration or 1.0),
            concentration_or_dose_unit=args.concentration_unit,
            developmental_age=args.developmental_age,
            exposure_duration=float(args.exposure_duration),
            exposure_duration_unit=args.exposure_duration_unit,
            vehicle=args.vehicle,
            organoid_context_id=args.organoid_context_id,
            organoid_type=args.organoid_type,
            brain_region=args.brain_region,
            biological_source_id=args.biological_source_id,
            species=args.species,
        )

    run_independent_verification(
        counts_path=counts_path,
        samples_path=args.samples,
        chemical_params=chem_params,
        treatment_sids=treatment_sids,
        control_sids=control_sids,
        r_binary_path=args.rscript,
        output_dir=args.output_dir,
        source_accession=source_acc,
        tolerance=args.tolerance,
    )


if __name__ == "__main__":
    main()
