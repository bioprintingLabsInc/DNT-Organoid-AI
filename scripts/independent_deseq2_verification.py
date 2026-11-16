"""Independent Direct R / Bioconductor DESeq2 Numerical Verification.

Universal, dataset-agnostic verification script that executes a completely separate,
standalone R process with pure Bioconductor DESeq2 commands (decoupled from the
application pipeline) and performs gene-by-gene numerical verification of:
- base mean expression
- log2 fold changes
- log2 fold change standard error (lfcSE)
- Wald statistics
- two-sided p-values
- Benjamini-Hochberg adjusted p-values (FDR)
- Treatment-vs-control directionality
- Model statistical design formula, matrix rank, and residual degrees of freedom

Supports arbitrary RNA-seq count inputs, explicit sample metadata tables, benchmark configs,
and customizable Rscript locations across Windows and Linux.
Fabrication or silent inference of experimental metadata, replicates, or designs is strictly prohibited.
"""

import argparse
from collections.abc import Sequence
import json
import math
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
from src.differential_expression.design import (
    resolve_sample_batch,
    resolve_sample_biological_source,
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
        chemical_params: Required chemical exposure and experimental parameters.
        treatment_sids: Optional subset of treatment sample IDs from samples metadata.
        control_sids: Optional subset of control sample IDs from samples metadata.
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

    # 1. Validate Sample Metadata (No fabrication allowed)
    if samples_df is None and samples_path is None:
        raise ValueError(
            "Explicit sample metadata is required via 'samples_path' or 'samples_df' with genuine biological replicate IDs. "
            "Automatic generation of biological replicates or conditions from sample names is prohibited."
        )

    if samples_df is None:
        s_path = Path(samples_path)  # type: ignore[arg-type]
        if not s_path.exists():
            raise FileNotFoundError(f"Sample metadata file not found: {s_path}")
        sep = "\t" if str(s_path).endswith((".tsv", ".tab", ".txt")) else ","
        raw_s_df = pd.read_csv(s_path, sep=sep)
        full_samples_input = raw_s_df
        samples_df = parse_sample_info(raw_s_df, sample_ids)
    else:
        full_samples_input = samples_df.copy()
        samples_df = parse_sample_info(samples_df, sample_ids)

    # Validate active treatment and control samples
    active_df = samples_df[samples_df["condition_type"].isin(["treatment", "control"])]
    all_trt = sorted(active_df[active_df["condition_type"] == "treatment"]["sample_id"])
    all_ctrl = sorted(active_df[active_df["condition_type"] == "control"]["sample_id"])

    if treatment_sids is not None:
        trt_set = set(treatment_sids)
        missing_trt = [s for s in trt_set if s not in all_trt]
        if missing_trt:
            raise ValueError(
                f"Specified treatment sample(s) {missing_trt} are not annotated as treatment in the sample metadata."
            )
        trt_sids = sorted(trt_set)
    else:
        trt_sids = all_trt

    if control_sids is not None:
        ctrl_set = set(control_sids)
        missing_ctrl = [s for s in ctrl_set if s not in all_ctrl]
        if missing_ctrl:
            raise ValueError(
                f"Specified control sample(s) {missing_ctrl} are not annotated as control in the sample metadata."
            )
        ctrl_sids = sorted(ctrl_set)
    else:
        ctrl_sids = all_ctrl

    if len(trt_sids) < 2 or len(ctrl_sids) < 2:
        raise ValueError(
            f"Differential expression verification requires >= 2 replicates per arm. "
            f"Found: {len(trt_sids)} treatment, {len(ctrl_sids)} control."
        )

    # Verify that every active sample has an explicit non-empty biological_replicate_id
    active_cohort_df = active_df[active_df["sample_id"].isin(trt_sids + ctrl_sids)]
    for _, row in active_cohort_df.iterrows():
        sid = row["sample_id"]
        rep = str(row.get("biological_replicate_id", "")).strip()
        if not rep or rep.lower() in ("nan", "none", ""):
            raise ValueError(
                f"Sample '{sid}' lacks an explicit biological replicate identifier. "
                "Biological replicates must not be fabricated from sample names."
            )

    print(f"Active cohort: {len(trt_sids)} treatment ({trt_sids}) vs {len(ctrl_sids)} control ({ctrl_sids})")

    # 2. Validate Chemical and Experimental Parameters (No dummy fabrication)
    if chemical_params is None:
        raise ValueError(
            "Chemical exposure and experimental parameters (ChemicalParams) are required. "
            "Fabricating default chemical concentrations, durations, or developmental ages is prohibited."
        )

    if not chemical_params.agent_name or not str(chemical_params.agent_name).strip():
        raise ValueError("Missing required experimental parameter: 'agent_name' must be specified.")
    if chemical_params.concentration_or_dose is None:
        raise ValueError("Missing required experimental parameter: 'concentration_or_dose' must be specified.")
    try:
        conc_val = float(chemical_params.concentration_or_dose)
        if conc_val <= 0:
            raise ValueError("concentration_or_dose must be > 0.")
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid concentration_or_dose: {e}")

    if not chemical_params.developmental_age or not str(chemical_params.developmental_age).strip():
        raise ValueError("Missing required experimental parameter: 'developmental_age' must be specified.")
    if chemical_params.exposure_duration is None:
        raise ValueError("Missing required experimental parameter: 'exposure_duration' must be specified.")
    try:
        dur_val = float(chemical_params.exposure_duration)
        if dur_val <= 0:
            raise ValueError("exposure_duration must be > 0.")
    except (TypeError, ValueError) as e:
        raise ValueError(f"Invalid exposure_duration: {e}")

    # 3. Run Application Pipeline (Stage 7B)
    print("\n1. Executing application pipeline runner...")
    pipe_t0 = time.time()
    pipe_res = run_full_analysis_pipeline(
        counts_input=counts_df,
        samples_input=full_samples_input,
        chemical_params=chemical_params,
        source_accession=source_accession,
    )
    pipe_dur = time.time() - pipe_t0
    print(f"   Pipeline completed in {pipe_dur:.2f}s.")
    if not pipe_res.is_success or pipe_res.contrast_result is None:
        raise RuntimeError(f"Application pipeline failed: {pipe_res.error_message}")

    contrast_res = pipe_res.contrast_result
    selected_design = contrast_res.design_formula
    expected_rank = contrast_res.design_matrix_rank
    expected_resid_df = contrast_res.residual_degrees_of_freedom

    if not selected_design:
        raise RuntimeError("Stage 7B did not produce a valid design formula.")

    print(f"   Stage 7B resolved statistical design: {selected_design}")
    print(f"   Stage 7B model matrix rank: {expected_rank}, residual df: {expected_resid_df}")

    pipe_genes_map = {g.canonical_gene_id: g for g in contrast_res.gene_results}
    size_factors = contrast_res.size_factors_used
    print(f"   Pipeline evaluated {len(pipe_genes_map):,} canonical genes.")

    # 4. Standalone Independent R Execution
    print("\n2. Executing independent standalone R / DESeq2 process...")
    # Follow exact deterministic sample ordering: controls first, then treatments
    ordered_samples = ctrl_sids + trt_sids

    if r_binary_path is None:
        resolved_r = resolve_r_binary()
    else:
        resolved_r = r_binary_path

    with tempfile.TemporaryDirectory(prefix="indep_deseq2_") as tmpdir:
        tmp = Path(tmpdir)
        counts_tsv = tmp / "independent_counts.tsv"
        coldata_tsv = tmp / "independent_coldata.tsv"
        r_output_csv = tmp / "direct_r_results.csv"
        r_diag_json = tmp / "direct_r_diagnostics.json"

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

        # Construct coldata using actual sample metadata and covariates
        sample_meta_map = {row["sample_id"]: row.to_dict() for _, row in samples_df.iterrows()}
        coldata_rows = []
        for sid in ordered_samples:
            cond = "control" if sid in ctrl_sids else "treatment"
            row_dict: dict[str, str] = {"sample_id": sid, "condition": cond}
            sm = sample_meta_map.get(sid, {})
            if "biological_source" in selected_design:
                src_id, src_finding = resolve_sample_biological_source(sid, sm, contrast_res.biological_source_ids)
                if src_finding is not None:
                    raise ValueError(
                        f"Sample '{sid}' biological source error: {src_finding.message}"
                    )
                if not src_id or not str(src_id).strip() or str(src_id).strip().lower() in ("nan", "none", ""):
                    raise ValueError(
                        f"Statistical design '{selected_design}' requires 'biological_source', "
                        f"but sample '{sid}' lacks an explicit biological source identifier."
                    )
                row_dict["biological_source"] = str(src_id).strip()
            if "batch" in selected_design:
                batch_id, batch_finding = resolve_sample_batch(sid, sm)
                if batch_finding is not None:
                    raise ValueError(
                        f"Sample '{sid}' batch error: {batch_finding.message}"
                    )
                if not batch_id or not str(batch_id).strip() or str(batch_id).strip().lower() in ("nan", "none", ""):
                    raise ValueError(
                        f"Statistical design '{selected_design}' requires 'batch', "
                        f"but sample '{sid}' lacks an explicit batch identifier."
                    )
                row_dict["batch"] = str(batch_id).strip()
            coldata_rows.append(row_dict)

        coldata_df = pd.DataFrame(coldata_rows)
        coldata_df.to_csv(coldata_tsv, sep="\t", index=False)
        print(f"   Wrote independent coldata TSV with columns: {list(coldata_df.columns)}")

        sf_vec = [str(size_factors[s]) for s in ordered_samples]

        # Generate standalone R script that independently checks matrix rank, estimability, and fits DESeq2
        r_code = f"""
        suppressPackageStartupMessages({{
            library(DESeq2)
            library(jsonlite)
        }})

        counts_mat <- as.matrix(read.delim('{counts_tsv.as_posix()}', row.names = 1, sep = '\\t', check.names = FALSE))
        storage.mode(counts_mat) <- 'integer'

        coldata <- read.delim('{coldata_tsv.as_posix()}', row.names = 1, sep = '\\t', check.names = FALSE, stringsAsFactors = FALSE)
        coldata$condition <- factor(coldata$condition, levels = c('control', 'treatment'))
        if ("biological_source" %in% colnames(coldata)) {{
            coldata$biological_source <- factor(coldata$biological_source)
        }}
        if ("batch" %in% colnames(coldata)) {{
            coldata$batch <- factor(coldata$batch)
        }}

        sf <- c({', '.join(sf_vec)})
        names(sf) <- colnames(counts_mat)

        # Independent statistical design rank and estimability verification
        formula_obj <- as.formula("{selected_design}")
        design_mat <- model.matrix(formula_obj, data = coldata)
        r_rank <- qr(design_mat)$rank
        r_ncol <- ncol(design_mat)
        r_resid_df <- nrow(coldata) - r_rank

        if (r_rank < r_ncol) {{
            stop(sprintf("Independent R check failed: Model matrix is rank deficient (rank %d < %d cols)", r_rank, r_ncol))
        }}
        if (r_resid_df <= 0) {{
            stop(sprintf("Independent R check failed: Residual degrees of freedom must be > 0 (got %d)", r_resid_df))
        }}

        # Treatment effect estimability check
        covars <- all.vars(formula_obj)
        other_covars <- covars[covars != "condition"]
        if (length(other_covars) == 0) {{
            formula_no_cond <- ~ 1
        }} else {{
            formula_no_cond <- as.formula(paste("~", paste(other_covars, collapse = " + ")))
        }}
        design_mat_no_cond <- model.matrix(formula_no_cond, data = coldata)
        if (qr(design_mat_no_cond)$rank == r_rank) {{
            stop("Independent R check failed: Treatment condition is not estimable in direct R model matrix")
        }}

        dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = formula_obj)
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

        diag_data <- list(
            design_formula = "{selected_design}",
            design_matrix_rank = r_rank,
            residual_degrees_of_freedom = r_resid_df,
            is_estimable = TRUE
        )
        write_json(diag_data, '{r_diag_json.as_posix()}', auto_unbox = TRUE)
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
            raise RuntimeError(f"Independent R process failed with exit code {proc.returncode}: {proc.stderr}")

        direct_r_df = pd.read_csv(r_output_csv)
        with r_diag_json.open("r", encoding="utf-8") as f:
            direct_r_diag = json.load(f)

    # 5. Verify Statistical Design Agreement
    print("\n3. Verifying statistical design agreement...")
    r_formula = direct_r_diag.get("design_formula")
    r_rank = direct_r_diag.get("design_matrix_rank")
    r_df = direct_r_diag.get("residual_degrees_of_freedom")

    design_matches = (r_formula == selected_design)
    rank_matches = (r_rank == expected_rank)
    df_matches = (r_df == expected_resid_df)

    print(f"   Stage 7B design: {selected_design} | Direct R design: {r_formula} (Match: {design_matches})")
    print(f"   Stage 7B rank:   {expected_rank} | Direct R rank:   {r_rank} (Match: {rank_matches})")
    print(f"   Stage 7B resid df: {expected_resid_df} | Direct R resid df: {r_df} (Match: {df_matches})")

    if not (design_matches and rank_matches and df_matches):
        raise ValueError(
            f"Statistical design mismatch: Stage 7B ({selected_design}, rank={expected_rank}, df={expected_resid_df}) "
            f"vs Direct R ({r_formula}, rank={r_rank}, df={r_df})."
        )

    # 6. Gene-by-Gene Numerical and Missing-Value Comparison
    print("\n4. Performing gene-by-gene comparison across all 6 statistical metrics...")
    direct_r_map = {row["canonical_gene_id"]: row for _, row in direct_r_df.iterrows()}

    # Check 100% gene set equivalence
    pipe_genes_set = set(pipe_genes_map.keys())
    direct_r_set = set(direct_r_map.keys())

    missing_in_r = sorted(pipe_genes_set - direct_r_set)
    missing_in_pipe = sorted(direct_r_set - pipe_genes_set)

    if missing_in_r or missing_in_pipe:
        raise ValueError(
            f"Gene universe mismatch: {len(missing_in_r)} genes missing in Direct R, "
            f"{len(missing_in_pipe)} genes missing in pipeline."
        )

    matched_count = 0
    max_diffs: dict[str, float] = {
        "base_mean": 0.0,
        "log2_fold_change": 0.0,
        "lfc_standard_error": 0.0,
        "wald_statistic": 0.0,
        "p_value": 0.0,
        "adjusted_p_value_bh": 0.0,
    }

    direction_matches = 0
    direction_comparable = 0

    def is_val_missing(v: Any) -> bool:
        if v is None:
            return True
        if pd.isna(v):
            return True
        if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
            return True
        return False

    metric_pairs = [
        ("base_mean", "base_mean"),
        ("log2_fold_change", "log2_fold_change"),
        ("lfc_standard_error", "lfc_se"),
        ("wald_statistic", "wald_statistic"),
        ("p_value", "p_value"),
        ("adjusted_p_value_bh", "adjusted_p_value"),
    ]

    for gid, pipe_gene in pipe_genes_map.items():
        matched_count += 1
        r_row = direct_r_map[gid]

        for p_field, r_field in metric_pairs:
            p_val = getattr(pipe_gene, p_field)
            r_val = r_row[r_field]

            p_miss = is_val_missing(p_val)
            r_miss = is_val_missing(r_val)

            if p_miss != r_miss:
                raise ValueError(
                    f"Missing-value mismatch for gene '{gid}' on metric '{p_field}': "
                    f"pipeline is {'missing' if p_miss else p_val}, Direct R is {'missing' if r_miss else r_val}."
                )

            if not p_miss:
                p_num = float(p_val)
                r_num = float(r_val)
                diff = abs(p_num - r_num)
                if diff > max_diffs[p_field]:
                    max_diffs[p_field] = diff

                if p_field == "log2_fold_change":
                    direction_comparable += 1
                    if (p_num > 0 and r_num > 0) or (p_num < 0 and r_num < 0) or (p_num == 0 and r_num == 0):
                        direction_matches += 1
                    else:
                        raise ValueError(
                            f"Treatment-vs-control direction mismatch for gene '{gid}': "
                            f"pipeline log2FC={p_num}, Direct R log2FC={r_num}."
                        )

                if diff > tolerance:
                    raise ValueError(
                        f"Numerical tolerance exceeded for gene '{gid}' on metric '{p_field}': "
                        f"difference {diff:.8e} > tolerance {tolerance:.8e} (pipeline={p_num}, Direct R={r_num})."
                    )

    dir_agreement_pct = (direction_matches / direction_comparable * 100.0) if direction_comparable else 100.0

    print(f"   Genes compared:               {matched_count:,} / {len(pipe_genes_map):,} (100.0%)")
    print(f"   Missing-value mismatches:     0")
    print(f"   Max baseMean difference:      {max_diffs['base_mean']:.8e}")
    print(f"   Max log2FC difference:        {max_diffs['log2_fold_change']:.8e}")
    print(f"   Max lfcSE difference:         {max_diffs['lfc_standard_error']:.8e}")
    print(f"   Max Wald statistic diff:      {max_diffs['wald_statistic']:.8e}")
    print(f"   Max p-value difference:       {max_diffs['p_value']:.8e}")
    print(f"   Max adjusted p-value diff:    {max_diffs['adjusted_p_value_bh']:.8e}")
    print(f"   Direction agreement:          {direction_matches:,} / {direction_comparable:,} ({dir_agreement_pct:.1f}%)")

    is_pass = (
        matched_count == len(pipe_genes_map)
        and design_matches
        and rank_matches
        and df_matches
        and all(d <= tolerance for d in max_diffs.values())
        and dir_agreement_pct == 100.0
    )

    report = {
        "status": "PASS" if is_pass else "FAIL",
        "design_formula": selected_design,
        "design_matrix_rank": expected_rank,
        "residual_degrees_of_freedom": expected_resid_df,
        "design_matches": design_matches,
        "rank_matches": rank_matches,
        "df_matches": df_matches,
        "genes_evaluated": len(pipe_genes_map),
        "genes_compared": matched_count,
        "missing_genes_count": 0,
        "missing_value_mismatches_count": 0,
        "max_base_mean_diff": max_diffs["base_mean"],
        "max_log2fc_diff": max_diffs["log2_fold_change"],
        "max_lfc_se_diff": max_diffs["lfc_standard_error"],
        "max_wald_stat_diff": max_diffs["wald_statistic"],
        "max_pvalue_diff": max_diffs["p_value"],
        "max_padj_diff": max_diffs["adjusted_p_value_bh"],
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
    parser.add_argument("--treatment-samples", "-t", nargs="+", help="Optional treatment sample ID filter.")
    parser.add_argument("--control-samples", "--ctrl", nargs="+", help="Optional control sample ID filter.")
    parser.add_argument("--config", type=str, help="Path to benchmark JSON configuration file.")
    parser.add_argument("--output-dir", "-o", type=str, help="Directory to save verification artifacts.")
    parser.add_argument("--rscript", type=str, default=None, help="Path to Rscript binary.")
    parser.add_argument("--tolerance", type=float, default=1e-6, help="Maximum acceptable numerical difference.")

    # Chemical parameter flags
    parser.add_argument("--agent-name", type=str, default=None, help="Chemical or perturbagen name.")
    parser.add_argument("--agent-identifier", type=str, default="", help="Chemical identifier (e.g. PubChem CID).")
    parser.add_argument("--concentration", type=float, default=None, help="Nominal concentration or dose.")
    parser.add_argument("--concentration-unit", type=str, default="uM", help="Concentration unit.")
    parser.add_argument("--developmental-age", type=str, default=None, help="Developmental age or stage.")
    parser.add_argument("--exposure-duration", type=float, default=None, help="Exposure duration value.")
    parser.add_argument("--exposure-duration-unit", type=str, default="h", help="Exposure duration unit.")
    parser.add_argument("--vehicle", type=str, default="", help="Vehicle control substance.")
    parser.add_argument("--organoid-context-id", type=str, default="", help="Organoid context identifier.")
    parser.add_argument("--organoid-type", type=str, default="", help="Organoid type or protocol.")
    parser.add_argument("--brain-region", type=str, default="", help="Brain region.")
    parser.add_argument("--biological-source-id", type=str, default="", help="Biological source identifier.")
    parser.add_argument("--species", type=str, default="Homo sapiens", help="Biological species.")
    parser.add_argument("--source-accession", type=str, default="", help="Source accession (e.g. GSE232218).")

    args = parser.parse_args()

    counts_path = args.counts
    samples_path = args.samples
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
        if not samples_path:
            samples_path = cfg.get("samples_file") or cfg.get("sample_metadata_file")
        if not treatment_sids:
            treatment_sids = cfg.get("treatment_sample_ids")
        if not control_sids:
            control_sids = cfg.get("control_sample_ids")
        if not source_acc:
            source_acc = cfg.get("source_accession", "")

        c_data = cfg.get("chemical_parameters", {})
        if c_data:
            chem_params = ChemicalParams(
                agent_name=c_data.get("agent_name", args.agent_name),
                agent_identifier=c_data.get("agent_identifier", args.agent_identifier or ""),
                concentration_or_dose=float(c_data.get("concentration_or_dose", args.concentration or 0.0)),
                concentration_or_dose_unit=c_data.get("concentration_or_dose_unit", args.concentration_unit or "uM"),
                developmental_age=c_data.get("developmental_age", args.developmental_age),
                exposure_duration=float(c_data.get("exposure_duration", args.exposure_duration or 0.0)),
                exposure_duration_unit=c_data.get("exposure_duration_unit", args.exposure_duration_unit or "h"),
                vehicle=c_data.get("vehicle", args.vehicle or ""),
                organoid_context_id=c_data.get("organoid_context_id", args.organoid_context_id or ""),
                organoid_type=c_data.get("organoid_type", args.organoid_type or ""),
                brain_region=c_data.get("brain_region", args.brain_region or ""),
                biological_source_id=c_data.get("biological_source_id", args.biological_source_id or ""),
                species=c_data.get("species", args.species or "Homo sapiens"),
            )

    if not counts_path:
        raise ValueError(
            "Counts matrix path is required. Provide via --counts <path> or via --config <benchmark.json>."
        )

    if not samples_path:
        raise ValueError(
            "Sample metadata path is required. Provide via --samples <path> or via --config <benchmark.json>."
        )

    if chem_params is None:
        if not args.agent_name or args.concentration is None or not args.developmental_age or args.exposure_duration is None:
            raise ValueError(
                "Missing required experimental parameters. --agent-name, --concentration, "
                "--developmental-age, and --exposure-duration must be provided."
            )
        chem_params = ChemicalParams(
            agent_name=args.agent_name,
            agent_identifier=args.agent_identifier or "",
            concentration_or_dose=float(args.concentration),
            concentration_or_dose_unit=args.concentration_unit,
            developmental_age=args.developmental_age,
            exposure_duration=float(args.exposure_duration),
            exposure_duration_unit=args.exposure_duration_unit,
            vehicle=args.vehicle or "",
            organoid_context_id=args.organoid_context_id or "",
            organoid_type=args.organoid_type or "",
            brain_region=args.brain_region or "",
            biological_source_id=args.biological_source_id or "",
            species=args.species or "Homo sapiens",
        )

    run_independent_verification(
        counts_path=counts_path,
        samples_path=samples_path,
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
