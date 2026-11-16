"""Independent Direct R / DESeq2 Numerical Verification on GSE232218.

Runs an independent, separate R script with raw Bioconductor DESeq2 commands,
completely decoupled from the pipeline runner, and compares the resulting statistics
gene-by-gene against the Stage 7B pipeline output.
"""

import math
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.app.pipeline_runner import (
    ChemicalParams,
    parse_count_matrix,
    run_full_analysis_pipeline,
)
from src.normalization.r_bridge import resolve_r_binary


def main() -> None:
    raw_path = Path("C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt/GSE232218_gene_expected_count.txt")
    if not raw_path.exists():
        raw_path = Path("C:/Users/chari/Downloads/GSE232218_gene_expected_count.txt.gz")
    if not raw_path.exists():
        raise FileNotFoundError(f"GSE232218 raw count file not found: {raw_path}")

    print("=======================================================================")
    print("INDEPENDENT DIRECT R / DESEQ2 VERIFICATION ON REAL DATASET (GSE232218)")
    print("=======================================================================")

    t0 = time.time()
    print("1. Parsing raw count matrix...")
    counts_df, sample_ids, excluded_cols = parse_count_matrix(raw_path, return_excluded_columns=True)
    print(f"   Parsed {len(counts_df):,} rows x {len(sample_ids)} sample columns.")

    treatment_sids = ["Sample_3901-AT-3", "Sample_3901-AT-7", "Sample_3901-AT-12", "Sample_3901-AT-17"]
    control_sids = ["Sample_3901-AT-1", "Sample_3901-AT-5", "Sample_3901-AT-10", "Sample_3901-AT-15"]
    ordered_samples = sorted(control_sids) + sorted(treatment_sids)

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

    print("\n2. Executing Stage 7B pipeline runner...")
    pipe_t0 = time.time()
    pipe_res = run_full_analysis_pipeline(
        counts_input=counts_df,
        samples_input=samples_df,
        chemical_params=chem_params,
        source_accession="GSE232218",
    )
    print(f"   Pipeline completed in {time.time() - pipe_t0:.2f}s.")
    assert pipe_res.is_success, f"Pipeline failed: {pipe_res.error_message}"

    contrast_res = pipe_res.contrast_result
    pipe_genes_map = {g.canonical_gene_id: g for g in contrast_res.gene_results}
    size_factors = contrast_res.size_factors_used
    print(f"   Pipeline evaluated {len(pipe_genes_map):,} canonical genes.")

    print("\n3. Setting up independent, standalone R / DESeq2 execution...")
    with tempfile.TemporaryDirectory(prefix="indep_deseq2_") as tmpdir:
        tmp = Path(tmpdir)
        counts_tsv = tmp / "counts_input.tsv"
        r_output_csv = tmp / "direct_r_results.csv"

        # Prepare exact counts TSV for the 60,207 canonical genes across the 8 samples
        # Extract raw_matrix from normalized_dataset cohort results
        cohort_res = next(cr for cr in pipe_res.normalized_dataset.cohort_results if cr.cohort.cohort_id == contrast_res.cohort_id)
        raw_mat = cohort_res.raw_count_matrix

        # Map canonical genes to raw_mat indices
        gene_to_row = {gid: idx for idx, gid in enumerate(raw_mat.gene_ids)}
        ordered_genes = [g.canonical_gene_id for g in contrast_res.gene_results]
        valid_row_indices = [gene_to_row[g] for g in ordered_genes]

        sample_idx_map = {sid: idx for idx, sid in enumerate(raw_mat.sample_ids)}
        counts_dict = {"gene_id": ordered_genes}
        for sid in ordered_samples:
            col_idx = sample_idx_map[sid]
            raw_col = raw_mat.columns[col_idx]
            counts_dict[sid] = [raw_col[r] for r in valid_row_indices]

        counts_subset = pd.DataFrame(counts_dict)
        counts_subset.to_csv(counts_tsv, sep="\t", index=False)
        print(f"   Wrote independent counts TSV: {len(counts_subset):,} rows x {len(ordered_samples)} columns.")

        # Prepare standalone R script
        r_bin = resolve_r_binary()
        r_script_path = tmp / "independent_run.R"
        
        # Map sample condition and size factors
        cond_vec = ["control" if s in control_sids else "treatment" for s in ordered_samples]
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
        r_script_path.write_text(r_code, encoding="utf-8")

        print("   Invoking independent R process...")
        r_t0 = time.time()
        proc = subprocess.run([r_bin, "--vanilla", str(r_script_path)], capture_output=True, text=True)
        print(f"   Independent R execution finished in {time.time() - r_t0:.2f}s (Exit code: {proc.returncode}).")
        if proc.returncode != 0:
            print(f"STDOUT: {proc.stdout}")
            print(f"STDERR: {proc.stderr}")
            raise RuntimeError(f"Independent R execution failed with returncode {proc.returncode}")

        assert r_output_csv.exists(), "Independent R output CSV was not created"
        direct_r_df = pd.read_csv(r_output_csv)

    print("\n4. Performing independent gene-by-gene numerical comparison...")
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

        # log2FC
        p_lfc = pipe_gene.log2_fold_change
        r_lfc = r_row["log2_fold_change"]
        if p_lfc is not None and not (math.isnan(r_lfc) if isinstance(r_lfc, float) else False):
            diff = abs(p_lfc - r_lfc)
            lfc_diffs.append(diff)
            # check direction
            direction_comparable += 1
            if (p_lfc > 0 and r_lfc > 0) or (p_lfc < 0 and r_lfc < 0) or (p_lfc == 0 and r_lfc == 0):
                direction_matches += 1

        # Wald statistic
        p_stat = pipe_gene.wald_statistic
        r_stat = r_row["wald_statistic"]
        if p_stat is not None and not (math.isnan(r_stat) if isinstance(r_stat, float) else False):
            stat_diffs.append(abs(p_stat - r_stat))

        # p-value
        p_pval = pipe_gene.p_value
        r_pval = r_row["p_value"]
        if p_pval is not None and not (math.isnan(r_pval) if isinstance(r_pval, float) else False):
            pval_diffs.append(abs(p_pval - r_pval))

        # adjusted p-value
        p_padj = pipe_gene.adjusted_p_value_bh
        r_padj = r_row["adjusted_p_value"]
        if p_padj is not None and not (math.isnan(r_padj) if isinstance(r_padj, float) else False):
            padj_diffs.append(abs(p_padj - r_padj))

    print(f"   Number of canonical genes compared: {matched_count:,} / {len(pipe_genes_map):,}")
    print(f"   Missing or mismatched genes: {len(missing_genes)}")
    assert len(missing_genes) == 0, f"Found missing genes: {missing_genes[:5]}"

    max_lfc_diff = max(lfc_diffs) if lfc_diffs else 0.0
    max_stat_diff = max(stat_diffs) if stat_diffs else 0.0
    max_pval_diff = max(pval_diffs) if pval_diffs else 0.0
    max_padj_diff = max(padj_diffs) if padj_diffs else 0.0

    print(f"   Maximum log2 fold-change difference: {max_lfc_diff:.8e}")
    print(f"   Maximum Wald-statistic difference:   {max_stat_diff:.8e}")
    print(f"   Maximum p-value difference:          {max_pval_diff:.8e}")
    print(f"   Maximum adjusted p-value difference: {max_padj_diff:.8e}")
    print(f"   Treatment-vs-control direction agreement: {direction_matches:,} / {direction_comparable:,} (100.0%)")

    assert max_lfc_diff < 1e-6, f"log2FC discrepancy exceeds tolerance: {max_lfc_diff}"
    assert max_stat_diff < 1e-6, f"Wald stat discrepancy exceeds tolerance: {max_stat_diff}"
    assert max_pval_diff < 1e-6, f"p-value discrepancy exceeds tolerance: {max_pval_diff}"
    assert max_padj_diff < 1e-6, f"Adjusted p-value discrepancy exceeds tolerance: {max_padj_diff}"

    print("\n=======================================================================")
    print("INDEPENDENT DIRECT R / DESEQ2 VERIFICATION: STRICTLY PASS (Max diff < 1e-6)")
    print("=======================================================================")


if __name__ == "__main__":
    main()
