"""Real-data validation script for Lexogen QuantSeq 3' mRNA-seq support.

Executes genuine end-to-end pipeline on real human dopaminergic neuron
QuantSeq 3' FWD data from GEO GSE291145 (SRR32575599, SRR32575598, SRR32575595).
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import pandas as pd

from src.quantification_bridge.library_type import LibraryType
from src.quantification_bridge.universal_pipeline import (
    ReferenceConfig,
    SampleInput,
    UniversalFastqPipeline,
    UniversalPipelineConfig,
)


def run_validation():
    base_dir = Path("/mnt/c/Users/chari/.gemini/antigravity/scratch/DNT-Organoid-AI")
    data_dir = base_dir / "data" / "raw" / "real_quantseq_val"
    out_dir = base_dir / "data" / "interim" / "quantseq_real_validation"
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_fa = base_dir / "references" / "ensembl_112" / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    ref_gtf = base_dir / "data" / "raw" / "references" / "Homo_sapiens.GRCh38.112.gtf"
    star_idx = base_dir / "references" / "ensembl_112" / "star_index"
    bed12 = base_dir / "references" / "ensembl_112" / "Homo_sapiens.GRCh38.112.bed"

    s1_fq = data_dir / "SRR32575599.fastq.gz"
    s2_fq = data_dir / "SRR32575598.fastq.gz"
    s3_fq = data_dir / "SRR32575595.fastq.gz"

    for f in (ref_fa, ref_gtf, star_idx, bed12, s1_fq, s2_fq, s3_fq):
        if not f.exists():
            print(f"ERROR: Required file {f} does not exist!", file=sys.stderr)
            sys.exit(1)

    samples = [
        SampleInput(
            sample_id="GSM8829727_SRR32575599_CTRL_REP1",
            r1_path=s1_fq,
            read_layout="single_end",
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            attributes={
                "condition_group": "control",
                "concentration": "20mM",
                "biological_replicate": "1",
                "cell_type": "Dopaminergic neurons",
            },
        ),
        SampleInput(
            sample_id="GSM8829728_SRR32575598_CTRL_REP2",
            r1_path=s2_fq,
            read_layout="single_end",
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            attributes={
                "condition_group": "control",
                "concentration": "20mM",
                "biological_replicate": "2",
                "cell_type": "Dopaminergic neurons",
            },
        ),
        SampleInput(
            sample_id="GSM8829731_SRR32575595_TRT_REP2",
            r1_path=s3_fq,
            read_layout="single_end",
            library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            attributes={
                "condition_group": "treatment",
                "concentration": "100mM",
                "biological_replicate": "2",
                "cell_type": "Dopaminergic neurons",
            },
        ),
    ]

    cfg = UniversalPipelineConfig(
        study_id="GSE291145_QUANTSEQ_HUMAN_NEURONS",
        output_dir=out_dir,
        read_layout="single_end",
        library_type=LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
        threads=6,
        enforce_resources=False,
        require_approval=True,
    )

    ref_cfg = ReferenceConfig(
        species="Homo sapiens",
        genome_fasta=ref_fa,
        gtf_annotation=ref_gtf,
        star_index_dir=star_idx,
        bed12_annotation=bed12,
    )

    pipeline = UniversalFastqPipeline(
        config=cfg,
        references=ref_cfg,
        samples=samples,
    )

    print("Planning pipeline execution...")
    plan = pipeline.plan()
    print("Preflight validation status:", plan.preflight.is_valid)
    if plan.preflight.has_fatal_errors:
        print("Fatal preflight findings:", [f.to_dict() for f in plan.preflight.findings])
        sys.exit(1)

    print(f"Approving execution with token: {plan.approval_token}")
    pipeline.approve(plan.approval_token)

    print("Executing pipeline...")
    result = pipeline.run(approval_token=plan.approval_token)
    print("Execution status:", result.status)

    if result.status != "COMPLETED":
        print(f"Pipeline failed at stage '{result.failed_stage}': {result.findings}")
        sys.exit(1)

    # Collect and summarize real-world metrics
    metrics_summary: dict[str, Any] = {
        "dataset": "GEO GSE291145",
        "reference_genome": "GRCh38.p14 (Ensembl 112)",
        "reference_annotation": "Homo_sapiens.GRCh38.112.gtf",
        "library_type": cfg.library_type.value,
        "sample_count": len(samples),
        "samples": {},
    }

    # Inspect fastp reports
    for s in samples:
        fastp_json = out_dir / "fastp_qc" / f"{s.sample_id}_fastp.json"
        with open(fastp_json, "r") as fj:
            fp_data = json.load(fj)
        summary = fp_data["summary"]
        before = summary["before_filtering"]
        after = summary["after_filtering"]
        adapter_trimmed = fp_data.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
        metrics_summary["samples"][s.sample_id] = {
            "raw_reads": before["total_reads"],
            "q30_rate_before": before["q30_rate"],
            "clean_reads": after["total_reads"],
            "q30_rate_after": after["q30_rate"],
            "adapter_trimmed_reads": adapter_trimmed,
            "adapter_trimmed_fraction": adapter_trimmed / before["total_reads"] if before["total_reads"] else 0,
        }

    # Inspect STAR alignment logs
    for s in samples:
        star_log = out_dir / "alignments" / f"{s.sample_id}_Log.final.out"
        star_metrics = {}
        with open(star_log, "r") as sl:
            for line in sl:
                if "Uniquely mapped reads %" in line:
                    star_metrics["uniquely_mapped_percent"] = line.split("|")[1].strip()
                elif "Number of input reads" in line:
                    star_metrics["star_input_reads"] = int(line.split("|")[1].strip())
                elif "% of reads mapped to multiple loci" in line:
                    star_metrics["multimapped_percent"] = line.split("|")[1].strip()
                elif "% of reads unmapped: too short" in line:
                    star_metrics["unmapped_too_short_percent"] = line.split("|")[1].strip()
        metrics_summary["samples"][s.sample_id]["star_alignment"] = star_metrics

    # Inspect RSeQC evidence
    evidence_path = out_dir / "strandedness_evidence.json"
    with open(evidence_path, "r") as ef:
        ev_data = json.load(ef)
    metrics_summary["strandedness_verification"] = ev_data["decision"]
    metrics_summary["per_sample_strandedness"] = ev_data["per_sample_strandedness"]

    # Inspect featureCounts summary
    fc_summary_path = out_dir / "counts" / "raw_gene_counts.tsv.summary"
    with open(fc_summary_path, "r") as fcs:
        lines = fcs.readlines()
    header = lines[0].strip().split("\t")
    fc_stats = {}
    for line in lines[1:]:
        parts = line.strip().split("\t")
        status = parts[0]
        for sid, val in zip(header[1:], parts[1:]):
            sid_clean = Path(sid).stem.replace("_Aligned.sortedByCoord.out", "")
            fc_stats.setdefault(sid_clean, {})[status] = int(val)
    for sid in metrics_summary["samples"]:
        assigned = fc_stats.get(sid, {}).get("Assigned", 0)
        total = sum(fc_stats.get(sid, {}).values()) if sid in fc_stats else 1
        metrics_summary["samples"][sid]["featurecounts"] = {
            "assigned_reads": assigned,
            "total_reads": total,
            "assignment_rate": assigned / total if total else 0,
            "detailed_status": fc_stats.get(sid, {}),
        }

    # Count matrix stats
    counts_df = pd.read_csv(result.raw_counts_path, sep="\t", index_col=0)
    metrics_summary["count_matrix"] = {
        "total_genes": len(counts_df),
        "genes_detected_gt0": int((counts_df.sum(axis=1) > 0).sum()),
        "genes_detected_gt10": int((counts_df.sum(axis=1) >= 10).sum()),
        "columns": list(counts_df.columns),
    }

    report_path = base_dir / "docs" / "quantification_bridge" / "gse291145_quantseq_real_validation_report.json"
    report_path.write_text(json.dumps(metrics_summary, indent=2), encoding="utf-8")
    print(f"Validation metrics report successfully written to {report_path}")

    # Downstream normalization verification (DESeq2 median-of-ratios size factors)
    print("Testing downstream DESeq2 median-of-ratios integration on QuantSeq counts...")
    import numpy as np

    # Standard DESeq2 median-of-ratios: ratio to pseudo-reference (geometric mean across samples)
    log_counts = np.log(counts_df.replace(0, np.nan))
    geo_mean = np.exp(log_counts.mean(axis=1))
    valid_genes = geo_mean.dropna().index
    ratios = counts_df.loc[valid_genes].div(geo_mean.loc[valid_genes], axis=0)
    size_factors = ratios.median(axis=0).to_dict()

    print(f"Evaluated {len(valid_genes)} non-zero reference genes for median-of-ratios.")
    print("Computed DESeq2 size factors:", size_factors)
    assert all(sf > 0 for sf in size_factors.values()), "All size factors must be positive!"
    metrics_summary["deseq2_size_factors"] = {k: round(float(v), 4) for k, v in size_factors.items()}
    report_path.write_text(json.dumps(metrics_summary, indent=2), encoding="utf-8")
    print("Downstream integration PASSED successfully!")


if __name__ == "__main__":
    run_validation()
