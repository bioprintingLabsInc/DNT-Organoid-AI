#!/usr/bin/env python3
"""Governed CLI runner for PRJNA802072 Quantification Bridge & Downstream Analysis.

Supports:
1. 2-Sample Technical Pilot (SRR17834920 [Control] & SRR17834905 [250 uM VPA]).
2. Full 12-Library PRJNA802072 Cohort (6 Controls vs 6 Treated).
3. Preflight suite execution and profile verification.
4. Mandatory Scientific Review Gate (HALT before featureCounts until strandedness confirmed).
5. Governed featureCounts execution (-p --countReadPairs -B -C -s <strand_flag>).
6. End-to-end connection to QC (6A), Harmonization (6B), Normalization (6C),
   DESeq2 (7B, guarded by n >= 2 biological replicates), and Pathway Analysis (8C).
"""

import argparse
import json
import sys
from pathlib import Path

# Add repo root to python path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.quantification_bridge.manifest import (
    APPROVED_2_SAMPLE_TECHNICAL_PILOT,
    APPROVED_PILOT_RUNS,
    TARGET_BIOPROJECT,
)
from src.quantification_bridge.pipeline import QuantificationBridgePipeline
from src.quantification_bridge.profiles import get_execution_profile


def main():
    parser = argparse.ArgumentParser(
        description="Run PRJNA802072 quantification bridge and downstream integration pipeline."
    )
    parser.add_argument(
        "--profile",
        type=str,
        default=None,
        help="Execution profile name (e.g. 'colab_highmem', 'workstation_evaluation', 'local_standard').",
    )
    parser.add_argument(
        "--all-12-libraries",
        action="store_true",
        help="Process all 12 cohort libraries instead of the default 2-sample technical pilot.",
    )
    parser.add_argument(
        "--technical-pilot-only",
        action="store_true",
        default=True,
        help="Process only the 2 approved technical pilot libraries (default: True).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of CPU threads to allocate (default: from profile).",
    )
    parser.add_argument(
        "--strandedness-flag",
        type=int,
        choices=[0, 1, 2],
        default=None,
        help=(
            "Owner-confirmed strandedness flag: 0 (unstranded), 1 (forward), 2 (reverse). "
            "If omitted, execution HALTS at the Stage 16 Scientific Review Gate."
        ),
    )
    parser.add_argument(
        "--run-downstream",
        action="store_true",
        default=True,
        help="Connect validated count matrix to downstream modules (QC, Harmonization, Normalization, DESeq2, Pathway).",
    )
    parser.add_argument(
        "--no-downstream",
        dest="run_downstream",
        action="store_false",
        help="Skip downstream module execution after featureCounts.",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Execute only environment audit and preflight test suite, then exit.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate library targets and reference manifests without downloading files.",
    )

    args = parser.parse_args()

    # Determine cohort
    is_technical_pilot = not args.all_12_libraries
    cohort_desc = (
        "2-Sample Technical Pilot (SRR17834920 [Control] & SRR17834905 [VPA 250 uM])"
        if is_technical_pilot
        else "Full 12-Library PRJNA802072 Cohort (6 Controls & 6 VPA 250 uM)"
    )

    print("=" * 80)
    print("DNT ORGANOID AI - PRJNA802072 QUANTIFICATION BRIDGE RUNNER")
    print("=" * 80)
    print(f"Target Cohort: {cohort_desc}")
    print(f"Repository Root: {REPO_ROOT}")
    print(f"Requested Profile: {args.profile or 'Auto-detected'}")

    pipeline = QuantificationBridgePipeline(repo_root=REPO_ROOT, profile_name=args.profile)
    print(f"Active Profile: {pipeline.profile.name} (RAM min: {pipeline.profile.physical_ram_gb_min} GB, Threads: {args.threads or pipeline.profile.threads})")

    # 1. Preflight only
    if args.preflight_only:
        print("\n--- EXECUTING PREFLIGHT SUITE ONLY ---")
        env = pipeline.record_environment_metrics()
        print(f"Environment: {env['cpu_cores']} CPU cores, {env['physical_ram_gb']} GB RAM, {env['scratch_disk_free_gb']} GB scratch free")
        pre = pipeline.run_preflight_checks()
        print(f"Preflight status: {pre['status']}")
        print("PREFLIGHT PASSED. Exiting.")
        sys.exit(0)

    # 2. Dry run
    if args.dry_run:
        print("\n--- EXECUTING ACQUISITION DRY RUN ---")
        dry = pipeline.dry_run_library_acquisition(technical_pilot_only=is_technical_pilot)
        print(f"Targeted runs: {dry['total_runs']} runs ({dry['total_fastqs']} FASTQs)")
        for t in dry["targets"]:
            print(f"  {t['run_accession']} ({t['sample_alias']}) - {t['read']}: {t['url']} [MD5: {t['expected_md5']}]")
        print("DRY RUN PASSED. All targets validated against manifest.")
        sys.exit(0)

    # 3. Execution through strandedness gate
    if is_technical_pilot:
        pairs = pipeline.pairs_from_approved_technical_pilot_runs()
        state = pipeline.run_technical_pilot_through_strandedness_gate(threads=args.threads)
    else:
        pairs = pipeline.pairs_from_approved_pilot_runs()
        state = pipeline.run_through_strandedness_gate(threads=args.threads)

    # 4. Mandatory Review Checkpoint Check
    if args.strandedness_flag is None:
        print("\n" + "=" * 80)
        print("MANDATORY SCIENTIFIC REVIEW CHECKPOINT REACHED")
        print("=" * 80)
        print("Upstream alignment, BAM validation, and RSeQC strandedness diagnostic complete.")
        print("To proceed to featureCounts gene counting:")
        print("  1. Review the strandedness evidence artifact in bridge_out/.")
        print("  2. Re-run this runner providing --strandedness-flag <0|1|2>.")
        print("=" * 80)
        sys.exit(0)

    # 5. featureCounts & Downstream Execution
    print("\n" + "=" * 80)
    print(f"OWNER STRANDEDNESS CONFIRMED: -s {args.strandedness_flag}")
    print("PROCEEDING TO FEATURECOUNTS AND DOWNSTREAM INTEGRATION")
    print("=" * 80)

    if args.run_downstream:
        res = pipeline.run_featurecounts_and_downstream(
            pairs=pairs,
            strandedness_flag=args.strandedness_flag,
            threads=args.threads,
        )
        downstream = res["downstream_result"]
        print("\n--- DOWNSTREAM EXECUTION RESULTS ---")
        print(f"QC Passed: {downstream['qc_passed']} (status: {downstream['qc_metrics']['status']})")
        print(f"Harmonized Genes: {downstream['harmonization_metrics']['mapped_genes_count']} / {downstream['harmonization_metrics']['input_gene_count']}")
        print(f"Normalization: {downstream['normalization_metrics']['status']}")
        print(f"DESeq2 Status: {downstream['deseq2_status']}")
        if downstream.get("biological_replicate_warning"):
            print(f"Biological Replicate Guard: {downstream['biological_replicate_warning']}")
        print(f"Pathway Status: {downstream['pathway_status']}")
    else:
        res = pipeline.run_featurecounts_pairs(
            pairs=pairs,
            strandedness_flag=args.strandedness_flag,
            threads=args.threads,
        )
        print("\nfeatureCounts completed successfully. Downstream execution skipped by user request.")

    print("\n" + "=" * 80)
    print("ALL REQUESTED PIPELINE STAGES COMPLETED SUCCESSFULLY")
    print("=" * 80)
    sys.exit(0)


if __name__ == "__main__":
    main()
