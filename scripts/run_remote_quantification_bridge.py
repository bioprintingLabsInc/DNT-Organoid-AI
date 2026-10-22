#!/usr/bin/env python3
"""Governed remote execution script for PRJNA802072 Quantification Bridge.

Executes only through the strandedness evidence gate:
1. Environment audit & hardware verification
2. Preflight test suite validation
3. Reference acquisition & cryptographic verification
4. Direct ENA acquisition of 12 approved libraries (24 FASTQ files)
5. FASTQ structural synchronization validation
6. Governed fastp treatment-blind preprocessing
7. STAR genome indexing
8. STAR alignment & coordinate-sorted BAM generation
9. SAMtools BAM indexing & alignment metrics collection
10. RSeQC infer_experiment.py strandedness diagnostic
11. Serialize authoritative prjna802072_strandedness_evidence.json
12. STOP and return AWAITING_OWNER_STRANDEDNESS_DECISION
"""

import argparse
import sys
from pathlib import Path

# Add repo root to python path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.quantification_bridge.pipeline import QuantificationBridgePipeline


def main():
    parser = argparse.ArgumentParser(
        description="Run governed PRJNA802072 quantification bridge through the strandedness evidence gate."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of CPU threads to allocate for fastp and STAR (default: 16)",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=REPO_ROOT,
        help="Path to repository root (default: auto-detected)",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Execute only environment audit and preflight test suite, then exit.",
    )

    args = parser.parse_args()

    print(f"=== PRJNA802072 QUANTIFICATION BRIDGE REMOTE RUNNER ===")
    print(f"Repo Root: {args.repo_root}")
    print(f"Allocated Threads: {args.threads}")

    pipeline = QuantificationBridgePipeline(repo_root=args.repo_root)

    if args.preflight_only:
        print("Executing preflight audit only...")
        env_metrics = pipeline.record_environment_metrics()
        print(f"Recorded Hardware: {env_metrics['cpu_cores']} cores | {env_metrics['physical_ram_gb']} GB RAM | {env_metrics['scratch_disk_free_gb']} GB free scratch")
        preflight_res = pipeline.run_preflight_checks()
        print(f"Preflight status: {preflight_res['status']}")
        sys.exit(0)

    # Execute through the strandedness gate
    state = pipeline.run_through_strandedness_gate(threads=args.threads)
    print(f"\nRemote execution package finished at governed state: {state}")
    sys.exit(0)


if __name__ == "__main__":
    main()
