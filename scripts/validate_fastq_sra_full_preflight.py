#!/usr/bin/env python3
"""Read-only entry gate for PRJNA802072 12-library FASTQ/SRA validation.

This command never downloads FASTQs, runs alignment, changes an index, or
authorizes featureCounts. Its purpose is to verify the exact pilot cohort
against live ENA metadata before the staged Colab execution starts.
"""

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from scripts.run_colab_reference_indexing import verify_index_artifacts
from src.quantification_bridge.manifest import APPROVED_PILOT_RUNS
from src.quantification_bridge.pipeline import QuantificationBridgePipeline


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--min-free-gb", type=float, default=70.0,
                        help="Minimum available local scratch required to start a single-library pilot")
    parser.add_argument("--report-path", type=Path, default=None,
                        help="Optional local JSON report path (not a Google Drive backup)")
    args = parser.parse_args()

    root = args.repo_root.expanduser().resolve()
    index_dir = root / "references" / "ensembl_112" / "star_index"
    errors = []
    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "repo_root": str(root),
        "index_dir": str(index_dir),
        "scope": "PRJNA802072: 12 runs, 24 paired FASTQ files",
        "operation": "READ_ONLY_PREFLIGHT_NO_FASTQ_DOWNLOAD",
    }

    print("=" * 78)
    print("DNT ORGANOID AI — FULL FASTQ/SRA PRE-EXECUTION VALIDATION")
    print("Read-only: no FASTQs downloaded; no alignment or index build")
    print("=" * 78)

    if len(APPROVED_PILOT_RUNS) != 12:
        errors.append(f"Expected exactly 12 approved runs, got {len(APPROVED_PILOT_RUNS)}.")
    runs = [x.run_accession for x in APPROVED_PILOT_RUNS]
    if len(set(runs)) != 12:
        errors.append("Pilot run IDs contain duplicates.")
    control = [x for x in APPROVED_PILOT_RUNS if x.condition_group == "control" and x.vpa_dose_um == 0]
    treated = [x for x in APPROVED_PILOT_RUNS if x.condition_group == "treatment" and x.vpa_dose_um == 250]
    if len(control) != 6 or len(treated) != 6:
        errors.append(f"Expected 6 controls and 6 VPA 250 uM libraries, got {len(control)} and {len(treated)}.")
    if any(x.developmental_day != 28 for x in APPROVED_PILOT_RUNS):
        errors.append("One or more pilot libraries have an unexpected collection day.")
    report["manifest"] = {
        "run_count": len(runs),
        "unique_run_count": len(set(runs)),
        "control_count": len(control),
        "treated_count": len(treated),
        "run_ids": runs,
    }

    try:
        pipeline = QuantificationBridgePipeline(repo_root=root, profile_name="colab_highmem")
        targets = pipeline.dry_run_library_acquisition()
        report["download_dry_run"] = {
            "status": targets["status"],
            "run_count": targets["total_runs"],
            "fastq_file_count": targets["total_fastqs"],
        }
        if targets["total_runs"] != 12 or targets["total_fastqs"] != 24:
            errors.append("The dry-run did not produce exactly 12 runs and 24 FASTQs.")
    except Exception as exc:
        errors.append(f"Internal 12-library dry-run failed: {exc}")

    if index_dir.is_dir():
        try:
            inventory = verify_index_artifacts(index_dir)
            report["index"] = {
                "is_complete": inventory["is_complete"],
                "required_files_present": len(inventory["file_details"]),
                "total_index_gb": inventory["total_index_gb"],
                "genome_parameters": inventory["genome_parameters"],
            }
            if not inventory["is_complete"]:
                errors.append(f"STAR index incomplete. Missing: {inventory['missing_files']}")
            gp = inventory["genome_parameters"]
            if gp.get("genomeSAsparseD") != "1" or gp.get("genomeSAindexNbases") != "14":
                errors.append("STAR index is not the expected dense sparseD=1, nbases=14 build.")
        except Exception as exc:
            errors.append(f"STAR index inventory failed: {exc}")
    else:
        errors.append("STAR index not found locally. Restore the verified archive first; do NOT rebuild by default.")

    try:
        usage = shutil.disk_usage(root)
        free_gb = round(usage.free / 1e9, 2)
        report["scratch"] = {
            "free_gb_decimal": free_gb,
            "total_gb_decimal": round(usage.total / 1e9, 2),
            "min_free_gb_decimal": args.min_free_gb,
            "full_12_library_unbatched_estimate_gb": 270.0,
        }
        print(f"Local scratch free: {free_gb} GB")
        if free_gb < args.min_free_gb:
            errors.append(f"Insufficient scratch to begin: {free_gb} GB free; need >= {args.min_free_gb} GB.")
        if free_gb < 270:
            print("[WARNING] Insufficient local storage for unbatched 12-library execution.")
            print("          Process one library at a time and persist verified outputs to Drive.")
    except Exception as exc:
        errors.append(f"Disk-space inspection failed: {exc}")

    print("Re-querying ENA for the exact 12 run IDs, metadata and FASTQ MD5 hashes...")
    check_script = root / "scripts" / "reverify_pilot_manifest.py"
    ena_proc = subprocess.run([sys.executable, str(check_script)],
                              cwd=str(root), capture_output=True, text=True, check=False)
    print(ena_proc.stdout[-14000:])
    if ena_proc.returncode != 0:
        errors.append("Live ENA re-verification did not pass; inspect output and do not download.")
        if ena_proc.stderr:
            print(ena_proc.stderr[-4000:])
    report["live_ena"] = {"passed": ena_proc.returncode == 0, "return_code": ena_proc.returncode}

    report["errors"] = errors
    report["status"] = "READY_FOR_GATED_TECHNICAL_PILOT" if not errors else "BLOCKED_REQUIRES_REVIEW"
    if args.report_path:
        output = args.report_path.expanduser()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"Local report: {output}")

    print("=" * 78)
    print(f"STATUS: {report['status']}")
    if errors:
        for issue in errors:
            print("BLOCKER:", issue)
    else:
        print("12 approved runs and 24 FASTQ targets checked against live ENA.")
        print("Reference and scratch checks passed. No FASTQ work performed.")
        print("Next: process 2 pilot runs; inspect mapping and strandedness before scaling to 12.")
    print("=" * 78)
    return 0 if not errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
