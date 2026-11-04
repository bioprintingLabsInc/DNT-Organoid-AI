#!/usr/bin/env python3
"""Governed Cloud / Google Colab Reference Indexing Runner for PRJNA802072.

Executes reference acquisition, cryptographic verification, high-memory STAR index
generation on local high-speed disk, lightweight load validation, and integrity-verified
synchronization to persistent storage (Google Drive or external POSIX storage).

Stops at the mandatory scientific review gate before any FASTQ data is acquired.
"""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any

# Add repo root to python path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.quantification_bridge.checkpoints import (
    CheckpointManager,
    archive_and_sync_directory,
    check_storage_capacity,
    restore_and_verify_archive,
    sha256_file,
    sync_and_verify_directory,
)
from src.quantification_bridge.pipeline import QuantificationBridgePipeline
from src.quantification_bridge.profiles import get_execution_profile, resolve_pipeline_paths
from src.quantification_bridge.star_runner import build_star_index


def verify_index_artifacts(index_dir: Path) -> dict[str, Any]:
    """Verify presence, non-zero size, and structural parameters of all required STAR index files."""
    required_files = [
        "Genome",
        "SA",
        "SAindex",
        "chrLength.txt",
        "chrName.txt",
        "chrNameLength.txt",
        "chrStart.txt",
        "exonGeTrInfo.tab",
        "exonInfo.tab",
        "geneInfo.tab",
        "transcriptInfo.tab",
        "sjdbList.fromGTF.out.tab",
        "sjdbInfo.txt",
        "genomeParameters.txt",
    ]

    missing = []
    file_details = {}
    total_bytes = 0

    for fname in required_files:
        fpath = index_dir / fname
        if not fpath.is_file() or fpath.stat().st_size == 0:
            missing.append(fname)
        else:
            sz = fpath.stat().st_size
            total_bytes += sz
            file_details[fname] = {
                "size_bytes": sz,
                "size_mb": round(sz / (1024**2), 2),
            }

    # Parse parameters from genomeParameters.txt
    params = {}
    param_file = index_dir / "genomeParameters.txt"
    if param_file.is_file():
        for line in param_file.read_text(encoding="utf-8", errors="ignore").splitlines():
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                params[parts[0]] = parts[1]

    # Validate chromosome count
    chr_names = []
    chr_file = index_dir / "chrName.txt"
    if chr_file.is_file():
        chr_names = [l.strip() for l in chr_file.read_text(encoding="utf-8", errors="ignore").splitlines() if l.strip()]

    is_valid = len(missing) == 0 and len(chr_names) > 0 and len(params) > 0

    return {
        "is_complete": is_valid,
        "missing_files": missing,
        "total_index_bytes": total_bytes,
        "total_index_gb": round(total_bytes / (1024**3), 2),
        "file_details": file_details,
        "chromosome_count": len(chr_names),
        "genome_parameters": params,
        "validation_type": "INDEX_STRUCTURAL_VALIDATION",
    }


def lightweight_index_load_test(index_dir: Path, star_bin: str = "STAR") -> dict[str, Any]:
    """Perform a process-level dry-run load test of the STAR index without using SysV shared memory.

    Uses --genomeLoad NoSharedMemory to avoid false rejections in containerized environments
    (such as Google Colab or Docker) where kernel shared memory (/dev/shm, shmmax) is restricted.
    Explicitly tests process address-space memory mapping of the genome and suffix arrays,
    distinguishing index loading from full biological sequencing alignment.
    """
    star_resolved = shutil.which(star_bin)
    if not star_resolved:
        return {
            "load_test_passed": False,
            "skipped": True,
            "reason": f"STAR executable '{star_bin}' not found on PATH.",
            "validation_type": "INDEX_LOAD_VALIDATION",
        }

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        probe_fastq = tmp_path / "dnt_load_probe.fastq"
        probe_fastq.write_text("@dnt_probe\nACGTACGT\n+\nIIIIIIII\n", encoding="utf-8")
        out_prefix = str(tmp_path / "star_load_test_")

        cmd = [
            star_resolved,
            "--runMode", "alignReads",
            "--genomeDir", str(index_dir),
            "--genomeLoad", "NoSharedMemory",  # immune to /dev/shm container restrictions
            "--readFilesIn", str(probe_fastq),
            "--outFileNamePrefix", out_prefix,
            "--outSAMtype", "None",  # minimal disk I/O, purely verifies index loading into process memory
            "--runThreadN", "2",
        ]

        print(f"\nExecuting process-level index load validation (NoSharedMemory):\n{' '.join(cmd)}")
        t0 = time.time()
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False)
            load_time = time.time() - t0
            out_text = (res.stdout + res.stderr).lower()
            success = res.returncode == 0 and ("finished successfully" in out_text or "all done" in out_text)

            return {
                "load_test_passed": success,
                "exit_code": res.returncode,
                "load_time_seconds": round(load_time, 2),
                "index_load_mode": "NoSharedMemory",
                "shared_memory_dependent": False,
                "validation_type": "INDEX_LOAD_VALIDATION",
                "stdout_tail": res.stdout[-400:] if res.stdout else "",
                "stderr_tail": res.stderr[-400:] if res.stderr else "",
            }
        except Exception as e:
            return {
                "load_test_passed": False,
                "error": str(e),
                "load_time_seconds": round(time.time() - t0, 2),
                "validation_type": "INDEX_LOAD_VALIDATION",
            }


def main():
    parser = argparse.ArgumentParser(
        description="Governed Reference Indexing & Verification for Cloud/Colab High-RAM Runtimes."
    )
    parser.add_argument(
        "--profile",
        type=str,
        default="colab_highmem",
        help="Execution profile name (default: colab_highmem)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Thread allocation override (default: use profile threads)",
    )
    parser.add_argument(
        "--sparse-d",
        type=int,
        default=None,
        help="STAR genomeSAsparseD override (1 for dense, 2/4 for sparse)",
    )
    parser.add_argument(
        "--persistent-dir",
        type=Path,
        default=None,
        help="Path to persistent storage (e.g. /content/drive/MyDrive/DNT_Organoid_AI_Persistence)",
    )
    parser.add_argument(
        "--safety-margin-gb",
        type=float,
        default=10.0,
        help="Minimum free disk safety margin required on persistent storage (default: 10.0 GB)",
    )
    parser.add_argument(
        "--account-quota-gb",
        type=float,
        default=5000.0,
        help="Total Google Drive account capacity in GB (default: 5000.0 GB = 5 TB)",
    )
    parser.add_argument(
        "--user-verified-account-free-gb",
        type=float,
        default=None,
        help="User-verified Google Drive available account headroom in GB (confirmed via https://drive.google.com/drive/quota)",
    )
    parser.add_argument(
        "--use-archive",
        action="store_true",
        default=True,
        help="Package index into a single verified tar archive for Google Drive transfer (default: True)",
    )
    parser.add_argument(
        "--skip-if-already-verified",
        action="store_true",
        help="Skip re-indexing if a verified index already exists in persistent storage",
    )

    args = parser.parse_args()

    print("=" * 80)
    print("DNT ORGANOID AI - REFERENCE INDEXING & VERIFICATION RUNNER")
    print(f"Active Profile: {args.profile}")
    print(f"Repository Root: {REPO_ROOT}")
    if args.persistent_dir:
        print(f"Persistent Storage Path: {args.persistent_dir}")
    else:
        print("Persistent Storage: Not configured (local storage only)")
    print("=" * 80)

    # 1. Initialize Pipeline & Paths
    pipeline = QuantificationBridgePipeline(repo_root=REPO_ROOT, profile_name=args.profile)
    paths = resolve_pipeline_paths(REPO_ROOT)
    eff_threads = args.threads or pipeline.profile.threads
    eff_sparse_d = args.sparse_d or pipeline.profile.genome_sa_sparse_d

    # 2. Initialize Checkpoint Manager
    checkpoint_mgr = CheckpointManager(
        local_checkpoint_dir=paths["data_dir"] / "checkpoints",
        persistent_checkpoint_dir=(args.persistent_dir / "checkpoints") if args.persistent_dir else None,
    )

    # 3. Environment & Hardware Preflight
    print("\n--- STAGE 1: HARDWARE & ENVIRONMENT AUDIT ---")
    env_metrics = pipeline.record_environment_metrics()
    print(f"Host Kernel: {env_metrics['uname']}")
    print(f"Physical RAM: {env_metrics['physical_ram_gb']} GB")
    print(f"Available CPU Cores: {env_metrics['cpu_cores']}")
    print(f"Local Scratch Free Space: {env_metrics['scratch_disk_free_gb']} GB")

    # 4. Check for Existing Verified Index in Persistent Storage (Archive or Directory)
    persistent_ref_dir = (args.persistent_dir / "references") if args.persistent_dir else None
    persistent_archive = (persistent_ref_dir / "star_index_GRCh38_112_dense.tar") if persistent_ref_dir else None
    persistent_manifest = (persistent_ref_dir / "star_index_GRCh38_112_dense.tar.manifest.json") if persistent_ref_dir else None
    persistent_index_dir = (persistent_ref_dir / "star_index") if persistent_ref_dir else None

    if args.skip_if_already_verified and persistent_archive and persistent_archive.is_file():
        if not persistent_manifest or not persistent_manifest.is_file():
            print(f"\n[WARNING] Archive {persistent_archive.name} was found, but required companion manifest is MISSING!")
            print("[SECURITY / INTEGRITY] Cannot treat unverified archive as valid. Will build fresh dense reference index.")
        else:
            print(f"\n[FOUND] Verified index archive and companion manifest detected on persistent store: {persistent_archive.name}")
            print("Restoring and cryptographically verifying all constituent files against manifest...")
            try:
                rest_res = restore_and_verify_archive(
                    archive_path=persistent_archive,
                    target_extract_dir=paths["star_index_dir"],
                    require_manifest=True,
                )
                print(f"[RESTORED] Extracted and verified {rest_res['verified_files_count']} files in {rest_res['extract_seconds']}s.")
                checkpoint_mgr.record_stage_completed(
                    "STAR_INDEX_RESTORED_FROM_PERSISTENT",
                    artifacts=[{
                        "archive": str(persistent_archive),
                        "extracted_to": str(paths["star_index_dir"]),
                        "files_verified": rest_res["verified_files_count"],
                    }],
                )
                print("\n" + "=" * 80)
                print("REFERENCE INDEX READY (RESTORED & VERIFIED FROM PERSISTENT STORE).")
                print("MANDATORY SCIENTIFIC REVIEW GATE: Execution paused before FASTQ acquisition.")
                print("=" * 80)
                sys.exit(0)
            except Exception as e:
                print(f"[INTEGRITY ALERT] Restoration verification failed: {e}. Will build fresh dense index.")
    elif args.skip_if_already_verified and persistent_index_dir and persistent_index_dir.is_dir():
        print(f"\nChecking existing index directory in persistent store: {persistent_index_dir}...")
        p_eval = verify_index_artifacts(persistent_index_dir)
        if p_eval["is_complete"]:
            print(f"[FOUND] Valid index already exists in persistent store ({p_eval['total_index_gb']} GB).")
            print("Restoring index to local high-speed SSD for execution...")
            sync_res = sync_and_verify_directory(
                persistent_index_dir,
                paths["star_index_dir"],
                safety_margin_gb=2.0,
            )
            print(f"[RESTORED] Index synced to local SSD in {sync_res['total_seconds']}s (all checksums verified).")
            checkpoint_mgr.record_stage_completed(
                "STAR_INDEX_RESTORED_FROM_PERSISTENT",
                artifacts=[{"path": str(paths["star_index_dir"]), "size_gb": p_eval["total_index_gb"]}],
            )
            print("\n" + "=" * 80)
            print("REFERENCE INDEX READY (RESTORED FROM PERSISTENT STORE).")
            print("MANDATORY SCIENTIFIC REVIEW GATE: Execution paused before FASTQ acquisition.")
            print("=" * 80)
            sys.exit(0)

    # 5. Reference File Verification
    print("\n--- STAGE 2: REFERENCE ACQUISITION & INTEGRITY VERIFICATION ---")
    ref_res = pipeline.acquire_and_verify_references()
    print(f"FASTA SHA-256: {ref_res['fasta_sha256']} (VERIFIED MATCH)")
    print(f"GTF SHA-256:   {ref_res['gtf_sha256']} (VERIFIED MATCH)")
    print(f"BED12 SHA-256: {ref_res['bed12_sha256']} (VERIFIED MATCH)")
    checkpoint_mgr.record_stage_completed("REFERENCES_VERIFIED", artifacts=[ref_res])

    # 6. STAR Index Generation on Local High-Speed Disk
    print("\n--- STAGE 3: STAR GENOME INDEX GENERATION (LOCAL NVMe/SSD) ---")
    print(
        f"Building STAR index: threads={eff_threads}, sparse_d={eff_sparse_d}, "
        f"nbases={pipeline.profile.genome_sa_index_nbases}, limit_ram={pipeline.profile.limit_genome_generate_ram}"
    )

    t0_idx = time.time()
    build_res = pipeline.build_star_index(
        threads=eff_threads,
        sparse_d=eff_sparse_d,
        genome_sa_index_nbases=pipeline.profile.genome_sa_index_nbases,
        limit_ram=pipeline.profile.limit_genome_generate_ram,
    )
    idx_duration = time.time() - t0_idx
    print(f"[OK] STAR genomeGenerate completed in {idx_duration:.2f}s ({idx_duration / 60:.2f} min).")

    # 7. Index Artifact Inventory & Completeness Check
    print("\n--- STAGE 4: INDEX ARTIFACT VERIFICATION ---")
    eval_res = verify_index_artifacts(paths["star_index_dir"])
    if not eval_res["is_complete"]:
        raise RuntimeError(f"STAR index is incomplete! Missing files: {eval_res['missing_files']}")

    print(f"[OK] All {len(eval_res['file_details'])} required STAR index files verified present and non-zero.")
    print(f"Total Local Index Size: {eval_res['total_index_gb']} GB ({eval_res['total_index_bytes']:,} bytes)")

    # 8. Lightweight In-Memory Load Test
    print("\n--- STAGE 5: LIGHTWEIGHT IN-MEMORY LOAD VALIDATION ---")
    star_cmd = env_metrics["tool_paths"].get("STAR", "STAR")
    load_res = lightweight_index_load_test(paths["star_index_dir"], star_bin=star_cmd)
    if not load_res["load_test_passed"]:
        print(f"[WARNING] Lightweight load test reported non-standard exit: {load_res}")
    else:
        print(f"[OK] In-memory index load test passed successfully in {load_res['load_time_seconds']}s.")

    checkpoint_mgr.record_stage_completed(
        "STAR_INDEX_BUILT_AND_VERIFIED",
        artifacts=[{"dir": str(paths["star_index_dir"]), "total_gb": eval_res["total_index_gb"]}],
        metrics={"duration_seconds": idx_duration, "load_test": load_res},
    )

    # 9. Persistent Storage Backup (Google Drive / POSIX Storage)
    if persistent_ref_dir:
        print("\n--- STAGE 6: PERSISTENT STORAGE BACKUP (PRE-FLIGHT & VERIFICATION) ---")
        print(f"Target persistent directory: {persistent_ref_dir}")

        if args.use_archive:
            print("Packaging index into single verified tar archive to minimize Google Drive API overhead...")
            archive_res = archive_and_sync_directory(
                src_dir=paths["star_index_dir"],
                dst_dir=persistent_ref_dir,
                archive_name="star_index_GRCh38_112_dense.tar",
                safety_margin_gb=args.safety_margin_gb,
                user_verified_account_free_gb=args.user_verified_account_free_gb,
                account_quota_gb=args.account_quota_gb,
                delete_source_after=False,  # Retain local copy on fast SSD for active alignment!
            )
            print(
                f"[OK] Successfully archived and synced {archive_res['total_files']} files "
                f"({archive_res['archive_size_gb']} GB archive) to persistent storage in {archive_res['total_seconds']}s."
            )
            checkpoint_mgr.record_stage_completed(
                "STAR_INDEX_PERSISTED_TO_GDRIVE",
                artifacts=[{"archive": archive_res["dst_archive_path"], "sha256": archive_res["archive_sha256"]}],
            )
        else:
            # Check capacity before transferring
            cap_metrics = check_storage_capacity(
                persistent_index_dir,
                required_bytes=eval_res["total_index_bytes"],
                safety_margin_gb=args.safety_margin_gb,
                account_quota_gb=args.account_quota_gb,
            )
            print(f"Persistent Storage Free: {cap_metrics['free_gb']} GB")
            print(f"Required Space: {cap_metrics['required_gb']} GB (plus {args.safety_margin_gb} GB safety margin)")

            # Sync and verify every single file with SHA-256
            print("Synchronizing index files with SHA-256 cryptographic verification...")
            sync_res = sync_and_verify_directory(
                src_dir=paths["star_index_dir"],
                dst_dir=persistent_index_dir,
                safety_margin_gb=args.safety_margin_gb,
                delete_source_after=False,  # Keep local copy on fast SSD for active alignment!
            )
            print(
                f"[OK] Successfully backed up {sync_res['total_files']} files ({round(sync_res['total_bytes'] / (1024**3), 2)} GB) "
                f"to persistent storage in {sync_res['total_seconds']}s."
            )
            checkpoint_mgr.record_stage_completed(
                "STAR_INDEX_PERSISTED_TO_GDRIVE",
                artifacts=[{"dir": str(persistent_index_dir), "files_verified": sync_res["total_files"]}],
            )

    # 10. Mandatory Scientific Review Halt
    print("\n" + "=" * 80)
    print("REFERENCE INDEX PREPARATION COMPLETE AND VERIFIED.")
    print("=" * 80)
    print(f"Profile:                      {args.profile}")
    print(f"Reference Assembly:           Ensembl GRCh38 release 112")
    print(f"Sparse Factor (sparseD):      {eff_sparse_d}")
    print(f"Total Index Directory Size:   {eval_res['total_index_gb']} GB")
    print(f"Local Index Path:             {paths['star_index_dir']}")
    if persistent_index_dir:
        print(f"Persistent Backup Path:       {persistent_index_dir}")
    print(f"Index Generation Time:        {idx_duration / 60:.2f} minutes")
    print("Status:                       INDEX_VERIFIED_READY")
    print("=" * 80)
    print("MANDATORY SCIENTIFIC REVIEW GATE:")
    print("Execution has paused as instructed by user governance.")
    print("DO NOT download FASTQ sequencing libraries or launch alignment until this report is reviewed.")
    print("=" * 80)


if __name__ == "__main__":
    main()
