# PRJNA802072 — Full FASTQ/SRA real-data validation (Colab)

**Review status:** Pre-execution verification ready; **no claim of FASTQ or 12-library execution**.  
**Branch:** `feature/fastq-sra-real-data-validation`  
**Reference index:** Ensembl GRCh38.112; STAR 2.7.11b; dense index `genomeSAsparseD=1`, `genomeSAindexNbases=14`.  
**Verified reference archive SHA-256 (owner audit):** `756ed0a395cd477db730f9f9a4199c77576da59fb5aab7169fa9a1a33dd56392`.

## Goal

Validate that **all 12 named paired-end libraries** (6 vehicle controls, 6 VPA 250 µM) can be acquired and processed with reproducible, verifiable outputs. Full FASTQ/SRA validation is **not** the same as completion of U3.2 or DNT prediction.

## Mandatory entry gate — no FASTQ download

In an already prepared Colab runtime that has run the reference notebook through Stage 9:

```bash
cd /content/DNT-Organoid-AI
git status --short
git pull --ff-only origin feature/fastq-sra-real-data-validation
python scripts/validate_fastq_sra_full_preflight.py --report-path data/interim/PRJNA802072_full_preflight.json
```

This checks the governed 12-run/24-FASTQ allowlist, repeats an **independent live ENA** metadata/MD5 check, verifies the **existing** dense index files, and checks actual remaining local disk space. It does not download or align data. The ENA file report is retrieved dynamically; do not invent or substitute missing read files. A failed check blocks downstream work.

**If the Colab session was reset:** restore the existing verified STAR archive using the Stage 7 restoration procedure; do **not** regenerate it by default. A checkpoint manifest alone is not a substitute for index files.

## Full-validation execution gates

| Gate | Input | Required evidence | Action |
|---|---|---|---|
| F0 | Existing index and 12-run manifest | Live ENA metadata, 24 FASTQ filenames + MD5, paired layout, space, index files | Stop on any mismatch |
| F1 | 2 technical pilot libraries, one per condition | R1/R2 MD5, gzip/FASTQ validity, fastp reports, STAR alignment stats, BAM/BAI integrity, per-library RSeQC | Review before proceeding to all 12 |
| F2 | Remaining 10 libraries | Same per-library validation, processing **sequentially** and checking persistence after each | Advance library-by-library |
| F3 | All 12 mapped BAMs and RSeQC | Per-library strand fractions, concordance, proportion ambiguous, mapping and read QC; saved provenance | **STOP** at `AWAITING_OWNER_STRANDEDNESS_DECISION` |
| F4 | Owner-approved strandedness setting | featureCounts v2.0.6 outputs, annotation version, BAM inputs, gene-count conservation, assignment QC | Requires **separate explicit strand decision**; quantification integration is not yet complete |
| F5 | All approved raw counts and verified biological-replicate metadata | Cohort-level QC, gene harmonization, DESeq2, pathway analysis | Not a validation claim until actually executed and reviewed |

## Why NOT run the existing one-command 12-library runner unchanged

The existing `QuantificationBridgePipeline.run_through_strandedness_gate()`:
1. Downloads **all 24 FASTQs** first, then trims **all 12 libraries**, and retains the intermediates.
2. Always calls `build_star_index()`, unnecessarily rebuilding the previously verified dense index.
3. Does not implement verified **per-library** Google Drive persistence and safe reclamation of local intermediates.

The repository execution profile estimates up to ~270 GB peak scratch for the full unbatched experiment. The previous Colab instance provided roughly 235 GB *total* scratch before index, references and other outputs. Do not run the unmodified all-at-once command on that instance.

## Safe processing policy

- Process one FASTQ pair / library at a time; keep a full manifest of statuses.
- Verify original ENA MD5 before reading data.
- Record read counts, fastp QC, STAR `Log.final.out`, BAM integrity/index, RSeQC raw output and interpretation.
- Upload each completed library's critical artifacts with cryptographic verification using an authorized Google Drive API transfer before deleting any local intermediates. Google Drive FUSE file copy alone is not sufficient for a >GB archive without verified checksums.
- Keep only data required for remaining analysis. Never delete the reference index. Do not delete unique unpersisted results.
- Record software versions, job settings, reference build SHA, timestamp and per-file checksums.
- Confirm scientific biological-replicate units before treating the six libraries per group as independent replicates.

## Key code

- `config/quantification_bridge/pilot_manifest.json`: governed 12-run allowlist.
- `scripts/reverify_pilot_manifest.py`: independent live ENA comparison.
- `scripts/validate_fastq_sra_full_preflight.py`: read-only entry gate.
- `src/quantification_bridge/pipeline.py`: existing acquisition and alignment functions.
- `src/quantification_bridge/drive_api_persistence.py`: verified archive persistence (already used for the reference index).
- `notebooks/colab_fastq_sra_bridge.ipynb`: reference preparation only; not yet the full 12-run processing notebook.

## Decision

**READY TO RUN THE READ-ONLY PREFLIGHT.** This does not certify live data processing, 12-run mapping, gene counting, differential expression, or a trained DNT model. These require actual runtime logs and validated results.
