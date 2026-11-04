# External-Laboratory Reproducibility Contract — DNT Organoid AI
**Status:** Implementation requirements and acceptance gates; **NOT** external validation or regulatory qualification.  
**Scope:** Raw FASTQ/SRA → bulk gene counts → differential expression → molecular features → independently evaluated DNT model.  
**Initial benchmark:** PRJNA802072 (12 paired-end sequencing libraries, 6 vehicle controls and 6 VPA 250 µM); replicate independence still needs verification.  
**Current code branch:** `feature/fastq-sra-real-data-validation`.

## 1. Non-negotiable architecture

**Scientific backend is independent of execution interface.** Colab is an execution environment, not a source-of-truth pipeline. The future website, a command-line tool, and an independent laboratory's compute cluster must call the *same versioned backend*, and use the same validation logic.

Separation of responsibilities:

| Layer | Responsibility | Current reference / boundary |
|---|---|---|
| Data ingestion | Recognize FASTQ/SRA, counts, DEG tables and metadata without guessing | `src/routing/`, `src/ingestion/`, `src/metadata/` |
| FASTQ/SRA quantification | Acquire/verify reads; fastp; STAR; RSeQC; strand review; featureCounts | `src/quantification_bridge/` — **partially implemented**, stops after RSeQC |
| Bulk statistics | Count QC, gene harmonization, normalization, matched contrasts, DESeq2 | `src/qc/`, `src/gene_harmonization/`, `src/normalization/`, `src/response_builder/`, `src/differential_expression/` |
| Molecular features | GO/KEGG, cross-study response, networks | `src/features/` — only parts implemented |
| DNT prediction | Curated reference labels, leak-safe splits, signatures, training, calibration, untouched test | **NOT implemented** |
| Operations | Job lifecycle, recovery, resource estimation, verified storage, logs and independent reruns | **To be completed** |
| Interfaces | CLI, optional notebook, future website API/UI | UI must not duplicate analysis/statistics |

**Scientific acceptance must not depend on a website or a Colab notebook.**

## 2. Required portable backend contract

Each job submission, whether from CLI or future website, must explicitly contain:
- `job_id`, `study_id`, `pipeline_version`, `requested_route`, `reference_version`
- Input asset references, format, source accession, file sizes and independent checksums
- Sample identifiers, exact assay and library details, biological replicate unit, control/treatment mapping
- Cell/organism, organoid developmental age at exposure and collection, dose, duration and units
- Requested operation scope and explicit user/scientific approvals for gated stages
- Execution resource profile (CPU/RAM/disk), storage adapter and output location
- Data provenance, source access permissions, versioned parameter set

Each job must provide machine-readable state:
`CREATED → PREFLIGHT_VERIFIED → INPUTS_VERIFIED → FASTP_COMPLETE → ALIGNMENT_COMPLETE → STRANDNESS_REVIEW_REQUIRED → COUNTS_COMPLETE → BULK_ANALYSIS_COMPLETE → SCIENTIFIC_REVIEW_COMPLETE`.
Failures are explicit terminal or resumable states, with diagnostic reasons. States are reported honestly, **never inferred from the mere presence of an output file**.

**API compatibility requirement:** The future website can query job state, inputs, validated stage results and report locations. No UI component is allowed to silently select biological replicates, strandedness modes, controls, DNT labels or statistical models.

## 3. Required reproducibility evidence per completed stage

Store (as versioned JSON/TSV/HTML outputs as appropriate):
1. Exact Git commit SHA / release tag, including confirmation of the working-tree state.
2. Operating system, Python/R/Bioconductor and executable/package versions; pinned dependencies/locks.
3. Exact tool commands, parameters, exit codes, and non-sensitive stdout/stderr logs.
4. Reference FASTA/GTF/STAR index versions, source URLs, file checksums and index-generation parameters.
5. Original input accession, FASTQ read-pair MD5, sample metadata, preprocessing policy and QC findings.
6. Per-library raw read QC, fastp reports, STAR `Log.final.out`, sorted BAM and BAI, `samtools quickcheck` output, RSeQC fractions, checksum and byte count of each produced artifact.
7. Approved featureCounts strandedness (0/1/2), exact GTF, assignment statistics, resulting integer count matrix, sample-to-column mapping and numerical conservation audit.
8. QC decisions, normalized size factors, contrasts, replicate model, DESeq2 version/params, results and excluded-feature ledger; downstream pathway-analysis parameters and retained gene universe.
9. Pipeline status, timestamps, approval decisions, aborted stages, rerun attempts, provenance lineage and material source-data differences.
10. Cryptographically verified, **accessible** persistent data/report artifacts outside temporary Colab scratch.

Do not commit raw FASTQs, BAMs, large reference indices, clinical or confidential sample data to Git. Git contains the **reproducible code, tests, pinned environments, SOPs, execution manifests, small fixtures and audit reports**; full datasets and outputs belong in approved controlled storage and are referenced by stable IDs, checksums and access instructions.

## 4. Execution, persistence and resource safety

- Run a two-library controlled pilot to check actual fastp/STAR/RSeQC performance before expanding to 12.
- Run 12 libraries sequentially or in bounded batches based on measured free disk and memory, not estimated total storage. A local peak estimate of ~270 GB exceeds the previously seen ~235 GB Colab scratch.
- Reuse the independently SHA-256-verified STAR index. Never rebuild it merely because another pipeline stage starts.
- Persist each library's required artifacts with remote byte-size + checksum verification, and record its manifest **before** removing regenerable local intermediates.
- Never delete the only valid copy of a BAM, QC report or count matrix.
- Do not proceed from strandedness to featureCounts without a source-specific approved strand decision; reject inconsistent/ambiguous libraries for review.
- Support local POSIX/HPC storage and cloud storage through explicit storage interfaces, rather than embedding Google Drive-only logic in scientific modules.
- Separate scientific approval from runtime-execution permission. The successful preflight is not proof of full analysis.

## 5. External laboratory / NIH review readiness — acceptance matrix

| Requirement | Evidence required | Status at contract creation |
|---|---|---|
| Same versioned program runs on different labs' machines | Documented CLI, pinned containers/environments, verified input fixtures | Partial |
| Public benchmark 12/12 libraries processed | Independently verified ENA reads, per-library reports, complete signed-off counts | **Not yet tested** |
| Resilience to Colab/compute interruption | Restart/replay with independently verified checkpoints and no lost outputs | Partial (reference archive only) |
| Correct strand and featureCounts | Per-sample RSeQC evidence, approved strand selection, matched count assignment | **Not implemented end to end** |
| Counts→DESeq2→pathways integration | Reproducible run with qualified biological replication and QC | Existing modules, integration pending |
| Independent second environment run | Different machine/lab using same versions/inputs; comparison report and expected variation | Not started |
| Cross-study synthesis and DNT model evaluation | Predefined biological comparability, independent labels, leakage-proof split, untouched tests | Not completed |
| Website readiness | Stable backend job schema, deterministic result models, authentication/permissions and API adapter | Architectural contract only |
| External NIH/lab scientific review | Reviewer-supplied protocol, independent findings, tracked remediation | **Not performed** |

Do **not** describe any item as NIH-validated, clinically validated, or a regulatory-approved method based only on internal tests.

## 6. Immediate implementation order (must be code-backed and tested)

1. **Safeguard full cohort:** run `scripts/validate_fastq_sra_full_preflight.py` in an existing verified Colab runtime; verify all 12 run metadata against live ENA before downloads.
2. **Build a resumable stage executor and machine-readable job manifest** in `src/quantification_bridge/`, not only notebook cells. Add unit tests for failed downloads, checksum mismatches, interrupted execution, insufficient disk, and duplicate samples.
3. **Two-run real-data pilot:** verify FASTQ → fastp → STAR sorted BAM → RSeQC, preserve artifacts and provenance. Do not claim success until actual run evidence exists.
4. **Full 12-library per-library validation:** verified upload/persistence and restart support before scaling.
5. **Add featureCounts backend** with explicit scientific strandedness approval and gene-level integer-count conservation; connect to the existing bulk runner.
6. **Reproduce bulk DESeq2 and GO/KEGG outputs** with verified biological replicates.
7. **Cross-platform external rerun** with controlled fixture and complete software/environment pins.
8. **Expose a stable job API for the future web UI** after the scientific backend passes integration acceptance.

**Code review / reproducibility evidence is necessary but does not by itself establish biological validity.** Scientific and model validation must be independently documented.
