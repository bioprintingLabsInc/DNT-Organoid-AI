# PRJNA802072 Day 28 — Scientific Review Gate

**Review date:** 2026-11-04  
**Repository:** `bioprintingLabsInc/DNT-Organoid-AI`  
**Review branch:** `feature/fastq-sra-real-data-validation`  
**Reported reference audit commit:** `ccbb524cd22650504c79c8b958c2736c976860af`  
**Decision:** **CONDITIONAL GO for a 2-library technical demonstration only; no execution authorization recorded here.**  
**Not approved:** differential-expression claims, model training, or expansion to all 12 libraries without separate scientific approval.

## Purpose

Determine whether the completed GRCh38/Ensembl 112 STAR reference-index work is ready to support a two-library pilot of raw FASTQ ingestion, preprocessing, alignment, and library strandedness assessment. No sequencing data were downloaded or processed during this review.

## Evidence reviewed

1. Owner-provided STAR reference reproducibility audit: STAR 2.7.11b; dense index `genomeSAsparseD=1`, `genomeSAindexNbases=14`, `sjdbOverhang=149`; 28.19 GB; 14 required artifacts; index-load verification passed; Google Drive archive SHA-256 `756ed0a395cd477db730f9f9a4199c77576da59fb5aab7169fa9a1a33dd56392`; Google Drive archive and checkpoint reported verified.
2. Repository `config/quantification_bridge/pilot_manifest.json`, `src/quantification_bridge/manifest.py`, `config/quantification_bridge/execution_profiles.yaml`, `src/quantification_bridge/pipeline.py`, `src/quantification_bridge/validator.py`, `src/quantification_bridge/rseqc_runner.py`, `src/quantification_bridge/fastp_runner.py`.
3. Zang et al., *Valproic acid exposure decreases neurogenic potential of outer radial glia in human brain organoids*, Frontiers in Molecular Neuroscience (2022), DOI: 10.3389/fnmol.2022.1023765, https://pmc.ncbi.nlm.nih.gov/articles/PMC9744776/.
4. ENA reference documentation for programmatic run/file reports, including run-specific FASTQ locations and checksums: https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html

**Evidence limitation:** The audit is owner-provided output, not an independent re-download of the archive. Live ENA run metadata/FASTQ MD5 could not be independently retrieved in this review; the two-run accession/file details below are from the repository, not separately confirmed by the reviewer.

## Study-level metadata

- Study: BioProject `PRJNA802072`, SRA study `SRP357550`.
- Biological material: **H9 human embryonic stem-cell-derived dorsal forebrain organoids**, not donor iPSC-derived organoids; bulk mRNA sequencing.
- Developmental age at collection: **Day 28**.
- Exposure: Valproic acid (VPA) **250 µM**, compared with **0 µM distilled-water vehicle control**, exposure from Day 16 to Day 28 (approximately **12 days**).
- Study reports six Day 28 samples per group for its biological analyses; **repository biological replicate unit remains `UNVERIFIED`**. Do not assume individual run equals independent organoid or independent differentiation batch until source metadata and sample preparation are checked.
- Organism: *Homo sapiens*; reference GRCh38 primary assembly with Ensembl release 112 annotations.

## Authorized technical-pilot scope proposed for owner approval

| Field | Control | VPA-treated |
|---|---|---|
| Run | `SRR17834920` | `SRR17834905` |
| Sample alias | `VPA_0_1_1` | `VPA_250_1_1` |
| Expected paired read count (manifest) | 21,899,521 | 22,634,745 |
| Dose | 0 µM | 250 µM |
| Collection | Day 28 | Day 28 |
| Replicate independence | Unverified | Unverified |

**Technical-pilot success definition:** verify exact input bytes and paired read structure; produce fastp QC reports; STAR-aligned sorted BAM/BAI and mapping metrics; run RSeQC strandedness; store provenance and halt at the strandedness decision. A one-run-per-condition contrast is **not a valid basis for differential-expression inference**.

## Conditions before launching FASTQ downloads

- [ ] Independently retrieve ENA `read_run` metadata for both runs: accession, layout, library strategy, instrument, expected FASTQ paths, file sizes and MD5s; compare to manifest and resolve discrepancies **before downloading**.
- [ ] Confirm BioSample/experiment metadata for cell source, organoid model, exposure window, RNA library protocol and replicate unit, preserving unknown fields explicitly.
- [ ] Confirm the verified archive is restorable from Google Drive in the current Colab runtime with SHA-256 and constituent-file verification; do not regenerate a working index unnecessarily.
- [ ] Confirm the two-run-only download function will be used; **do not invoke `run_through_strandedness_gate()`**, which currently acquires all 12 libraries and may rebuild an index.
- [ ] Recheck local scratch, memory, Google Drive quotas, active binary versions, and per-step output persistence before transfer.
- [ ] Obtain explicit owner authorization for the two-run pilot acquisition and alignment.

## Planned method and QC policy

1. Acquire just four paired FASTQ.gz files for the two named runs, verifying MD5 against independently refreshed ENA file reports.
2. Validate gzip integrity, four-line FASTQ structure, total mate-pair counts, and paired-read identifiers. The current validator performs intensive per-record read-pair checks only for the first 100,000 pairs; extended read-ID verification may be warranted.
3. Run fastp 0.23.4 with documented adapter/poly-G trimming; quality and length filters are intentionally disabled in the current wrapper. Report Q30, adapter burden, and surviving read fraction without silently imposing new thresholds.
4. Align to the verified STAR 2.7.11b GRCh38/Ensembl 112 dense index; retain STAR mapping statistics, BAM indexes, tool/parameter provenance.
5. Run RSeQC `infer_experiment.py` on each BAM (configured sample size 200,000, MAPQ minimum 10). Interpret strand direction from each sample and compare them; do **not** select featureCounts strand mode prematurely.
6. Halt at `AWAITING_OWNER_STRANDEDNESS_DECISION` before featureCounts, comparative expression analysis, or model ingestion.

The execution profiles record **advisory** warnings including Q30 fraction below 0.85, read survival below 0.90, unique alignment below 0.75, multimapping above 0.15, and featureCounts assignment below 0.60. These are review flags, not proof of failure, and do not replace inspection of assay-specific distributions. Missing/invalid FASTQs, zero mapped reads, and conflicting strandedness warrant stopping and investigating.

## Scientific-review conclusion

The reported reference indexing and index-load checks are sufficient technical groundwork for a **small, gated pipeline pilot**. Study identity and Day 28 VPA exposure design are supported by the publication and repository. **Remaining pre-execution checks** are independent ENA run metadata verification, replicate interpretation, and confirmation of storage/persistence and exact two-run execution.

**Decision remains:** `AWAITING_OWNER_PILOT_FASTQ_AUTHORIZATION`. A user statement agreeing to scientific review is **not** recorded as permission to download FASTQs or start alignment.
