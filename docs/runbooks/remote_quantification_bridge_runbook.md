# PRJNA802072 Remote Compute & HPC Execution Runbook

## 1. Scope & Objective
This runbook governs the offload and execution of the **PRJNA802072 Raw-Read Quantification Bridge** on an external server, HPC cluster, or cloud compute node.

Due to local hardware memory constraints (16 GB physical host RAM vs. ~32 GB RAM required for human GRCh38 STAR indexing) and disk storage limits (186.91 GB free vs. ~207.54 GB peak concurrent footprint), this pipeline is executed on an appropriate high-memory Linux node.

The execution workflow covers:
$$\text{ENA FASTQs} \longrightarrow \text{fastp QC} \longrightarrow \text{STAR Alignment} \longrightarrow \text{RSeQC Diagnostics} \longrightarrow \mathbf{AWAITING\_OWNER\_STRANDEDNESS\_DECISION}$$

> [!IMPORTANT]
> **HARD STOP AT STRANDEDNESS GATE:**
> The automated runbook terminates after generating `prjna802072_strandedness_evidence.json`. It **MUST NOT** proceed to `featureCounts` or frozen downstream pipeline stages without explicit Project Owner authorization.

---

## 2. Remote Node Minimum Specifications

The target execution node must satisfy the governed specifications in [`config/quantification_bridge/remote_compute_spec.yaml`](../../config/quantification_bridge/remote_compute_spec.yaml):

- **Architecture:** Linux x86_64 (Ubuntu 22.04/24.04 LTS or Rocky Linux 9 / RHEL 9).
- **Physical RAM:** **$\ge 64\text{ GB}$ physical memory**.
  - *Policy:* Virtual swap memory is strictly prohibited as a substitute for physical RAM.
- **CPU:** **$\ge 16\text{ physical cores}$** (dedicated compute).
- **Scratch Storage:** **$\ge 500\text{ GB}$ high-speed local NVMe SSD**.
  - Filesystem must support large files (>50 GB BAMs) with POSIX file locking.
- **Network:** Unrestricted outbound HTTPS to:
  - `ftp.sra.ebi.ac.uk` / `www.ebi.ac.uk` (ENA/SRA downloads)
  - `ftp.ensembl.org` (Ensembl references)
  - `conda.anaconda.org` (bioinformatics package channels)

---

## 3. Environment Provisioning on Remote Node

### Option A: Explicit Lockfile Reproduction (Recommended)
Reproduce the exact bit-for-bit validated bioinformatics toolchain:
```bash
# Clone the governed feature branch
git clone -b feature/prjna802072-quantification-bridge <REPO_URL> DNT-Organoid-AI
cd DNT-Organoid-AI

# Create conda environment from explicit lockfile
conda create --name dnt_bridge --file config/quantification_bridge/environment_bridge.lock.txt -y
conda activate dnt_bridge
```

### Option B: Declarative YAML Environment
```bash
conda env create -f config/quantification_bridge/environment_bridge.yml
conda activate dnt_bridge
```

### Pinned Tool Versions
- Python: `3.11.9`
- fastp: `0.23.4`
- STAR: `2.7.11b`
- Subread / featureCounts: `2.0.6`
- samtools: `1.19.2` (htslib 1.21)
- RSeQC: `5.0.3` (`infer_experiment.py 5.0.2/5.0.3`)
- UCSC tools: `gtfToGenePred` (v469), `genePredToBed` (v469)

---

## 4. Remote Execution Workflow

### Step 1: Preflight Verification
Verify node hardware and execute the complete synthetic preflight suite:
```bash
python scripts/run_remote_quantification_bridge.py --preflight-only
```
*Verification requirement:* All 5 test modules covering Tests A–G must pass (100% green).

### Step 2: Full Upstream Bridge Execution Through Strandedness Gate
```bash
chmod +x scripts/run_remote_quantification_bridge.sh
./scripts/run_remote_quantification_bridge.sh 16
```
*(Or invoke via python: `python scripts/run_remote_quantification_bridge.py --threads 16`)*

### Sequential Stages Executed:
1. **Environment Audit:** Hardware metrics, OS, and tool binaries logged.
2. **Preflight Test Suite:** Synthetic BAM, adapter conservation, and immutability verified.
3. **Reference Acquisition & Verification:**
   - FASTA: SHA-256 verified against `d8c3af0094a7bba6125763bad779ec18a81483c739c6ed122094bdf86c187b92`.
   - GTF: SHA-256 verified against `8c87436bab973c871447743c2f5cbb2b557be0523fe6efe3ff33307a4c836002`.
   - BED12: Derived deterministically and verified against `4d1231bb5fcb25731dac508592dfbf78ef04db041034b006d15a53209c309a06`.
4. **Library Acquisition Directly from ENA:**
   - 24 FASTQ files downloaded via HTTPS with retry.
   - Published ENA MD5s verified; local SHA-256 and byte sizes recorded.
5. **FASTQ Structural Synchronization:**
   - Decompressibility, header matching, quality string length equality confirmed.
6. **Treatment-Blind fastp Preprocessing:**
   - Poly-G and adapter trimming executed. Quality and length filters disabled per governance.
7. **STAR Genome Index Generation:**
   - Human GRCh38 primary assembly index generated (`sjdbOverhang 149`).
8. **STAR Alignment & BAM Indexing:**
   - 12 libraries aligned to coordinate-sorted BAMs with governed flags.
   - SAMtools indices (`.bai`) generated.
9. **RSeQC Strandedness Diagnostic:**
   - `infer_experiment.py` executed across all 12 BAMs against Ensembl 112 BED12.
10. **Evidence Serialization:**
    - Output written to `data/interim/prjna802072_bridge/prjna802072_strandedness_evidence.json`.
11. **HALT:** Execution stops and reports `AWAITING_OWNER_STRANDEDNESS_DECISION`.

---

## 5. Storage Preservation Policy
All primary data and intermediate artifacts must be preserved on the $\ge 500\text{ GB}$ scratch storage throughout the Owner strandedness decision:
- `data/raw/prjna802072_fastqs/` (24 raw FASTQs, ~36 GB)
- `data/interim/prjna802072_trimmed/` (24 trimmed FASTQs, ~30 GB)
- `references/ensembl_112/` (FASTA, GTF, BED12, ~5 GB)
- `references/ensembl_112/star_index/` (STAR index, ~29 GB)
- `data/interim/prjna802072_alignments/` (12 coordinate-sorted BAMs + BAIs, ~60 GB)
- `data/interim/prjna802072_qc_reports/` (fastp JSON/HTML, RSeQC logs, STAR logs, ~1 GB)
- `data/interim/prjna802072_bridge/prjna802072_strandedness_evidence.json`

Do not delete or prune any intermediate files until the Project Owner reviews the strandedness evidence and authorizes downstream featureCounts execution.
