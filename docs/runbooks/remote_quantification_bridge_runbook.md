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

## 2. Governed Non-Contiguous 12-Run Pilot Manifest

The pilot manifest is governed strictly in [`config/quantification_bridge/pilot_manifest.json`](../../config/quantification_bridge/pilot_manifest.json) and [`src/quantification_bridge/manifest.py`](../../src/quantification_bridge/manifest.py).

> [!CAUTION]
> **PROHIBITION OF ACCESSION-RANGE GENERATION:**
> The 12-run pilot manifest is **non-contiguous**. It must NEVER be generated or derived via `range(...)`, numeric accession intervals, contiguous assumptions, or sorting assumptions. Deriving or displaying the pilot as contiguous ranges (e.g. `SRR17834920–SRR17834925` or `SRR17834914–SRR17834919`) is strictly prohibited.

### Authoritative 12-Sample Allowlist
| Run Accession | Sample Alias | Condition Group | VPA Dose | Dev Day | BioSample | Experiment | Read Pairs |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **SRR17834920** | VPA_0_1_1 | control | 0 µM | Day 28 | SAMN25510021 | SRX13996121 | 21,899,521 |
| **SRR17834908** | VPA_0_2_1 | control | 0 µM | Day 28 | SAMN25510023 | SRX13996133 | 23,152,105 |
| **SRR17834897** | VPA_0_3_1 | control | 0 µM | Day 28 | SAMN25510025 | SRX13996144 | 22,237,894 |
| **SRR17834886** | VPA_0_4_1 | control | 0 µM | Day 28 | SAMN25510027 | SRX13996155 | 22,600,697 |
| **SRR17834875** | VPA_0_5_1 | control | 0 µM | Day 28 | SAMN25510029 | SRX13996166 | 22,052,958 |
| **SRR17834919** | VPA_0_6_1 | control | 0 µM | Day 28 | SAMN25510031 | SRX13996122 | 23,085,831 |
| **SRR17834905** | VPA_250_1_1 | treatment | 250 µM | Day 28 | SAMN25510045 | SRX13996136 | 22,634,745 |
| **SRR17834904** | VPA_250_2_1 | treatment | 250 µM | Day 28 | SAMN25510047 | SRX13996137 | 20,783,074 |
| **SRR17834903** | VPA_250_3_1 | treatment | 250 µM | Day 28 | SAMN25510049 | SRX13996138 | 21,073,965 |
| **SRR17834902** | VPA_250_4_1 | treatment | 250 µM | Day 28 | SAMN25510051 | SRX13996139 | 21,110,376 |
| **SRR17834901** | VPA_250_5_1 | treatment | 250 µM | Day 28 | SAMN25510053 | SRX13996140 | 20,605,954 |
| **SRR17834900** | VPA_250_6_1 | treatment | 250 µM | Day 28 | SAMN25510055 | SRX13996141 | 22,878,898 |

*Biological Replicate Unit:* `UNVERIFIED` across all samples.

---

## 3. Remote Node Minimum Specifications

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

## 4. Environment Provisioning on Remote Node

### Option A: Explicit Lockfile Reproduction (Recommended)
Reproduce the exact bit-for-bit validated bioinformatics toolchain:
```bash
git clone -b feature/prjna802072-quantification-bridge <REPO_URL> DNT-Organoid-AI
cd DNT-Organoid-AI

conda create --name dnt_bridge --file config/quantification_bridge/environment_bridge.lock.txt -y
conda activate dnt_bridge
```

### Option B: Declarative YAML Environment
```bash
conda env create -f config/quantification_bridge/environment_bridge.yml
conda activate dnt_bridge
```

### Formally Governed Toolchain Versions
- Python: `3.11.9`
- fastp: `0.23.4`
- STAR: `2.7.11b`
- Subread / featureCounts: `2.0.6`
- samtools: `1.19.2` (htslib 1.21)
- RSeQC: `5.0.3` (`infer_experiment.py 5.0.2/5.0.3`)
- UCSC tools: `gtfToGenePred` (release 469), `genePredToBed` (release 469)

---

## 5. Remote Execution Workflow

### Step 1: Preflight Hardware & Unit Test Verification
Verify node hardware and execute the complete synthetic preflight suite:
```bash
python scripts/run_remote_quantification_bridge.py --preflight-only
```
*Verification requirement:* All preflight tests (Tests A–G and manifest regression) must pass (100% green).

### Step 2: Acquisition Dry-Run Verification (No Downloads)
Verify that the remote script targets exactly the 24 FASTQ files for the 12 approved runs:
```bash
python scripts/run_remote_quantification_bridge.py --dry-run-download
```
*Requirement:* Asserts 0 missing, 0 unexpected against `pilot_manifest.json`.

### Step 3: Full Upstream Bridge Execution Through Strandedness Gate
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

## 6. Expected Output Artifact Inventory

Upon reaching the strandedness decision gate, the following primary evidence artifacts will exist:

```text
data/interim/prjna802072_bridge/
└── prjna802072_strandedness_evidence.json    <-- AUTHORITATIVE EVIDENCE GATE ARTIFACT

data/interim/prjna802072_alignments/
├── SRR17834920_Aligned.sortedByCoord.out.bam
├── SRR17834920_Aligned.sortedByCoord.out.bam.bai
├── SRR17834908_Aligned.sortedByCoord.out.bam
├── SRR17834908_Aligned.sortedByCoord.out.bam.bai
├── SRR17834897_Aligned.sortedByCoord.out.bam
├── SRR17834897_Aligned.sortedByCoord.out.bam.bai
├── SRR17834886_Aligned.sortedByCoord.out.bam
├── SRR17834886_Aligned.sortedByCoord.out.bam.bai
├── SRR17834875_Aligned.sortedByCoord.out.bam
├── SRR17834875_Aligned.sortedByCoord.out.bam.bai
├── SRR17834919_Aligned.sortedByCoord.out.bam
├── SRR17834919_Aligned.sortedByCoord.out.bam.bai
├── SRR17834905_Aligned.sortedByCoord.out.bam
├── SRR17834905_Aligned.sortedByCoord.out.bam.bai
├── SRR17834904_Aligned.sortedByCoord.out.bam
├── SRR17834904_Aligned.sortedByCoord.out.bam.bai
├── SRR17834903_Aligned.sortedByCoord.out.bam
├── SRR17834903_Aligned.sortedByCoord.out.bam.bai
├── SRR17834902_Aligned.sortedByCoord.out.bam
├── SRR17834902_Aligned.sortedByCoord.out.bam.bai
├── SRR17834901_Aligned.sortedByCoord.out.bam
├── SRR17834901_Aligned.sortedByCoord.out.bam.bai
├── SRR17834900_Aligned.sortedByCoord.out.bam
└── SRR17834900_Aligned.sortedByCoord.out.bam.bai

data/interim/prjna802072_qc_reports/
├── SRR17834920_fastp.html / json
├── ... (individual fastp reports for all 12 runs)
└── SRR17834920_infer_experiment.txt ... (individual RSeQC reports for all 12 runs)
```

---

## 7. Storage Preservation Policy
All primary data and intermediate artifacts must be preserved on the $\ge 500\text{ GB}$ scratch storage throughout the Owner strandedness decision:
- `data/raw/prjna802072_fastqs/` (24 raw FASTQs, ~36 GB)
- `data/interim/prjna802072_trimmed/` (24 trimmed FASTQs, ~30 GB)
- `references/ensembl_112/` (FASTA, GTF, BED12, ~5 GB)
- `references/ensembl_112/star_index/` (STAR index, ~29 GB)
- `data/interim/prjna802072_alignments/` (12 coordinate-sorted BAMs + BAIs, ~60 GB)
- `data/interim/prjna802072_qc_reports/` (fastp JSON/HTML, RSeQC logs, STAR logs, ~1 GB)
- `data/interim/prjna802072_bridge/prjna802072_strandedness_evidence.json`

Do not delete or prune any intermediate files until the Project Owner reviews the strandedness evidence and authorizes downstream featureCounts execution.
