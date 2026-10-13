# Project Overview: DNT Organoid AI Platform

## 1. Scientific Mission & Motivation

**Developmental Neurotoxicity (DNT)** represents one of the most critical challenges in chemical safety assessment. Traditional regulatory testing relies heavily on rodent developmental studies, which suffer from substantial interspecies physiological differences, high costs, low throughput, and ethical limitations.

Human pluripotent stem-cell-derived brain organoids (e.g., human cerebral cortex organoids, SOSR-COs) provide physiologically relevant, 3D multicellular neurodevelopmental models with human genetic backgrounds. However, systematic computational pipelines are needed to translate genome-wide transcriptomic profiling of chemically exposed organoids into robust, generalizable toxicological insights.

The **DNT Organoid AI** platform is engineered to:
1. Intake diverse public and in-house RNA-seq datasets across disparate labs and protocols;
2. Apply strict quality control and canonical gene harmonization without losing low-expression genes;
3. Construct scientifically rigorous treatment-versus-matched-control contrasts;
4. Estimate differential molecular responses using authoritative Bioconductor statistics (DESeq2);
5. Extract unbiased pathway- and signaling-level feature representations;
6. Build machine learning models for DNT risk prediction using leak-free, compound-safe validation.

---

## 2. Core Architectural Principles

### Principle 1: Zero Information Leakage
* **DNT outcome labels are never used during feature construction**: Pathway perturbations, gene rankings, and differential expression statistics reflect the cellular response to chemical exposure alone. DNT toxicity classifications are introduced only in downstream machine learning stages.
* **Compound-safe validation**: Chemical compounds evaluated during training must be held out in test splits to verify generalizability to prospective, unseen chemicals.

### Principle 2: Experimental Context & Matched Controls
* Treatment organoids are strictly compared against **matched vehicle controls** derived from the same study, experiment, developmental stage/age, and organoid protocol.
* Pooled controls across distinct experiments or disparate developmental ages are prohibited.
* Sample metadata (dose, duration, replicate identity, cell line) must be explicitly verified from primary sources. Automatic guessing or synthetic biological pooling is rejected.

### Principle 3: Gene Universe Conservation
* The pipeline evaluates the full canonical human genome from **Ensembl Release 112** (60,207 genes) rather than pre-filtering to arbitrary subsets.
* Statistically valid genes lacking external database keys (e.g., 7,122 valid genes without NCBI Entrez IDs) are **preserved in the molecular feature universe** rather than discarded as QC failures.

### Principle 4: Reproducibility & Dual-Stack Architecture
* Computational tasks leverage a hybrid architecture:
  * **Python 3.12**: Data parsing, schema validation, web UI, orchestration, feature matrices, machine learning.
  * **R 4.4.3 / Bioconductor 3.20**: Official gold-standard bioinformatics engines (`DESeq2 1.46.0`, `clusterProfiler 4.14.6`, `fgsea 1.32.4`).
* Headless subprocess bridges enforce strict runtime version locking and capture execution provenance.

---

## 3. Biological & Experimental Scope

* **Biological Systems**: Human neural organoids, cerebral cortex organoids (SOSR-CO), neural progenitor spheroids derived from WA09 (H9) or other approved human pluripotent stem cell lines.
* **Exposure Regimens**: Chemical exposures spanning acute (24–48 h) to sub-chronic windows during early organoid neurogenesis (e.g., Day 6, Day 35).
* **Primary Pilot Compound**: **Valproic Acid (VPA)**, CAS 99-66-1, CID:3121. Known human teratogen and benchmark DNT reference compound evaluated at 200 µM in public study **GSE232218**.

---

## 4. Intended Use and Regulatory Disclaimer

> [!IMPORTANT]
> This platform is a research software system designed for computational toxicology, biomarker discovery, and methodological development. It is **not** a clinical diagnostic system or regulatory safety determination engine. Its outputs must not be used for medical decisions or standalone chemical licensing without independent experimental verification.
