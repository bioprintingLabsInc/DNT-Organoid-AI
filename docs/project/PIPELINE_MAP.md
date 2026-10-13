# Scientific Pipeline Architecture Map

This document defines the high-level architecture of the DNT Organoid AI platform.

---

## 1. High-Level Flowchart

```mermaid
flowchart TD
    subgraph Intake["1. Data Intake, QC & Standardization (IMPLEMENTED)"]
        A1["Ingestion & Format Detection\n[VALIDATED]"]
        A2["Metadata Validation\n[VALIDATED]"]
        A3["Quality Control & Preprocessing\n[VALIDATED]"]
        A4["Gene Harmonization\n[VALIDATED]"]
        A5["Normalization\n[VALIDATED]"]
        A1 --> A2 --> A3 --> A4 --> A5
    end

    subgraph Contrast["2. Treatment-Control Molecular Response (IMPLEMENTED)"]
        B1["Treatment-Control Response Construction\n[VALIDATED - Step 7A]"]
        B2["Differential Molecular Response\n[VALIDATED - Step 7B]"]
        B1 --> B2
    end

    subgraph Features["3. Molecular Feature Engineering"]
        C1["Pathway & Signaling Features\n[FROZEN v1 - Stage 8C]"]
        C2["Cross-Study Gene Response Features\n[PLANNED - Stage 8A]"]
        C3["Gene Coresponse & Networks\n[PLANNED - Stage 8B]"]
        C4["Compound-Gene-Pathway Relationships\n[PLANNED - Stage 8D]"]
    end

    subgraph Labels["4. Evidence & Evaluation"]
        D1["Independent DNT Label Curation\n[PLANNED - Stage 9]"]
        D2["Compound-Safe Dataset Splitting\n[PLANNED - Stage 10]"]
    end

    subgraph Modeling["5. Machine Learning & Calibration"]
        E1["DNT Signature Discovery\n[PLANNED - Stage 11]"]
        E2["Model Training, Calibration & Validation\n[PLANNED - Stage 12]"]
    end

    A5 --> B1
    B2 --> C1
    B2 --> C2
    B2 --> C3
    B2 --> C4

    C1 --> E1
    C2 --> E1
    C3 --> E1
    C4 --> E1

    D1 --> E2
    D2 --> E2
    E1 --> E2

    classDef frozen fill:#d4edda,stroke:#28a745,stroke-width:2px;
    classDef validated fill:#cce5ff,stroke:#004085,stroke-width:2px;
    classDef planned fill:#f8f9fa,stroke:#6c757d,stroke-width:2px,stroke-dasharray: 5 5;

    class C1 frozen;
    class A1,A2,A3,A4,A5,B1,B2 validated;
    class C2,C3,C4,D1,D2,E1,E2 planned;
```

---

## 2. Stage Breakdown & Data Transition Contracts

### Stage 1: Data Intake, QC & Standardization (Implemented & Validated)
* **Goal**: Convert heterogeneous public and in-house RNA-seq files into canonically identified, technically qualified, normalized expression matrices.
* **Internal Layers**:
  1. *Ingestion & Format Detection*: Identifies format (`FASTQ`, `raw_counts`, `UMI_counts`, `H5AD`) and assigns modality route (`BULK_RAW_COUNT_PIPELINE`, etc.).
  2. *Metadata Validation*: Validates study, experiment, assay, sample, and exposure records against the canonical schema.
  3. *Quality Control & Preprocessing*: Checks matrix integrity (integer counts, library size, zero fraction, replicate correlation).
  4. *Gene Harmonization*: Maps raw gene identifiers to canonical Ensembl Release 112 IDs while retaining unmapped genes.
  5. *Normalization*: Estimates official DESeq2 size factors using the median-of-ratios method within metadata-defined cohorts.
* **Transition Contract**: Produces a `NormalizedDataset` containing size factors and normalized diagnostic counts.

---

### Stage 2: Molecular Response Foundations (Implemented & Validated)
* **Goal**: Construct rigorous treatment-versus-matched-control contrasts and compute genome-wide differential expression statistics.
* **Internal Layers**:
  1. *Treatment-Control Response* (`src/response_builder`): Pairs treated organoid conditions with matched vehicle controls sharing the same study, experiment, developmental age, and organoid protocol ($\ge 2$ true biological replicates per arm).
  2. *Differential Molecular Response* (`src/differential_expression`): Executes negative binomial generalized linear model fitting (DESeq2 Wald test) across all 60,207 canonical genes with gene-universe conservation auditing.
* **Transition Contract**: Produces a `DifferentialExpressionContrastResult` containing log2 fold changes, standard errors, Wald statistics, p-values, and BH-adjusted FDR for all evaluated genes.

---

### Stage 3: Molecular Feature Engineering (Partially Implemented & Frozen)
* **Goal**: Extract biologically meaningful, standardized feature representations from the unbiased molecular response.
* **Components**:
  1. *Pathway & Signaling Features* (**FROZEN v1 - tag: `stage8c-v1-frozen`**): Preranked GSEA using clusterProfiler for GO BP and KEGG based on finite Wald statistics, with strict preservation of the 32,140 statistically eligible gene universe.
  2. *Cross-Study Gene Response Features* (*Planned*): Standardizes baseline gene-level sensitivity and effect sizes across distinct organoid lines and exposure regimens.
  3. *Gene Coresponse & Networks* (*Planned*): Models co-regulated gene modules and functional network topology.
  4. *Compound-Gene-Pathway Relationships* (*Planned*): Integrates known chemical targets and bioactivity data.

---

### Stage 4: Independent Labels & Splitting (Planned)
* **Goal**: Curate ground-truth reference evidence and enforce leakage-safe train/test splits.
* **Components**:
  1. *DNT Label Curation*: Curates high-confidence reference positive and negative toxicant evidence independently of gene expression data.
  2. *Compound-Safe Splitting*: Partitions datasets so that entire chemical classes and studies are held out from training sets.

---

### Stage 5: Machine Learning & Calibration (Planned)
* **Goal**: Train predictive models, discover sparse DNT signatures, and calibrate prospective risk scores.
* **Components**:
  1. *DNT Signature Discovery*: Identifies minimal, robust gene/pathway signatures indicative of organoid developmental impairment.
  2. *Model Training & Calibration*: Fits supervised classifiers, evaluates cross-study generalizability, and calibrates probabilistic confidence.
