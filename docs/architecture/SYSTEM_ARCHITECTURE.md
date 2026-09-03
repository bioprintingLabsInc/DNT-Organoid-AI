# System Architecture

## Purpose and scope

DNT Organoid AI is intended to support reproducible development and validation of AI models that assess developmental-neurotoxicity-associated molecular responses in human brain organoids. The platform must accept public datasets and future in-house experiments without assuming that different studies, RNA-seq modalities, or source formats are directly comparable.

This document defines high-level responsibilities and scientific safeguards only. It does not prescribe QC thresholds, normalization methods, biomarkers, DNT classification rules, or model algorithms.

## Intended workflow

> Data sources → ingestion → format detection → metadata validation → format-specific QC/preprocessing → gene harmonization → normalization → treatment-versus-matched-control response construction → standardized features → model-ready dataset → AI training → calibration/validation → locked model → prospective prediction

Each transition produces a versioned, traceable artifact or an explicit validation record. Prospective prediction must use preprocessing compatible with the preprocessing used to develop the locked model.

## Input routes

The routing layer must identify source format and RNA-seq modality before preprocessing. Supported routes include:

- bulk RNA-seq;
- single-cell RNA sequencing (scRNA-seq);
- single-nucleus RNA sequencing (snRNA-seq);
- raw FASTQ or SRA inputs; and
- processed expression matrices when their provenance, metadata, processing history, and scientific usability are sufficient.

Raw read inputs and processed matrices enter through different format-specific routes. Bulk, single-cell, and single-nucleus data retain their modality identity and receive modality-appropriate QC and preprocessing. Modalities must not be blindly combined or normalized together. Any later integration requires an explicit, scientifically reviewed decision and a traceable transformation.

## Architectural principles

### Data integrity and provenance

- Raw source data are immutable.
- Original source provenance is always preserved, including public accession or in-house source identity where applicable.
- Every derived artifact is traceable to its source artifacts, processing version, configuration, and relevant software versions.
- Reprocessing creates a new derived artifact; it does not replace source data or erase prior provenance.
- QC failures, exclusions, and their reasons are explicit and auditable rather than silently removed.

### Experimental context

- Treatment samples remain linked to scientifically appropriate matched controls.
- Required experimental context includes compound, concentration and unit, exposure duration, organoid type, developmental age or stage, cell line or source, treatment/control status, and biological replicate identity.
- Developmental-stage information remains available so early-versus-later DNT prediction can be evaluated scientifically.
- Missing, ambiguous, or conflicting metadata is surfaced for resolution rather than inferred silently.

### Scientific governance

- Scientific thresholds are configurable, versioned, and scientifically approved rather than invented in code.
- RNA-seq modalities retain distinct processing routes unless a documented and approved integration step is applied.
- DNT reference labels are independent evidence and are not inferred automatically from gene-expression changes.
- Positive, negative/reference, and unresolved/uncertain reference states are representable without forcing uncertain evidence into a binary label.
- Training, validation, and test splitting prevents study-level and compound-level leakage.
- Model scores are not described automatically as absolute DNT probabilities unless calibration and external validation support that interpretation.
- A locked model includes the preprocessing contract required for compatible prospective prediction.

## Repository responsibilities

### Top-level directories

- `config/`: Version-controlled configuration inputs and references to scientifically approved settings. Configuration records choices; it does not embed unreviewed scientific assumptions.
- `data/`: Data lifecycle storage. `manifests/` identifies source assets and relationships; `metadata/` holds metadata templates or controlled records; `raw/` holds immutable source data; `interim/` holds traceable intermediate artifacts; `processed/` holds completed modality-specific processing outputs; and `model_ready/` holds versioned datasets approved for modeling.
- `src/`: Pipeline and platform implementation, organized by responsibility. Scientific logic is not yet implemented.
- `tests/`: Unit and integration checks corresponding to pipeline responsibilities, with emphasis on contracts, traceability, and prevention of silent data loss.
- `scripts/`: Minimal operational entry points for reproducible workflows and maintenance tasks.
- `docs/`: Architecture, scientific decisions, data dictionaries, and validation documentation.
- `outputs/`: Generated QC, processing, and model reports. These outputs are reproducible artifacts rather than source data.

### Major `src/` modules

- `ingestion/`: Registers incoming public or in-house assets without changing raw source content and captures source identifiers and checksums.
- `metadata/`: Parses and validates required experimental metadata while preserving missingness, ambiguity, and source values.
- `routing/`: Detects supported formats and modalities, then selects the appropriate processing route without treating modalities as interchangeable.
- `qc/`: Records modality-appropriate QC results, failures, exclusions, reasons, and approved configuration versions.
- `normalization/`: Applies a scientifically approved, versioned normalization route appropriate to the modality and intended comparison.
- `gene_harmonization/`: Maps gene identifiers to a documented reference while retaining original identifiers and mapping provenance.
- `response_builder/`: Constructs treatment-versus-scientifically-matched-control response records while preserving treatment, control, and replicate links.
- `features/`: Produces standardized, versioned feature representations from approved response records without defining DNT evidence labels.
- `labels/`: Manages independent DNT reference evidence and supports positive, negative/reference, and unresolved/uncertain states.
- `provenance/`: Connects every derived artifact to source artifacts, configurations, processing steps, versions, and decisions.
- `validation/`: Enforces data contracts and evaluates dataset, preprocessing, calibration, and model validation requirements, including leakage controls.
- `models/`: Contains future training, evaluation, calibration, locking, and prediction interfaces; no model implementation is currently defined.
- `reporting/`: Produces auditable summaries of QC, processing, exclusions, datasets, validation, and model results.

## Traceability and artifact lifecycle

Raw assets enter `data/raw/` and remain unchanged. Manifests and metadata records associate each asset with its source, modality, study, experiment, sample, replicate, treatment, and matched-control context. Derived artifacts move through interim, processed, and model-ready stages only with provenance links and recorded validation outcomes.

Exclusion does not mean deletion: excluded samples and failed artifacts remain represented in audit records with reasons and the applicable configuration version. A model-ready dataset must identify its inputs, feature definition, independent reference-label version, split assignments, and processing lineage.

## Model development and prospective use

Dataset splitting must isolate studies and compounds as required to prevent leakage across training, validation, and test partitions. Calibration and validation results determine how scores may be interpreted; absent adequate evidence, scores remain model outputs rather than absolute DNT probabilities.

A locked model is paired with its feature contract, gene reference, modality scope, preprocessing configuration, and validation record. Future in-house samples may be used for prospective prediction only when their processing is compatible with that contract. Incompatible inputs must be rejected or routed through a separately validated workflow rather than transformed silently.
