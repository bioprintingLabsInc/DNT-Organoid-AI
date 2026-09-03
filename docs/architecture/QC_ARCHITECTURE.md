# Quality-Control Architecture

## Purpose and scope

The QC architecture places modality-specific technical assessment between Step 5 ingestion/routing and later scientific processing. Each route is evaluated by a QC layer that understands the structure and limitations of that data type while sharing project-wide rules for metadata authority, provenance, reproducibility, and non-destructive review.

This document describes the architectural boundary. Step 6A is implemented and accepted. Step 6B and Step 6C are planned boundaries only and are not implemented by this documentation change.

## Routing to modality-specific QC

```text
Data sources
    │
    ▼
Step 4: canonical metadata validation
    │
    ▼
Step 5: asset registration, bounded inspection,
        format detection, consistency checking, routing
    │
    ├── BULK_RAW_COUNT_PIPELINE
    │       └── Step 6A: bulk RNA-seq raw-count QC
    │
    ├── SINGLE_CELL_COUNT_PIPELINE
    │       └── Step 6B: scRNA-seq/snRNA-seq count QC (planned)
    │
    ├── SINGLE_CELL_OBJECT_PIPELINE
    │       └── Step 6B: scRNA-seq/snRNA-seq object QC (planned)
    │
    ├── PROCESSED_EXPRESSION_PIPELINE
    │       └── Step 6C: processed-expression QC (planned)
    │
    └── MANUAL_REVIEW
            └── no automatic QC route until ambiguity is resolved
```

Routing establishes the appropriate processing destination; it does not perform QC or silently reinterpret biological metadata. A route-specific QC validator must confirm that the asset, assay, and canonical metadata remain compatible with that route.

## Why QC is modality-specific

### Bulk raw counts: Step 6A

A bulk raw-count column represents an assay-level sample profile across genes. QC therefore evaluates matrix integrity, per-sample library characteristics, and similarity among independent biological replicates within a condition. Raw counts remain unnormalized, and biological replication occurs at the sample level.

### Single-cell and single-nucleus data: Step 6B

Single-cell and single-nucleus inputs contain many cell or nucleus observations nested within biological samples. They require barcode integrity, cell-to-sample mapping, sparse/object-format handling, cell-level metrics, and explicit protection against treating cells as independent biological replicates. Object containers such as H5AD also require safe inspection of count-layer semantics. These needs cannot be represented correctly by the bulk matrix QC model.

Step 6B is not implemented yet.

### Processed expression: Step 6C

Processed matrices may contain normalized, transformed, scaled, summarized, or otherwise derived values. Raw-count invariants such as non-negative integer-like values may not apply, and the preprocessing history becomes part of scientific usability. Processed data therefore require separate provenance and compatibility checks instead of being forced through raw-count QC.

Step 6C is not implemented yet.

Combining these routes into one generic validator would risk applying scientifically invalid checks, normalization assumptions, or replicate definitions across incompatible representations.

## Shared QC principles

### Canonical metadata authority

Study, experiment, condition, sample, assay, source asset, modality, replicate, treatment/control, agent, concentration, exposure duration, and developmental-stage identities come from validated canonical metadata. QC does not infer them from filenames, sample labels, barcode prefixes, or expression patterns.

### Immutable raw data

QC reads source assets without rewriting them. It does not repair counts, overwrite metadata, rename identifiers, aggregate observations, or replace source files. Derived QC results are separate artifacts linked back to their inputs.

### Deterministic results

The same inputs, canonical metadata, configuration, and implementation version must produce the same measurements, finding order, status, and configuration checksum. Determinism makes review artifacts comparable and supports reproducible audit logs.

### Version-controlled rules

Structural invariants and scientific thresholds live in explicit versioned configurations. Scientifically contingent thresholds remain `null` until approved; Python code must not invent them.

### Explicit provenance

Every result links to the source asset and assay and preserves the applicable study, experiment, condition, sample, and replicate identifiers. It also records the ruleset version, semantic configuration checksum, calculation implementation/version, and method identifiers needed to reproduce the assessment.

### No silent repair or exclusion

QC produces metrics and findings. It does not silently delete, merge, impute, repair, aggregate, or exclude data. If a future approved rule supports exclusion, the preferred product is an auditable disposition or mask while the original data remain unchanged.

### Ambiguity requires review

Ambiguous format, biological mapping, replicate handling, or processing semantics produce a review finding and `NEEDS_REVIEW`. The system does not select the most convenient interpretation.

### Biological replicate integrity

Biological replication is defined by canonical experimental metadata. Technical replicates, files, columns, cells, and nuclei do not automatically increase biological *n*. Modality-specific QC must preserve this distinction and prevent pseudoreplication.

### One framework for public and in-house data

Public and future in-house datasets use the same metadata, routing, QC result, provenance, and audit principles. Source-specific information may differ, but no source receives an undocumented scientific shortcut.

### Scientific thresholds are separate from implementation

Code calculates measurements and applies only configured rules. Threshold selection, sensitivity testing, and validation are scientific governance activities. Separating them from calculation logic permits review and versioning without changing the meaning of the underlying metric implementation.

## Shared finding and status model

QC findings use four severities:

- `ERROR`: deterministic failure that prevents valid QC interpretation;
- `REVIEW`: unresolved ambiguity or decision requiring review;
- `WARNING`: an active, approved warning rule was triggered;
- `INFO`: descriptive observation that does not alter disposition.

Status precedence is consistent across modality-specific QC where practical:

```text
ERROR present   → FAIL
else REVIEW     → NEEDS_REVIEW
else WARNING    → PASS_WITH_WARNINGS
else            → PASS
```

`INFO` findings do not change status. A QC status is scoped to data quality under a particular ruleset and is never a DNT classification.

## Step 6A component flow

The accepted bulk raw-count implementation separates responsibilities:

```text
config/qc_rules.yaml
        │
        ▼
src/qc/config.py ────── rules validation and semantic checksum
        │
        ▼
src/qc/validator.py ─── orchestration
        ├── matrix.py ───── parsing and structural checks
        ├── metadata.py ─── canonical sample/assay/asset reconciliation
        ├── metrics.py ──── deterministic sample metrics
        ├── rules.py ────── opt-in configured numerical warnings
        ├── replicates.py ─ condition-level Spearman summaries
        └── result.py ───── immutable results, findings, ordering, status
```

The validator returns dataset-, sample-, and condition-level records. Structural errors prevent metric construction from invalid values. Valid matrices remain complete even when genes, samples, or replicate comparisons produce review findings.

## Downstream boundary

Modality-specific QC ends with measurements, findings, status, and provenance. It does not perform normalization, transformation, harmonization, treatment-control response construction, or modeling.

Downstream stages must consume an explicit reviewed QC disposition and remain compatible with the route, modality, metadata, and preprocessing history. A `PASS` does not authorize a later scientific method that has not itself been approved.

## Extension requirements for Step 6B and Step 6C

Future QC modules should:

- reuse stable shared concepts such as severity, deterministic ordering, status precedence, and canonical configuration checksums;
- retain separate modality-specific input, metric, and result models where their scientific meanings differ;
- avoid changing the locked Step 6A behavior;
- add route-compatibility tests from Step 5;
- preserve canonical metadata and raw assets;
- document every newly activated threshold and scientific decision;
- treat ambiguous input semantics as review rather than inference.

This architecture intentionally does not define Step 6B or Step 6C algorithms, thresholds, or scientific policies.
