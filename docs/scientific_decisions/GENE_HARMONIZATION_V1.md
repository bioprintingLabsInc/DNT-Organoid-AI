# Gene Harmonization v1: Scientific Decision and Architecture

## Purpose and Scope

Gene Harmonization v1 provides a reproducible, provenance-preserving gene-identifier harmonization layer for human bulk RNA-seq expression data that has already passed Step 6A bulk raw-count quality control (QC). It standardizes heterogeneous gene identifiers (such as Ensembl Gene IDs, gene symbols, and Entrez Gene IDs) across different public datasets and future in-house studies into a common canonical human gene representation based on stable Ensembl Gene IDs, while preserving the approved HGNC gene symbol and original source identifier for every feature.

This layer operates strictly as a preprocessing checkpoint between QC and downstream normalization. It does not perform normalization, filter genes, alter expression counts, or assign toxicological classifications.

## Position in the Pipeline

```text
metadata contract
→ metadata validation
→ ingestion / format routing
→ bulk raw-count QC
→ gene harmonization (Step 6B / Gene Harmonization v1)
→ normalization
→ treatment-vs-control response construction
→ standardized features
→ model-ready dataset
→ AI training / calibration / validation
```

Gene harmonization precedes normalization because cross-study gene space alignment and gene-level feature attributes (such as GC content, gene length, and canonical biotypes) must be established before normalization algorithms and cross-study comparisons can operate on consistent gene definitions.

## Canonical Gene Reference

Gene Harmonization v1 implements a truly version-locked dual-source reference architecture:

1. **Canonical Gene Universe Authority (Ensembl Release 112):**
   - **Source Authority:** EMBL-EBI / Ensembl
   - **Release:** `112`
   - **Genome Assembly:** `GRCh38.p14`
   - **Release Date:** `2024-05`
   - **Canonical Source Asset:** `Homo_sapiens.GRCh38.112.gtf.gz`
   - **Canonical Source URL:** `https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz`
   - **Source SHA-256 Checksum:** `8c87436bab973c871447743c2f5cbb2b557be0523fe6efe3ff33307a4c836002`
   - **Total Canonical Genes:** 63,140

2. **Nomenclature and Cross-Reference Layer (HGNC Pinned Snapshot):**
   - **Source Authority:** HUGO Gene Nomenclature Committee (HGNC)
   - **Snapshot Date:** `2026-09-24`
   - **Live Source URL:** `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt`
   - **Uncompressed Snapshot SHA-256:** `69bb5722d5a42bb355580deb2c9f197ce3d7f7d13807173b9674a7db65e52191`
   - **Archived In-Repo Asset:** [`data/references/sources/hgnc_complete_set_20260924.tsv.gz`](file:///C:/Users/chari/.gemini/antigravity/scratch/DNT-Organoid-AI/data/references/sources/hgnc_complete_set_20260924.tsv.gz) (Gzip SHA-256: `86e67e6fa725d1db6b42442f45238ebf9fcd810b720aea4f65a23988020d1719`)
   - **Ensembl-Linked Records:** 42,396
   - **Reproducibility Guarantee:** Even if the external live HGNC URL changes or retires, the exact retrieved snapshot is preserved in repository storage for full offline recovery.

3. **Generated Production Reference Asset:**
   - **File Path:** [`data/references/ensembl_human_genes_v112.tsv`](file:///C:/Users/chari/.gemini/antigravity/scratch/DNT-Organoid-AI/data/references/ensembl_human_genes_v112.tsv)
   - **SHA-256 Checksum:** `3054aa734c15bb0a9cc64b67a8ead2a023e7dfd2cc0d0d6b7706372f4b13940e`
   - **Record Count:** 63,140
   - **Provenance Metadata:** [`data/references/ensembl_human_genes_v112_provenance.json`](file:///C:/Users/chari/.gemini/antigravity/scratch/DNT-Organoid-AI/data/references/ensembl_human_genes_v112_provenance.json)
   - **Automated Reference Builder:** [`scripts/build_gene_reference.py`](file:///C:/Users/chari/.gemini/antigravity/scratch/DNT-Organoid-AI/scripts/build_gene_reference.py) (v2.0.0)
   - **Hardcoding Policy:** Zero gene rows, mappings, synonyms, or expected mappings are hardcoded.

### Reference Catalog Coverage Statistics (Ensembl 112 Canonical Human Gene Universe)

> [!NOTE]
> The 63,140 canonical records define the **Ensembl Release 112 human gene universe** at the gene level (covering all protein-coding genes, long non-coding RNAs, pseudogenes, and small RNA loci annotated in GRCh38.p14). This gene-level catalog is not described as the complete human transcriptome, which comprises alternative splice isoforms and transcript-level variants.

- **Total Canonical Genes:** 63,140 (100.0%)
- **Genes with HGNC Cross-Reference:** 42,076 (66.64%)
- **Genes without HGNC Cross-Reference (Ensembl-Only):** 21,064 (33.36%)
- **Genes with Approved HGNC Symbols:** 42,076 (66.64%)
- **Genes with HGNC IDs:** 42,076 (66.64%)
- **Genes with Entrez IDs:** 41,578 (65.85%)
- **Genes with Synonyms / Previous Symbols:** 25,060 (39.69%)
- **Top Gene Biotypes Represented:**
  - `protein_coding`: 20,089 (31.81%)
  - `lncRNA`: 19,258 (30.50%)
  - `processed_pseudogene`: 10,144 (16.07%)
  - `unprocessed_pseudogene`: 2,602 (4.12%)
  - `misc_RNA`: 2,217 (3.51%)
  - `snRNA`: 1,910 (3.03%)
  - `miRNA`: 1,879 (2.98%)
  - `TEC` (to be experimentally confirmed): 1,052 (1.67%)
  - `transcribed_unprocessed_pseudogene`: 962 (1.52%)
  - `snoRNA`: 942 (1.49%)

### Preservation of Ensembl Genes Lacking HGNC Nomenclature

The reference catalog explicitly preserves all 21,064 Ensembl Release 112 genes that lack HGNC nomenclature cross-references (such as novel lncRNAs, antisense transcripts, and pseudogenes).
- When queried by their canonical Ensembl Gene ID (e.g. `ENSG00000228037`), these genes resolve to `GeneMappingStatus.UNIQUELY_MAPPED` with `canonical_gene_id = "ENSG00000228037"` and `approved_symbol = None`.
- They are **never treated as unmapped** simply because HGNC has not assigned a symbol or HGNC ID.

## Upstream Quality Control (QC) Gating

Gene harmonization is strictly gated by the upstream bulk raw-count QC stage (Step 6A) through an auditable disposition contract:

1. **`QCStatus.FAIL`:**
   - Downstream harmonization is **completely blocked**.
   - Produces finding `qc_status_failure` (`Severity.ERROR`).
   - Dataset status is `Status.FAIL`, and `matrix` is set to `None`.
2. **`QCStatus.NEEDS_REVIEW`:**
   - Requires an explicit, auditable **`QCReviewDisposition`** (specifying non-empty `decision` and scientific `reason`/`rationale`, with optional `reviewer`, `review_date`, `affected_qc_status`, `affected_asset_id`, and `provenance`).
   - **Rejection of Bare Booleans:** A bare production Boolean such as `allow_needs_review_qc=True` is **strictly rejected**. If passed without a valid `QCReviewDisposition`, processing fails immediately with `unauthorized_upstream_qc_needs_review` (`Severity.ERROR`), dataset status `Status.FAIL`, and `matrix = None`.
   - **Authorized Processing:** When a valid `QCReviewDisposition` is provided, harmonization proceeds with `upstream_qc_needs_review_authorized` (`Severity.REVIEW`), dataset status `Status.NEEDS_REVIEW`, and the full disposition object preserved on `HarmonizedDataset.qc_review_disposition`.
3. **`QCStatus.PASS_WITH_WARNINGS`:**
   - Harmonization proceeds with `upstream_qc_pass_with_warnings` (`Severity.WARNING`).
4. **`QCStatus.PASS`:**
   - Harmonization proceeds cleanly without QC warnings.

## Supported Identifier Types and Handling

The harmonization layer supports four principal identifier types:

1. **Ensembl Gene IDs (`ensembl_gene_id`):**
   - Plain stable IDs: `ENSGxxxxxxxxxxx`
   - Versioned IDs: `ENSGxxxxxxxxxxx.15`
   - Version suffixes are recognized and stripped solely for canonical lookup to ensure compatibility across Ensembl annotation releases.
   - The original identifier string (including any version suffix) is preserved verbatim in `HarmonizedGene.original_gene_id`.
2. **Gene Symbols (`gene_symbol`):**
   - Case-insensitive resolution against approved HGNC symbols and curated synonyms/previous symbols.
   - Symbols matching multiple canonical genes (shared ambiguous synonyms) are flagged as `AMBIGUOUS` and never arbitrarily mapped.
3. **Entrez / NCBI Gene IDs (`entrez_gene_id`):**
   - Discrete integer IDs resolved against NCBI-to-Ensembl cross-references.
4. **HGNC IDs (`hgnc_id`):**
   - Stable HGNC identifiers (e.g. `HGNC:11998`) resolved to canonical Ensembl IDs.

### Identifier Type Authority and Non-Guessing Policy

The canonical metadata field `gene_identifier_type` on `input_data_assets` is authoritative. When declared:
- The declared type is normalized against configured aliases.
- Bounded syntax checks verify whether matrix identifiers conform to the declared pattern.
- If an incompatibility is observed, an `identifier_type_mismatch` review finding is produced.

When `gene_identifier_type` is missing or unknown:
- Bounded inspection is applied to the matrix identifiers.
- If identifiers unambiguously match the Ensembl ID pattern (`^ENSG\d{11}(\.\d+)?$`), they are recognized without ambiguity.
- If identifiers are ambiguous, heterogeneous, or purely numeric without metadata context (which could represent row indices or Entrez IDs), the system **refuses to silently guess**. Instead, it generates an `ambiguous_identifier_type` finding and assigns an `AMBIGUOUS` mapping status to the affected features.

## Explicit Mapping Status Semantics

Every source gene row in the input matrix receives an explicit, deterministic mapping status:

| Mapping Status | Definition | Result Fields Populated |
|---|---|---|
| `UNIQUELY_MAPPED` | Resolved unambiguously to exactly one canonical Ensembl Gene ID. | `canonical_gene_id`, `approved_symbol` (when available), `candidate_canonical_ids = (canonical_id,)` |
| `AMBIGUOUS` | Identifier matches multiple canonical genes (e.g. paralogs, multi-locus repeats, shared synonyms) or identifier type cannot be determined without guessing. | `canonical_gene_id = None`, `approved_symbol = None`, `candidate_canonical_ids = (cand1, cand2, ...)` |
| `UNMAPPED` | Validly formed identifier of the declared/detected type that is not present in the active reference version. | `canonical_gene_id = None`, `approved_symbol = None`, `candidate_canonical_ids = ()` |
| `INVALID` | Blank, whitespace, or malformed non-gene token (e.g. header leftovers or comment strings). | `canonical_gene_id = None`, `approved_symbol = None`, `candidate_canonical_ids = ()` |

**Scientific Safeguard:** Ambiguous mappings are never silently converted into unique mappings. Unmapped genes are never silently dropped.

### Unmapped Genes and Mapping Coverage Statistics

Real-world bulk RNA-seq count matrices routinely contain non-coding features, novel transcripts, or deprecated loci that do not map to the current canonical reference catalog:
- Individual unmapped identifiers produce `unmapped_gene_identifier` findings at `Severity.INFO` by default. They do not trigger a fatal abort or arbitrarily block processing of the entire matrix.
- Summary metrics track mapping efficiency via `HarmonizationSummary`:
  - `percentage_uniquely_mapped`: Percentage of input rows uniquely mapped to canonical Ensembl IDs.
  - `percentage_unresolved`: Percentage of input rows that are unmapped, ambiguous, or invalid.
- An optional mapping coverage threshold can be enforced via `coverage_thresholds.minimum_mapping_coverage_percentage` in `config/gene_harmonization.yaml`. If configured and the dataset fails to meet the threshold, a `low_mapping_coverage` finding is generated (severity configurable; defaults to `REVIEW`).

## Non-Destructive Invariant and Collision Policy

Public transcriptomic datasets occasionally contain:
1. Duplicate source identifiers (e.g. multiple rows with identical gene symbols or identical Ensembl IDs).
2. Multiple distinct source identifiers that map to the same canonical Ensembl Gene ID (e.g. `ENSG00000141510.15` and `ENSG00000141510.16`, or an approved symbol and an obsolete alias).

### Harmonization Collision Rules

- **Zero Expression Loss:** All original matrix rows and columns remain intact. No counts are modified, zeroed, or removed.
- **No Automatic Aggregation:** The harmonization layer **never** sums, averages, selects, or merges count values for colliding or duplicate genes.
- **Explicit Collision Surfacing:** All collision groups are recorded in `HarmonizedDataset.collisions` with their canonical ID, original IDs, and source row indices.
- **Review Precedence:** Any collision or duplicate source ID generates a `Severity.REVIEW` finding (`canonical_gene_collision` or `duplicate_source_gene_identifiers`), setting dataset status to `NEEDS_REVIEW`.
- **Downstream Deferral:** Scientific decisions on how to resolve collisions (e.g. summing raw counts, selecting highest variance, or isoform-level quantification) are deferred to future approved downstream stages.

## Status Model

Harmonization status follows strict severity precedence:

```text
any ERROR present  → FAIL
else any REVIEW    → NEEDS_REVIEW
else any WARNING   → PASS_WITH_WARNINGS
else               → PASS
```

- `FAIL`: Severe failure (e.g. input asset failed raw-count QC, unauthorized upstream `NEEDS_REVIEW`, matrix is unreadable, empty identifiers, or unsupported declared identifier type).
- `NEEDS_REVIEW`: Matrix is readable and harmonized, but contains ambiguous genes, canonical collisions, duplicate source identifiers, low mapping coverage, or authorized upstream `NEEDS_REVIEW` requiring scientific curation.
- `PASS_WITH_WARNINGS`: Harmonized with non-critical warnings under active rules (e.g. input passed upstream QC with warnings).
- `PASS`: 100% of genes uniquely resolved to canonical Ensembl IDs with zero collisions, ambiguities, or unmapped features.

## Explicit Boundaries: What Gene Harmonization v1 Does and Does NOT Do

### What It Does
- Maps bulk RNA-seq gene identifiers to canonical human Ensembl Gene IDs.
- Retains approved HGNC gene symbols and original source identifiers.
- Preserves input matrix expression counts in an immutable `CountMatrix`.
- Surfaces unmapped, ambiguous, invalid, and colliding identifiers with explicit reasons.
- Records reference identity, release version, and SHA-256 provenance checksums.
- Computes dataset-level summary and coverage metrics.

### What It Does NOT Do
- Does **not** perform count normalization (e.g. CPM, TPM, DESeq2 size factors, or RUVSeq).
- Does **not** perform differential gene expression analysis.
- Does **not** aggregate technical or biological replicates.
- Does **not** construct treatment-versus-matched-control responses.
- Does **not** filter or discard low-expression or unmapped genes.
- Does **not** select features or discover DNT biomarker panels.
- Does **not** train AI models, calibrate predictions, or assign toxicological labels.

## Known Limitations and Unresolved Scientific Decisions

1. **Collision Resolution Policy:** An approved scientific policy for resolving canonical collisions before normalization (e.g. whether to sum raw counts for duplicated features, select highest-variance rows, or maintain an unresolved feature mask) is deferred to the downstream stage.
2. **Species Scope:** Current reference catalog is restricted to *Homo sapiens* (Taxonomy ID 9606, GRCh38). Non-human cross-species mappings are out of scope.
3. **Ambiguous Historical Synonyms:** Certain historical gene aliases (e.g. `FLIP` matching both `CFLAR` and `NDUFA13`) represent multi-gene mapping candidates. Gene Harmonization v1 flags these as `AMBIGUOUS` with candidate sets preserved. Resolving which gene was intended in historical datasets requires study-specific assay curation.
4. **Microarray Probe Sets & Transcript Isoforms:** Microarray probe identifiers (e.g. Affymetrix / Illumina probe IDs) and transcript-level quantifications (ENST IDs) are not natively mapped in v1 and require dedicated translation layers or transcript-level quantification tools.
