# Stage 8C v1: Pathway and Signaling Feature Foundation Specification (Hardened)

## 1. Scientific Purpose
The scientific purpose of Stage 8C is to convert the unbiased genome-wide treatment-versus-matched-control molecular response from bulk RNA-seq (produced in Stage 7B) into interpretable, standardized pathway- and signaling-level feature representations.

Crucially, Stage 8C operates **without using any developmental neurotoxicity (DNT) outcome labels, reference toxicant classifications, or expected biological outcomes**. It strictly captures the cellular perturbation phenotype of the chemical exposure in human neural organoids.

> **Foundational Principle**: Stage 8C identifies treatment-associated pathway perturbations. It does not establish that a pathway is a DNT biomarker. DNT association is evaluated only after independent labels and leakage-safe model development in later stages.

---

## 2. Input Contract & Immutable Boundaries
Stage 8C consumes the authoritative output of Stage 7B:
- **Primary Input**: `DifferentialExpressionContrastResult` (and its contained `gene_results: tuple[GeneDifferentialExpressionResult, ...]`).
- **Gene Universe Reference**: Governed Ensembl Release 112 human gene catalog (`data/references/ensembl_human_genes_v112.tsv`).

Stage 8C does **not**:
- Modify any completed Stage 1–7 statistical logic;
- Alter DESeq2 calculations, dispersion estimates, or normalization size factors;
- Filter input genes to only FDR-significant DEGs ($p_{\text{adj}} < 0.05$);
- Make live internet queries to external gene mapping web services;
- Hard-code expected pathways or use expected biology as an acceptance gate;
- Discard statistically valid Stage 7 genes merely because they lack an Entrez ID.

---

## 3. Preserved Molecular Universe vs. Pathway Database Eligibility

### Architectural Imperative:
\[
\boxed{\text{No statistically valid gene is discarded just because it lacks an Entrez ID}}
\]

Stage 8C establishes a strict architectural distinction:
1. **`STATISTICALLY_ELIGIBLE_MOLECULAR_UNIVERSE`**:
   The complete set of genes from Stage 7 that pass governed statistical quality standards:
   - Non-zero count distribution (`baseMean > 0`);
   - Mathematically finite Wald statistic ($\text{Wald} \in \mathbb{R}$, non-NaN, non-Inf);
   - Governed Stage 7 status `result_status == "OK"`;
   - Passed Cook's distance outlier filtering.
   For the GSE232218 VPA 200 µM pilot, this contains exactly **32,140 genes**.

2. **`PATHWAY_DATABASE_ANNOTATION_ELIGIBILITY`**:
   Membership in external pathway knowledge bases (such as GO or KEGG) requires standard database identifiers (e.g., NCBI Entrez Gene ID).
   - An Entrez ID is an annotation property, **not** a measure of whether an RNA-seq result is statistically or scientifically valid.
   - Genes lacking Entrez IDs are **not** QC failures; they are valid cellular perturbation responses preserved in `retained_non_entrez_gene_features`.

```
                    Stage 7 Canonical Evaluated Genes (60,207)
                                      │
                     ┌────────────────┴────────────────┐
                     ▼                                 ▼
             All-Zero Counts (27,980)      Genes with Valid Wald (32,227)
                                                       │
                                      ┌────────────────┴────────────────┐
                                      ▼                                 ▼
                             Cook's Outliers (87)    STATISTICALLY ELIGIBLE UNIVERSE (32,140)
                                                                        │
                                      ┌─────────────────────────────────┴─────────────────┐
                                      ▼                                                   ▼
                       Entrez-Mapped Eligible (25,018)                   Retained Non-Entrez Genes (7,122)
                                      │                                                   │
                       ┌──────────────┴──────────────┐                           (lncRNAs, pseudogenes, TEC;
                       ▼                             ▼                            preserved for Stage 8A/8B/8D)
             Deduplicated (25,016)         Duplicates Dropped (2)
                       │
         ┌─────────────┴─────────────┐
         ▼                           ▼
     GO BP GSEA                  KEGG GSEA
```

---

## 4. Reconciled GSE232218 Gene Accounting

Every gene in the input is accounted for deterministically without silent drops:

| Accounting Stage | Gene Count | Percentage of Evaluated | Scientific Category / Reason |
| :--- | :--- | :--- | :--- |
| **Stage 7 Evaluated Canonical Genes** | **60,207** | 100.00% | Full canonical gene universe from Stage 2 |
| Less: All-zero count genes | -27,980 | 46.47% | `baseMean == 0.0`; no read support across samples |
| **Genes with Valid Wald Statistic** | **32,227** | 53.53% | Non-zero counts; Wald statistic calculated |
| Less: Cook's distance outliers | -87 | 0.14% | Flagged by DESeq2 as extreme sample-specific outliers |
| **Statistically Eligible Molecular Universe** | **32,140** | 53.38% | Reliable, QC-qualified Stage 7 molecular response |
| Less: Unmapped to NCBI Entrez ID in Ensembl 112 | -7,122 | 11.83% | Valid genes lacking Entrez ID; preserved in retained table |
| **Entrez-Mapped Eligible Genes** | **25,018** | 41.55% | Eligible to enter Entrez-based pathway analysis |
| Less: Duplicate Entrez entries resolved | -2 | 0.003% | 4 genes map to 2 Entrez IDs; resolved via `MAX_ABS_WALD` |
| **Final Genes Entering GSEA (GO & KEGG)** | **25,016** | 41.55% | Unique Entrez-ranked input for preranked GSEA |

---

## 5. Audit of the 7,122 Entrez-Unmapped Genes

The 7,122 statistically eligible genes lacking Entrez IDs in Ensembl Release 112 were empirically classified using governed local metadata:

### Biotype Breakdown:
- **lncRNA**: 4,933 (69.26%)
- **processed_pseudogene**: 821 (11.53%)
- **TEC (To be Experimentally Confirmed)**: 560 (7.86%)
- **protein_coding**: 390 (5.48%)
- **unprocessed_pseudogene**: 144 (2.02%)
- **transcribed_unprocessed_pseudogene**: 99 (1.39%)
- **transcribed_processed_pseudogene**: 91 (1.28%)
- **misc_RNA**: 43 (0.60%)
- **artifact**: 16 (0.22%)
- **transcribed_unitary_pseudogene**: 11 (0.15%)
- **snoRNA / rRNA / snRNA / pseudogene**: 14 (0.20%)

### Identifier & Annotation Availability:
- **Has Canonical Ensembl ID**: 7,122 (100.0%)
- **Has HGNC ID**: 241 (3.38%)
- **Has Reference Gene Symbol**: 241 (3.38%)
- **Recognized in `org.Hs.eg.db` SQLite keys**: 879 (12.34%)
- **Has GO BP Annotation in `org.Hs.eg.db`**: 31 (0.44%)
- **Eligible for KEGG Pathway Membership**: 0 (0.00%)

### Fate and Governance:
These 7,122 genes represent real transcriptional perturbations (primarily non-coding RNAs and pseudogenes). They are **not** classified as failed QC. Instead, they are exported to `data/processed/GSE232218_VPA_200uM_retained_non_entrez_genes.csv` with `future_feature_eligibility = True` for cross-study gene feature aggregation in Stages 8A, 8B, and 8D.

---

## 6. Governed Identifier Decisions for GO and KEGG

### GO Identifier Evaluation (Ensembl vs. Entrez):
We empirically tested direct Ensembl identification (`keyType = "ENSEMBL"` in `clusterProfiler::gseGO`) against the governed Ensembl $\to$ Entrez route:

1. **Gene Coverage**:
   - `org.Hs.eg.db` is an Entrez-centric database (`eg` = Entrez Gene). Its internal Ensembl xref table contains only 40,839 Ensembl IDs.
   - Out of the 32,140 statistically eligible genes, direct Ensembl recognizes only **22,111 genes (68.80%)**, missing 10,029 genes.
   - In contrast, our governed Ensembl 112 reference maps **25,018 genes (77.84%)** to Entrez IDs. Direct Ensembl loses **2,907 valid genes**.
   - Among the 7,122 unmapped genes, direct Ensembl recovers only **31 genes** with GO BP annotations.
2. **Mapping Ambiguity & Ranking Ties**:
   - Direct Ensembl generates 1:many mappings (`'select()' returned 1:many mapping between keys and columns`).
   - fgsea warns of extensive preranked statistic ties (**19.76% of list**), resulting in arbitrary gene ordering.
   - Ensembl $\to$ Entrez resolves duplicates deterministically via `MAX_ABS_WALD` with zero arbitrary tie-breaking.
3. **Concordance**:
   - For shared pathways ($N=4,890$), the Pearson correlation of Normalized Enrichment Scores (NES) is **$r = 0.9947$**.

**Governed Decision**: Direct Ensembl does not provide superior coverage (it reduces coverage by 2,907 genes) and introduces ambiguity. Therefore, the **Ensembl $\to$ Entrez $\to$ GO** route is retained as the governed Stage 8C standard, while all 32,140 genes remain fully preserved in the molecular universe.

### KEGG Identifier Decision:
KEGG is strictly Entrez-based. We do not manufacture synthetic mappings. The KEGG route legitimately consumes the 25,016 Entrez-compatible genes.

---

## 7. Explicit GSEA Parameters & Reproducibility Lock

To eliminate reliance on implicit defaults, all GSEA parameters are explicitly governed:
- `exponent / gseaParam`: `1.0` (standard weighted Kolmogorov-Smirnov statistic)
- `scoreType`: `"std"`
- `minGSSize`: `15`
- `maxGSSize`: `500`
- `eps`: `0.0` (exact permutation p-value calculation via fgsea multilevel)
- `pvalueCutoff`: `1.0` (preserves the full evaluated pathway universe)
- `pAdjustMethod`: `"BH"` (Benjamini-Hochberg FDR)
- `seed`: `42`
- `nPermSimple`: `1000`

### Governed Reference Provenance & Hashes:
- **KEGG Offline Snapshot**:
  - `data/references/kegg_hsa_pathways.rds`: SHA-256 `A54298B77CE24D0059FB58CE0B3E22EB6837589D02FA4B2437E7545A02967C20`
  - `data/references/kegg_hsa_pathways.json`: SHA-256 `B454F803A2935333FEA52076BAF5E7D98C4FDEFE29EFA8EB75903BD5BBBA3DB8`
  - Date: `2026-10-13` (372 pathways, 39,671 associations)
- **GO Annotation Databases**:
  - `org.Hs.eg.db`: v3.20.0 (extdata SQLite SHA-256: `7b88d1cf970e8acef5da60b0535bbc1b53768381940d2b4fd538524f9b612d18`)
  - `GO.db`: v3.20.0 (extdata SQLite SHA-256: `a20386d939b63e94e50cfcb5c2e154ed17d7996611482b54d96807258cacc7ef`)
  - R v4.4.3 / Bioconductor 3.20 locked in `scripts/r/environment-r.yaml` and `scripts/r/renv.lock`.

---

## 8. Scientific Wording Standards

Pathway enrichment metrics reflect empirical statistical distributions of transcript changes, not definitive causal mechanisms. User-facing outputs and reports adhere strictly to standard scientific terminology:
- **Approved**: *"positively enriched in VPA vs vehicle"*
- **Approved**: *"negatively enriched in VPA vs vehicle"*
- **Approved**: *"statistically significant positive enrichment (FDR < 0.05)"*
- **Approved**: *"statistically significant negative enrichment (FDR < 0.05)"*
- **Approved**: *"directionally concordant with published findings"*
- **Prohibited**: Causal claims that chemical exposure "activates" or "shuts down" an entire cellular pathway without biochemical verification.
- **Hypotheses**: Any mechanistic speculation is segregated into explicitly labeled sections:
  `POST-HOC BIOLOGICAL INTERPRETATION`.
