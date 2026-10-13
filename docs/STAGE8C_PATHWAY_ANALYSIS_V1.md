# Stage 8C v1: Pathway and Signaling Feature Foundation Specification

## 1. Scientific Purpose
The scientific purpose of Stage 8C is to convert the unbiased genome-wide treatment-versus-matched-control molecular response from bulk RNA-seq (produced in Stage 7B) into interpretable, standardized pathway- and signaling-level feature representations.

Crucially, Stage 8C operates **without using any developmental neurotoxicity (DNT) outcome labels, reference toxicant classifications, or expected biological outcomes**. It strictly captures the cellular perturbation phenotype of the chemical exposure in human neural organoids.

> **Foundational Principle**: Stage 8C identifies treatment-associated pathway perturbations. It does not establish that a pathway is a DNT biomarker. DNT association is evaluated only after independent labels and leakage-safe model development in later stages.

---

## 2. Input Contract & Immutable Boundaries
Stage 8C consumes the authoritative output of Stage 7B:
- Primary Input: `DifferentialExpressionContrastResult` (and its contained `gene_results: tuple[GeneDifferentialExpressionResult, ...]`).
- Gene Universe Reference: Governed Ensembl Release 112 human gene catalog (`data/references/ensembl_human_genes_v112.tsv`).

Stage 8C does **not**:
- Modify any completed Stage 1–7 statistical logic;
- Alter DESeq2 calculations, dispersion estimates, or normalization size factors;
- Filter input genes to only FDR-significant DEGs ($p_{\text{adj}} < 0.05$);
- Make live internet queries to external gene mapping web services;
- Hard-code expected pathways or use expected biology as an acceptance gate.

---

## 3. Why the Full Gene Response is Preserved
Standard over-representation analysis (ORA, e.g. hypergeometric test) requires an arbitrary statistical cutoff (e.g. $p_{\text{adj}} < 0.05$ and $|\log_2\text{FC}| > 1$), discarding 99%+ of the measured genome. In contrast, **Preranked Gene Set Enrichment Analysis (GSEA)** evaluates the coordinated distribution of an entire gene set across the complete spectrum of differential expression.

Preserving the full evaluated canonical gene response ($N=60,207$ genes in the validated GSE232218 pilot):
1. Captures subtle, coordinated shifts in multi-subunit complexes and signaling cascades that fail individual single-gene significance thresholds;
2. Prevents threshold-dependent selection bias;
3. Accurately estimates pathway enrichment scores ($ES$) and normalized enrichment scores ($NES$) against the empirical genome-wide background.

---

## 4. Ranking Formula & Contract
The primary approved ranking metric is the **DESeq2 Wald statistic**:

\[
\text{Wald}_g = \frac{\log_2 \text{FC}_g}{\text{SE}(\log_2 \text{FC}_g)}
\]

### Rationale:
1. **Continuous & Signed**: Positive values represent chemical induction; negative values represent transcriptional repression.
2. **Variance-Stabilized**: Normalizes effect size ($\log_2\text{FC}$) by its asymptotic standard error ($\text{SE}$), naturally down-weighting noisy, low-count transcripts.
3. **Continuous Distribution**: Avoids the p-value truncation cliff ($p < 10^{-16}$) common in standard GLM tests.

### Deterministic Tie-Breaking:
If multiple genes share identical Wald statistics, ordering is resolved deterministically by ascending Ensembl canonical gene ID (`canonical_gene_id`).

### Duplicate Entrez Resolution (`MAX_ABS_WALD`):
When multiple canonical Ensembl genes map to the identical NCBI Entrez Gene ID:
- **Governed Rule**: The representative gene with the maximum absolute Wald statistic ($\max |\text{Wald}|$) is selected.
- **Scientific Justification**: Preranked GSEA tests pathway locus enrichment. Retaining the representative transcript/locus with the highest magnitude perturbation captures the strongest biological signal while preventing artificial duplication of gene weights in the Kolmogorov-Smirnov running sum. Excluded duplicates are explicitly logged in the audit trail.

---

## 5. Governed Gene-to-Pathway Mapping Strategy
Mapping proceeds strictly through a local, governed adapter:
```
canonical Ensembl ID (ENSG...)
       ↓
governed local Ensembl 112 catalog (data/references/ensembl_human_genes_v112.tsv)
       ↓
NCBI Entrez Gene ID
       ↓
Pathway Database (GO BP / KEGG)
```
- **Local Integrity**: Fully reproducible offline; zero live HTTP/REST queries.
- **Audit Accounting**: Identifies and quantifies mapped genes, unmapped genes, duplicate Entrez targets, and one-to-many mappings.
- **No Silent Drops**: Every input gene row is explicitly accounted for in `PathwayGeneMappingResult`.

---

## 6. Pathway Databases & Analysis Engines
- **Gene Ontology (GO) Biological Process (BP)**:
  - Source: Bioconductor `org.Hs.eg.db` (v3.20.0) / `GO.db` (v3.20.0).
  - Scope: All valid human Biological Process terms with size $15 \le N \le 500$.
- **KEGG Pathways (Homo sapiens, hsa)**:
  - Source: Governed KEGG pathway catalog (`data/references/kegg_hsa_pathways.rds` and `.json`), generated via `clusterProfiler::download_KEGG('hsa')`.
  - Scope: All valid human KEGG metabolic and signaling pathways with size $15 \le N \le 500$.
- **Statistical Engine**: Headless R bridge invoking `clusterProfiler` (v4.14.6) and `fgsea` (v1.32.4) with adaptive multilevel permutation estimation (`eps = 0`, `seed = 42`).
- **Complete Preservation**: Retains all evaluated pathways ($p_{\text{value}} \le 1.0$), with an explicit boolean flag `is_significant = (adjusted_p_value < 0.05)`.

---

## 7. Quality Control (QC) & Provenance
Every run generates a comprehensive `PathwayAnalysisQC` report answering:
- Total Stage 7 evaluated genes;
- Valid ranking statistic genes (count & %);
- Invalid statistical result genes (count & %);
- Mapped Entrez genes (count & %);
- Unmapped Ensembl genes (count & %);
- One-to-many and duplicate Entrez mappings;
- Genes entering GSEA;
- Pathways evaluated and significant at FDR < 0.05.

Complete provenance is stored in `PathwayAnalysisProvenance`:
```
Pathway Feature
  ↓
Leading-Edge Genes & Entrez IDs
  ↓
DESeq2 Wald Statistics & Log2FC
  ↓
Stage 7 Contrast Specification
  ↓
Experimental Metadata (Compound, Dose, Duration, Developmental Stage)
  ↓
Replicate Sample IDs & Raw Matrix
  ↓
GEO Study Series
```

---

## 8. Leakage Prevention & Audit
Stage 8C code is strictly verified against outcome leakage:
- No toxicant vs. non-toxicant labels are accepted, required, or stored.
- No machine learning classifiers or predictive models are implemented.
- Pathway evaluation is blind to chemical class or expected neurodevelopmental outcome.

---

## 9. Limitations & Future Relationship to Stages 8A/8B/8D
- **Single-Contrast Limitation**: A single contrast (e.g. GSE232218 VPA 200 µM) reflects the cellular response of one compound at one concentration in one organoid context. It cannot establish a generalizable DNT signature.
- **Relationship to Stage 8A**: Stage 8A will aggregate multiple contrasts across diverse chemicals and studies into a unified cross-study gene response matrix.
- **Relationship to Stage 8B**: Stage 8B will construct gene-gene co-perturbation networks across multiple chemical treatments.
- **Relationship to Stage 8D**: Stage 8D will map chemical structures and exposure concentrations directly to the pathway activation profiles computed in Stage 8C.
