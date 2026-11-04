# Track B — First source-checked pilot (6 of 1,666 GSE records)

Status: **scientific characterization pilot**, not a finished 1,666-study census or completed sample-level metadata validation. Source check: 2026-11-04.

## Scope
- Existing Phase 3 accession inventory is retained unchanged as history.
- Characterization is based on original public GEO Series records and a linked peer-reviewed full text for GSE101577 and GSE169751. A GEO Series record is not proof of detailed supplementary-material review.
- One study can contain many chemical agents and experimental conditions; study-level sample counts do **not** certify independent biological replication.
- The primary collection is chemically exposed human brain organoids and related 3D neural models; secondary includes chemical-exposed other human neural cells.
- These classifications concern collection scope, **not** chemical-level DNT outcome labels.

## Initial findings
1. **GSE119861** — Cerebral organoids exposed to methamphetamine (5 µM, 7 days); 8 mock and 6 treated organoids analyzed using scRNA-seq. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119861)
2. **GSE169751** — Human differentiating cortical neural cultures exposed to methylmercury at 0.1 or 1 µM in two temporal paradigms. This is a **secondary neural model**, not evidence of a 3D brain organoid. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169751) · [Publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC8761390/)
3. **GSE150005** — LUHMES dopaminergic-neuron bulk RNA-seq, 7 toxicants, triplicate vehicle and treatment cultures; one GEO Series has **multiple chemical identities**. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150005)
4. **GSE166297** — Human neural-progenitor differentiation and neuron–astrocyte cultures, 5 chemicals at 2 study-defined dose levels. GEO reports N=1 and n=8; these symbols require experimental-design scrutiny before accepting independent biological replication. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166297)
5. **GSE147633** — THAP1/SNCA genetic perturbation in neuronal models; no administered compound verified in GEO; **exclude from chemical-exposure collection**. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147633)
6. **GSE101577** — Cerebral-organoid brain-tumor study with engineered genetic alterations; **not** a normal unperturbed developmental reference, despite previous Phase 3 labeling. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101577) · [Publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC6071863/)

## Quality-control rules for scaling to all 1,666
- Use actual GSE summary and study objective, GSM treatment fields, source publications and supplement where available. Do not copy old Phase 3 eligibility as factual study characterization.
- Track **study topic**, **evidence for administered chemical**, **model**, **assay**, **DNT relevance**, **replicate independence**, and **source tier** separately.
- Preserve overlapping themes, but supply one primary study category for reconciliation.
- For any published paper not checked, status must remain GEO_SERIES_CONFIRMED, not PUBLICATION_VERIFIED.
- No pipeline code changes and no DNT-positive/negative outcome labeling.
- Chemical-exposure candidates require sample and study-level evidence checks before inclusion in downstream analysis.

## Next Track B tasks
Extend the same source-first methodology across the remaining 1,660 accession records, prioritize suspicious prior inclusions and the 'REFERENCE_DATASET' category, then revisit uncertain/excluded candidates. Document scientific counts only when derived from confirmed per-study records.
