# Track B — Independent evidence audit checkpoint (2026-11-04)

## Coverage and limitations
- **All 1,666 existing GSE accessions** were independently checked for basic inventory uniqueness, title-derived 3D organoid, chemical-exposure, and genetic/disease signals. Output: `data/metadata/track_b_census/all_1666_accession_title_flags.csv`.
- Title-level screening is **not** publication verification; all 1,666 census rows correctly retain `source_verification_completed=0`. Existing Phase 3 INCLUDE/EXCLUDE designations are stored only as historical labels.
- **12 individual GEO records** were checked against original public GEO accession descriptions in two pilot batches. See `data/metadata/track_b_scientific_pilot_6.csv` and `data/metadata/track_b_direct_geo_checks_batch02.csv`.
- **Title flags across the 1,666 accession census:** 369 titles contain one of the literal 3D keywords (organoid/spheroid/assembloid/brainsphere), 117 contain one of the prespecified chemical/exposure title keywords, and 303 contain one of the genetic/disease keyword patterns. These are **overlapping indicators**, not eligible-study counts, and absence of words from a title does not prove absence of chemical exposure.
- Full GEO/GSM/publication/supplementary verification of the remaining records has **not** been completed. A scientifically valid final inclusion count cannot be derived from title flags.

## Direct evidence correcting the original inventory
1. GSE147633 is THAP1/SNCA genetic research, not a proven protein-aggregate chemical exposure: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147633
2. GSE101577 concerns genetically engineered cerebral-organoid brain tumors, not ordinary unexposed developmental maturation: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101577
3. GSE231638 and GSE231639 study APPV717I/AD genetic perturbations, linked under GSE231640 SuperSeries; the three accessions are not three independent chemical experiments: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231640
4. GSE289218 is a cerebral-organoid irradiation experiment (0 versus 24 Gy), not a chemical exposure: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289218
5. GSE186814 is a genuine VPA-exposed forebrain-organoid scRNA-seq series, but groups of three organoids were pooled into sequencing samples; pooled organoids do not establish independent biological sequencing replicates: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186814
6. GSE282273 carries an organoid/DNT title, but the deposited overall design describes chemical exposure in undifferentiated hiPSCs from day 2 and collection on day 4; original publication and model context need examination before assigning a mature brain-organoid model: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282273

## Evidence tiers
- `TITLE_ONLY_NOT_SOURCE_VERIFIED`: full 1,666-row first-pass accession census; no scientific inclusion decision.
- `GEO_SERIES_CONFIRMED`: source-reviewed pilot examples; verify GSM-level design and papers where missing.
- `PUBLICATION_AND_SUPPLEMENT_CONFIRMED`: only after exact original document sections and sample relationships are examined.

## Why previous model results must not be assumed correct
- A GSE can be a SuperSeries combining subseries: not an independent experiment by default.
- Biological agent, mutation, radiation, and chemical treatments are different intervention classes.
- Each study must preserve distinct condition-level doses/durations while counting chemicals by canonical identity.
- A GSM count, pooled organoid count, and biological replicate count are distinct quantities.

## Status
**Complete:** independent title-based census of 1,666 records, direct GEO checks on 12 records, stored on a separate GitHub branch.

**Not complete:** source-level scientific verification of all 1,666 records. This document must not be cited as proof of full independent validation or AI-training readiness.
