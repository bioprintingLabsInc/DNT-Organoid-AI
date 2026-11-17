# DNT AI — GEO search protocol v1.0

**Status:** Registered search plan; expanded search NOT yet executed.  
**Database:** NCBI Gene Expression Omnibus (GEO), using NCBI E-utilities (`db=gds`) and Series (`GSE`) filtering.  
**Research objective:** Identify candidate human neural/organoid chemical-exposure transcriptomics studies for sample-level eligibility review, data collection, and eventually leakage-safe DNT model development.

## Scientific distinction

1. **Original discovery (historical notebook):** Earlier 11 human organoid search phrases reportedly yielded 590 query matches and 516 deduplicated GEO Series accession IDs. These are **historical notebook-reported counts**, not outputs of the expanded search.
2. **Expanded discovery (this protocol):** 45 controlled search phrases across four groups, executed as three documented Boolean queries. Query counts, date and final unique results must be filled **only from an actual API run**.
3. **Eligibility:** A search hit is not a verified DNT dataset. Each study and relevant individual sample must undergo biological, exposure, control, replication, and data-availability verification.
4. **Training readiness:** Counts of independent chemicals, DNT-positive/negative evidence, studies, controls, and biological replicates must be evaluated separately from accession totals.

## Controlled keyword vocabulary

`config/dnt_geo_search_keywords_v1.json` is the authoritative machine-readable list of **45 terms**: 16 organoid/3D terms, 10 neural-cell terms, 12 exposure terms, and 7 transcriptomics terms. The original 11 terms are labeled separately within this file.

## Exact query construction

`scripts/search_geo_dnt.py` constructs and logs complete queries from this frozen JSON configuration.

- **Q1_BIO_EXPOSURE:** (all organoid/3D OR neural-cell terms) AND (all exposure terms).
- **Q2_BIO_TRANSCRIPTOMICS:** (all organoid/3D OR neural-cell terms) AND (all transcriptomics terms).
- **Q3_EXPOSURE_TRANSCRIPTOMICS:** (all exposure terms) AND (all transcriptomics terms).
- All three append `"Homo sapiens"[Organism]` AND `gse[Entry Type]`.
- Each phrase is quoted and queried in `[All Fields]`; phrase matching and GEO indexing may affect sensitivity and should be reviewed. Searches are intentionally overlapping, and their counts must not be summed as unique studies.
- Q3 deliberately prioritizes sensitivity; it may find non-neural studies requiring exclusion. Q2 may find non-exposure work. Never label every query hit eligible.

## How to run (requires network, not run by this commit)

From the repository root with Python 3:

```bash
export NCBI_EMAIL="researcher@example.org"  # replace with an actual contact
python scripts/search_geo_dnt.py --output-dir outputs/geo_search/2026-11-04_run01
```

On Windows PowerShell use `$env:NCBI_EMAIL="researcher@example.org"` instead of `export`. Choose a **new, unique** directory for each run; do not overwrite an existing run. The script uses only Python standard library modules. It requests paginated GEO search IDs and maps them to GSE accessions using ESummary; requests are rate-limited and retried. Do not place API keys in code or log files. Follow NCBI usage policies.

## Required evidence from each actual search

The run directory records:
- `run_manifest.json`: query IDs, exact full Boolean query strings, UTC start/end times, ESearch hit counts, unique accession count, code environment and SHA-256 hash of the keyword configuration.
- `Q*_accessions.csv`: the GSE IDs retrieved per query.
- `all_unique_geo_series.csv`: unique accessions and which query IDs matched each accession.

Additionally preserve Git commit SHA, original notebook source, search engine/service, accession-level retrieval outcomes, API failures, and any query modifications. Never use results from an incomplete run to claim a final count. Once executed, archive a checksum-verified copy of each output and record the run's Git SHA.

## Next stages, after search

1. Compare newly discovered GSE IDs with the original 516 by accession. Preserve original source identity and distinguish historical vs newly discovered IDs.
2. Reinspect all original 516 **independently**, regardless of historical classification.
3. At GEO Series and GEO Sample level verify species, neural model, organoid differentiation age, exposure timing, compound identity and concentration, controls, biological replicates, sequencing modality, available counts/FASTQ, and source links.
4. Record one source/evidence reference for each extracted value and log inclusion, exclusion, or review reason. Missing metadata is not automatically an exclusion.
5. Transfer eligible sample-level records into `data/metadata/dnt_sample_metadata_v1.xlsx` and matching CSV; do not confuse a series-level search catalog with the sample-level 22-field sheet.
6. Summarize deduplication, screening exclusions and final eligibility in a reproducible selection diagram (with counts supported by run logs).

## Reporting limits

The exact date, full executable queries and search engine responses for the **historical** 516 need separate reconstruction from the original notebook and/or its saved run logs. A present-day rerun cannot prove what a historical search returned on an unknown date. The expanded search's result counts must not be reported until executed and logged.

The attached method supports a *prospective* protocol; it does not establish diagnostic performance, prediction validity, biological comparability, or eligible training size.
