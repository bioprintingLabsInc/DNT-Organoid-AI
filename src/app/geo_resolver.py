"""NCBI GEO Accession Metadata Resolver.

Resolves public series (GSE) and sample (GSM) metadata to eliminate manual,
error-prone matching of cryptic sample identifiers (e.g. Sample_3901-AT-X).

Provides:
1. `fetch_geo_metadata`: Fetches and parses GEO Series records (GSM brief format) with local disk caching.
2. `resolve_geo_samples`: Maps count matrix columns to official GSM descriptions/titles,
   proposing explicit condition assignments ('control', 'treatment', or 'exclude') and biological replicate IDs.
3. `suggest_experiment_parameters`: Extracts inferred biological source, developmental stage,
   and vehicle information from verified GEO characteristics.
"""

from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any
import urllib.error
import urllib.request

import pandas as pd


@dataclass(frozen=True, slots=True)
class GeoSampleRecord:
    """Structured metadata record for a single GEO Sample (GSM)."""

    gsm_accession: str
    title: str
    description: str
    source_name: str
    characteristics: list[str]
    organism: str = "Homo sapiens"
    instrument_model: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def get_default_geo_cache_dir() -> Path:
    """Return cache directory for downloaded GEO metadata."""
    cache_dir = Path(tempfile.gettempdir()) / "dnt_geo_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def parse_gsm_brief_text(text: str, series_accession: str) -> dict[str, Any]:
    """Parse text returned by NCBI GEO acc.cgi with targ=gsm&form=text&view=brief."""
    samples: dict[str, dict[str, Any]] = {}
    blocks = text.split("^SAMPLE = ")

    for block in blocks[1:]:
        lines = block.splitlines()
        if not lines:
            continue
        gsm_id = lines[0].strip()

        title = ""
        desc = ""
        source_name = ""
        organism = "Homo sapiens"
        instrument = ""
        characteristics: list[str] = []

        for line in lines[1:]:
            line_str = line.strip()
            if line_str.startswith("!Sample_title = "):
                title = line_str.split(" = ", 1)[1].strip()
            elif line_str.startswith("!Sample_description = "):
                desc = line_str.split(" = ", 1)[1].strip()
            elif line_str.startswith("!Sample_source_name_ch1 = "):
                source_name = line_str.split(" = ", 1)[1].strip()
            elif line_str.startswith("!Sample_organism_ch1 = "):
                organism = line_str.split(" = ", 1)[1].strip()
            elif line_str.startswith("!Sample_instrument_model = "):
                instrument = line_str.split(" = ", 1)[1].strip()
            elif line_str.startswith("!Sample_characteristics_ch1 = "):
                char_val = line_str.split(" = ", 1)[1].strip()
                characteristics.append(char_val)

        samples[gsm_id] = {
            "gsm_accession": gsm_id,
            "title": title,
            "description": desc,
            "source_name": source_name,
            "characteristics": characteristics,
            "organism": organism,
            "instrument_model": instrument,
        }

    return {
        "series_accession": series_accession,
        "sample_count": len(samples),
        "samples": samples,
    }


def fetch_geo_metadata(
    accession: str,
    cache_dir: Path | None = None,
    timeout: int = 15,
) -> dict[str, Any]:
    """Fetch and parse GEO series sample metadata from NCBI, caching locally on disk.

    Args:
        accession: GEO series accession (e.g. 'GSE232218').
        cache_dir: Directory to cache metadata JSON files.
        timeout: Network timeout in seconds.

    Returns:
        Dictionary containing series accession and parsed sample records.
    """
    clean_acc = accession.strip().upper()
    if not clean_acc.startswith("GSE"):
        raise ValueError(f"Invalid GEO series accession '{accession}'. Must start with 'GSE' (e.g., 'GSE232218').")

    if cache_dir is None:
        cache_dir = get_default_geo_cache_dir()

    cache_file = cache_dir / f"{clean_acc}_metadata.json"
    if cache_file.exists():
        try:
            with cache_file.open("r", encoding="utf-8") as f:
                data = json.load(f)
                if isinstance(data, dict) and data.get("samples"):
                    return data
        except Exception:
            pass  # Re-fetch if cache is corrupt

    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={clean_acc}&targ=gsm&form=text&view=brief"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "DNT-Organoid-AI/1.0 (Public Research Metadata Resolver)"},
    )

    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            text = response.read().decode("utf-8", errors="replace")
    except urllib.error.URLError as e:
        # Check if local fallback cache exists
        if cache_file.exists():
            with cache_file.open("r", encoding="utf-8") as f:
                return json.load(f)
        raise RuntimeError(f"Failed to connect to NCBI GEO for accession '{clean_acc}': {e}") from e

    parsed = parse_gsm_brief_text(text, clean_acc)
    if not parsed["samples"]:
        raise ValueError(f"No sample records found for GEO accession '{clean_acc}'. Please check the accession ID.")

    try:
        with cache_file.open("w", encoding="utf-8") as f:
            json.dump(parsed, f, indent=2)
    except Exception:
        pass  # Non-fatal if cache write fails

    return parsed


def extract_replicate_id(title: str, description: str, fallback_id: str) -> str:
    """Extract a biological replicate identifier from sample title or description."""
    # Pattern: rep1, rep_1, rep.1, replicate 1, replicate_1
    m = re.search(r"rep(?:licate)?[_.\s-]?(\d+)", title, re.IGNORECASE)
    if m:
        return f"rep{m.group(1)}"

    m_desc = re.search(r"rep(?:licate)?[_.\s-]?(\d+)", description, re.IGNORECASE)
    if m_desc:
        return f"rep{m_desc.group(1)}"

    return f"rep_{fallback_id}"


def extract_developmental_age(characteristics: list[str], source_name: str) -> str:
    """Extract developmental stage string from characteristics or source name (e.g. 'Day 6')."""
    text_to_search = " ".join([source_name] + characteristics)
    m = re.search(r"day\s*(\d+)", text_to_search, re.IGNORECASE)
    if m:
        return f"Day {m.group(1)}"
    m_wk = re.search(r"week\s*(\d+)", text_to_search, re.IGNORECASE)
    if m_wk:
        return f"Week {m_wk.group(1)}"
    return ""


def extract_cell_line(characteristics: list[str]) -> str:
    """Extract cell line from characteristics (e.g. 'iPSC (AICS-0023)')."""
    for c in characteristics:
        if c.lower().startswith("cell line:"):
            return c.split(":", 1)[1].strip()
    return ""


def extract_exposure_duration(
    characteristics: list[str],
    title: str,
    description: str,
) -> tuple[float | None, str]:
    """Extract explicit exposure duration from sample characteristics, title, or description.

    Returns:
        (duration_float, duration_unit) or (None, "") if not explicitly documented.
    """
    text_to_search = " ".join([title, description] + characteristics)

    # Pattern 1: e.g. "treated for 48 h", "duration: 48h", "exposure: 24 hours"
    m = re.search(
        r"(?:duration|exposure|treated for|incubated for)[:\s]+(\d+(?:\.\d+)?)\s*(h|hr|hrs|hours?|days?|d|weeks?|wk)\b",
        text_to_search,
        re.IGNORECASE,
    )
    if m:
        val = float(m.group(1))
        unit_raw = m.group(2).lower()
        unit = "h" if unit_raw.startswith("h") else ("days" if unit_raw.startswith("d") else "weeks")
        return val, unit

    # Pattern 2: e.g. "48h treatment", "24-hour exposure"
    m2 = re.search(r"\b(\d+(?:\.\d+)?)\s*(?:-| )(?:hour|hr)s?\b", text_to_search, re.IGNORECASE)
    if m2:
        return float(m2.group(1)), "h"

    m3 = re.search(r"\b(\d+(?:\.\d+)?)\s*(?:-| )(?:day)s?\b", text_to_search, re.IGNORECASE)
    if m3:
        return float(m3.group(1)), "days"

    return None, ""


# Verified peer-reviewed publication evidence for public series where GEO sample records omit duration.
# Grounded strictly in official literature citations (DOI / PMCID / Section evidence).
VERIFIED_PUBLIC_METADATA: dict[str, dict[str, Any]] = {
    "GSE232218": {
        "publication": "Cells 2023, 12(13), 1697",
        "doi": "10.3390/cells12131697",
        "pmcid": "PMC10340169",
        "title": "A Shared Pathogenic Mechanism for Valproic Acid and SHROOM3 Knockout in a Brain Organoid Model of Neural Tube Defects",
        "section": "Section 2.9 (RNA Sequencing Experiment)",
        "exposure_duration": 48.0,
        "exposure_duration_unit": "h",
        "evidence_quote": (
            "Day 6 brain organoids were generated from a male iPSC line (AICS-0023) and treated for 48 h "
            "with either a vehicle or one of the following 3 compounds: ... (2) the teratogenic anti-seizure "
            "medication VPA at a concentration of either 200 µM or 400 µM ... All the treatments were performed in 4 replicates."
        ),
        "organoid_context": "Cerebral Cortex Organoid (SOSR-CO)",
        "organoid_type": "cerebral_organoid",
        "brain_region": "cerebral_cortex",
        "biological_source": "iPSC (AICS-0023)",
        "species": "Homo sapiens",
        "developmental_age": "Day 6",
        "vehicle": "vehicle control",
    }
}


def resolve_geo_samples(
    matrix_sample_ids: list[str],
    geo_data: dict[str, Any],
    target_agent: str,
    target_dose: float | None = None,
    target_dose_unit: str = "uM",
    vehicle_keyword: str = "vehicle",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Map matrix sample IDs to official GEO GSM records and propose condition/replicate assignments.

    Args:
        matrix_sample_ids: Sample column headers from uploaded count matrix.
        geo_data: Parsed GEO series metadata dictionary.
        target_agent: Chemical / stressor name (e.g. 'VPA' or 'Valproic Acid').
        target_dose: Target concentration/dose (e.g. 200).
        target_dose_unit: Concentration unit (e.g. 'uM').
        vehicle_keyword: Vehicle string (e.g. 'vehicle' or 'control').

    Returns:
        (assignments_df, metadata_suggestions)
        where assignments_df has columns:
        ['sample_id', 'condition_type', 'biological_replicate_id', 'gsm_accession', 'geo_title', 'geo_treatment']
    """
    gsm_records: dict[str, dict[str, Any]] = geo_data.get("samples", {})
    if not gsm_records:
        raise ValueError("GEO metadata contains no sample records.")

    # Build lookup maps by description, title, and GSM accession
    desc_to_gsm: dict[str, dict[str, Any]] = {}
    title_to_gsm: dict[str, dict[str, Any]] = {}
    acc_to_gsm: dict[str, dict[str, Any]] = {}

    for gsm_id, rec in gsm_records.items():
        acc_to_gsm[gsm_id.lower()] = rec
        desc = rec.get("description", "").strip()
        if desc:
            desc_to_gsm[desc.lower()] = rec
        title = rec.get("title", "").strip()
        if title:
            title_to_gsm[title.lower()] = rec

    # Terms for matching
    clean_agent = target_agent.strip().lower()
    clean_veh = vehicle_keyword.strip().lower() or "vehicle"
    dose_str = str(int(target_dose)) if (target_dose is not None and target_dose.is_integer()) else (str(target_dose) if target_dose is not None else "")

    rows: list[dict[str, Any]] = []
    matched_ages: list[str] = []
    matched_cell_lines: list[str] = []
    matched_durations: list[tuple[float, str]] = []

    for sid in matrix_sample_ids:
        sid_clean = str(sid).strip()
        sid_lower = sid_clean.lower()

        # Find matching GSM record
        matched_rec = (
            desc_to_gsm.get(sid_lower)
            or acc_to_gsm.get(sid_lower)
            or title_to_gsm.get(sid_lower)
        )

        if not matched_rec:
            # Partial search if exact match fails
            for d_key, rec in desc_to_gsm.items():
                if sid_lower == d_key or sid_lower in d_key or d_key in sid_lower:
                    matched_rec = rec
                    break

        if not matched_rec:
            rows.append(
                {
                    "sample_id": sid_clean,
                    "condition_type": "",
                    "biological_replicate_id": "",
                    "gsm_accession": "",
                    "geo_title": "",
                    "geo_treatment": "",
                }
            )
            continue

        # Extract treatment characteristic
        chars = matched_rec.get("characteristics", [])
        trt_str = ""
        for c in chars:
            c_low = c.lower()
            if c_low.startswith("treatment:") or c_low.startswith("agent:") or c_low.startswith("compound:"):
                trt_str = c.split(":", 1)[1].strip()
                break
        if not trt_str:
            trt_str = matched_rec.get("title", "")

        trt_lower = trt_str.lower()
        title_lower = matched_rec.get("title", "").lower()
        all_trt_context = f"{trt_lower} {title_lower}"

        # Determine condition assignment
        proposed_condition = ""
        is_control = any(k in all_trt_context for k in [clean_veh, "control", "ctrl", "untreated", "vehicle control"])

        # Check agent match
        is_agent_match = False
        if clean_agent:
            # Match agent name or acronym (e.g. VPA for Valproic Acid)
            if clean_agent in all_trt_context:
                is_agent_match = True
            elif clean_agent == "valproic acid" and "vpa" in all_trt_context:
                is_agent_match = True
            elif clean_agent == "vpa" and "valproic" in all_trt_context:
                is_agent_match = True

        if is_control:
            proposed_condition = "control"
        elif is_agent_match:
            # Check dose if specified
            if dose_str:
                # If dose matches this sample
                if dose_str in all_trt_context:
                    proposed_condition = "treatment"
                else:
                    # Different dose of same agent -> exclude from this specific contrast
                    proposed_condition = "exclude"
            else:
                proposed_condition = "treatment"
        else:
            # Different chemical/condition -> exclude
            proposed_condition = "exclude"

        rep_id = extract_replicate_id(
            matched_rec.get("title", ""),
            matched_rec.get("description", ""),
            matched_rec["gsm_accession"],
        )

        age = extract_developmental_age(chars, matched_rec.get("source_name", ""))
        if age and age not in matched_ages:
            matched_ages.append(age)

        cline = extract_cell_line(chars)
        if cline and cline not in matched_cell_lines:
            matched_cell_lines.append(cline)

        dur_val, dur_unit = extract_exposure_duration(
            chars,
            matched_rec.get("title", ""),
            matched_rec.get("description", ""),
        )
        if dur_val is not None and (dur_val, dur_unit) not in matched_durations:
            matched_durations.append((dur_val, dur_unit))

        rows.append(
            {
                "sample_id": sid_clean,
                "condition_type": proposed_condition,
                "biological_replicate_id": rep_id,
                "gsm_accession": matched_rec["gsm_accession"],
                "geo_title": matched_rec.get("title", ""),
                "geo_treatment": trt_str,
            }
        )

    assignments_df = pd.DataFrame(rows)

    series_acc = geo_data.get("series_accession", "").strip().upper()
    verified_entry = VERIFIED_PUBLIC_METADATA.get(series_acc)

    # Resolve exposure duration strictly from explicit source evidence
    if verified_entry and "exposure_duration" in verified_entry:
        exp_duration = float(verified_entry["exposure_duration"])
        exp_dur_unit = str(verified_entry.get("exposure_duration_unit", "h"))
        duration_status = "verified_publication"
        duration_citation = (
            f"{verified_entry['publication']} (PMCID: {verified_entry['pmcid']}, {verified_entry['section']}): "
            f"\"{verified_entry['evidence_quote']}\""
        )
    elif matched_durations:
        exp_duration = float(matched_durations[0][0])
        exp_dur_unit = str(matched_durations[0][1])
        duration_status = "verified_geo_sample"
        duration_citation = "Extracted directly from GEO sample characteristics."
    else:
        exp_duration = None
        exp_dur_unit = "h"
        duration_status = "unresolved_requires_review"
        duration_citation = "Unresolved: duration not declared in GEO sample records or verified publication metadata."

    dev_age = matched_ages[0] if matched_ages else (verified_entry.get("developmental_age", "") if verified_entry else "")
    bio_src = matched_cell_lines[0] if matched_cell_lines else (verified_entry.get("biological_source", "") if verified_entry else "")
    org_ctx = verified_entry.get("organoid_context", "Cerebral Cortex Organoid") if verified_entry else "Cerebral Cortex Organoid"
    org_type = verified_entry.get("organoid_type", "cerebral_organoid") if verified_entry else "cerebral_organoid"
    brain_reg = verified_entry.get("brain_region", "cerebral_cortex") if verified_entry else "cerebral_cortex"
    veh_val = verified_entry.get("vehicle", "vehicle control") if verified_entry else "vehicle control"

    suggestions = {
        "developmental_age": dev_age,
        "biological_source": bio_src,
        "species": "Homo sapiens",
        "organoid_context": org_ctx,
        "organoid_type": org_type,
        "brain_region": brain_reg,
        "vehicle": veh_val,
        "source_accession": geo_data.get("series_accession", ""),
        "exposure_duration": exp_duration,
        "exposure_duration_unit": exp_dur_unit,
        "duration_status": duration_status,
        "duration_citation": duration_citation,
        "verified_evidence": verified_entry,
    }

    return assignments_df, suggestions
