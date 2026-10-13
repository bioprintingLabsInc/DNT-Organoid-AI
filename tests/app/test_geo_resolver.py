"""Unit tests for GEO accession resolver and sample metadata mapper."""

import json
from pathlib import Path
import tempfile
import pytest

from src.app.geo_resolver import (
    extract_cell_line,
    extract_developmental_age,
    extract_exposure_duration,
    extract_replicate_id,
    fetch_geo_metadata,
    parse_gsm_brief_text,
    resolve_geo_samples,
)

SAMPLE_GEO_TEXT = """
^SAMPLE = GSM7321024
!Sample_title = brain organoids, Vehicle; rep1
!Sample_geo_accession = GSM7321024
!Sample_source_name_ch1 = brain organoids (day 6)
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: brain organoids (day 6)
!Sample_characteristics_ch1 = cell line: iPSC (AICS-0023)
!Sample_characteristics_ch1 = Sex: male
!Sample_characteristics_ch1 = treatment: vehicle control
!Sample_description = Sample_3901-AT-1

^SAMPLE = GSM7321026
!Sample_title = brain organoids, VPA200; rep1
!Sample_geo_accession = GSM7321026
!Sample_source_name_ch1 = brain organoids (day 6)
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: brain organoids (day 6)
!Sample_characteristics_ch1 = cell line: iPSC (AICS-0023)
!Sample_characteristics_ch1 = Sex: male
!Sample_characteristics_ch1 = treatment: VPA; 200 uM
!Sample_description = Sample_3901-AT-12

^SAMPLE = GSM7321027
!Sample_title = brain organoids, VPA400; rep1
!Sample_geo_accession = GSM7321027
!Sample_source_name_ch1 = brain organoids (day 6)
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: brain organoids (day 6)
!Sample_characteristics_ch1 = cell line: iPSC (AICS-0023)
!Sample_characteristics_ch1 = Sex: male
!Sample_characteristics_ch1 = treatment: VPA; 400 uM
!Sample_description = Sample_3901-AT-13

^SAMPLE = GSM7321028
!Sample_title = brain organoids, CHIR99021; rep1
!Sample_geo_accession = GSM7321028
!Sample_source_name_ch1 = brain organoids (day 6)
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: brain organoids (day 6)
!Sample_characteristics_ch1 = cell line: iPSC (AICS-0023)
!Sample_characteristics_ch1 = Sex: male
!Sample_characteristics_ch1 = treatment: CHIR99021; 10uM
!Sample_description = Sample_3901-AT-14
"""


def test_parse_gsm_brief_text() -> None:
    """Verify parsing of NCBI GEO brief format blocks."""
    parsed = parse_gsm_brief_text(SAMPLE_GEO_TEXT, "GSE232218")
    assert parsed["series_accession"] == "GSE232218"
    assert parsed["sample_count"] == 4

    gsm1 = parsed["samples"]["GSM7321024"]
    assert gsm1["title"] == "brain organoids, Vehicle; rep1"
    assert gsm1["description"] == "Sample_3901-AT-1"
    assert gsm1["source_name"] == "brain organoids (day 6)"
    assert "treatment: vehicle control" in gsm1["characteristics"]


def test_extract_replicate_id() -> None:
    """Verify regex-based biological replicate identification."""
    assert extract_replicate_id("brain organoids, Vehicle; rep1", "", "GSM1") == "rep1"
    assert extract_replicate_id("Sample rep_2", "", "GSM2") == "rep2"
    assert extract_replicate_id("", "Sample_replicate-3", "GSM3") == "rep3"
    assert extract_replicate_id("no replicate mentioned", "", "GSM4") == "rep_GSM4"


def test_extract_developmental_age() -> None:
    """Verify developmental stage extraction."""
    chars = ["tissue: brain organoids (day 6)", "cell line: iPSC"]
    assert extract_developmental_age(chars, "brain organoids (day 6)") == "Day 6"

    chars_wk = ["tissue: brain organoids (week 4)"]
    assert extract_developmental_age(chars_wk, "") == "Week 4"

    assert extract_developmental_age([], "") == ""


def test_extract_cell_line() -> None:
    """Verify cell line string extraction."""
    chars = ["tissue: brain organoids", "cell line: iPSC (AICS-0023)", "Sex: male"]
    assert extract_cell_line(chars) == "iPSC (AICS-0023)"
    assert extract_cell_line(["tissue: only"]) == ""


def test_resolve_geo_samples() -> None:
    """Verify mapping matrix sample headers to conditions, replicates, and exclusions."""
    geo_data = parse_gsm_brief_text(SAMPLE_GEO_TEXT, "GSE232218")
    matrix_sids = [
        "Sample_3901-AT-1",   # Vehicle
        "Sample_3901-AT-12",  # VPA 200
        "Sample_3901-AT-13",  # VPA 400 (other dose -> exclude)
        "Sample_3901-AT-14",  # CHIR99021 (other agent -> exclude)
    ]

    df, suggestions = resolve_geo_samples(
        matrix_sample_ids=matrix_sids,
        geo_data=geo_data,
        target_agent="VPA",
        target_dose=200.0,
        target_dose_unit="uM",
        vehicle_keyword="vehicle",
    )

    assert len(df) == 4
    cond_map = dict(zip(df["sample_id"], df["condition_type"]))
    rep_map = dict(zip(df["sample_id"], df["biological_replicate_id"]))

    assert cond_map["Sample_3901-AT-1"] == "control"
    assert rep_map["Sample_3901-AT-1"] == "rep1"

    assert cond_map["Sample_3901-AT-12"] == "treatment"
    assert rep_map["Sample_3901-AT-12"] == "rep1"

    assert cond_map["Sample_3901-AT-13"] == "exclude"
    assert cond_map["Sample_3901-AT-14"] == "exclude"

    assert suggestions["developmental_age"] == "Day 6"
    assert suggestions["biological_source"] == "iPSC (AICS-0023)"
    assert suggestions["source_accession"] == "GSE232218"
    assert suggestions["exposure_duration"] == 48.0
    assert suggestions["exposure_duration_unit"] == "h"
    assert suggestions["duration_status"] == "verified_publication"
    assert "PMC10340169" in suggestions["duration_citation"]


def test_extract_exposure_duration() -> None:
    """Verify explicit duration extraction from sample metadata strings."""
    assert extract_exposure_duration(["treatment: VPA", "duration: 48h"], "", "") == (48.0, "h")
    assert extract_exposure_duration([], "brain organoids, 24-hour treatment", "") == (24.0, "h")
    assert extract_exposure_duration(["exposure: 7 days"], "", "") == (7.0, "days")
    assert extract_exposure_duration(["treatment: vehicle control"], "sample 1", "") == (None, "")


def test_unresolved_duration_requires_review() -> None:
    """When duration is absent from GEO sample records and publication metadata, it is marked unresolved."""
    geo_data = {
        "series_accession": "GSE999999",
        "sample_count": 2,
        "samples": {
            "GSM1": {
                "gsm_accession": "GSM1",
                "title": "sample 1 ctrl",
                "description": "s1",
                "source_name": "organoid",
                "characteristics": ["treatment: vehicle control"],
            },
            "GSM2": {
                "gsm_accession": "GSM2",
                "title": "sample 2 trt",
                "description": "s2",
                "source_name": "organoid",
                "characteristics": ["treatment: compound_x 10 uM"],
            },
        },
    }
    _, suggestions = resolve_geo_samples(["s1", "s2"], geo_data, target_agent="compound_x", target_dose=10.0)
    assert suggestions["exposure_duration"] is None
    assert suggestions["duration_status"] == "unresolved_requires_review"


def test_fetch_geo_metadata_caching(tmp_path: Path) -> None:
    """Verify that cached metadata JSON avoids network requests."""
    fake_data = {
        "series_accession": "GSE999999",
        "sample_count": 1,
        "samples": {
            "GSM1": {
                "gsm_accession": "GSM1",
                "title": "fake title",
                "description": "s1",
                "source_name": "organoid",
                "characteristics": ["treatment: control"],
                "organism": "Homo sapiens",
                "instrument_model": "Illumina",
            }
        },
    }
    cache_file = tmp_path / "GSE999999_metadata.json"
    cache_file.write_text(json.dumps(fake_data), encoding="utf-8")

    # Fetch with cache_dir pointing to tmp_path
    loaded = fetch_geo_metadata("GSE999999", cache_dir=tmp_path)
    assert loaded["series_accession"] == "GSE999999"
    assert "GSM1" in loaded["samples"]


def test_fetch_geo_metadata_invalid_accession() -> None:
    """Invalid accession format must raise ValueError."""
    with pytest.raises(ValueError, match="Invalid GEO series accession"):
        fetch_geo_metadata("NOT_A_GEO_ID")
