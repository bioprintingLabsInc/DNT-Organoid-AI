"""NCBI GEO Series and GSM sample metadata adapter.

Parses GEO brief text format, series matrix text, or pre-fetched GEO JSON records,
preserving 100% of sample titles, descriptions, and characteristics for provenance tracking.
"""

from collections.abc import Mapping
import json
from pathlib import Path
from typing import Any

from .base import MappingFinding, MetadataSourceAdapter, RawSampleRecord, StudyContext


class GeoMetadataAdapter(MetadataSourceAdapter):
    """Adapter for NCBI GEO series and sample records."""

    def can_read(self, source: Any) -> bool:
        if isinstance(source, dict):
            return "series_accession" in source or "samples" in source
        if isinstance(source, (str, Path)):
            p = Path(source)
            name = p.name.lower()
            return name.endswith(("_metadata.json", ".soft", ".soft.gz", ".txt")) or "geo" in name
        return False

    def parse_study(
        self,
        source: Any,
        study_id: str | None = None,
        experiment_id: str = "exp_01",
        **kwargs: Any,
    ) -> tuple[StudyContext, list[RawSampleRecord], list[MappingFinding]]:
        """Parse GEO metadata dictionary or file into canonical intermediate records."""
        findings: list[MappingFinding] = []
        source_name = str(source) if isinstance(source, (str, Path)) else "geo_metadata_dict"

        data: dict[str, Any]
        if isinstance(source, dict):
            data = source
        else:
            p = Path(source)
            if not p.exists():
                findings.append(
                    MappingFinding(
                        severity="ERROR",
                        rule="file_not_found",
                        field_path=str(source),
                        message=f"GEO metadata file not found: {source}",
                    )
                )
                return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

            text = p.read_text(encoding="utf-8-sig", errors="replace")
            if p.name.lower().endswith(".json"):
                try:
                    data = json.loads(text)
                except Exception as e:
                    findings.append(
                        MappingFinding(
                            severity="ERROR",
                            rule="json_parse_error",
                            field_path=str(p),
                            message=f"Failed to parse JSON: {e}",
                        )
                    )
                    return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings
            else:
                # Text format: parse using brief text format parser
                from src.app.geo_resolver import parse_gsm_brief_text
                series_acc = p.stem.split("_")[0]
                data = parse_gsm_brief_text(text, series_acc)

        samples_map = data.get("samples", {})
        series_acc = data.get("series_accession", study_id or "GSE_UNKNOWN")

        raw_records: list[RawSampleRecord] = []
        for gsm_id, rec in samples_map.items():
            if not isinstance(rec, dict):
                continue

            attrs: dict[str, Any] = {
                "gsm_accession": gsm_id,
                "title": rec.get("title", ""),
                "description": rec.get("description", ""),
                "source_name": rec.get("source_name", ""),
                "organism": rec.get("organism", "Homo sapiens"),
                "instrument_model": rec.get("instrument_model", ""),
            }
            key_prov: dict[str, str] = {
                "gsm_accession": "!Sample_geo_accession",
                "title": "!Sample_title",
                "description": "!Sample_description",
                "source_name": "!Sample_source_name_ch1",
                "organism": "!Sample_organism_ch1",
                "instrument_model": "!Sample_instrument_model",
            }

            # Parse characteristics tag-value pairs
            chars = rec.get("characteristics", [])
            for c_idx, c_str in enumerate(chars):
                if ":" in c_str:
                    tag, val = c_str.split(":", 1)
                    clean_tag = tag.strip().lower().replace(" ", "_")
                    attrs[clean_tag] = val.strip()
                    key_prov[clean_tag] = f"!Sample_characteristics_ch1[{c_idx}]"
                else:
                    attrs[f"characteristic_{c_idx}"] = c_str.strip()
                    key_prov[f"characteristic_{c_idx}"] = f"!Sample_characteristics_ch1[{c_idx}]"

            raw_records.append(
                RawSampleRecord(
                    sample_id=gsm_id,
                    attributes=attrs,
                    source_file=source_name,
                    source_locator=f"GSM: {gsm_id}",
                    source_key_provenance=key_prov,
                )
            )

        pub_accessions = (
            ({"namespace": "GSE", "identifier": series_acc},)
            if series_acc and series_acc != "GSE_UNKNOWN"
            else ()
        )

        context = StudyContext(
            study_id=study_id or f"study_{series_acc}",
            study_title=data.get("title", f"GEO Series {series_acc}"),
            source_type="public",
            public_accessions=pub_accessions,
            experiment_id=experiment_id,
            attributes={"series_accession": series_acc, "sample_count": len(raw_records)},
            findings=tuple(findings),
        )

        return context, raw_records, findings
