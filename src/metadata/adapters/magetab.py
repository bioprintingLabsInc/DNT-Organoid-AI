"""MAGE-TAB metadata adapter for ArrayExpress / BioStudies datasets (IDF and SDRF).

Supports:
1. Investigation Description Format (.idf.txt): study title, experimental factors, protocols, accessions.
2. Sample and Data Relationship Format (.sdrf.txt):
   - Sample identifiers (Source Name)
   - Biological characteristics (Characteristics[...])
   - Experimental factors (Factor Value[...]) and associated units (Unit[...])
   - Library strategy, layout, and sequencing platform
   - Associated raw data files and run accessions (Comment[ENA_RUN], Comment[FASTQ_URI])
   - Grouping multiple paired-read rows belonging to the same biological sample.
"""

from collections import defaultdict
import io
from pathlib import Path
import re
from typing import Any

from .base import MappingFinding, MetadataSourceAdapter, RawSampleRecord, StudyContext


class MageTabAdapter(MetadataSourceAdapter):
    """Adapter for ArrayExpress and BioStudies MAGE-TAB formats."""

    def can_read(self, source: Any) -> bool:
        if isinstance(source, (str, Path)):
            p = Path(source)
            name = p.name.lower()
            return name.endswith((".sdrf.txt", ".sdrf", ".idf.txt", ".idf"))
        return False

    def parse_idf(self, idf_path: Path) -> dict[str, Any]:
        """Parse an Investigation Description Format (.idf.txt) file."""
        info: dict[str, Any] = {}
        if not idf_path.exists():
            return info

        text = idf_path.read_text(encoding="utf-8-sig", errors="replace")
        for line in text.splitlines():
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            parts = line.split("\t")
            if not parts:
                continue
            key = parts[0].strip()
            vals = [p.strip() for p in parts[1:] if p.strip()]
            val_str = vals[0] if len(vals) == 1 else (vals if len(vals) > 1 else "")
            info[key] = val_str

        return info

    def parse_study(
        self,
        source: Any,
        idf_source: Any | None = None,
        study_id: str | None = None,
        experiment_id: str = "exp_01",
        **kwargs: Any,
    ) -> tuple[StudyContext, list[RawSampleRecord], list[MappingFinding]]:
        """Parse MAGE-TAB SDRF and optional IDF file into canonical intermediate records."""
        findings: list[MappingFinding] = []
        sdrf_path = Path(source)

        if not sdrf_path.exists():
            findings.append(
                MappingFinding(
                    severity="ERROR",
                    rule="file_not_found",
                    field_path=str(source),
                    message=f"SDRF file not found: {source}",
                )
            )
            return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

        # Check for companion IDF file if not supplied explicitly
        idf_info: dict[str, Any] = {}
        if idf_source:
            idf_info = self.parse_idf(Path(idf_source))
        else:
            # Look for companion .idf.txt next to .sdrf.txt
            candidate_idf = sdrf_path.with_name(sdrf_path.name.replace(".sdrf.txt", ".idf.txt").replace(".sdrf", ".idf.txt"))
            if candidate_idf.exists():
                idf_info = self.parse_idf(candidate_idf)

        # Parse SDRF tabular contents
        text = sdrf_path.read_text(encoding="utf-8-sig", errors="replace")
        lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
        if not lines:
            findings.append(
                MappingFinding(
                    severity="ERROR",
                    rule="empty_sdrf",
                    field_path=str(sdrf_path),
                    message="SDRF file is empty.",
                )
            )
            return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

        headers = lines[0].split("\t")
        col_indices = {h.strip(): idx for idx, h in enumerate(headers)}

        # Find Source Name column
        source_name_idx = col_indices.get("Source Name")
        if source_name_idx is None:
            # Try case-insensitive lookup
            for h, idx in col_indices.items():
                if h.lower() == "source name":
                    source_name_idx = idx
                    break

        if source_name_idx is None:
            findings.append(
                MappingFinding(
                    severity="ERROR",
                    rule="missing_source_name",
                    field_path="sdrf.header",
                    message="MAGE-TAB SDRF requires a 'Source Name' column.",
                )
            )
            return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

        # Group rows by biological Source Name (handling paired-end multi-row entries)
        sample_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row_idx, line in enumerate(lines[1:], start=1):
            parts = line.split("\t")
            if len(parts) <= source_name_idx:
                continue
            s_name = parts[source_name_idx].strip()
            if not s_name:
                continue

            row_dict: dict[str, str] = {}
            for col_name, idx in col_indices.items():
                val = parts[idx].strip() if idx < len(parts) else ""
                row_dict[col_name] = val
            sample_rows[s_name].append(row_dict)

        raw_records: list[RawSampleRecord] = []
        for s_name, rows in sample_rows.items():
            primary_row = rows[0]
            consolidated_attrs: dict[str, Any] = {}
            key_prov: dict[str, str] = {}

            # Collect FASTQ URIs and file names across paired rows
            all_files: list[str] = []
            all_runs: list[str] = []

            for r in rows:
                for k, v in r.items():
                    if not v:
                        continue
                    if "submitted_file_name" in k.lower() or "scan name" in k.lower() or "fastq_uri" in k.lower():
                        if v not in all_files:
                            all_files.append(v)
                    elif "ena_run" in k.lower():
                        if v not in all_runs:
                            all_runs.append(v)

            # Consolidate standard attributes from primary row
            for k, v in primary_row.items():
                consolidated_attrs[k] = v
                key_prov[k] = k

            if all_files:
                consolidated_attrs["associated_files"] = all_files
            if all_runs:
                consolidated_attrs["ena_runs"] = all_runs

            raw_records.append(
                RawSampleRecord(
                    sample_id=s_name,
                    attributes=consolidated_attrs,
                    source_file=str(sdrf_path),
                    source_locator=f"Source Name: {s_name} ({len(rows)} rows)",
                    source_key_provenance=key_prov,
                )
            )

        # Build StudyContext
        accession = (
            idf_info.get("Comment[ArrayExpressAccession]")
            or idf_info.get("Comment[SecondaryAccession]")
            or sdrf_path.stem.split(".")[0]
        )
        study_title = idf_info.get("Investigation Title") or f"Study {accession}"

        pub_accessions: list[dict[str, str]] = []
        if accession:
            pub_accessions.append({"namespace": "ArrayExpress", "identifier": str(accession)})
        sec_acc = idf_info.get("Comment[SecondaryAccession]")
        if sec_acc and sec_acc != accession:
            pub_accessions.append({"namespace": "ENA", "identifier": str(sec_acc)})

        context = StudyContext(
            study_id=study_id or f"study_{accession}",
            study_title=study_title,
            source_type="public",
            public_accessions=tuple(pub_accessions),
            experiment_id=experiment_id,
            attributes={"idf": idf_info, "sdrf_file": str(sdrf_path)},
            findings=tuple(findings),
        )

        return context, raw_records, findings
