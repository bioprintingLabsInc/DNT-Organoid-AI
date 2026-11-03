"""Tabular source metadata adapter for CSV, TSV, and XLSX sample tables.

Extracts sample annotations, exposure metadata, and biological parameters while
preserving 100% of the original column headers and cell values for provenance tracking.
Zero external dependency: reads Excel (.xlsx) workbooks using standard library zipfile and xml.etree.
"""

from collections.abc import Iterable
import csv
import io
from pathlib import Path
import re
from typing import Any
import xml.etree.ElementTree as ET
import zipfile

import pandas as pd

from src.ingestion.xlsx_reader import (
    NS_MAIN,
    NS_PKG_REL,
    NS_REL,
    parse_cell_ref,
    read_shared_strings,
)
from .base import MappingFinding, MetadataSourceAdapter, RawSampleRecord, StudyContext

# Common column aliases for sample identifier
SAMPLE_ID_ALIASES = {
    "sample_id",
    "sample",
    "sample_name",
    "sample_alias",
    "samplename",
    "sampleid",
    "gsm",
    "gsm_accession",
    "source_name",
    "source name",
    "run_accession",
    "run",
    "id",
}


def _detect_delimiter(sample_text: str) -> str:
    """Detect comma, tab, or semicolon delimiter safely."""
    first_line = sample_text.splitlines()[0] if sample_text else ""
    counts = {
        "\t": first_line.count("\t"),
        ",": first_line.count(","),
        ";": first_line.count(";"),
    }
    sorted_delims = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    if sorted_delims[0][1] > 0:
        return sorted_delims[0][0]
    return ","


def read_xlsx_sheet_to_df(path: Path, sheet_name_or_index: str | int = 0) -> pd.DataFrame:
    """Pure Python standard library XLSX worksheet reader (zero external dependencies)."""
    with zipfile.ZipFile(path, "r") as z:
        shared = read_shared_strings(z)
        wb_root = ET.fromstring(z.read("xl/workbook.xml"))
        sheets_elem = wb_root.find(f"{{{NS_MAIN}}}sheets")
        if sheets_elem is None:
            raise ValueError("No sheets found in workbook.")

        sheet_nodes = sheets_elem.findall(f"{{{NS_MAIN}}}sheet")
        if not sheet_nodes:
            raise ValueError("No sheets found in workbook.")

        target_node = None
        if isinstance(sheet_name_or_index, int):
            if sheet_name_or_index < len(sheet_nodes):
                target_node = sheet_nodes[sheet_name_or_index]
            else:
                raise ValueError(f"Sheet index {sheet_name_or_index} out of range ({len(sheet_nodes)} sheets).")
        else:
            for node in sheet_nodes:
                if node.attrib.get("name") == sheet_name_or_index:
                    target_node = node
                    break
            if target_node is None:
                raise ValueError(f"Sheet '{sheet_name_or_index}' not found in workbook.")

        r_id = target_node.attrib.get(f"{{{NS_REL}}}id")
        rels_root = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        target_path = None
        for rel in rels_root.findall(f"{{{NS_PKG_REL}}}Relationship"):
            if rel.attrib.get("Id") == r_id:
                target_path = "xl/" + rel.attrib.get("Target").lstrip("/")
                break

        if not target_path or target_path not in z.namelist():
            idx = sheet_name_or_index if isinstance(sheet_name_or_index, int) else sheet_nodes.index(target_node)
            target_path = f"xl/worksheets/sheet{idx + 1}.xml"

        rows_data: dict[int, dict[int, Any]] = {}
        with z.open(target_path) as stream:
            for event, elem in ET.iterparse(stream, events=("end",)):
                if elem.tag == f"{{{NS_MAIN}}}row":
                    r_attr = elem.attrib.get("r")
                    row_idx = int(r_attr) - 1 if r_attr else len(rows_data)
                    row_cells = {}
                    for c in elem.findall(f"{{{NS_MAIN}}}c"):
                        ref = c.attrib.get("r", "")
                        c_type = c.attrib.get("t", "")
                        _, col_idx = parse_cell_ref(ref)

                        v_elem = c.find(f"{{{NS_MAIN}}}v")
                        val_text = v_elem.text if v_elem is not None else None

                        if c_type == "s" and val_text is not None:
                            try:
                                s_idx = int(val_text)
                                cell_val = shared[s_idx] if s_idx < len(shared) else val_text
                            except ValueError:
                                cell_val = val_text
                        elif c_type == "b" and val_text is not None:
                            cell_val = val_text == "1"
                        else:
                            cell_val = val_text

                        row_cells[col_idx] = cell_val
                    rows_data[row_idx] = row_cells
                    elem.clear()

    if not rows_data:
        return pd.DataFrame()

    sorted_rows = sorted(rows_data.keys())
    min_row = sorted_rows[0]
    max_row = sorted_rows[-1]
    all_col_indices = [col_idx for cells in rows_data.values() for col_idx in cells.keys()]
    if not all_col_indices:
        return pd.DataFrame()
    min_col = min(all_col_indices)
    max_col = max(all_col_indices)

    header_cells = rows_data.get(min_row, {})
    headers = [str(header_cells.get(c, f"col_{c}") or f"col_{c}").strip() for c in range(min_col, max_col + 1)]

    data_rows = []
    for r in sorted_rows[1:]:
        row_dict = rows_data.get(r, {})
        row_vals = [row_dict.get(c) for c in range(min_col, max_col + 1)]
        if any(v is not None for v in row_vals):
            data_rows.append(row_vals)

    return pd.DataFrame(data_rows, columns=headers)


class TabularMetadataAdapter(MetadataSourceAdapter):
    """Adapter for structured sample metadata tables in CSV, TSV, or XLSX format."""

    def can_read(self, source: Any) -> bool:
        """Check if source is a tabular file (CSV, TSV, XLSX) or DataFrame."""
        if isinstance(source, pd.DataFrame):
            return True
        if isinstance(source, (str, Path)):
            p = Path(source)
            if not p.exists() or not p.is_file():
                return False
            name = p.name.lower()
            return name.endswith((".csv", ".tsv", ".txt", ".tab", ".xlsx", ".csv.gz", ".tsv.gz"))
        return False

    def parse_study(
        self,
        source: Any,
        study_id: str | None = None,
        experiment_id: str = "exp_01",
        sheet_name: str | int | None = 0,
        sample_id_column: str | None = None,
        **kwargs: Any,
    ) -> tuple[StudyContext, list[RawSampleRecord], list[MappingFinding]]:
        """Parse tabular sample metadata table."""
        findings: list[MappingFinding] = []
        source_name = str(source) if isinstance(source, (str, Path)) else "in_memory_table"

        df: pd.DataFrame
        if isinstance(source, pd.DataFrame):
            df = source.copy()
        else:
            p = Path(source)
            filename = p.name.lower()
            if filename.endswith(".xlsx"):
                try:
                    df = read_xlsx_sheet_to_df(p, sheet_name_or_index=sheet_name or 0)
                except Exception as e:
                    findings.append(
                        MappingFinding(
                            severity="ERROR",
                            rule="tabular_read_error",
                            field_path=source_name,
                            message=f"Failed to read Excel workbook: {e}",
                        )
                    )
                    return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings
            else:
                # CSV or TSV
                try:
                    raw_bytes = p.read_bytes()
                    if raw_bytes.startswith(b"\x1f\x8b"):
                        import gzip
                        raw_bytes = gzip.decompress(raw_bytes)
                    text = raw_bytes.decode("utf-8-sig", errors="replace")
                    sep = _detect_delimiter(text)
                    df = pd.read_csv(io.StringIO(text), sep=sep)
                except Exception as e:
                    findings.append(
                        MappingFinding(
                            severity="ERROR",
                            rule="tabular_read_error",
                            field_path=source_name,
                            message=f"Failed to read delimited text table: {e}",
                        )
                    )
                    return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

        if df.empty:
            findings.append(
                MappingFinding(
                    severity="ERROR",
                    rule="empty_table",
                    field_path=source_name,
                    message="Metadata table is empty.",
                )
            )
            return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

        # Discover sample ID column
        col_map_lower = {str(c).strip().lower(): c for c in df.columns}
        chosen_sample_col: Any = None
        if sample_id_column and sample_id_column in df.columns:
            chosen_sample_col = sample_id_column
        elif sample_id_column and sample_id_column.lower() in col_map_lower:
            chosen_sample_col = col_map_lower[sample_id_column.lower()]
        else:
            for alias in SAMPLE_ID_ALIASES:
                if alias in col_map_lower:
                    chosen_sample_col = col_map_lower[alias]
                    break

        if not chosen_sample_col:
            # Fallback to first column if it looks like identifiers
            chosen_sample_col = df.columns[0]
            findings.append(
                MappingFinding(
                    severity="WARNING",
                    rule="sample_id_column_inferred",
                    field_path=str(chosen_sample_col),
                    message=f"No explicit sample identifier column found; defaulting to first column '{chosen_sample_col}'.",
                )
            )

        # Build raw sample records
        raw_records: list[RawSampleRecord] = []
        seen_sids: set[str] = set()

        for idx, row in df.iterrows():
            raw_sid = row[chosen_sample_col]
            sid_str = str(raw_sid).strip() if pd.notna(raw_sid) else ""
            if not sid_str:
                findings.append(
                    MappingFinding(
                        severity="ERROR",
                        rule="missing_sample_id",
                        field_path=f"row_{idx}.{chosen_sample_col}",
                        message=f"Row {idx} is missing a sample identifier.",
                    )
                )
                continue

            if sid_str in seen_sids:
                findings.append(
                    MappingFinding(
                        severity="ERROR",
                        rule="duplicate_sample_id",
                        field_path=f"row_{idx}.{chosen_sample_col}",
                        message=f"Duplicate sample identifier '{sid_str}' encountered at row {idx}.",
                    )
                )
            seen_sids.add(sid_str)

            attrs: dict[str, Any] = {}
            key_prov: dict[str, str] = {}
            for col in df.columns:
                val = row[col]
                clean_val = None if pd.isna(val) else val
                col_str = str(col).strip()
                attrs[col_str] = clean_val
                key_prov[col_str] = str(col)

            raw_records.append(
                RawSampleRecord(
                    sample_id=sid_str,
                    attributes=attrs,
                    source_file=source_name,
                    source_locator=f"row_{idx}",
                    source_key_provenance=key_prov,
                )
            )

        inferred_study_id = study_id or f"study_{Path(source_name).stem}"
        context = StudyContext(
            study_id=inferred_study_id,
            experiment_id=experiment_id,
            attributes={"source_table": source_name, "total_rows": len(raw_records)},
            findings=tuple(findings),
        )

        return context, raw_records, findings
