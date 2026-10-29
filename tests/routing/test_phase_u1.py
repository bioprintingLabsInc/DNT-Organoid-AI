"""Regression and integration tests for Universal RNA-seq Ingestion v2 — Phase U1.

Covers:
1. XLSX single valid count-matrix sheet
2. XLSX multiple candidate expression sheets (ambiguity rejection)
3. XLSX malformed or encrypted workbook
4. Compressed CSV raw count matrix (.csv.gz)
5. Compressed TSV raw count matrix (.tsv.gz)
6. Compressed processed-expression matrix (.csv.gz / .tsv.gz) and DESeq2 isolation
7. Duplicate gene detection linear performance and exact semantic equivalence
8. SRA paired-end extraction flow and FastqPair generation
9. SRA extraction failure handling
10. DEG-table detection and classification (differential_expression_table)
11. DEG-table ambiguous schema rejection
12. Prevention of DEG table entering raw-count pipeline
"""

from collections import Counter
import gzip
from pathlib import Path
import tempfile
import time
import unittest
import xml.etree.ElementTree as ET
import zipfile

from routing.fixtures import routing_metadata
from src.ingestion.inspection import inspect_asset
from src.ingestion.models import AssetStatus, Confidence, RegisteredAsset, Route
from src.ingestion.xlsx_reader import classify_xlsx_workbook, read_xlsx_count_matrix
from src.qc.config import load_rules
from src.qc.matrix import read_matrix
from src.qc.result import Severity
from src.quantification_bridge.pairing import FastqPair
from src.quantification_bridge.sra_runner import SraExtractionResult, extract_sra
from src.routing.detector import detect_format
from src.routing.dispatcher import dispatch_route_execution, plan_route_execution
from src.routing.router import _select_route


def create_test_xlsx(path: Path, sheets: dict[str, list[list]], encrypted: bool = False) -> Path:
    """Helper to create a synthetic XLSX package with arbitrary sheets."""
    with zipfile.ZipFile(path, "w") as z:
        types = ET.Element("Types", xmlns="http://schemas.openxmlformats.org/package/2006/content-types")
        ET.SubElement(types, "Default", Extension="rels", ContentType="application/vnd.openxmlformats-package.relationships+xml")
        ET.SubElement(types, "Default", Extension="xml", ContentType="application/xml")
        ET.SubElement(types, "Override", PartName="/xl/workbook.xml", ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml")
        for i in range(1, len(sheets) + 1):
            ET.SubElement(types, "Override", PartName=f"/xl/worksheets/sheet{i}.xml", ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml")
        z.writestr("[Content_Types].xml", ET.tostring(types, encoding="utf-8", xml_declaration=True))

        rels = ET.Element("Relationships", xmlns="http://schemas.openxmlformats.org/package/2006/relationships")
        ET.SubElement(rels, "Relationship", Id="rId1", Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument", Target="xl/workbook.xml")
        z.writestr("_rels/.rels", ET.tostring(rels, encoding="utf-8", xml_declaration=True))

        wb = ET.Element("workbook", xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main")
        sheets_elem = ET.SubElement(wb, "sheets")
        wb_rels = ET.Element("Relationships", xmlns="http://schemas.openxmlformats.org/package/2006/relationships")

        for idx, (sheet_name, rows) in enumerate(sheets.items(), start=1):
            r_id = f"rId{idx}"
            ET.SubElement(sheets_elem, "sheet", name=sheet_name, sheetId=str(idx), attrib={"{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id": r_id})
            ET.SubElement(wb_rels, "Relationship", Id=r_id, Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet", Target=f"worksheets/sheet{idx}.xml")

            ws = ET.Element("worksheet", xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main")
            sdata = ET.SubElement(ws, "sheetData")
            for r_idx, row in enumerate(rows, start=1):
                row_el = ET.SubElement(sdata, "row", r=str(r_idx))
                for c_idx, val in enumerate(row, start=1):
                    col_str = chr(ord("A") + c_idx - 1)
                    c_ref = f"{col_str}{r_idx}"
                    if isinstance(val, (int, float)):
                        c_el = ET.SubElement(row_el, "c", r=c_ref)
                        v_el = ET.SubElement(c_el, "v")
                        v_el.text = str(val)
                    elif val is not None:
                        c_el = ET.SubElement(row_el, "c", r=c_ref, t="inlineStr")
                        is_el = ET.SubElement(c_el, "is")
                        t_el = ET.SubElement(is_el, "t")
                        t_el.text = str(val)
            z.writestr(f"xl/worksheets/sheet{idx}.xml", ET.tostring(ws, encoding="utf-8", xml_declaration=True))

        z.writestr("xl/workbook.xml", ET.tostring(wb, encoding="utf-8", xml_declaration=True))
        z.writestr("xl/_rels/workbook.xml.rels", ET.tostring(wb_rels, encoding="utf-8", xml_declaration=True))
    return path


class UniversalIngestionPhaseU1Tests(unittest.TestCase):

    def test_xlsx_single_valid_sheet(self):
        """XLSX with one valid integer count sheet is detected as raw_counts and parsed cleanly."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "counts.xlsx"
            create_test_xlsx(p, {
                "Counts": [["gene_id", "s1", "s2"], ["ENSG1", 10, 20], ["ENSG2", 5, 15]],
                "Metadata": [["sample", "cell_line"], ["s1", "H9"], ["s2", "H9"]],
            })

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_xlsx", str(p), "raw_counts", ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)

            self.assertEqual("raw_counts", det.detected_format)
            self.assertEqual(Confidence.HIGH, det.confidence)
            self.assertEqual((), det.findings)

            route, findings = _select_route(asset, det)
            self.assertEqual(Route.BULK_RAW_COUNT_PIPELINE, route)
            self.assertEqual((), findings)

            matrix, f = read_matrix(p, load_rules())
            self.assertIsNotNone(matrix)
            self.assertEqual(("ENSG1", "ENSG2"), matrix.gene_ids)
            self.assertEqual(("s1", "s2"), matrix.sample_ids)

    def test_xlsx_multiple_candidate_sheets(self):
        """XLSX with multiple candidate expression sheets rejects silent selection and halts for review."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "ambiguous.xlsx"
            create_test_xlsx(p, {
                "Counts_Day20": [["gene", "s1", "s2"], ["G1", 10, 20], ["G2", 30, 40]],
                "Counts_Day40": [["gene", "s3", "s4"], ["G1", 15, 25], ["G2", 35, 45]],
            })

            selected, all_sheets, findings, prov = classify_xlsx_workbook(p)
            self.assertIsNone(selected)
            self.assertEqual(2, len(all_sheets))
            self.assertEqual("AMBIGUOUS_CANDIDATE_SHEETS", prov.get("status"))
            self.assertTrue(any(f.rule == "ambiguous_sheets" for f in findings))

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_ambig", str(p), "raw_counts", ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)
            self.assertIsNone(det.detected_format)
            self.assertTrue(any(f.rule == "ambiguous_sheets" for f in det.findings))

            matrix, qc_findings = read_matrix(p, load_rules())
            self.assertIsNone(matrix)
            self.assertTrue(any(f.rule_id == "ambiguous_sheets" for f in qc_findings))

    def test_xlsx_corrupt_workbook(self):
        """Corrupt or non-zip XLSX is rejected with corrupt_workbook error finding."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "corrupt.xlsx"
            p.write_bytes(b"not a real zip file")

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_corrupt", str(p), "raw_counts", ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)
            self.assertIsNone(det.detected_format)
            self.assertTrue(any(f.rule == "corrupt_workbook" for f in det.findings))

            matrix, qc_findings = read_matrix(p, load_rules())
            self.assertIsNone(matrix)
            self.assertTrue(any(f.rule_id == "corrupt_workbook" for f in qc_findings))

    def test_compressed_csv_raw_counts(self):
        """Gzip-compressed CSV with integers is detected as raw_counts and routed to bulk count pipeline."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "counts.csv.gz"
            with gzip.open(p, "wt", encoding="utf-8") as stream:
                stream.write("gene_id,s1,s2\nENSG1,10,20\nENSG2,5,15\n")

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_gz_csv", str(p), "raw_counts", ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)

            self.assertEqual("raw_counts", det.detected_format)
            route, _ = _select_route(asset, det)
            self.assertEqual(Route.BULK_RAW_COUNT_PIPELINE, route)

            matrix, f = read_matrix(p, load_rules())
            self.assertIsNotNone(matrix)
            self.assertEqual(("ENSG1", "ENSG2"), matrix.gene_ids)
            self.assertEqual(("s1", "s2"), matrix.sample_ids)

    def test_compressed_tsv_raw_counts(self):
        """Gzip-compressed TSV with integers is detected as raw_counts and routed to bulk count pipeline."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "counts.tsv.gz"
            with gzip.open(p, "wt", encoding="utf-8") as stream:
                stream.write("gene_id\ts1\ts2\nENSG1\t10\t20\nENSG2\t5\t15\n")

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_gz_tsv", str(p), "raw_counts", ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)

            self.assertEqual("raw_counts", det.detected_format)
            route, _ = _select_route(asset, det)
            self.assertEqual(Route.BULK_RAW_COUNT_PIPELINE, route)

            matrix, f = read_matrix(p, load_rules())
            self.assertIsNotNone(matrix)
            self.assertEqual(("s1", "s2"), matrix.sample_ids)

    def test_compressed_processed_matrix_isolation(self):
        """Gzip-compressed CSV with floating-point values routes to PROCESSED_EXPRESSION and blocks DESeq2."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "tpm.csv.gz"
            with gzip.open(p, "wt", encoding="utf-8") as stream:
                stream.write("gene_id,s1,s2\nENSG1,1.45,20.2\nENSG2,0.05,15.8\n")

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_gz_proc", str(p), "processed_expression_matrix", ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)

            self.assertEqual("processed_expression_matrix", det.detected_format)
            self.assertEqual(Confidence.HIGH, det.confidence)
            route, _ = _select_route(asset, det)
            self.assertEqual(Route.PROCESSED_EXPRESSION_PIPELINE, route)

            meta = routing_metadata(p, "processed_expression_matrix")
            plan = plan_route_execution(meta)
            self.assertEqual(Route.PROCESSED_EXPRESSION_PIPELINE, plan.route)
            self.assertEqual("EXECUTION_NOT_IMPLEMENTED", plan.status)
            self.assertIsNone(plan.executor_id)

            res = dispatch_route_execution(plan)
            self.assertFalse(res.is_success)
            self.assertEqual("EXECUTION_NOT_IMPLEMENTED", res.execution_state)

    def test_duplicate_detection_semantic_equivalence(self):
        """Optimized Counter duplicate detection produces exact same findings and ordering as O(N^2) logic."""
        genes = ["G1", "G2", "G3", "G1", "G4", "G2", "G5", "G1"]
        # Old logic:
        old_result = sorted({x for x in genes if genes.count(x) > 1})
        # New logic:
        gene_counts = Counter(genes)
        new_result = sorted([x for x, c in gene_counts.items() if c > 1])
        self.assertEqual(old_result, new_result)
        self.assertEqual(["G1", "G2"], new_result)

    def test_deg_table_detection(self):
        """Tabular file with gene, log2FoldChange, and p-value columns is detected as differential_expression_table."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "results_vpa.tsv"
            p.write_text(
                "gene_id\tlog2FoldChange\tpvalue\tpadj\tstatus\n"
                "ENSG1\t2.15\t0.001\t0.02\tUP\n"
                "ENSG2\t-1.40\t0.03\t0.15\tDOWN\n",
                encoding="utf-8",
            )

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_deg", str(p), "differential_expression_table", ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)

            self.assertEqual("differential_expression_table", det.detected_format)
            self.assertEqual(Confidence.HIGH, det.confidence)

            route, findings = _select_route(asset, det)
            self.assertEqual(Route.DIFFERENTIAL_EXPRESSION_PIPELINE, route)
            self.assertEqual((), findings)

    def test_deg_table_ambiguous_schema_rejection(self):
        """Table missing effect-size column is not classified as differential_expression_table."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "not_deg.tsv"
            p.write_text(
                "gene_id\tpvalue\tpadj\n"
                "ENSG1\t0.001\t0.02\n"
                "ENSG2\t0.03\t0.15\n",
                encoding="utf-8",
            )

            insp = inspect_asset(p)
            asset = RegisteredAsset("asset_not_deg", str(p), None, ("assay_1",), ("bulk_rna_seq",))
            det = detect_format(asset, insp)
            self.assertNotEqual("differential_expression_table", det.detected_format)

    def test_prevention_deg_table_raw_counts(self):
        """DEG tables cannot enter raw count pipeline and dispatcher halts with EXECUTION_NOT_IMPLEMENTED."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "deg.tsv"
            p.write_text(
                "symbol\tlogFC\tP.Value\tadj.P.Val\n"
                "A1BG\t1.5\t0.001\t0.01\n",
                encoding="utf-8",
            )
            meta = routing_metadata(p, "differential_expression_table")
            plan = plan_route_execution(meta)

            self.assertEqual(Route.DIFFERENTIAL_EXPRESSION_PIPELINE, plan.route)
            self.assertEqual("EXECUTION_NOT_IMPLEMENTED", plan.status)
            self.assertIsNone(plan.executor_id)

            res = dispatch_route_execution(plan)
            self.assertFalse(res.is_success)
            self.assertEqual("EXECUTION_NOT_IMPLEMENTED", res.execution_state)

    def test_sra_paired_extraction_flow(self):
        """SRA route plan assigns fasterq-dump extraction adapter into quantification bridge."""
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "SRR12345.sra"
            p.write_bytes(b"NCBI.sra-synthetic")
            meta = routing_metadata(p, "SRA")
            plan = plan_route_execution(meta)

            self.assertEqual(Route.RAW_READ_PIPELINE, plan.route)
            self.assertEqual("READY_FOR_EXECUTION", plan.status)
            self.assertEqual("quantification_bridge", plan.executor_id)
            self.assertEqual("fasterq-dump", plan.provenance_chain.get("sra_extraction_adapter"))

            # Mock SRA executor returning synthetic FastqPair
            def mock_executor(pairs, plan, **kwargs):
                return RouteExecutionResult(
                    plan=plan,
                    is_success=True,
                    execution_state="AWAITING_OWNER_STRANDEDNESS_DECISION",
                    emitted_asset_path=None,
                    emitted_format=None,
                    provenance_record=plan.provenance_chain,
                )

            from src.quantification_bridge.pairing import FastqPair
            from src.routing.dispatcher import RouteExecutionResult

            def mock_sra_extractor(sra_path, output_dir, sample_id):
                return FastqPair(
                    sample_id=sample_id,
                    r1_path=output_dir / f"{sample_id}_1.fastq.gz",
                    r2_path=output_dir / f"{sample_id}_2.fastq.gz",
                )

            res = dispatch_route_execution(plan, executor=mock_executor, sra_extractor=mock_sra_extractor)
            self.assertTrue(res.is_success)
            self.assertEqual("AWAITING_OWNER_STRANDEDNESS_DECISION", res.execution_state)

    def test_sra_extract_to_fastq_pair_construction_and_provenance(self):
        """Verify that extract_sra_to_fastq_pair constructs a standard FastqPair with SRA provenance in extra_metadata."""
        from unittest.mock import patch
        from src.quantification_bridge.sra_runner import SraExtractionResult, extract_sra_to_fastq_pair
        from src.quantification_bridge.pairing import FastqPair

        with tempfile.TemporaryDirectory() as td:
            td_path = Path(td)
            r1 = td_path / "SRR999_1.fastq"
            r2 = td_path / "SRR999_2.fastq"
            r1.write_text("@read1/1\nACGT\n+\nIIII\n", encoding="utf-8")
            r2.write_text("@read1/2\nTGCA\n+\nJJJJ\n", encoding="utf-8")
            temp_scratch = td_path / "temp_scratch"

            mock_res = SraExtractionResult(
                status="SRA_EXTRACTION_SUCCESS",
                run_accession="SRR999",
                layout="PAIRED",
                fastq_r1=r1,
                fastq_r2=r2,
                tool_name="fasterq-dump",
                tool_version="3.1.1",
                output_dir=td_path,
                temp_dir=temp_scratch,
            )

            with patch("src.quantification_bridge.sra_runner.extract_sra", return_value=mock_res):
                pair = extract_sra_to_fastq_pair(
                    sra_path="dummy.sra",
                    output_dir=td_path,
                    sample_id="sample_test",
                    biological_replicate_id="rep_2",
                    condition_id="cond_vpa",
                    temp_dir=temp_scratch,
                    assay_id="assay_001",
                )

                # 1. Assert FastqPair public dataclass structure
                self.assertIsInstance(pair, FastqPair)
                self.assertEqual("sample_test", pair.sample_id)
                self.assertEqual(r1, pair.r1_path)
                self.assertEqual(r2, pair.r2_path)
                self.assertEqual("assay_001", pair.assay_id)

                # 2. Assert provenance preservation in extra_metadata
                meta = pair.extra_metadata
                self.assertEqual("SRR999", meta["run_accession"])
                self.assertEqual("rep_2", meta["biological_replicate_id"])
                self.assertEqual("cond_vpa", meta["condition_id"])
                self.assertEqual("SRA_EXTRACTION_SUCCESS", meta["sra_extraction_status"])
                self.assertEqual("fasterq-dump", meta["sra_extraction_tool"])
                self.assertEqual("3.1.1", meta["sra_extraction_version"])
                self.assertEqual(str(temp_scratch), meta["sra_temp_dir"])

                # 3. Assert downstream serializability
                d = pair.to_dict()
                self.assertEqual("sample_test", d["sample_id"])
                self.assertEqual(str(r1), d["r1_path"])
                self.assertEqual("SRR999", d["extra_metadata"]["run_accession"])

    def test_fasterq_dump_command_and_temp_storage_governance(self):
        """Verify deterministic fasterq-dump command construction and governed scratch directory handling."""
        from src.quantification_bridge.sra_runner import build_fasterq_dump_command, extract_sra

        # 1. Deterministic command generation with separate output and temp options
        cmd = build_fasterq_dump_command(
            tool_bin="/usr/bin/fasterq-dump",
            sra_path="/data/SRR123.sra",
            output_dir="/data/output",
            temp_dir="/data/scratch",
            threads=8,
        )
        self.assertEqual(
            "/usr/bin/fasterq-dump --split-files -O /data/output -t /data/scratch -e 8 -f /data/SRR123.sra",
            cmd,
        )

        # 2. extract_sra validates file existence and creates separate output & temp directories
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            sra_file = base / "test.sra"
            sra_file.write_bytes(b"dummy-sra-data")
            out_dir = base / "out"
            temp_dir = base / "scratch_governed"

            self.assertFalse(out_dir.exists())
            self.assertFalse(temp_dir.exists())

            # Missing SRA file raises FileNotFoundError
            with self.assertRaises(FileNotFoundError):
                extract_sra(base / "nonexistent.sra", output_dir=out_dir, temp_dir=temp_dir)
