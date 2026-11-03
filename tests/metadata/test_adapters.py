"""Unit tests for source metadata adapters (Tabular, MAGE-TAB, GEO, JSON)."""

import io
from pathlib import Path
import tempfile
import pytest

import pandas as pd

from src.metadata.adapters import (
    GeoMetadataAdapter,
    JsonMetadataAdapter,
    MageTabAdapter,
    TabularMetadataAdapter,
)


class TestTabularMetadataAdapter:
    def test_csv_metadata_read(self, tmp_path):
        csv_file = tmp_path / "sample_sheet.csv"
        csv_file.write_text(
            "sample_id,condition,dose,dose_unit,duration,replicate\n"
            "s1,control,0,uM,48,rep1\n"
            "s2,control,0,uM,48,rep2\n"
            "s3,treatment,200,uM,48,rep1\n"
            "s4,treatment,200,uM,48,rep2\n",
            encoding="utf-8",
        )
        adapter = TabularMetadataAdapter()
        assert adapter.can_read(csv_file)

        ctx, recs, findings = adapter.parse_study(csv_file)
        assert len(recs) == 4
        assert recs[0].sample_id == "s1"
        assert recs[2].attributes["dose"] == 200
        assert recs[2].attributes["condition"] == "treatment"
        assert not any(f.severity == "ERROR" for f in findings)

    def test_tsv_metadata_read(self, tmp_path):
        tsv_file = tmp_path / "sample_sheet.tsv"
        tsv_file.write_text(
            "sample_name\ttreatment\tdose\ttime\n"
            "sample_A\tvehicle\t0\t24h\n"
            "sample_B\tvalproic acid\t100\t24h\n",
            encoding="utf-8",
        )
        adapter = TabularMetadataAdapter()
        assert adapter.can_read(tsv_file)

        ctx, recs, findings = adapter.parse_study(tsv_file)
        assert len(recs) == 2
        assert recs[0].sample_id == "sample_A"
        assert recs[1].attributes["treatment"] == "valproic acid"

    def test_xlsx_metadata_read(self):
        fixture_xlsx = Path(__file__).resolve().parents[1] / "fixtures" / "Table_S1.xlsx"
        assert fixture_xlsx.exists(), "Table_S1.xlsx missing"

        adapter = TabularMetadataAdapter()
        assert adapter.can_read(fixture_xlsx)

        ctx, recs, findings = adapter.parse_study(fixture_xlsx)
        assert len(recs) > 0
        assert recs[0].sample_id != ""

    def test_empty_table_fails_safely(self, tmp_path):
        empty_file = tmp_path / "empty.csv"
        empty_file.write_text("", encoding="utf-8")
        adapter = TabularMetadataAdapter()
        ctx, recs, findings = adapter.parse_study(empty_file)
        assert len(recs) == 0
        assert any(f.rule == "empty_table" or f.rule == "tabular_read_error" for f in findings)

    def test_duplicate_sample_id_detected(self, tmp_path):
        dup_file = tmp_path / "dup.csv"
        dup_file.write_text(
            "sample_id,condition\n"
            "s1,control\n"
            "s1,treatment\n",
            encoding="utf-8",
        )
        adapter = TabularMetadataAdapter()
        ctx, recs, findings = adapter.parse_study(dup_file)
        assert any(f.rule == "duplicate_sample_id" for f in findings)


class TestMageTabAdapter:
    def test_parse_sdrf_and_idf(self, tmp_path):
        idf_file = tmp_path / "E-TEST-01.idf.txt"
        idf_file.write_text(
            "Investigation Title\tDNT Test Study on Human Cerebral Organoids\n"
            "Comment[ArrayExpressAccession]\tE-TEST-01\n"
            "Comment[SecondaryAccession]\tERP999999\n",
            encoding="utf-8",
        )

        sdrf_file = tmp_path / "E-TEST-01.sdrf.txt"
        sdrf_file.write_text(
            "Source Name\tCharacteristics[organism]\tCharacteristics[organism part]\tFactor Value[compound]\tFactor Value[dose]\tUnit[concentration unit]\tComment[ENA_RUN]\n"
            "Sample 1\tHomo sapiens\tcerebral organoid\tvalproic acid\t200\tmicromolar\tERR10001\n"
            "Sample 1\tHomo sapiens\tcerebral organoid\tvalproic acid\t200\tmicromolar\tERR10002\n"
            "Sample 2\tHomo sapiens\tcerebral organoid\tnone\t0\tmicromolar\tERR10003\n",
            encoding="utf-8",
        )

        adapter = MageTabAdapter()
        assert adapter.can_read(sdrf_file)
        assert adapter.can_read(idf_file)

        ctx, recs, findings = adapter.parse_study(sdrf_file, idf_source=idf_file)
        # Sample 1 had 2 rows (e.g. paired reads), so consolidated records = 2
        assert len(recs) == 2
        assert recs[0].sample_id == "Sample 1"
        assert recs[0].attributes["Factor Value[compound]"] == "valproic acid"
        assert recs[0].attributes["Factor Value[dose]"] == "200"
        assert recs[1].sample_id == "Sample 2"
        assert ctx.study_id == "study_E-TEST-01"
        assert any(a["identifier"] == "E-TEST-01" for a in ctx.public_accessions)


class TestGeoMetadataAdapter:
    def test_parse_geo_brief_text(self, tmp_path):
        geo_file = tmp_path / "GSE99999_brief.txt"
        geo_file.write_text(
            "^SAMPLE = GSM100001\n"
            "!Sample_title = Cortical organoid DMSO rep 1\n"
            "!Sample_description = Vehicle control\n"
            "!Sample_source_name_ch1 = Cortical organoid\n"
            "!Sample_characteristics_ch1 = treatment: vehicle control\n"
            "!Sample_characteristics_ch1 = cell line: WTC11\n"
            "!Sample_characteristics_ch1 = developmental stage: day 45\n"
            "^SAMPLE = GSM100002\n"
            "!Sample_title = Cortical organoid VPA 200uM rep 1\n"
            "!Sample_description = VPA treated\n"
            "!Sample_source_name_ch1 = Cortical organoid\n"
            "!Sample_characteristics_ch1 = treatment: valproic acid\n"
            "!Sample_characteristics_ch1 = dose: 200 uM\n"
            "!Sample_characteristics_ch1 = cell line: WTC11\n"
            "!Sample_characteristics_ch1 = developmental stage: day 45\n",
            encoding="utf-8",
        )

        adapter = GeoMetadataAdapter()
        assert adapter.can_read(geo_file)

        ctx, recs, findings = adapter.parse_study(geo_file, study_id="study_GSE99999")
        assert len(recs) == 2
        assert recs[0].sample_id == "GSM100001"
        assert recs[0].attributes["treatment"] == "vehicle control"
        assert recs[0].attributes["cell_line"] == "WTC11"
        assert recs[1].sample_id == "GSM100002"
        assert recs[1].attributes["dose"] == "200 uM"


class TestJsonMetadataAdapter:
    def test_parse_study_manifest(self, tmp_path):
        json_file = tmp_path / "study_manifest.json"
        json_file.write_text(
            '{\n'
            '  "bioproject": "PRJNA802072",\n'
            '  "sra_study": "SRP357550",\n'
            '  "contrast": "VPA vs Control",\n'
            '  "samples": [\n'
            '    {"sample_id": "VPA_0_1", "group": "control", "dose": 0, "day": 28},\n'
            '    {"sample_id": "VPA_250_1", "group": "treatment", "dose": 250, "day": 28}\n'
            '  ]\n'
            '}',
            encoding="utf-8",
        )

        adapter = JsonMetadataAdapter()
        assert adapter.can_read(json_file)

        ctx, recs, findings = adapter.parse_study(json_file)
        assert len(recs) == 2
        assert recs[0].sample_id == "VPA_0_1"
        assert recs[1].attributes["dose"] == 250
        assert ctx.study_id == "PRJNA802072"
