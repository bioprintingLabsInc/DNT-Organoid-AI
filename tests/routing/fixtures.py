"""Small synthetic file and canonical metadata fixtures for Step 5."""

import gzip
from pathlib import Path

from metadata.fixtures import canonical_bulk_study, reported, with_provenance
from src.routing import HDF5_MAGIC, SRA_MAGIC


FASTQ_TEXT = "@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n"


def write_asset(path: Path, kind: str) -> Path:
    if kind == "fastq":
        path.write_text(FASTQ_TEXT, encoding="utf-8")
    elif kind == "fastq_gz":
        with gzip.open(path, "wt", encoding="utf-8") as stream:
            stream.write(FASTQ_TEXT)
    elif kind == "malformed_fastq":
        path.write_text("@read1\nACGT\nnot-plus\nIII\n", encoding="utf-8")
    elif kind == "integer_matrix":
        path.write_text("gene\tc1\tc2\nGENE1\t1\t0\nGENE2\t2\t3\n", encoding="utf-8")
    elif kind == "processed_matrix":
        path.write_text("gene\ts1\ts2\nGENE1\t1.25\t0.5\nGENE2\t-0.2\t3.7\n", encoding="utf-8")
    elif kind == "empty_matrix":
        path.write_text("gene\ts1\n", encoding="utf-8")
    elif kind == "h5ad":
        path.write_bytes(HDF5_MAGIC + b"synthetic")
    elif kind == "sra":
        path.write_bytes(SRA_MAGIC + b"synthetic")
    elif kind == "rds":
        path.write_bytes(b"X\nsynthetic-r-serialization")
    elif kind == "unknown":
        path.write_bytes(b"unrecognized")
    elif kind == "empty":
        path.write_bytes(b"")
    else:
        raise ValueError(kind)
    return path


def routing_metadata(path: Path, declared_format: str, modality: str = "bulk_rna_seq"):
    document = canonical_bulk_study()
    document["sequencing_assays"][0]["rna_seq_modality"] = reported(modality)
    asset = document["input_data_assets"][0]
    asset["input_data_format"] = reported(declared_format)
    asset["source_file_identifier_or_reference"] = reported(str(path))
    if declared_format not in {"raw_counts", "UMI_counts", "H5AD", "RDS", "processed_expression_matrix"}:
        asset.pop("gene_identifier_type", None)
    return with_provenance(document)


def add_asset(document, path: Path, asset_id: str):
    second = dict(document["input_data_assets"][0])
    second["asset_id"] = asset_id
    second["source_file_identifier_or_reference"] = reported(str(path))
    document["input_data_assets"].append(second)
    return with_provenance(document)

