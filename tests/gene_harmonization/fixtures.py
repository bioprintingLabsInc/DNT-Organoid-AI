"""Reusable test fixtures and matrix helpers for gene harmonization tests."""

from copy import deepcopy
from pathlib import Path
from typing import Any

from tests.metadata.fixtures import canonical_bulk_study, reported, with_provenance


def write_test_matrix(path: Path, text: str | None = None) -> Path:
    """Write a synthetic bulk raw-count matrix for harmonization testing."""
    default_text = (
        "gene\tsample_c1\tsample_c2\tsample_c3\tsample_t1\tsample_t2\tsample_t3\n"
        "ENSG00000141510\t10\t20\t30\t40\t50\t60\n"
        "ENSG00000111640\t100\t110\t120\t130\t140\t150\n"
        "ENSG00000075624\t500\t520\t510\t530\t540\t550\n"
    )
    path.write_text(text or default_text, encoding="utf-8")
    return path


def harmonization_metadata(
    path: Path,
    gene_identifier_type: str = "Ensembl Gene ID",
    asset_id: str = "counts_asset",
    assay_id: str = "bulk_assay",
) -> dict[str, list[dict[str, Any]]]:
    """Create a validated canonical metadata document configured for the test asset."""
    doc = canonical_bulk_study()
    asset = doc["input_data_assets"][0]
    asset["asset_id"] = asset_id
    asset["source_file_identifier_or_reference"] = reported(str(path))
    asset["gene_identifier_type"] = reported(gene_identifier_type)
    return with_provenance(doc)
