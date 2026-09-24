"""Deterministic, non-guessing identification and validation of gene identifier types."""

import re
from typing import Any

from .models import Finding, GeneIdentifierType, Severity

_ENSEMBL_PATTERN = re.compile(r"^ENSG\d{11}(\.\d+)?$")
_ENTREZ_PATTERN = re.compile(r"^\d+$")
_HGNC_PATTERN = re.compile(r"^HGNC:\d+$", re.IGNORECASE)
_SYMBOL_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def normalize_declared_type(
    declared: str,
    supported_config: dict[str, Any],
) -> GeneIdentifierType | None:
    """Normalize a declared type string against configured aliases."""
    cleaned = declared.strip().lower()
    for type_key, type_spec in supported_config.items():
        if cleaned == type_key.lower():
            return GeneIdentifierType(type_key)
        aliases = [a.lower() for a in type_spec.get("aliases", [])]
        if cleaned in aliases:
            return GeneIdentifierType(type_key)
    return None


def resolve_identifier_type(
    gene_ids: tuple[str, ...],
    declared_type_raw: str | None,
    supported_config: dict[str, Any],
    asset_id: str,
) -> tuple[GeneIdentifierType | None, list[Finding]]:
    """Determine the gene identifier type using canonical metadata, without guessing when ambiguous."""
    findings: list[Finding] = []

    # Filter non-empty gene IDs for bounded syntax checking
    valid_ids = [gid.strip() for gid in gene_ids if gid and gid.strip()]
    if not valid_ids:
        findings.append(
            Finding(
                Severity.ERROR,
                "empty_gene_identifiers",
                "asset",
                asset_id,
                "matrix.gene_ids",
                "Matrix contains no non-empty gene identifiers.",
            )
        )
        return None, findings

    if declared_type_raw is not None and declared_type_raw.strip():
        resolved = normalize_declared_type(declared_type_raw, supported_config)
        if resolved is None:
            findings.append(
                Finding(
                    Severity.ERROR,
                    "unsupported_identifier_type",
                    "asset",
                    asset_id,
                    "gene_identifier_type",
                    f"Declared gene identifier type '{declared_type_raw}' is not supported by the active configuration.",
                )
            )
            return None, findings

        # Syntax sanity check against declared type
        if resolved == GeneIdentifierType.ENSEMBL_GENE_ID:
            ensembl_matches = sum(1 for gid in valid_ids if _ENSEMBL_PATTERN.match(gid))
            if ensembl_matches == 0:
                findings.append(
                    Finding(
                        Severity.REVIEW,
                        "identifier_type_mismatch",
                        "asset",
                        asset_id,
                        "gene_identifier_type",
                        f"Declared type '{declared_type_raw}' conflicts with matrix identifiers: none match Ensembl Gene ID pattern.",
                    )
                )
        elif resolved == GeneIdentifierType.ENTREZ_GENE_ID:
            entrez_matches = sum(1 for gid in valid_ids if _ENTREZ_PATTERN.match(gid))
            if entrez_matches == 0:
                findings.append(
                    Finding(
                        Severity.REVIEW,
                        "identifier_type_mismatch",
                        "asset",
                        asset_id,
                        "gene_identifier_type",
                        f"Declared type '{declared_type_raw}' conflicts with matrix identifiers: none match Entrez numeric pattern.",
                    )
                )
        elif resolved == GeneIdentifierType.HGNC_ID:
            hgnc_matches = sum(1 for gid in valid_ids if _HGNC_PATTERN.match(gid))
            if hgnc_matches == 0:
                findings.append(
                    Finding(
                        Severity.REVIEW,
                        "identifier_type_mismatch",
                        "asset",
                        asset_id,
                        "gene_identifier_type",
                        f"Declared type '{declared_type_raw}' conflicts with matrix identifiers: none match HGNC ID pattern.",
                    )
                )

        return resolved, findings

    # When undeclared, check if identifiers are unambiguously identifiable
    sample_size = min(len(valid_ids), 100)
    sample_ids = valid_ids[:sample_size]

    ensembl_count = sum(1 for gid in sample_ids if _ENSEMBL_PATTERN.match(gid))
    if ensembl_count == len(sample_ids):
        # All sample IDs unambiguously match Ensembl pattern
        return GeneIdentifierType.ENSEMBL_GENE_ID, findings

    entrez_count = sum(1 for gid in sample_ids if _ENTREZ_PATTERN.match(gid))
    if entrez_count == len(sample_ids):
        # Pure numeric strings without metadata declaration could be Entrez IDs or row indices
        # Per rule: Do not silently guess an identifier type when ambiguous
        findings.append(
            Finding(
                Severity.REVIEW,
                "ambiguous_identifier_type",
                "asset",
                asset_id,
                "gene_identifier_type",
                "Gene identifier type is undeclared; numeric identifiers could represent Entrez IDs or indices.",
            )
        )
        return None, findings

    # Mixed or symbol-like without declaration
    symbol_count = sum(1 for gid in sample_ids if _SYMBOL_PATTERN.match(gid))
    if symbol_count == len(sample_ids) and ensembl_count == 0 and entrez_count == 0:
        # Standard alphabetic symbols
        return GeneIdentifierType.GENE_SYMBOL, findings

    # Truly ambiguous / mixed
    findings.append(
        Finding(
            Severity.REVIEW,
            "ambiguous_identifier_type",
            "asset",
            asset_id,
            "gene_identifier_type",
            "Gene identifier type is not declared in metadata and contains ambiguous or heterogeneous identifiers.",
        )
    )
    return None, findings
