"""Orchestration of Gene Harmonization v1 for human bulk RNA-seq expression data."""

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from src.ingestion.registry import reported_value
from src.qc.matrix import CountMatrix, read_matrix
from src.qc.result import DatasetResult, Status as QCStatus

from .config import config_checksum, load_config, validate_config
from .models import (
    CanonicalCollision,
    Finding,
    GeneIdentifierType,
    GeneMappingStatus,
    HarmonizationSummary,
    HarmonizedDataset,
    HarmonizedGene,
    QCReviewDisposition,
    Severity,
    Status,
    ordered,
    status_for,
)
from .reference import GeneReference
from .type_detector import resolve_identifier_type

HARMONIZATION_VERSION = "bulk_gene_harmonization_v1"


def harmonize_bulk_counts(
    matrix: CountMatrix | str | Path,
    metadata: dict[str, Any],
    asset_id: str,
    assay_id: str,
    qc_result: DatasetResult | None = None,
    qc_disposition: QCReviewDisposition | dict[str, Any] | str | None = None,
    allow_needs_review_qc: bool = False,
    config_path: str | Path | None = None,
    reference: GeneReference | None = None,
    config: dict[str, Any] | None = None,
) -> HarmonizedDataset:
    """Harmonize gene identifiers from a QC-passed bulk expression matrix into canonical Ensembl IDs.

    Preserves original identifiers, expression values, and mapping provenance.
    Never aggregates duplicate or colliding genes.
    """
    findings: list[Finding] = []

    # 1. Load and validate configuration
    if config is None:
        config = load_config(config_path)
    else:
        validate_config(config)
    cfg_checksum = config_checksum(config)
    structural_rules = config.get("structural_rules", {})
    coverage_thresholds = config.get("coverage_thresholds", {})

    upstream_qc_status_str: str | None = None
    review_disposition: QCReviewDisposition | None = None

    # 2. Strict upstream QC gate
    if qc_result is not None:
        upstream_qc_status_str = str(qc_result.status)
        if qc_result.status == QCStatus.FAIL:
            findings.append(
                Finding(
                    Severity.ERROR,
                    "qc_status_failure",
                    "asset",
                    asset_id,
                    "qc_result.status",
                    f"Input asset failed bulk raw-count QC with status '{qc_result.status}'; downstream harmonization cannot proceed.",
                )
            )
            empty_summary = HarmonizationSummary(0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
            return HarmonizedDataset(
                status=Status.FAIL,
                findings=ordered(findings),
                summary=empty_summary,
                genes=(),
                collisions=(),
                matrix=None,
                source_asset_id=asset_id,
                assay_id=assay_id,
                sample_ids=(),
                reference_identity=config["reference"]["reference_id"],
                reference_version=config["reference"]["reference_version"],
                reference_checksum="",
                config_version=config["config_version"],
                config_checksum=cfg_checksum,
                harmonization_version=HARMONIZATION_VERSION,
                upstream_qc_status=upstream_qc_status_str,
                qc_review_disposition=None,
            )
        elif qc_result.status == QCStatus.NEEDS_REVIEW:
            review_disposition = QCReviewDisposition.from_any(
                qc_disposition,
                asset_id=asset_id,
                qc_status=upstream_qc_status_str,
            )
            if review_disposition is None:
                findings.append(
                    Finding(
                        Severity.ERROR,
                        "unauthorized_upstream_qc_needs_review",
                        "asset",
                        asset_id,
                        "qc_result.status",
                        "Input asset has upstream QC status 'NEEDS_REVIEW' without an explicit, auditable QCReviewDisposition (decision and scientific rationale required). A bare Boolean is insufficient for scientific authorization.",
                    )
                )
                empty_summary = HarmonizationSummary(0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
                return HarmonizedDataset(
                    status=Status.FAIL,
                    findings=ordered(findings),
                    summary=empty_summary,
                    genes=(),
                    collisions=(),
                    matrix=None,
                    source_asset_id=asset_id,
                    assay_id=assay_id,
                    sample_ids=(),
                    reference_identity=config["reference"]["reference_id"],
                    reference_version=config["reference"]["reference_version"],
                    reference_checksum="",
                    config_version=config["config_version"],
                    config_checksum=cfg_checksum,
                    harmonization_version=HARMONIZATION_VERSION,
                    upstream_qc_status=upstream_qc_status_str,
                    qc_review_disposition=None,
                )
            else:
                findings.append(
                    Finding(
                        Severity.REVIEW,
                        "upstream_qc_needs_review_authorized",
                        "asset",
                        asset_id,
                        "qc_result.status",
                        f"Input asset has upstream QC status 'NEEDS_REVIEW' and was explicitly authorized by reviewed disposition (decision: '{review_disposition.decision}', reason: '{review_disposition.reason}').",
                    )
                )
        elif qc_result.status == QCStatus.PASS_WITH_WARNINGS:
            findings.append(
                Finding(
                    Severity.WARNING,
                    "upstream_qc_pass_with_warnings",
                    "asset",
                    asset_id,
                    "qc_result.status",
                    "Input asset passed upstream QC with warnings.",
                )
            )

    # 3. Load or obtain gene reference
    if reference is None:
        project_root = Path(config_path).resolve().parents[1] if config_path else None
        reference = GeneReference.load_from_config(config, project_root=project_root)

    # 4. Resolve count matrix
    count_matrix: CountMatrix | None = None
    if isinstance(matrix, CountMatrix):
        count_matrix = matrix
    else:
        count_matrix, matrix_findings = read_matrix(matrix, {})
        for mf in matrix_findings:
            if mf.severity.value == "ERROR":
                findings.append(
                    Finding(
                        Severity.ERROR,
                        mf.rule_id,
                        mf.entity_type,
                        mf.entity_id,
                        mf.path,
                        mf.message,
                    )
                )

    if count_matrix is None or not count_matrix.gene_ids:
        findings.append(
            Finding(
                Severity.ERROR,
                "empty_or_unreadable_matrix",
                "asset",
                asset_id,
                "matrix",
                "Expression matrix could not be read or contains no genes.",
            )
        )
        empty_summary = HarmonizationSummary(0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)
        return HarmonizedDataset(
            status=Status.FAIL,
            findings=ordered(findings),
            summary=empty_summary,
            genes=(),
            collisions=(),
            matrix=None,
            source_asset_id=asset_id,
            assay_id=assay_id,
            sample_ids=(),
            reference_identity=reference.reference_id,
            reference_version=reference.reference_version,
            reference_checksum=reference.checksum,
            config_version=config["config_version"],
            config_checksum=cfg_checksum,
            harmonization_version=HARMONIZATION_VERSION,
            upstream_qc_status=upstream_qc_status_str,
        )

    # 5. Extract metadata declarations
    asset_record = next(
        (a for a in metadata.get("input_data_assets", []) if isinstance(a, dict) and a.get("asset_id") == asset_id),
        None,
    )
    declared_type_raw = reported_value(asset_record.get("gene_identifier_type")) if asset_record else None

    # 6. Resolve identifier type
    id_type, type_findings = resolve_identifier_type(
        count_matrix.gene_ids,
        declared_type_raw,
        config["supported_identifier_types"],
        asset_id,
    )
    findings.extend(type_findings)

    # If fatal error resolving identifier type
    if any(f.severity == Severity.ERROR for f in findings):
        empty_summary = HarmonizationSummary(
            total_genes=len(count_matrix.gene_ids),
            uniquely_mapped_genes=0,
            unmapped_genes=0,
            ambiguous_genes=0,
            invalid_genes=len(count_matrix.gene_ids),
            canonical_collision_genes=0,
            colliding_source_gene_count=0,
            percentage_uniquely_mapped=0.0,
            percentage_unresolved=100.0,
        )
        return HarmonizedDataset(
            status=Status.FAIL,
            findings=ordered(findings),
            summary=empty_summary,
            genes=(),
            collisions=(),
            matrix=count_matrix,
            source_asset_id=asset_id,
            assay_id=assay_id,
            sample_ids=count_matrix.sample_ids,
            reference_identity=reference.reference_id,
            reference_version=reference.reference_version,
            reference_checksum=reference.checksum,
            config_version=config["config_version"],
            config_checksum=cfg_checksum,
            harmonization_version=HARMONIZATION_VERSION,
            upstream_qc_status=upstream_qc_status_str,
        )

    # 7. Map each gene row
    unmapped_severity = Severity(structural_rules.get("unmapped_gene_identifier", "INFO"))
    ambiguous_severity = Severity(structural_rules.get("ambiguous_gene_identifier", "REVIEW"))

    harmonized_genes: list[HarmonizedGene] = []
    for idx, raw_gid in enumerate(count_matrix.gene_ids):
        stripped = raw_gid.strip()
        if not stripped:
            findings.append(
                Finding(
                    Severity.ERROR,
                    "empty_gene_identifiers",
                    "gene",
                    None,
                    f"matrix.gene_ids[{idx}]",
                    "Gene identifier is empty or whitespace.",
                )
            )
            harmonized_genes.append(
                HarmonizedGene(
                    source_index=idx,
                    original_gene_id=raw_gid,
                    identifier_type="unknown",
                    canonical_gene_id=None,
                    approved_symbol=None,
                    mapping_status=GeneMappingStatus.INVALID,
                    mapping_reason="Identifier is empty or whitespace.",
                    candidate_canonical_ids=(),
                    reference_id=reference.reference_id,
                    reference_version=reference.reference_version,
                    source_asset_id=asset_id,
                )
            )
            continue

        if id_type is None:
            harmonized_genes.append(
                HarmonizedGene(
                    source_index=idx,
                    original_gene_id=raw_gid,
                    identifier_type="ambiguous",
                    canonical_gene_id=None,
                    approved_symbol=None,
                    mapping_status=GeneMappingStatus.AMBIGUOUS,
                    mapping_reason="Gene identifier type could not be determined unambiguously without guessing.",
                    candidate_canonical_ids=(),
                    reference_id=reference.reference_id,
                    reference_version=reference.reference_version,
                    source_asset_id=asset_id,
                )
            )
            continue

        lookup_res = reference.lookup(raw_gid, id_type)

        if lookup_res.status == GeneMappingStatus.UNMAPPED:
            findings.append(
                Finding(
                    unmapped_severity,
                    "unmapped_gene_identifier",
                    "gene",
                    raw_gid,
                    f"matrix.gene_ids[{idx}]",
                    lookup_res.reason or f"Identifier '{raw_gid}' is unmapped.",
                )
            )
        elif lookup_res.status == GeneMappingStatus.AMBIGUOUS:
            findings.append(
                Finding(
                    ambiguous_severity,
                    "ambiguous_gene_identifier",
                    "gene",
                    raw_gid,
                    f"matrix.gene_ids[{idx}]",
                    lookup_res.reason or f"Identifier '{raw_gid}' is ambiguous.",
                )
            )
        elif lookup_res.status == GeneMappingStatus.INVALID:
            findings.append(
                Finding(
                    Severity.ERROR,
                    "invalid_gene_identifier",
                    "gene",
                    raw_gid,
                    f"matrix.gene_ids[{idx}]",
                    lookup_res.reason or f"Identifier '{raw_gid}' is invalid.",
                )
            )

        harmonized_genes.append(
            HarmonizedGene(
                source_index=idx,
                original_gene_id=raw_gid,
                identifier_type=id_type.value,
                canonical_gene_id=lookup_res.canonical_gene_id,
                approved_symbol=lookup_res.approved_symbol,
                mapping_status=lookup_res.status,
                mapping_reason=lookup_res.reason,
                candidate_canonical_ids=lookup_res.candidate_canonical_ids,
                reference_id=reference.reference_id,
                reference_version=reference.reference_version,
                source_asset_id=asset_id,
            )
        )

    # 8. Check for duplicate source identifiers in input matrix
    source_counts = Counter(count_matrix.gene_ids)
    for orig_id, count in sorted(source_counts.items(), key=lambda x: x[0]):
        if count > 1:
            findings.append(
                Finding(
                    Severity.REVIEW,
                    "duplicate_source_gene_identifiers",
                    "gene",
                    orig_id,
                    "matrix.gene_ids",
                    f"Duplicate source gene identifier '{orig_id}' appears {count} times in input matrix; no aggregation performed.",
                )
            )

    # 9. Detect canonical gene collisions (multiple rows mapping to the same canonical Ensembl ID)
    canonical_groups: dict[str, list[HarmonizedGene]] = defaultdict(list)
    for g in harmonized_genes:
        if g.canonical_gene_id:
            canonical_groups[g.canonical_gene_id].append(g)

    collisions: list[CanonicalCollision] = []
    for canon_id, group in sorted(canonical_groups.items(), key=lambda x: x[0]):
        if len(group) > 1:
            orig_ids = tuple(g.original_gene_id for g in group)
            src_indices = tuple(g.source_index for g in group)
            collisions.append(
                CanonicalCollision(
                    canonical_gene_id=canon_id,
                    original_gene_ids=orig_ids,
                    source_indices=src_indices,
                    approved_symbol=group[0].approved_symbol,
                )
            )
            findings.append(
                Finding(
                    Severity.REVIEW,
                    "canonical_gene_collision",
                    "gene",
                    canon_id,
                    "canonical_gene_id",
                    f"Multiple source identifiers ({', '.join(orig_ids)}) map to canonical gene '{canon_id}'; no aggregation performed.",
                )
            )

    # 10. Compute dataset summary and coverage metrics
    total_genes = len(harmonized_genes)
    uniquely_mapped = sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.UNIQUELY_MAPPED)
    unmapped = sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.UNMAPPED)
    ambiguous = sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.AMBIGUOUS)
    invalid = sum(1 for g in harmonized_genes if g.mapping_status == GeneMappingStatus.INVALID)
    collision_genes = len(collisions)
    colliding_source_count = sum(len(c.source_indices) for c in collisions)

    pct_uniquely_mapped = round((uniquely_mapped / total_genes * 100.0), 4) if total_genes > 0 else 0.0
    pct_unresolved = round(((unmapped + ambiguous + invalid) / total_genes * 100.0), 4) if total_genes > 0 else 0.0

    # 11. Evaluate mapping coverage threshold if configured
    min_cov = coverage_thresholds.get("minimum_mapping_coverage_percentage")
    if min_cov is not None and pct_uniquely_mapped < float(min_cov):
        findings.append(
            Finding(
                Severity(structural_rules.get("low_mapping_coverage", "REVIEW")),
                "low_mapping_coverage",
                "asset",
                asset_id,
                "summary.percentage_uniquely_mapped",
                f"Uniquely mapped coverage ({pct_uniquely_mapped:.2f}%) is below configured threshold ({min_cov:.2f}%).",
            )
        )

    summary = HarmonizationSummary(
        total_genes=total_genes,
        uniquely_mapped_genes=uniquely_mapped,
        unmapped_genes=unmapped,
        ambiguous_genes=ambiguous,
        invalid_genes=invalid,
        canonical_collision_genes=collision_genes,
        colliding_source_gene_count=colliding_source_count,
        percentage_uniquely_mapped=pct_uniquely_mapped,
        percentage_unresolved=pct_unresolved,
    )

    # 12. Final findings ordering and dataset status
    ordered_findings = ordered(findings)
    final_status = status_for(ordered_findings)

    return HarmonizedDataset(
        status=final_status,
        findings=ordered_findings,
        summary=summary,
        genes=tuple(harmonized_genes),
        collisions=tuple(collisions),
        matrix=count_matrix,
        source_asset_id=asset_id,
        assay_id=assay_id,
        sample_ids=count_matrix.sample_ids,
        reference_identity=reference.reference_id,
        reference_version=reference.reference_version,
        reference_checksum=reference.checksum,
        config_version=config["config_version"],
        config_checksum=cfg_checksum,
        harmonization_version=HARMONIZATION_VERSION,
        upstream_qc_status=upstream_qc_status_str,
        qc_review_disposition=review_disposition,
    )
