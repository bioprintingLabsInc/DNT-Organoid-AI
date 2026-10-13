"""Data models and immutable result contracts for Stage 8C Pathway and Signaling Analysis."""

from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any


class MappingStatus(str, Enum):
    """Mapping status of a canonical gene to an Entrez Gene ID."""

    MAPPED = "MAPPED"
    UNMAPPED = "UNMAPPED"
    DUPLICATE_ENTREZ = "DUPLICATE_ENTREZ"
    ONE_TO_MANY = "ONE_TO_MANY"
    INVALID_STAT = "INVALID_STAT"


@dataclass(frozen=True, slots=True)
class PathwayGeneMappingRecord:
    """Mapping audit record linking an Ensembl canonical gene to its Entrez ID."""

    canonical_gene_id: str
    approved_symbol: str | None
    entrez_gene_id: str | None
    mapping_status: str
    source_reference_version: str
    gene_biotype: str | None = None
    notes: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PathwayGeneMappingResult:
    """Full mapping audit outcome across all evaluated genes in a contrast."""

    mappings: tuple[PathwayGeneMappingRecord, ...]
    total_input_genes: int
    mapped_count: int
    unmapped_count: int
    duplicate_entrez_count: int
    one_to_many_count: int
    invalid_stat_count: int
    reference_version: str
    reference_checksum: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "total_input_genes": self.total_input_genes,
            "mapped_count": self.mapped_count,
            "unmapped_count": self.unmapped_count,
            "duplicate_entrez_count": self.duplicate_entrez_count,
            "one_to_many_count": self.one_to_many_count,
            "invalid_stat_count": self.invalid_stat_count,
            "reference_version": self.reference_version,
            "reference_checksum": self.reference_checksum,
            "mappings": [m.to_dict() for m in self.mappings],
        }


@dataclass(frozen=True, slots=True)
class StatisticallyEligibleGeneFeature:
    """A statistically valid Stage 7 gene in the preserved molecular universe."""

    canonical_gene_id: str
    approved_symbol: str | None
    gene_biotype: str | None
    entrez_gene_id: str | None
    base_mean: float
    log2_fold_change: float | None
    lfc_standard_error: float | None
    wald_statistic: float | None
    p_value: float | None
    adjusted_p_value_bh: float | None
    result_status: str
    annotation_status: str  # "ENTREZ_MAPPED" or "ENTREZ_UNMAPPED"
    statistical_eligibility: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class RetainedNonEntrezGeneFeature:
    """A statistically valid gene retained for feature construction despite lacking Entrez ID."""

    canonical_gene_id: str
    approved_symbol: str | None
    gene_biotype: str | None
    wald_statistic: float | None
    log2_fold_change: float | None
    lfc_standard_error: float | None
    p_value: float | None
    adjusted_p_value_bh: float | None
    reason_not_in_kegg: str  # "NO_ENTREZ_ID_IN_ENSEMBL_112"
    go_annotation_status: str  # "NOT_ANNOTATED_IN_ENTREZ_GO"
    future_feature_eligibility: bool = True

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class RankedGeneFeature:
    """Single gene feature evaluated and ordered for preranked GSEA."""

    contrast_id: str
    study_id: str
    canonical_gene_id: str
    approved_symbol: str | None
    entrez_gene_id: str | None
    log2_fold_change: float | None
    lfc_standard_error: float | None
    wald_statistic: float | None
    p_value: float | None
    adjusted_p_value_bh: float | None
    result_status: str
    ranking_metric_name: str
    ranking_value: float
    mapping_status: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class RankedGeneSet:
    """Full sorted collection of gene features entering preranked GSEA."""

    contrast_id: str
    study_id: str
    ranking_metric_name: str
    ranked_genes: tuple[RankedGeneFeature, ...]
    total_eligible_genes: int
    excluded_genes_count: int
    duplicate_resolution_rule: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "contrast_id": self.contrast_id,
            "study_id": self.study_id,
            "ranking_metric_name": self.ranking_metric_name,
            "total_eligible_genes": self.total_eligible_genes,
            "excluded_genes_count": self.excluded_genes_count,
            "duplicate_resolution_rule": self.duplicate_resolution_rule,
            "ranked_genes": [g.to_dict() for g in self.ranked_genes],
        }


@dataclass(frozen=True, slots=True)
class PathwayEnrichmentResult:
    """Enrichment statistics for a single pathway/gene set."""

    contrast_id: str
    source_database: str  # "GO_BP" or "KEGG"
    pathway_id: str
    pathway_name: str
    enrichment_score: float
    normalized_enrichment_score: float
    p_value: float
    adjusted_p_value: float
    pathway_size: int
    direction: str  # "POSITIVE" or "NEGATIVE"
    leading_edge_gene_ids: tuple[str, ...]
    leading_edge_symbols: tuple[str, ...]
    ranking_metric: str
    analysis_version: str
    is_significant: bool  # adjusted_p_value < 0.05

    def to_dict(self) -> dict[str, Any]:
        return {
            "contrast_id": self.contrast_id,
            "source_database": self.source_database,
            "pathway_id": self.pathway_id,
            "pathway_name": self.pathway_name,
            "enrichment_score": self.enrichment_score,
            "normalized_enrichment_score": self.normalized_enrichment_score,
            "p_value": self.p_value,
            "adjusted_p_value": self.adjusted_p_value,
            "pathway_size": self.pathway_size,
            "direction": self.direction,
            "leading_edge_gene_ids": list(self.leading_edge_gene_ids),
            "leading_edge_symbols": list(self.leading_edge_symbols),
            "ranking_metric": self.ranking_metric,
            "analysis_version": self.analysis_version,
            "is_significant": self.is_significant,
        }


@dataclass(frozen=True, slots=True)
class PathwayAnalysisQC:
    """Comprehensive quality control and accounting report for Stage 8C."""

    stage7_evaluated_genes: int
    zero_count_genes: int
    genes_with_valid_wald: int
    cooks_outlier_genes: int
    statistically_eligible_genes: int
    entrez_mapped_eligible_genes: int
    entrez_unmapped_eligible_genes: int
    duplicate_entrez_entries_removed: int
    final_kegg_ranked_genes: int
    final_go_ranked_genes: int
    go_gene_sets_tested: int
    kegg_gene_sets_tested: int
    go_significant_pathways_count: int
    go_significant_pathways_pct: float
    kegg_significant_pathways_count: int
    kegg_significant_pathways_pct: float
    qc_status: str  # "PASS", "WARNING", "FAIL"
    qc_notes: tuple[str, ...]

    # Compatibility properties for Stage 8C v1 callers
    @property
    def valid_ranking_statistic_genes(self) -> int:
        return self.statistically_eligible_genes

    @property
    def valid_ranking_statistic_pct(self) -> float:
        return round((self.statistically_eligible_genes / self.stage7_evaluated_genes) * 100.0, 2) if self.stage7_evaluated_genes > 0 else 0.0

    @property
    def invalid_statistical_result_genes(self) -> int:
        return self.zero_count_genes + self.cooks_outlier_genes

    @property
    def invalid_statistical_result_pct(self) -> float:
        return round((self.invalid_statistical_result_genes / self.stage7_evaluated_genes) * 100.0, 2) if self.stage7_evaluated_genes > 0 else 0.0

    @property
    def mapped_entrez_genes(self) -> int:
        return self.entrez_mapped_eligible_genes

    @property
    def mapped_entrez_pct(self) -> float:
        return round((self.entrez_mapped_eligible_genes / self.stage7_evaluated_genes) * 100.0, 2) if self.stage7_evaluated_genes > 0 else 0.0

    @property
    def unmapped_ensembl_genes(self) -> int:
        return self.entrez_unmapped_eligible_genes

    @property
    def unmapped_ensembl_pct(self) -> float:
        return round((self.entrez_unmapped_eligible_genes / self.stage7_evaluated_genes) * 100.0, 2) if self.stage7_evaluated_genes > 0 else 0.0

    @property
    def one_to_many_mappings(self) -> int:
        return 0

    @property
    def duplicate_entrez_mappings(self) -> int:
        return self.duplicate_entrez_entries_removed

    @property
    def genes_entering_go_gsea(self) -> int:
        return self.final_go_ranked_genes

    @property
    def genes_entering_kegg_gsea(self) -> int:
        return self.final_kegg_ranked_genes

    def to_dict(self) -> dict[str, Any]:
        return {
            "stage7_evaluated_genes": self.stage7_evaluated_genes,
            "zero_count_genes": self.zero_count_genes,
            "genes_with_valid_wald": self.genes_with_valid_wald,
            "cooks_outlier_genes": self.cooks_outlier_genes,
            "statistically_eligible_genes": self.statistically_eligible_genes,
            "entrez_mapped_eligible_genes": self.entrez_mapped_eligible_genes,
            "entrez_unmapped_eligible_genes": self.entrez_unmapped_eligible_genes,
            "duplicate_entrez_entries_removed": self.duplicate_entrez_entries_removed,
            "final_kegg_ranked_genes": self.final_kegg_ranked_genes,
            "final_go_ranked_genes": self.final_go_ranked_genes,
            "go_gene_sets_tested": self.go_gene_sets_tested,
            "kegg_gene_sets_tested": self.kegg_gene_sets_tested,
            "go_significant_pathways_count": self.go_significant_pathways_count,
            "go_significant_pathways_pct": self.go_significant_pathways_pct,
            "kegg_significant_pathways_count": self.kegg_significant_pathways_count,
            "kegg_significant_pathways_pct": self.kegg_significant_pathways_pct,
            "qc_status": self.qc_status,
            "qc_notes": list(self.qc_notes),
            # Compatibility fields
            "valid_ranking_statistic_genes": self.valid_ranking_statistic_genes,
            "valid_ranking_statistic_pct": self.valid_ranking_statistic_pct,
            "invalid_statistical_result_genes": self.invalid_statistical_result_genes,
            "invalid_statistical_result_pct": self.invalid_statistical_result_pct,
            "mapped_entrez_genes": self.mapped_entrez_genes,
            "mapped_entrez_pct": self.mapped_entrez_pct,
            "unmapped_ensembl_genes": self.unmapped_ensembl_genes,
            "unmapped_ensembl_pct": self.unmapped_ensembl_pct,
            "one_to_many_mappings": self.one_to_many_mappings,
            "duplicate_entrez_mappings": self.duplicate_entrez_mappings,
            "genes_entering_go_gsea": self.genes_entering_go_gsea,
            "genes_entering_kegg_gsea": self.genes_entering_kegg_gsea,
        }


@dataclass(frozen=True, slots=True)
class PathwayAnalysisProvenance:
    """Full upstream and downstream metadata traceability record."""

    study_id: str
    experiment_id: str
    contrast_id: str
    agent_name: str
    agent_identifier: str | None
    concentration: float
    concentration_unit: str
    exposure_duration: float
    duration_unit: str
    developmental_stage: str
    treatment_sample_ids: tuple[str, ...]
    control_sample_ids: tuple[str, ...]
    gene_reference_version: str
    ranking_metric: str
    pathway_databases: tuple[str, ...]
    pathway_database_version: str | None
    software_packages: dict[str, str]
    stage7_source_commit: str
    stage8c_source_commit: str
    statistical_eligibility_criteria: str
    pathway_annotation_eligibility_criteria: dict[str, str]
    go_provenance: dict[str, Any]
    kegg_provenance: dict[str, Any]
    gsea_parameters: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PathwayEnrichmentDataset:
    """Governed container for all Stage 8C pathway enrichment artifacts."""

    contrast_id: str
    study_id: str
    go_results: tuple[PathwayEnrichmentResult, ...]
    kegg_results: tuple[PathwayEnrichmentResult, ...]
    ranked_genes: RankedGeneSet
    statistically_eligible_universe: tuple[StatisticallyEligibleGeneFeature, ...]
    retained_non_entrez_genes: tuple[RetainedNonEntrezGeneFeature, ...]
    qc: PathwayAnalysisQC
    provenance: PathwayAnalysisProvenance
    mapping_result: PathwayGeneMappingResult

    def to_dict(self) -> dict[str, Any]:
        return {
            "contrast_id": self.contrast_id,
            "study_id": self.study_id,
            "go_results": [r.to_dict() for r in self.go_results],
            "kegg_results": [r.to_dict() for r in self.kegg_results],
            "ranked_genes": self.ranked_genes.to_dict(),
            "statistically_eligible_universe": [g.to_dict() for g in self.statistically_eligible_universe],
            "retained_non_entrez_genes": [g.to_dict() for g in self.retained_non_entrez_genes],
            "qc": self.qc.to_dict(),
            "provenance": self.provenance.to_dict(),
            "mapping_result": self.mapping_result.to_dict(),
        }
