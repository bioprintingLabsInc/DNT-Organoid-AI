"""Governed gene ranking builder and contract implementation for Stage 8C preranked GSEA."""

from dataclasses import dataclass
import math
from typing import Sequence

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)

from .mapping import EnsemblToEntrezAdapter
from .models import (
    MappingStatus,
    PathwayGeneMappingResult,
    RankedGeneFeature,
    RankedGeneSet,
    RetainedNonEntrezGeneFeature,
    StatisticallyEligibleGeneFeature,
)


@dataclass(frozen=True, slots=True)
class ExcludedRankingGene:
    """Record of a gene excluded from the final preranked gene list."""

    canonical_gene_id: str
    approved_symbol: str | None
    entrez_gene_id: str | None
    result_status: str
    wald_statistic: float | None
    exclusion_reason: str


class RankingResult(tuple):
    """Immutable outcome of ranking construction.

    Unpacks as a 3-tuple (ranked_gene_set, mapping_result, excluded_records) for
    strict backwards-compatibility with existing callers, while exposing
    the full statistically eligible universe and retained non-Entrez genes.
    """

    ranked_gene_set: RankedGeneSet
    mapping_result: PathwayGeneMappingResult
    excluded_records: tuple[ExcludedRankingGene, ...]
    statistically_eligible_universe: tuple[StatisticallyEligibleGeneFeature, ...]
    retained_non_entrez_genes: tuple[RetainedNonEntrezGeneFeature, ...]

    def __new__(
        cls,
        ranked_gene_set: RankedGeneSet,
        mapping_result: PathwayGeneMappingResult,
        excluded_records: tuple[ExcludedRankingGene, ...],
        statistically_eligible_universe: tuple[StatisticallyEligibleGeneFeature, ...],
        retained_non_entrez_genes: tuple[RetainedNonEntrezGeneFeature, ...],
    ):
        instance = super().__new__(
            cls, (ranked_gene_set, mapping_result, excluded_records)
        )
        instance.ranked_gene_set = ranked_gene_set
        instance.mapping_result = mapping_result
        instance.excluded_records = excluded_records
        instance.statistically_eligible_universe = statistically_eligible_universe
        instance.retained_non_entrez_genes = retained_non_entrez_genes
        return instance


class GeneRankingBuilder:
    """Builds deterministically ordered, statistically qualified gene rankings for preranked GSEA."""

    DUPLICATE_RULE_MAX_ABS_WALD = "MAX_ABS_WALD: select representative with highest absolute Wald statistic (|wald_statistic|), tie-broken by canonical_gene_id"

    def __init__(
        self,
        mapping_adapter: EnsemblToEntrezAdapter | None = None,
        ranking_metric_name: str = "wald_statistic",
    ) -> None:
        self.mapping_adapter = mapping_adapter or EnsemblToEntrezAdapter()
        self.ranking_metric_name = ranking_metric_name

    def build_ranking(
        self,
        contrast: DifferentialExpressionContrastResult,
        deduplicate_entrez: bool = True,
    ) -> RankingResult:
        """Construct ranked gene set and preserve the full statistically eligible universe."""
        gene_results = contrast.gene_results
        mapping_result = self.mapping_adapter.map_genes(gene_results)

        # Build lookup from canonical_gene_id to mapping record
        mapping_by_id = {m.canonical_gene_id: m for m in mapping_result.mappings}

        eligible_features: list[RankedGeneFeature] = []
        excluded_records: list[ExcludedRankingGene] = []
        stat_eligible_universe: list[StatisticallyEligibleGeneFeature] = []
        retained_non_entrez: list[RetainedNonEntrezGeneFeature] = []

        for gene in gene_results:
            map_rec = mapping_by_id.get(gene.canonical_gene_id)
            entrez_id = map_rec.entrez_gene_id if map_rec else None
            biotype = map_rec.gene_biotype if map_rec else None
            symbol = map_rec.approved_symbol if map_rec else gene.approved_symbol

            # Check validity of Wald statistic
            wald = gene.wald_statistic
            is_valid_wald = (
                wald is not None
                and not math.isnan(wald)
                and not math.isinf(wald)
            )

            if not is_valid_wald:
                excluded_records.append(
                    ExcludedRankingGene(
                        canonical_gene_id=gene.canonical_gene_id,
                        approved_symbol=gene.approved_symbol,
                        entrez_gene_id=entrez_id,
                        result_status=gene.result_status,
                        wald_statistic=wald,
                        exclusion_reason=f"Invalid wald_statistic ({wald}) with status {gene.result_status}",
                    )
                )
                continue

            if gene.result_status != "OK":
                excluded_records.append(
                    ExcludedRankingGene(
                        canonical_gene_id=gene.canonical_gene_id,
                        approved_symbol=gene.approved_symbol,
                        entrez_gene_id=entrez_id,
                        result_status=gene.result_status,
                        wald_statistic=wald,
                        exclusion_reason=f"Non-OK Stage 7 result_status: {gene.result_status}",
                    )
                )
                continue

            # Gene is statistically eligible!
            stat_feat = StatisticallyEligibleGeneFeature(
                canonical_gene_id=gene.canonical_gene_id,
                approved_symbol=symbol,
                gene_biotype=biotype,
                entrez_gene_id=entrez_id,
                base_mean=gene.base_mean,
                log2_fold_change=gene.log2_fold_change,
                lfc_standard_error=gene.lfc_standard_error,
                wald_statistic=wald,
                p_value=gene.p_value,
                adjusted_p_value_bh=gene.adjusted_p_value_bh,
                result_status=gene.result_status,
                annotation_status="ENTREZ_MAPPED" if entrez_id else "ENTREZ_UNMAPPED",
                statistical_eligibility=True,
            )
            stat_eligible_universe.append(stat_feat)

            if not entrez_id:
                # Gene is statistically valid but cannot enter Entrez-based pathway GSEA.
                # Must NOT be classified as QC failure; retained for feature construction.
                retained_feat = RetainedNonEntrezGeneFeature(
                    canonical_gene_id=gene.canonical_gene_id,
                    approved_symbol=symbol,
                    gene_biotype=biotype,
                    wald_statistic=wald,
                    log2_fold_change=gene.log2_fold_change,
                    lfc_standard_error=gene.lfc_standard_error,
                    p_value=gene.p_value,
                    adjusted_p_value_bh=gene.adjusted_p_value_bh,
                    reason_not_in_kegg="NO_ENTREZ_ID_IN_ENSEMBL_112",
                    go_annotation_status="NOT_ANNOTATED_IN_ENTREZ_GO",
                    future_feature_eligibility=True,
                )
                retained_non_entrez.append(retained_feat)
                excluded_records.append(
                    ExcludedRankingGene(
                        canonical_gene_id=gene.canonical_gene_id,
                        approved_symbol=gene.approved_symbol,
                        entrez_gene_id=None,
                        result_status=gene.result_status,
                        wald_statistic=wald,
                        exclusion_reason="Unmapped to Entrez Gene ID in Ensembl 112 reference (retained in non-Entrez feature universe)",
                    )
                )
                continue

            feat = RankedGeneFeature(
                contrast_id=contrast.contrast_id,
                study_id=contrast.study_id,
                canonical_gene_id=gene.canonical_gene_id,
                approved_symbol=symbol,
                entrez_gene_id=entrez_id,
                log2_fold_change=gene.log2_fold_change,
                lfc_standard_error=gene.lfc_standard_error,
                wald_statistic=wald,
                p_value=gene.p_value,
                adjusted_p_value_bh=gene.adjusted_p_value_bh,
                result_status=gene.result_status,
                ranking_metric_name=self.ranking_metric_name,
                ranking_value=wald,
                mapping_status=map_rec.mapping_status if map_rec else MappingStatus.MAPPED.value,
            )
            eligible_features.append(feat)

        # Handle duplicate Entrez mappings if requested
        final_features: list[RankedGeneFeature]
        if deduplicate_entrez:
            by_entrez: dict[str, list[RankedGeneFeature]] = {}
            for f in eligible_features:
                by_entrez.setdefault(f.entrez_gene_id, []).append(f)

            final_features = []
            for e_id, feats in by_entrez.items():
                if len(feats) == 1:
                    final_features.append(feats[0])
                else:
                    # Deterministic selection: max abs(wald), tie-breaker: positive wald, then canonical_id
                    sorted_dups = sorted(
                        feats,
                        key=lambda x: (abs(x.ranking_value), x.ranking_value, -1 if x.canonical_gene_id else 0),
                        reverse=True,
                    )
                    chosen = sorted_dups[0]
                    final_features.append(chosen)
                    for dropped in sorted_dups[1:]:
                        excluded_records.append(
                            ExcludedRankingGene(
                                canonical_gene_id=dropped.canonical_gene_id,
                                approved_symbol=dropped.approved_symbol,
                                entrez_gene_id=dropped.entrez_gene_id,
                                result_status=dropped.result_status,
                                wald_statistic=dropped.wald_statistic,
                                exclusion_reason=f"Duplicate Entrez {e_id} resolved: smaller |wald| ({dropped.ranking_value}) than chosen {chosen.canonical_gene_id} ({chosen.ranking_value})",
                            )
                        )
        else:
            final_features = eligible_features

        # Deterministic sorting: ranking_value descending, then canonical_gene_id ascending as tie-breaker
        final_features.sort(
            key=lambda x: (-x.ranking_value, x.canonical_gene_id)
        )

        ranked_set = RankedGeneSet(
            contrast_id=contrast.contrast_id,
            study_id=contrast.study_id,
            ranking_metric_name=self.ranking_metric_name,
            ranked_genes=tuple(final_features),
            total_eligible_genes=len(final_features),
            excluded_genes_count=len(excluded_records),
            duplicate_resolution_rule=self.DUPLICATE_RULE_MAX_ABS_WALD if deduplicate_entrez else "NONE",
        )

        return RankingResult(
            ranked_gene_set=ranked_set,
            mapping_result=mapping_result,
            excluded_records=tuple(excluded_records),
            statistically_eligible_universe=tuple(stat_eligible_universe),
            retained_non_entrez_genes=tuple(retained_non_entrez),
        )
