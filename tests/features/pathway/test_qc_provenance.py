"""Unit tests for Stage 8C QC reporting and provenance tracking."""

from src.differential_expression.models import DifferentialExpressionContrastResult
from src.features.pathway.models import (
    PathwayEnrichmentResult,
    PathwayGeneMappingRecord,
    PathwayGeneMappingResult,
    RankedGeneFeature,
    RankedGeneSet,
)
from src.features.pathway.provenance import build_pathway_analysis_provenance
from src.features.pathway.qc import compute_pathway_analysis_qc


def test_qc_metrics_and_percentages_calculation() -> None:
    mapping_recs = (
        PathwayGeneMappingRecord("ENSG01", "G1", "101", "MAPPED", "v112"),
        PathwayGeneMappingRecord("ENSG02", "G2", None, "UNMAPPED", "v112"),
        PathwayGeneMappingRecord("ENSG03", "G3", "103", "INVALID_STAT", "v112"),
    )
    mapping_res = PathwayGeneMappingResult(
        mappings=mapping_recs,
        total_input_genes=3,
        mapped_count=2,
        unmapped_count=1,
        duplicate_entrez_count=0,
        one_to_many_count=0,
        invalid_stat_count=1,
        reference_version="Ensembl_112",
        reference_checksum="chk1",
    )
    ranked_feats = (
        RankedGeneFeature(
            contrast_id="c1",
            study_id="s1",
            canonical_gene_id="ENSG01",
            approved_symbol="G1",
            entrez_gene_id="101",
            log2_fold_change=1.0,
            lfc_standard_error=0.2,
            wald_statistic=5.0,
            p_value=1e-5,
            adjusted_p_value_bh=1e-4,
            result_status="OK",
            ranking_metric_name="wald_statistic",
            ranking_value=5.0,
            mapping_status="MAPPED",
        ),
    )
    ranked_set = RankedGeneSet(
        contrast_id="c1",
        study_id="s1",
        ranking_metric_name="wald_statistic",
        ranked_genes=ranked_feats,
        total_eligible_genes=1,
        excluded_genes_count=2,
        duplicate_resolution_rule="MAX_ABS_WALD",
    )
    go_res = (
        PathwayEnrichmentResult(
            contrast_id="c1",
            source_database="GO_BP",
            pathway_id="GO:0001",
            pathway_name="Pathway 1",
            enrichment_score=0.8,
            normalized_enrichment_score=2.1,
            p_value=0.001,
            adjusted_p_value=0.01,
            pathway_size=25,
            direction="POSITIVE",
            leading_edge_gene_ids=("101",),
            leading_edge_symbols=("G1",),
            ranking_metric="wald_statistic",
            analysis_version="clusterProfiler_4.14.6",
            is_significant=True,
        ),
        PathwayEnrichmentResult(
            contrast_id="c1",
            source_database="GO_BP",
            pathway_id="GO:0002",
            pathway_name="Pathway 2",
            enrichment_score=-0.4,
            normalized_enrichment_score=-1.2,
            p_value=0.15,
            adjusted_p_value=0.35,
            pathway_size=40,
            direction="NEGATIVE",
            leading_edge_gene_ids=(),
            leading_edge_symbols=(),
            ranking_metric="wald_statistic",
            analysis_version="clusterProfiler_4.14.6",
            is_significant=False,
        ),
    )
    kegg_res = ()

    qc = compute_pathway_analysis_qc(
        stage7_gene_count=3,
        mapping_result=mapping_res,
        ranked_set=ranked_set,
        excluded_records=(),
        go_results=go_res,
        kegg_results=kegg_res,
    )

    assert qc.stage7_evaluated_genes == 3
    assert qc.valid_ranking_statistic_genes == 2
    assert qc.valid_ranking_statistic_pct == 66.67
    assert qc.invalid_statistical_result_genes == 1
    assert qc.invalid_statistical_result_pct == 33.33
    assert qc.mapped_entrez_genes == 2
    assert qc.mapped_entrez_pct == 66.67
    assert qc.unmapped_ensembl_genes == 1
    assert qc.unmapped_ensembl_pct == 33.33
    assert qc.genes_entering_go_gsea == 1
    assert qc.go_gene_sets_tested == 2
    assert qc.go_significant_pathways_count == 1
    assert qc.go_significant_pathways_pct == 50.0
    assert qc.kegg_gene_sets_tested == 0


def test_provenance_preserves_experimental_metadata() -> None:
    contrast = DifferentialExpressionContrastResult(
        contrast_id="c_vpa_01",
        study_id="GSE232218",
        experiment_id="EXP_VPA",
        cohort_id="COHORT_01",
        relationship_id="REL_01",
        treatment_condition_id="TRT_VPA_200",
        matched_control_condition_ids=("CTL_VEH",),
        treatment_sample_ids=("S_TRT1", "S_TRT2"),
        control_sample_ids=("S_CTL1", "S_CTL2"),
        treatment_biological_replicate_ids=("R1", "R2"),
        control_biological_replicate_ids=("R3", "R4"),
        treatment_replicate_count=2,
        control_replicate_count=2,
        agent_name="valproic acid",
        agent_identifier="CID:3121",
        concentration_or_dose=200.0,
        concentration_or_dose_unit="uM",
        exposure_duration=48.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 6",
        design_formula="~ condition",
        residual_degrees_of_freedom=2,
        size_factors_used={},
        total_input_genes_count=2,
        gene_results=(),
    )

    prov = build_pathway_analysis_provenance(
        contrast=contrast,
        ranking_metric="wald_statistic",
        gene_reference_version="Ensembl_112",
        pathway_databases=("GO_BP", "KEGG"),
    )

    assert prov.study_id == "GSE232218"
    assert prov.contrast_id == "c_vpa_01"
    assert prov.agent_name == "valproic acid"
    assert prov.agent_identifier == "CID:3121"
    assert prov.concentration == 200.0
    assert prov.concentration_unit == "uM"
    assert prov.exposure_duration == 48.0
    assert prov.duration_unit == "h"
    assert prov.developmental_stage == "Day 6"
    assert prov.treatment_sample_ids == ("S_TRT1", "S_TRT2")
    assert prov.control_sample_ids == ("S_CTL1", "S_CTL2")
    assert prov.ranking_metric == "wald_statistic"
    assert "clusterProfiler" in prov.software_packages
    assert prov.stage7_source_commit == "ca26d34e6c0263e743425cc7d6752d0e5773a9fe"
