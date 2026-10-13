"""Unit tests for Stage 8C QC reporting and provenance tracking."""

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.pathway.models import (
    PathwayEnrichmentResult,
    PathwayGeneMappingRecord,
    PathwayGeneMappingResult,
    RankedGeneFeature,
    RankedGeneSet,
    RetainedNonEntrezGeneFeature,
    StatisticallyEligibleGeneFeature,
)
from src.features.pathway.provenance import build_pathway_analysis_provenance
from src.features.pathway.qc import compute_pathway_analysis_qc
from src.features.pathway.ranking import ExcludedRankingGene


def test_qc_metrics_and_percentages_calculation() -> None:
    # 5 test genes:
    # Gene 1: Valid Wald, Entrez mapped (in ranking)
    # Gene 2: Valid Wald, unmapped to Entrez (retained non-Entrez)
    # Gene 3: Valid Wald, Entrez duplicate dropped
    # Gene 4: Zero counts
    # Gene 5: Outlier
    gene_results = (
        GeneDifferentialExpressionResult("ENSG01", "G1", 10.0, 1.0, 0.2, 5.0, 1e-5, 1e-4, "OK"),
        GeneDifferentialExpressionResult("ENSG02", "G2", 20.0, 0.5, 0.1, 2.5, 0.01, 0.05, "OK"),
        GeneDifferentialExpressionResult("ENSG03", "G3", 15.0, 0.8, 0.2, 4.0, 1e-4, 1e-3, "OK"),
        GeneDifferentialExpressionResult("ENSG04", "G4", 0.0, 0.0, None, None, None, None, "ZERO_COUNTS"),
        GeneDifferentialExpressionResult("ENSG05", "G5", 100.0, 2.0, 0.3, 6.0, 1e-6, 1e-5, "OUTLIER"),
    )

    contrast = DifferentialExpressionContrastResult(
        contrast_id="c1",
        study_id="s1",
        experiment_id="e1",
        cohort_id="c1",
        relationship_id="r1",
        treatment_condition_id="t",
        matched_control_condition_ids=("c",),
        treatment_sample_ids=("s1",),
        control_sample_ids=("s2",),
        treatment_biological_replicate_ids=("r1",),
        control_biological_replicate_ids=("r2",),
        treatment_replicate_count=1,
        control_replicate_count=1,
        agent_name="compound",
        agent_identifier=None,
        concentration_or_dose=1.0,
        concentration_or_dose_unit="uM",
        exposure_duration=1.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="D1",
        design_formula="~ condition",
        residual_degrees_of_freedom=1,
        size_factors_used={},
        total_input_genes_count=5,
        gene_results=gene_results,
    )

    stat_universe = (
        StatisticallyEligibleGeneFeature("ENSG01", "G1", "protein_coding", "101", 10.0, 1.0, 0.2, 5.0, 1e-5, 1e-4, "OK", "ENTREZ_MAPPED"),
        StatisticallyEligibleGeneFeature("ENSG02", "G2", "lncRNA", None, 20.0, 0.5, 0.1, 2.5, 0.01, 0.05, "OK", "ENTREZ_UNMAPPED"),
        StatisticallyEligibleGeneFeature("ENSG03", "G3", "protein_coding", "101", 15.0, 0.8, 0.2, 4.0, 1e-4, 1e-3, "OK", "ENTREZ_MAPPED"),
    )
    retained_non_entrez = (
        RetainedNonEntrezGeneFeature("ENSG02", "G2", "lncRNA", 2.5, 0.5, 0.1, 0.01, 0.05, "NO_ENTREZ_ID_IN_ENSEMBL_112", "NOT_ANNOTATED_IN_ENTREZ_GO"),
    )

    ranked_feats = (
        RankedGeneFeature("c1", "s1", "ENSG01", "G1", "101", 1.0, 0.2, 5.0, 1e-5, 1e-4, "OK", "wald_statistic", 5.0, "MAPPED"),
    )
    ranked_set = RankedGeneSet("c1", "s1", "wald_statistic", ranked_feats, 1, 4, "MAX_ABS_WALD")

    mapping_res = PathwayGeneMappingResult(
        mappings=(),
        total_input_genes=5,
        mapped_count=2,
        unmapped_count=3,
        duplicate_entrez_count=1,
        one_to_many_count=0,
        invalid_stat_count=2,
        reference_version="Ensembl_112",
        reference_checksum="chk1",
    )

    excluded_recs = (
        ExcludedRankingGene("ENSG04", "G4", None, "ZERO_COUNTS", None, "Invalid wald_statistic"),
        ExcludedRankingGene("ENSG05", "G5", None, "OUTLIER", 6.0, "Non-OK Stage 7 result_status: OUTLIER"),
        ExcludedRankingGene("ENSG02", "G2", None, "OK", 2.5, "Unmapped to Entrez Gene ID in Ensembl 112 reference"),
        ExcludedRankingGene("ENSG03", "G3", "101", "OK", 4.0, "Duplicate Entrez 101 resolved: smaller |wald|"),
    )

    go_res = (
        PathwayEnrichmentResult("c1", "GO_BP", "GO:0001", "Pathway 1", 0.8, 2.1, 0.001, 0.01, 25, "POSITIVE", ("101",), ("G1",), "wald_statistic", "clusterProfiler_4.14.6", True),
        PathwayEnrichmentResult("c1", "GO_BP", "GO:0002", "Pathway 2", -0.4, -1.2, 0.15, 0.35, 40, "NEGATIVE", (), (), "wald_statistic", "clusterProfiler_4.14.6", False),
    )
    kegg_res = ()

    qc = compute_pathway_analysis_qc(
        stage7_gene_count=5,
        mapping_result=mapping_res,
        ranked_set=ranked_set,
        excluded_records=excluded_recs,
        go_results=go_res,
        kegg_results=kegg_res,
        contrast=contrast,
        statistically_eligible_universe=stat_universe,
        retained_non_entrez_genes=retained_non_entrez,
    )

    # Exact accounting verification
    assert qc.stage7_evaluated_genes == 5
    assert qc.zero_count_genes == 1
    assert qc.genes_with_valid_wald == 4
    assert qc.cooks_outlier_genes == 1
    assert qc.statistically_eligible_genes == 3
    assert qc.entrez_mapped_eligible_genes == 2
    assert qc.entrez_unmapped_eligible_genes == 1
    assert qc.duplicate_entrez_entries_removed == 1
    assert qc.final_kegg_ranked_genes == 1
    assert qc.final_go_ranked_genes == 1
    assert qc.go_gene_sets_tested == 2
    assert qc.go_significant_pathways_count == 1
    assert qc.go_significant_pathways_pct == 50.0
    assert qc.kegg_gene_sets_tested == 0
    assert qc.qc_status == "PASS"

    # Mathematical reconciliation assertions
    assert qc.stage7_evaluated_genes == qc.zero_count_genes + qc.genes_with_valid_wald
    assert qc.statistically_eligible_genes == qc.genes_with_valid_wald - qc.cooks_outlier_genes
    assert qc.statistically_eligible_genes == qc.entrez_mapped_eligible_genes + qc.entrez_unmapped_eligible_genes
    assert qc.final_kegg_ranked_genes == qc.entrez_mapped_eligible_genes - qc.duplicate_entrez_entries_removed


def test_provenance_preserves_experimental_metadata_and_hashes() -> None:
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

    # Provenance hardening assertions
    assert "baseMean > 0" in prov.statistical_eligibility_criteria
    assert "GO_BP" in prov.pathway_annotation_eligibility_criteria
    assert "RETAINED_NON_ENTREZ" in prov.pathway_annotation_eligibility_criteria
    assert prov.go_provenance["key_type"] == "ENTREZID"
    assert "org_hs_eg_sqlite_sha256" in prov.go_provenance
    assert prov.kegg_provenance["organism"] == "hsa"
    assert prov.kegg_provenance["pathway_count"] == 372
    assert prov.kegg_provenance["sha256_rds"] == "A54298B77CE24D0059FB58CE0B3E22EB6837589D02FA4B2437E7545A02967C20"
    assert prov.gsea_parameters["gsea_param"] == 1.0
    assert prov.gsea_parameters["score_type"] == "std"
    assert prov.gsea_parameters["nperm_simple"] == 1000
