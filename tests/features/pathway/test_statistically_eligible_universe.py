"""Comprehensive tests for Stage 8C statistically eligible molecular universe preservation."""

from pathlib import Path
import tempfile
import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    GeneDifferentialExpressionResult,
)
from src.features.pathway.mapping import EnsemblToEntrezAdapter
from src.features.pathway.models import (
    RetainedNonEntrezGeneFeature,
    StatisticallyEligibleGeneFeature,
)
from src.features.pathway.pipeline import (
    export_retained_non_entrez_genes_to_csv,
    export_statistically_eligible_universe_to_csv,
)
from src.features.pathway.qc import compute_pathway_analysis_qc
from src.features.pathway.ranking import GeneRankingBuilder
from src.gene_harmonization.reference import GeneReference, ReferenceRecord


@pytest.fixture
def synthetic_multibiotype_contrast() -> tuple[DifferentialExpressionContrastResult, GeneReference]:
    """Fixture providing a mock contrast with protein-coding, lncRNAs, zero-counts, and outliers."""
    records = (
        ReferenceRecord("ENSG_PC01", "BRCA1", "672", "HGNC:1100", "protein_coding", (), "17", ""),
        ReferenceRecord("ENSG_PC02", "TP53", "7157", "HGNC:11998", "protein_coding", (), "17", ""),
        ReferenceRecord("ENSG_LNC01", "HOTAIR", None, "HGNC:33144", "lncRNA", (), "12", ""),
        ReferenceRecord("ENSG_LNC02", "NEAT1", None, "HGNC:33145", "lncRNA", (), "11", ""),
        ReferenceRecord("ENSG_TEC01", "TEC_GENE", None, "", "TEC", (), "1", ""),
        ReferenceRecord("ENSG_OUT01", "OUT_GENE", "9999", "HGNC:9999", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG_ZERO1", "ZERO_GENE", "8888", "HGNC:8888", "protein_coding", (), "1", ""),
        ReferenceRecord("ENSG_DUP01", "TP53_DUP", "7157", "HGNC:11999", "protein_coding", (), "17", ""),
    )
    ref = GeneReference("ref_test", "v112", "GRCh38", "Ensembl", records, "chk_test")

    gene_results = (
        # 1. BRCA1: Valid Wald, Entrez mapped -> GSEA eligible
        GeneDifferentialExpressionResult("ENSG_PC01", "BRCA1", 100.0, 1.2, 0.2, 6.0, 1e-8, 1e-7, "OK"),
        # 2. TP53: Valid Wald, Entrez mapped -> GSEA eligible (higher wald than DUP)
        GeneDifferentialExpressionResult("ENSG_PC02", "TP53", 150.0, 2.0, 0.3, 6.67, 1e-9, 1e-8, "OK"),
        # 3. HOTAIR: Valid Wald, lncRNA without Entrez -> Retained non-Entrez eligible
        GeneDifferentialExpressionResult("ENSG_LNC01", "HOTAIR", 40.0, -1.5, 0.3, -5.0, 1e-6, 1e-5, "OK"),
        # 4. NEAT1: Valid Wald, lncRNA without Entrez -> Retained non-Entrez eligible
        GeneDifferentialExpressionResult("ENSG_LNC02", "NEAT1", 80.0, 0.8, 0.2, 4.0, 1e-4, 1e-3, "OK"),
        # 5. TEC_GENE: Valid Wald, TEC biotype without Entrez -> Retained non-Entrez eligible
        GeneDifferentialExpressionResult("ENSG_TEC01", "TEC_GENE", 25.0, 0.4, 0.2, 2.0, 0.05, 0.1, "OK"),
        # 6. OUT_GENE: Cook's distance outlier -> Excluded from statistical universe
        GeneDifferentialExpressionResult("ENSG_OUT01", "OUT_GENE", 500.0, 4.0, 0.5, 8.0, 1e-12, 1e-11, "OUTLIER"),
        # 7. ZERO_GENE: All zero counts -> Excluded from statistical universe
        GeneDifferentialExpressionResult("ENSG_ZERO1", "ZERO_GENE", 0.0, 0.0, None, None, None, None, "ZERO_COUNTS"),
        # 8. TP53_DUP: Duplicate Entrez 7157 with lower |wald| (4.0 vs 6.67) -> Dropped by MAX_ABS_WALD
        GeneDifferentialExpressionResult("ENSG_DUP01", "TP53_DUP", 90.0, 1.2, 0.3, 4.0, 1e-4, 1e-3, "OK"),
    )

    contrast = DifferentialExpressionContrastResult(
        contrast_id="contrast_synth_01",
        study_id="STUDY_SYNTH",
        experiment_id="EXP_SYNTH",
        cohort_id="COHORT_01",
        relationship_id="REL_01",
        treatment_condition_id="TRT",
        matched_control_condition_ids=("CTL",),
        treatment_sample_ids=("S1", "S2"),
        control_sample_ids=("S3", "S4"),
        treatment_biological_replicate_ids=("R1", "R2"),
        control_biological_replicate_ids=("R3", "R4"),
        treatment_replicate_count=2,
        control_replicate_count=2,
        agent_name="Test Chemical",
        agent_identifier="CID:0000",
        concentration_or_dose=50.0,
        concentration_or_dose_unit="uM",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        developmental_age_or_stage_at_exposure="Day 10",
        design_formula="~ condition",
        residual_degrees_of_freedom=2,
        size_factors_used={},
        total_input_genes_count=8,
        gene_results=gene_results,
    )
    return contrast, ref


def test_statistically_eligible_universe_preserves_all_valid_genes(
    synthetic_multibiotype_contrast: tuple[DifferentialExpressionContrastResult, GeneReference],
) -> None:
    contrast, ref = synthetic_multibiotype_contrast
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    res = builder.build_ranking(contrast, deduplicate_entrez=True)

    # 1. Total input genes = 8
    assert len(contrast.gene_results) == 8

    # 2. Statistically eligible universe must contain exactly 6 genes
    # (BRCA1, TP53, HOTAIR, NEAT1, TEC_GENE, TP53_DUP)
    # OUT_GENE and ZERO_GENE are excluded
    universe = res.statistically_eligible_universe
    assert len(universe) == 6
    universe_ids = {g.canonical_gene_id for g in universe}
    assert universe_ids == {"ENSG_PC01", "ENSG_PC02", "ENSG_LNC01", "ENSG_LNC02", "ENSG_TEC01", "ENSG_DUP01"}
    assert "ENSG_OUT01" not in universe_ids
    assert "ENSG_ZERO1" not in universe_ids

    # Every gene in universe must have statistical_eligibility == True
    for g in universe:
        assert g.statistical_eligibility is True
        assert g.result_status == "OK"
        assert g.wald_statistic is not None


def test_unmapped_entrez_genes_are_not_classified_as_qc_failures(
    synthetic_multibiotype_contrast: tuple[DifferentialExpressionContrastResult, GeneReference],
) -> None:
    contrast, ref = synthetic_multibiotype_contrast
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    res = builder.build_ranking(contrast, deduplicate_entrez=True)

    retained = res.retained_non_entrez_genes
    # 3 genes lack Entrez: HOTAIR, NEAT1, TEC_GENE
    assert len(retained) == 3
    retained_ids = {g.canonical_gene_id for g in retained}
    assert retained_ids == {"ENSG_LNC01", "ENSG_LNC02", "ENSG_TEC01"}

    for g in retained:
        assert g.future_feature_eligibility is True
        assert g.reason_not_in_kegg == "NO_ENTREZ_ID_IN_ENSEMBL_112"
        assert g.go_annotation_status == "NOT_ANNOTATED_IN_ENTREZ_GO"


def test_kegg_entrez_limitation_and_duplicate_handling_are_explicit(
    synthetic_multibiotype_contrast: tuple[DifferentialExpressionContrastResult, GeneReference],
) -> None:
    contrast, ref = synthetic_multibiotype_contrast
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    res = builder.build_ranking(contrast, deduplicate_entrez=True)

    # Entrez-ranked genes entering KEGG/GO GSEA:
    # 3 Entrez-mapped eligible genes: BRCA1 (672), TP53 (7157), TP53_DUP (7157)
    # Deduplication drops TP53_DUP (|wald|=4.0 < |wald|=6.67)
    # Final ranked genes = 2: TP53 (6.67) and BRCA1 (6.0)
    ranked_set = res.ranked_gene_set
    assert len(ranked_set.ranked_genes) == 2
    assert ranked_set.ranked_genes[0].canonical_gene_id == "ENSG_PC02"  # TP53 wald=6.67
    assert ranked_set.ranked_genes[1].canonical_gene_id == "ENSG_PC01"  # BRCA1 wald=6.00


def test_mathematical_reconciliation_of_qc_counts(
    synthetic_multibiotype_contrast: tuple[DifferentialExpressionContrastResult, GeneReference],
) -> None:
    contrast, ref = synthetic_multibiotype_contrast
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    res = builder.build_ranking(contrast, deduplicate_entrez=True)

    qc = compute_pathway_analysis_qc(
        stage7_gene_count=len(contrast.gene_results),
        mapping_result=res.mapping_result,
        ranked_set=res.ranked_gene_set,
        excluded_records=res.excluded_records,
        go_results=(),
        kegg_results=(),
        contrast=contrast,
        statistically_eligible_universe=res.statistically_eligible_universe,
        retained_non_entrez_genes=res.retained_non_entrez_genes,
    )

    # Check exact counts:
    assert qc.stage7_evaluated_genes == 8
    assert qc.zero_count_genes == 1
    assert qc.genes_with_valid_wald == 7
    assert qc.cooks_outlier_genes == 1
    assert qc.statistically_eligible_genes == 6
    assert qc.entrez_mapped_eligible_genes == 3
    assert qc.entrez_unmapped_eligible_genes == 3
    assert qc.duplicate_entrez_entries_removed == 1
    assert qc.final_kegg_ranked_genes == 2
    assert qc.final_go_ranked_genes == 2

    # Invariants:
    assert qc.stage7_evaluated_genes == qc.zero_count_genes + qc.genes_with_valid_wald
    assert qc.statistically_eligible_genes == qc.genes_with_valid_wald - qc.cooks_outlier_genes
    assert qc.statistically_eligible_genes == qc.entrez_mapped_eligible_genes + qc.entrez_unmapped_eligible_genes
    assert qc.final_kegg_ranked_genes == qc.entrez_mapped_eligible_genes - qc.duplicate_entrez_entries_removed


def test_export_helpers_produce_valid_tabular_outputs(
    synthetic_multibiotype_contrast: tuple[DifferentialExpressionContrastResult, GeneReference],
) -> None:
    contrast, ref = synthetic_multibiotype_contrast
    adapter = EnsemblToEntrezAdapter(gene_reference=ref)
    builder = GeneRankingBuilder(mapping_adapter=adapter)

    res = builder.build_ranking(contrast, deduplicate_entrez=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        u_csv = tmp_path / "universe.csv"
        r_csv = tmp_path / "retained.csv"

        export_statistically_eligible_universe_to_csv(res.statistically_eligible_universe, u_csv)
        export_retained_non_entrez_genes_to_csv(res.retained_non_entrez_genes, r_csv)

        assert u_csv.exists()
        assert r_csv.exists()

        u_lines = u_csv.read_text(encoding="utf-8").strip().splitlines()
        r_lines = r_csv.read_text(encoding="utf-8").strip().splitlines()

        # Header + rows
        assert len(u_lines) == 1 + 6
        assert len(r_lines) == 1 + 3

        assert "canonical_gene_id,approved_symbol" in u_lines[0]
        assert "future_feature_eligibility" in r_lines[0]
        assert "True" in r_lines[1]
