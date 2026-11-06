"""Comprehensive unit, integration, and regression test suite for Stage 8D.

Tests:
1. Chemical-Gene relationship construction and context preservation.
2. Law of Conservation of the Gene Universe (0 genes dropped).
3. Pathway evidence linking for available Stage 8C outputs.
4. Graceful handling of unavailable pathway evidence (GSE187006 / unanalyzed).
5. Causality and direct binding safeguards.
6. Zero DNT label leakage / circularity prevention.
7. Generalization across non-VPA compounds (e.g. Chlorpyrifos, Rotenone) and diverse regimens.
8. Real artifact integration with GSE232218 (with GO/KEGG).
9. Real artifact integration with GSE187006 (U1M & U2F distinct arms, pathway unavailable).
"""

from pathlib import Path
import pytest
import pandas as pd

from src.features.compound_pathway.builder import (
    build_chemical_exposure_context,
    build_chemical_gene_records,
    build_chemical_pathway_records,
    link_chemical_gene_pathway_evidence,
)
from src.features.compound_pathway.models import (
    ChemicalExposureContext,
    ChemicalGenePathwayEvidenceDataset,
    ChemicalGenePathwayEvidenceRecord,
    ChemicalGeneRecord,
    ChemicalPathwayRecord,
    EvidenceType,
    GenePathwayMembershipType,
    PathwayEvidenceStatus,
)
from src.features.compound_pathway.qc import (
    EvidenceQCSummary,
    audit_evidence_dataset_qc,
)
from src.features.compound_pathway.runner import (
    build_stage8d_evidence_pipeline,
    load_pathway_results_from_csv,
    load_registry_from_csv,
)
from src.features.compound_pathway.serialization import (
    CANONICAL_EVIDENCE_COLUMNS,
    export_evidence_to_csv,
    export_evidence_to_dataframe,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)
from src.features.pathway.models import PathwayEnrichmentResult


def _make_dummy_response_record(
    gene_id: str = "ENSG00000141510",
    symbol: str = "TP53",
    compound_id: str = "CID:3121",
    compound_name: str = "Valproic Acid",
    dose_val: float = 200.0,
    dose_unit: str = "uM",
    duration: float = 48.0,
    stage: str = "Day 6",
    model: str = "Cortical Organoid",
    log2_fc: float = 1.5,
    p_val: float = 0.01,
    direction: ResponseDirection = ResponseDirection.UP,
    contrast_id: str = "contrast_1",
    study_id: str = "study_1",
) -> ConditionGeneResponseRecord:
    return ConditionGeneResponseRecord(
        canonical_ensembl_gene_id=gene_id,
        gene_symbol=symbol,
        study_id=study_id,
        study_family_id=study_id,
        independent_study_id=study_id,
        experiment_id="exp_1",
        contrast_id=contrast_id,
        treatment_condition_id="cond_trt",
        matched_control_condition_ids=("cond_ctrl",),
        compound_id=compound_id,
        compound_name=compound_name,
        dose_value=dose_val,
        dose_unit=dose_unit,
        dose_original=f"{dose_val} {dose_unit}",
        exposure_duration=duration,
        exposure_duration_unit="h",
        developmental_stage=stage,
        model_type=model,
        biological_source_ids=("donor_1",),
        organoid_context_ids=("ctx_1",),
        treatment_replicate_count=3,
        control_replicate_count=3,
        treatment_sample_ids=("s1", "s2", "s3"),
        control_sample_ids=("s4", "s5", "s6"),
        treatment_biological_replicate_ids=("r1", "r2", "r3"),
        control_biological_replicate_ids=("r4", "r5", "r6"),
        base_mean=100.0,
        log2_fold_change=log2_fc,
        log2_fold_change_se=0.2,
        wald_statistic=log2_fc / 0.2,
        p_value=p_val,
        adjusted_p_value=p_val,
        response_direction=direction,
        estimate_status=EstimateStatus.OK,
        missing_reason=None,
        cross_study_status=CrossStudyStatus.SINGLE_STUDY_ONLY,
    )


def _make_dummy_pathway_result(
    contrast_id: str = "contrast_1",
    pathway_id: str = "GO:0007399",
    name: str = "nervous system development",
    nes: float = 1.85,
    p_val: float = 0.001,
    leading_edge_ids: tuple[str, ...] = ("ENSG00000141510",),
    leading_edge_syms: tuple[str, ...] = ("TP53",),
) -> PathwayEnrichmentResult:
    return PathwayEnrichmentResult(
        contrast_id=contrast_id,
        source_database="GO_BP",
        pathway_id=pathway_id,
        pathway_name=name,
        enrichment_score=0.65,
        normalized_enrichment_score=nes,
        p_value=p_val,
        adjusted_p_value=p_val,
        pathway_size=120,
        direction="POSITIVE",
        leading_edge_gene_ids=leading_edge_ids,
        leading_edge_symbols=leading_edge_syms,
        ranking_metric="wald_statistic",
        analysis_version="clusterProfiler_4.14.6",
        is_significant=True,
    )


def test_chemical_gene_record_construction() -> None:
    """Test ChemicalGeneRecord accurately preserves Stage 8A.1 fields and disclaims direct binding."""
    r = _make_dummy_response_record()
    reg = ConditionGeneResponseRegistry(records=(r,), registry_id="test_reg")
    cg_records = build_chemical_gene_records(reg)

    assert len(cg_records) == 1
    cg = cg_records[0]
    assert cg.canonical_ensembl_gene_id == "ENSG00000141510"
    assert cg.gene_symbol == "TP53"
    assert cg.log2_fold_change == 1.5
    assert cg.p_value == 0.01
    assert cg.response_direction == ResponseDirection.UP
    assert cg.is_direct_target_interaction is False
    assert cg.evidence_type == EvidenceType.OBSERVED_TRANSCRIPTOMIC_PERTURBATION
    assert "Does NOT constitute verified direct binding" in cg.scientific_limitation
    assert cg.exposure.compound_name == "Valproic Acid"
    assert cg.exposure.dose_value == 200.0


def test_gene_universe_conservation_and_pathway_linking() -> None:
    """Test that all genes are preserved and leading edge genes are properly linked to pathways."""
    r1 = _make_dummy_response_record(gene_id="ENSG00000141510", symbol="TP53")
    r2 = _make_dummy_response_record(gene_id="ENSG00000139618", symbol="BRCA2", log2_fc=-1.2, direction=ResponseDirection.DOWN)
    reg = ConditionGeneResponseRegistry(records=(r1, r2), registry_id="test_reg")

    pw = _make_dummy_pathway_result(leading_edge_ids=("ENSG00000141510",), leading_edge_syms=("TP53",))

    dataset, qc = build_stage8d_evidence_pipeline(
        registry=reg,
        pathway_results=[pw],
        dataset_id="test_ds",
    )

    assert qc.is_valid is True
    assert qc.input_genes_count == 2
    assert qc.unique_genes_count == 2
    assert qc.pathway_annotated_records_count == 1
    assert qc.unannotated_records_count == 1

    # Check TP53 record (linked to pathway)
    tp53_recs = dataset.get_records_for_gene("ENSG00000141510")
    assert len(tp53_recs) == 1
    assert tp53_recs[0].pathway_id == "GO:0007399"
    assert tp53_recs[0].pathway_evidence_status == PathwayEvidenceStatus.AVAILABLE
    assert tp53_recs[0].is_leading_edge is True
    assert tp53_recs[0].gene_membership_type == GenePathwayMembershipType.LEADING_EDGE_MEMBER

    # Check BRCA2 record (not in leading edge, conserved as baseline)
    brca2_recs = dataset.get_records_for_gene("ENSG00000139618")
    assert len(brca2_recs) == 1
    assert brca2_recs[0].pathway_id is None
    assert brca2_recs[0].pathway_evidence_status == PathwayEvidenceStatus.AVAILABLE
    assert brca2_recs[0].is_leading_edge is False
    assert brca2_recs[0].gene_membership_type == GenePathwayMembershipType.UNANNOTATED


def test_pathway_evidence_unavailable_handling() -> None:
    """Test that when pathway evidence is absent, records are marked NOT_AVAILABLE without failing."""
    r1 = _make_dummy_response_record(gene_id="ENSG00000141510", symbol="TP53")
    reg = ConditionGeneResponseRegistry(records=(r1,), registry_id="test_reg")

    dataset, qc = build_stage8d_evidence_pipeline(
        registry=reg,
        pathway_results=None,
        pathway_status=PathwayEvidenceStatus.NOT_AVAILABLE,
        dataset_id="test_no_pw",
    )

    assert qc.is_valid is True
    assert qc.unique_genes_count == 1
    assert qc.pathway_annotated_records_count == 0
    assert qc.unannotated_records_count == 1

    rec = dataset.records[0]
    assert rec.pathway_id is None
    assert rec.pathway_evidence_status == PathwayEvidenceStatus.NOT_AVAILABLE
    assert "Pathway analysis not performed or unavailable" in rec.scientific_limitations


def test_causality_and_direction_safeguards() -> None:
    """Test that individual gene log2FC sign is decoupled from pathway direction and binding claims."""
    r_down = _make_dummy_response_record(gene_id="ENSG00000141510", symbol="TP53", log2_fc=-2.0, direction=ResponseDirection.DOWN)
    reg = ConditionGeneResponseRegistry(records=(r_down,), registry_id="test_reg")

    # Pathway has POSITIVE NES despite gene being DOWN
    pw_pos = _make_dummy_pathway_result(nes=2.1, leading_edge_ids=("ENSG00000141510",), leading_edge_syms=("TP53",))

    dataset, qc = build_stage8d_evidence_pipeline(registry=reg, pathway_results=[pw_pos])
    assert qc.is_valid is True
    rec = dataset.records[0]

    assert rec.gene_response_direction == ResponseDirection.DOWN
    assert rec.pathway_direction == "POSITIVE"
    assert rec.is_direct_target_interaction is False
    assert "Does NOT imply direct chemical-gene binding or directional pathway causality" in rec.scientific_limitations


def test_leakage_and_circularity_prevention() -> None:
    """Test that QC auditor detects and rejects DNT supervisory outcome labels."""
    r = _make_dummy_response_record()
    reg = ConditionGeneResponseRegistry(records=(r,), registry_id="test_reg")
    dataset, qc = build_stage8d_evidence_pipeline(registry=reg, pathway_results=[])
    assert qc.is_valid is True
    assert qc.leakage_violations_count == 0

    # Inject artificial leakage
    corrupted_rec = ChemicalGenePathwayEvidenceRecord(
        evidence_id=dataset.records[0].evidence_id,
        contrast_id=dataset.records[0].contrast_id,
        study_id=dataset.records[0].study_id,
        independent_study_id=dataset.records[0].independent_study_id,
        compound_id=dataset.records[0].compound_id,
        compound_name=dataset.records[0].compound_name,
        dose_value=dataset.records[0].dose_value,
        dose_unit=dataset.records[0].dose_unit,
        exposure_duration=dataset.records[0].exposure_duration,
        exposure_duration_unit=dataset.records[0].exposure_duration_unit,
        developmental_stage=dataset.records[0].developmental_stage,
        model_type=dataset.records[0].model_type,
        canonical_ensembl_gene_id=dataset.records[0].canonical_ensembl_gene_id,
        gene_symbol=dataset.records[0].gene_symbol,
        log2_fold_change=dataset.records[0].log2_fold_change,
        log2_fold_change_se=dataset.records[0].log2_fold_change_se,
        wald_statistic=dataset.records[0].wald_statistic,
        p_value=dataset.records[0].p_value,
        adjusted_p_value=dataset.records[0].adjusted_p_value,
        gene_response_direction=dataset.records[0].gene_response_direction,
        gene_estimate_status=dataset.records[0].gene_estimate_status,
        pathway_evidence_status=dataset.records[0].pathway_evidence_status,
        pathway_id=dataset.records[0].pathway_id,
        pathway_name=dataset.records[0].pathway_name,
        source_database=dataset.records[0].source_database,
        database_version=dataset.records[0].database_version,
        pathway_normalized_enrichment_score=dataset.records[0].pathway_normalized_enrichment_score,
        pathway_p_value=dataset.records[0].pathway_p_value,
        pathway_adjusted_p_value=dataset.records[0].pathway_adjusted_p_value,
        pathway_direction=dataset.records[0].pathway_direction,
        is_pathway_significant=dataset.records[0].is_pathway_significant,
        gene_membership_type=dataset.records[0].gene_membership_type,
        is_leading_edge=dataset.records[0].is_leading_edge,
        evidence_type=dataset.records[0].evidence_type,
        is_direct_target_interaction=dataset.records[0].is_direct_target_interaction,
        scientific_limitations="Classified as dnt_positive training_class",  # LEAKAGE VIOLATION
        provenance_notes=dataset.records[0].provenance_notes,
    )
    corrupted_dataset = ChemicalGenePathwayEvidenceDataset(
        dataset_id="corrupted",
        records=(corrupted_rec,),
        chemical_gene_records=dataset.chemical_gene_records,
        chemical_pathway_records=dataset.chemical_pathway_records,
        created_at=dataset.created_at,
    )
    corrupted_qc = audit_evidence_dataset_qc(corrupted_dataset)
    assert corrupted_qc.is_valid is False
    assert corrupted_qc.leakage_violations_count > 0


def test_generalization_non_vpa_synthetic_dataset() -> None:
    """Test that Stage 8D functions seamlessly on arbitrary non-VPA compounds (e.g. Chlorpyrifos) without hardcoding."""
    # Synthetic test fixture: Chlorpyrifos in human dopaminergic organoids
    r_cpf1 = _make_dummy_response_record(
        gene_id="ENSG00000105383",
        symbol="CDKN1A",
        compound_id="CID:2730",
        compound_name="Chlorpyrifos",
        dose_val=10.0,
        dose_unit="uM",
        duration=24.0,
        stage="Day 21",
        model="Human Midbrain Dopaminergic Organoid",
        log2_fc=2.45,
        p_val=1e-5,
        contrast_id="CPF_10uM_24h_vs_DMSO",
        study_id="STUDY_CPF_PILOT",
    )
    r_cpf2 = _make_dummy_response_record(
        gene_id="ENSG00000146648",
        symbol="EGFR",
        compound_id="CID:2730",
        compound_name="Chlorpyrifos",
        dose_val=10.0,
        dose_unit="uM",
        duration=24.0,
        stage="Day 21",
        model="Human Midbrain Dopaminergic Organoid",
        log2_fc=-1.85,
        direction=ResponseDirection.DOWN,
        p_val=2e-4,
        contrast_id="CPF_10uM_24h_vs_DMSO",
        study_id="STUDY_CPF_PILOT",
    )
    reg_cpf = ConditionGeneResponseRegistry(records=(r_cpf1, r_cpf2), registry_id="reg_cpf")

    pw_cpf = _make_dummy_pathway_result(
        contrast_id="CPF_10uM_24h_vs_DMSO",
        pathway_id="hsa04020",
        name="Calcium signaling pathway",
        nes=-2.05,
        leading_edge_ids=("ENSG00000146648",),
        leading_edge_syms=("EGFR",),
    )

    dataset, qc = build_stage8d_evidence_pipeline(
        registry=reg_cpf,
        pathway_results=[pw_cpf],
        dataset_id="cpf_evidence_ds",
    )

    assert qc.is_valid is True
    assert qc.unique_compounds_count == 1
    assert "CID:2730" in dataset.unique_compounds
    assert len(dataset.records) == 2
    assert qc.leakage_violations_count == 0

    df = export_evidence_to_dataframe(dataset)
    assert len(df) == 2
    assert list(df.columns) == CANONICAL_EVIDENCE_COLUMNS
    assert (df["compound_name"] == "Chlorpyrifos").all()
    assert (df["dose_value"] == 10.0).all()
    assert (df["model_type"] == "Human Midbrain Dopaminergic Organoid").all()


def test_real_gse232218_evidence_synthesis() -> None:
    """Integration test: Real GSE232218 condition-level gene responses linked to Stage 8C GO/KEGG outputs."""
    reg_csv = Path("data/processed/GSE232218_VPA_200uM_condition_gene_response_registry.csv")
    go_csv = Path("data/processed/GSE232218_VPA_200uM_GO_BP_results.csv")
    kegg_csv = Path("data/processed/GSE232218_VPA_200uM_KEGG_results.csv")

    if not reg_csv.exists() or not go_csv.exists() or not kegg_csv.exists():
        pytest.skip("Required GSE232218 real artifacts not found.")

    # Load registry
    from scripts.run_gse187006_and_comparability import load_gse232218_registry
    reg = load_gse232218_registry(reg_csv)

    # Load pathway results
    go_results = load_pathway_results_from_csv(go_csv)
    kegg_results = load_pathway_results_from_csv(kegg_csv)
    all_pathways = go_results + kegg_results

    out_csv = Path("data/processed/GSE232218_VPA_200uM_chemical_gene_pathway_evidence.csv")
    dataset, qc = build_stage8d_evidence_pipeline(
        registry=reg,
        pathway_results=all_pathways,
        dataset_id="evds_GSE232218_VPA_200uM",
        output_csv=out_csv,
    )

    assert qc.is_valid is True
    assert qc.leakage_violations_count == 0
    assert qc.numerical_mismatches_count == 0
    assert qc.input_genes_count == 60207
    assert qc.unique_genes_count == 60207
    assert qc.pathway_annotated_records_count > 0
    assert out_csv.exists()


def test_real_gse187006_u1m_and_u2f_evidence_synthesis() -> None:
    """Integration test: Real GSE187006 U1M and U2F registries synthesized with pathway status NOT_AVAILABLE."""
    u1m_csv = Path("data/processed/GSE187006_U1M_VPA_1mM_condition_gene_response_registry.csv")
    u2f_csv = Path("data/processed/GSE187006_U2F_VPA_1mM_condition_gene_response_registry.csv")

    if not u1m_csv.exists() or not u2f_csv.exists():
        pytest.skip("Required GSE187006 registry artifacts not found.")

    # Load using standard registry loader
    reg_u1m = load_registry_from_csv(
        u1m_csv,
        registry_id="reg_u1m",
        default_treatment_condition="U1M_VPA",
        default_control_condition="U1M_NC",
    )

    out_u1m_csv = Path("data/processed/GSE187006_U1M_VPA_1mM_chemical_gene_pathway_evidence.csv")
    ds_u1m, qc_u1m = build_stage8d_evidence_pipeline(
        registry=reg_u1m,
        pathway_results=None,
        pathway_status=PathwayEvidenceStatus.NOT_AVAILABLE,
        dataset_id="evds_GSE187006_U1M",
        output_csv=out_u1m_csv,
    )

    assert qc_u1m.is_valid is True
    assert qc_u1m.unique_genes_count == 15561
    assert qc_u1m.pathway_annotated_records_count == 0
    assert qc_u1m.unannotated_records_count == 15561
    assert out_u1m_csv.exists()

    # Repeat for U2F
    reg_u2f = load_registry_from_csv(
        u2f_csv,
        registry_id="reg_u2f",
        default_treatment_condition="U2F_VPA",
        default_control_condition="U2F_NC",
    )

    out_u2f_csv = Path("data/processed/GSE187006_U2F_VPA_1mM_chemical_gene_pathway_evidence.csv")
    ds_u2f, qc_u2f = build_stage8d_evidence_pipeline(
        registry=reg_u2f,
        pathway_results=None,
        pathway_status=PathwayEvidenceStatus.NOT_AVAILABLE,
        dataset_id="evds_GSE187006_U2F",
        output_csv=out_u2f_csv,
    )

    assert qc_u2f.is_valid is True
    assert qc_u2f.unique_genes_count == 15561
    assert qc_u2f.pathway_annotated_records_count == 0
    assert qc_u2f.unannotated_records_count == 15561
    assert out_u2f_csv.exists()
