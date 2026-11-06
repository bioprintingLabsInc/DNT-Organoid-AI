"""Tests for Stage 8B Partitioning and Data Leakage Safeguards.

Verifies:
1. Strict filtering of ConditionGeneResponseRegistry by training compounds and study families.
2. Prevention of held-out validation/test compounds from leaking into data-adaptive prediction networks.
3. Strict detection and rejection of forbidden supervisory DNT terms (e.g. dnt_positive, gold_standard_positive).
4. Full isolation validation across reference, descriptive, and training partition networks.
"""

import pytest
from src.features.coresponse.models import (
    BiologicalNetwork,
    EvidenceCategory,
    GeneNetworkEdge,
    InteractionNature,
    NetworkRelationshipType,
    PartitionAssignment,
    derive_edge_id,
)
from src.features.coresponse.partitioning import (
    FORBIDDEN_LEAKAGE_TERMS,
    filter_registry_by_training_partition,
    verify_partition_isolation,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    CrossStudyStatus,
    EstimateStatus,
    ResponseDirection,
)


def _make_dummy_response_record(
    gene_id: str,
    contrast_id: str,
    compound_id: str,
    compound_name: str,
    study_family_id: str = "study_fam_1",
) -> ConditionGeneResponseRecord:
    return ConditionGeneResponseRecord(
        canonical_ensembl_gene_id=gene_id,
        gene_symbol="DUMMY",
        study_id="study_1",
        study_family_id=study_family_id,
        independent_study_id=f"indep_{study_family_id}",
        experiment_id="exp_1",
        contrast_id=contrast_id,
        treatment_condition_id=f"{contrast_id}_treatment",
        matched_control_condition_ids=(f"{contrast_id}_control",),
        compound_id=compound_id,
        compound_name=compound_name,
        dose_value=10.0,
        dose_unit="uM",
        dose_original="10 uM",
        exposure_duration=48.0,
        exposure_duration_unit="h",
        developmental_stage="Day 35",
        model_type="forebrain_organoid",
        biological_source_ids=("source_1",),
        organoid_context_ids=("ctx_1",),
        treatment_replicate_count=3,
        control_replicate_count=3,
        treatment_sample_ids=("s1", "s2", "s3"),
        control_sample_ids=("c1", "c2", "c3"),
        treatment_biological_replicate_ids=("b1", "b2", "b3"),
        control_biological_replicate_ids=("cb1", "cb2", "cb3"),
        base_mean=100.0,
        log2_fold_change=1.5,
        log2_fold_change_se=0.2,
        wald_statistic=7.5,
        p_value=1e-5,
        adjusted_p_value=1e-4,
        response_direction=ResponseDirection.UP,
        estimate_status=EstimateStatus.OK,
        missing_reason=None,
        cross_study_status=CrossStudyStatus.SINGLE_STUDY_ONLY,
        source_stage="stage7",
        source_contrast_id=contrast_id,
        de_version="1.0.0",
        registry_version="1.0.0",
        input_hashes={},
    )


def _make_edge(
    gene_a: str,
    gene_b: str,
    relationship_type: NetworkRelationshipType,
    compounds: tuple[str, ...] = (),
    notes: str = "",
    limitations: str = "Observational",
) -> GeneNetworkEdge:
    return GeneNetworkEdge(
        edge_id=derive_edge_id(gene_a, gene_b, str(relationship_type), "Stage8B_Inference"),
        gene_a_canonical_id=gene_a,
        gene_b_canonical_id=gene_b,
        relationship_type=relationship_type,
        interaction_nature=InteractionNature.STATISTICAL_CO_RESPONSE,
        evidence_source="Stage8B_Inference",
        evidence_version="1.0.0",
        evidence_category=EvidenceCategory.EMPIRICAL_CORESPONSE,
        statistical_metric_name="spearman_rho",
        statistical_estimate=0.85,
        p_value=0.001,
        adjusted_p_value=0.01,
        compound_context=compounds,
        scientific_limitations=limitations,
        provenance_notes=notes,
    )


def test_filter_registry_by_training_partition_strictly_isolates_records() -> None:
    """Ensure records from non-training chemicals or families are strictly excluded."""
    r1 = _make_dummy_response_record("ENSG01", "c1", "CHEM_TRAIN_A", "TrainA", "fam_1")
    r2 = _make_dummy_response_record("ENSG02", "c2", "CHEM_TRAIN_B", "TrainB", "fam_1")
    r3 = _make_dummy_response_record("ENSG03", "c3", "CHEM_TEST_HELD_OUT", "TestX", "fam_1")
    r4 = _make_dummy_response_record("ENSG04", "c4", "CHEM_TRAIN_A", "TrainA", "fam_unapproved")

    reg = ConditionGeneResponseRegistry(records=(r1, r2, r3, r4), registry_id="mixed_reg")

    # Filter with training chemical set and approved families
    train_reg = filter_registry_by_training_partition(
        registry=reg,
        training_chemical_ids={"CHEM_TRAIN_A", "CHEM_TRAIN_B"},
        training_study_family_ids={"fam_1"},
    )

    assert len(train_reg.records) == 2
    compounds_in_train = {r.compound_id for r in train_reg.records}
    assert compounds_in_train == {"CHEM_TRAIN_A", "CHEM_TRAIN_B"}
    assert "CHEM_TEST_HELD_OUT" not in compounds_in_train
    assert all(r.study_family_id == "fam_1" for r in train_reg.records)


def test_verify_partition_isolation_flags_held_out_compound_leakage() -> None:
    """Ensure prediction features containing held-out chemicals are caught and reported."""
    edge1 = _make_edge("ENSG01", "ENSG02", NetworkRelationshipType.DATA_ADAPTIVE_PREDICTION_FEATURE, compounds=("CHEM_A", "CHEM_B"))
    edge2 = _make_edge("ENSG02", "ENSG03", NetworkRelationshipType.DATA_ADAPTIVE_PREDICTION_FEATURE, compounds=("CHEM_B", "CHEM_TEST_HELD_OUT"))

    net = BiologicalNetwork(
        network_id="pred_net",
        network_name="Prediction Feature Network",
        relationship_type=NetworkRelationshipType.DATA_ADAPTIVE_PREDICTION_FEATURE,
        edges=(edge1, edge2),
        source_database="DataAdaptivePipeline",
    )

    is_isolated, violations = verify_partition_isolation(
        network=net,
        held_out_chemicals={"CHEM_TEST_HELD_OUT", "CHEM_VAL_ANOTHER"},
    )

    assert not is_isolated
    assert len(violations) == 1
    assert "CHEM_TEST_HELD_OUT" in violations[0]
    assert "leaked into training-fitted prediction network" in violations[0]


def test_verify_partition_isolation_detects_forbidden_supervisory_terms() -> None:
    """Ensure labels or supervisory outcomes like dnt_positive trigger strict leakage violations."""
    edge_clean = _make_edge(
        "ENSG01", "ENSG02",
        NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION,
        notes="Computed from unexposed control replicates",
    )
    edge_leaky = _make_edge(
        "ENSG03", "ENSG04",
        NetworkRelationshipType.DATA_ADAPTIVE_PREDICTION_FEATURE,
        notes="Fitted using gold_standard_positive chemical responses",
    )

    net = BiologicalNetwork(
        network_id="supervisory_leak_net",
        network_name="Supervisory Leak Network",
        relationship_type=NetworkRelationshipType.DATA_ADAPTIVE_PREDICTION_FEATURE,
        edges=(edge_clean, edge_leaky),
        source_database="Adaptive",
    )

    is_isolated, violations = verify_partition_isolation(
        network=net,
        held_out_chemicals=set(),
    )

    assert not is_isolated
    assert any("gold_standard_positive" in v for v in violations)


def test_reference_network_unrestricted_by_held_out_chemicals() -> None:
    """Established reference networks (e.g. STRING, OmniPath) are external and do not depend on chemical partitions."""
    ref_edge = GeneNetworkEdge(
        edge_id=derive_edge_id("ENSG000001", "ENSG000002", str(NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION), "STRING"),
        gene_a_canonical_id="ENSG000001",
        gene_b_canonical_id="ENSG000002",
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        interaction_nature=InteractionNature.PHYSICAL_INTERACTION,
        evidence_source="STRING",
        evidence_version="v12.0",
        evidence_category=EvidenceCategory.DATABASE_CURATED,
        confidence_score=0.95,
        partition_assignment=PartitionAssignment.UNRESTRICTED_REFERENCE,
        is_chemical_specific=False,
        scientific_limitations="Literature curated reference interaction.",
    )

    net = BiologicalNetwork(
        network_id="string_v12",
        network_name="STRING PPI",
        relationship_type=NetworkRelationshipType.ESTABLISHED_REFERENCE_INTERACTION,
        edges=(ref_edge,),
        source_database="STRING",
    )

    # Even if testing a held-out set, reference networks have no compound context and pass isolation
    is_isolated, violations = verify_partition_isolation(
        network=net,
        held_out_chemicals={"VALPROIC_ACID", "ROTENONE"},
    )

    assert is_isolated
    assert len(violations) == 0
