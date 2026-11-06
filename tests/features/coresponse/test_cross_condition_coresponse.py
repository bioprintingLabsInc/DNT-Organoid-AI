"""Unit and regression tests for cross-condition gene co-response analysis in Stage 8B."""

from pathlib import Path
import pytest

from src.features.compound_pathway.runner import load_registry_from_csv
from src.features.coresponse.cross_condition import (
    analyze_cross_condition_coresponse,
    validate_cross_condition_eligibility,
)
from src.features.coresponse.models import (
    CoresponseQCStatus,
    InteractionNature,
    NetworkRelationshipType,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
    EstimateStatus,
    ResponseDirection,
)


def _make_dummy_response_record(
    gene_id: str,
    contrast_id: str,
    compound_id: str,
    compound_name: str,
    log2_fc: float,
    dose_val: float = 10.0,
    duration: float = 24.0,
    study_id: str = "study_synthetic",
) -> ConditionGeneResponseRecord:
    return ConditionGeneResponseRecord(
        canonical_ensembl_gene_id=gene_id,
        gene_symbol=f"SYM_{gene_id[-4:]}",
        study_id=study_id,
        study_family_id=study_id,
        independent_study_id=study_id,
        experiment_id=f"exp_{contrast_id}",
        contrast_id=contrast_id,
        treatment_condition_id=f"trt_{contrast_id}",
        matched_control_condition_ids=(f"ctrl_{contrast_id}",),
        compound_id=compound_id,
        compound_name=compound_name,
        dose_value=dose_val,
        dose_unit="uM",
        dose_original=f"{dose_val} uM",
        exposure_duration=duration,
        exposure_duration_unit="h",
        developmental_stage="Day 30",
        model_type="Cortical Organoid",
        biological_source_ids=("WA09",),
        organoid_context_ids=("ctx_01",),
        treatment_replicate_count=3,
        control_replicate_count=3,
        treatment_sample_ids=("s1", "s2", "s3"),
        control_sample_ids=("c1", "c2", "c3"),
        treatment_biological_replicate_ids=("b1", "b2", "b3"),
        control_biological_replicate_ids=("cb1", "cb2", "cb3"),
        base_mean=50.0,
        log2_fold_change=log2_fc,
        log2_fold_change_se=0.2,
        wald_statistic=log2_fc / 0.2,
        p_value=0.001,
        adjusted_p_value=0.01,
        response_direction=ResponseDirection.UP if log2_fc > 0 else ResponseDirection.DOWN,
        estimate_status=EstimateStatus.OK,
        missing_reason=None,
    )


def test_real_gse232218_gating_honestly_reports_insufficient_conditions() -> None:
    """Test that real GSE232218 (1 contrast) honestly emits INSUFFICIENT_CONDITIONS."""
    csv_path = Path("data/processed/GSE232218_VPA_200uM_condition_gene_response_registry.csv")
    if not csv_path.exists():
        pytest.skip("GSE232218 registry CSV not found.")

    reg = load_registry_from_csv(csv_path)
    net, summary = analyze_cross_condition_coresponse(
        registry=reg,
        min_conditions=8,
        min_chemicals=4,
    )

    # CRITICAL: 1 contrast cannot compute cross-condition covariance!
    assert summary.is_eligible is False
    assert summary.status == CoresponseQCStatus.INSUFFICIENT_CONDITIONS
    assert summary.evaluated_conditions_count == 1
    assert summary.evaluated_chemicals_count == 1
    assert summary.inferred_edges_count == 0
    assert net.edge_count == 0
    assert "Insufficient independent conditions: found 1 contrast" in summary.reasons[0]


def test_real_gse187006_u1m_and_u2f_honestly_reports_insufficient_chemicals() -> None:
    """Test that combined GSE187006 U1M/U2F (2 contrasts, 1 chemical) emits INSUFFICIENT_CHEMICALS."""
    csv_u1m = Path("data/processed/GSE187006_U1M_VPA_1mM_condition_gene_response_registry.csv")
    csv_u2f = Path("data/processed/GSE187006_U2F_VPA_1mM_condition_gene_response_registry.csv")
    if not csv_u1m.exists() or not csv_u2f.exists():
        pytest.skip("GSE187006 registries not found.")

    reg1 = load_registry_from_csv(csv_u1m)
    reg2 = load_registry_from_csv(csv_u2f)
    combined_reg = ConditionGeneResponseRegistry(
        records=reg1.records + reg2.records,
        registry_id="reg_gse187006_combined",
    )

    net, summary = analyze_cross_condition_coresponse(
        registry=combined_reg,
        min_conditions=4,
        min_chemicals=4,
    )

    assert summary.is_eligible is False
    assert summary.status in {CoresponseQCStatus.INSUFFICIENT_CONDITIONS, CoresponseQCStatus.INSUFFICIENT_CHEMICALS}
    assert summary.inferred_edges_count == 0
    assert net.edge_count == 0


def test_non_independent_conditions_of_same_chemical_rejected() -> None:
    """Test that multiple concentrations of the SAME chemical are not accepted as independent chemical evidence."""
    # 10 dose conditions of the same chemical (VPA)
    records = []
    for dose_idx in range(1, 11):
        for gene_idx in range(1, 6):
            r = _make_dummy_response_record(
                gene_id=f"ENSG000{gene_idx}",
                contrast_id=f"VPA_dose_{dose_idx}",
                compound_id="CID:3121",
                compound_name="Valproic Acid",
                log2_fc=0.5 * dose_idx,
            )
            records.append(r)

    reg = ConditionGeneResponseRegistry(records=tuple(records), registry_id="reg_dose_response")
    is_eligible, n_cond, n_chem, reason = validate_cross_condition_eligibility(
        registry=reg,
        min_conditions=8,
        min_chemicals=4,
    )

    assert is_eligible is False
    assert n_cond == 10
    assert n_chem == 1  # Only 1 unique chemical!
    assert "Insufficient chemical diversity" in reason


def test_cross_condition_coresponse_inference_on_diverse_synthetic_cohort() -> None:
    """Test co-response network inference on an eligible multi-chemical cohort (12 conditions, 6 chemicals)."""
    # 6 diverse chemicals: VPA, Chlorpyrifos, Rotenone, Paraquat, Lead, Methylmercury
    chemicals = [
        ("CID:3121", "Valproic Acid"),
        ("CID:2730", "Chlorpyrifos"),
        ("CID:6758", "Rotenone"),
        ("CID:15938", "Paraquat"),
        ("CID:5352425", "Lead"),
        ("CID:24085", "Methylmercury"),
    ]

    # 12 conditions (2 doses per chemical)
    contrasts = []
    for cid, cname in chemicals:
        contrasts.append((cid, cname, f"{cname}_LowDose", 10.0))
        contrasts.append((cid, cname, f"{cname}_HighDose", 50.0))

    # Gene 1 and Gene 2 co-respond strongly across all conditions (r > 0.85)
    records = []
    for c_idx, (chem_id, chem_name, cont_id, dose_val) in enumerate(contrasts):
        base_resp = float(c_idx - 5.5) * 0.4  # Response shifts across conditions

        # Gene 1 and Gene 2 track each other
        g1_fc = base_resp + 0.05
        g2_fc = base_resp - 0.05
        # Gene 3 is random
        g3_fc = float((c_idx % 3) - 1.0) * 0.2
        # Gene 4 is anti-correlated
        g4_fc = -base_resp

        records.append(_make_dummy_response_record("ENSG000001", cont_id, chem_id, chem_name, g1_fc, dose_val))
        records.append(_make_dummy_response_record("ENSG000002", cont_id, chem_id, chem_name, g2_fc, dose_val))
        records.append(_make_dummy_response_record("ENSG000003", cont_id, chem_id, chem_name, g3_fc, dose_val))
        records.append(_make_dummy_response_record("ENSG000004", cont_id, chem_id, chem_name, g4_fc, dose_val))

    reg = ConditionGeneResponseRegistry(records=tuple(records), registry_id="reg_multi_chemical_eligible")

    net, summary = analyze_cross_condition_coresponse(
        registry=reg,
        min_conditions=8,
        min_chemicals=4,
        method="spearman",
        min_abs_correlation=0.75,
        max_fdr=0.05,
    )

    assert summary.is_eligible is True
    assert summary.status == CoresponseQCStatus.VALID
    assert net.edge_count >= 1

    # Check ENSG000001 - ENSG000002 co-response edge
    g1_edges = net.get_edges_for_gene("ENSG000001")
    g2_partner = [e for e in g1_edges if "ENSG000002" in (e.gene_a_canonical_id, e.gene_b_canonical_id)][0]
    assert g2_partner.relationship_type == NetworkRelationshipType.CROSS_CONDITION_CO_RESPONSE
    assert g2_partner.interaction_nature == InteractionNature.STATISTICAL_CO_RESPONSE
    assert g2_partner.is_chemical_specific is True
    assert g2_partner.statistical_estimate is not None
    assert g2_partner.statistical_estimate > 0.85
    assert g2_partner.condition_count == 12
    assert len(g2_partner.compound_context) == 6
