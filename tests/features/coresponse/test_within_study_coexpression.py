"""Unit and regression tests for within-study gene co-expression analysis in Stage 8B."""

import numpy as np
import pandas as pd
import pytest

from src.features.coresponse.models import (
    CoresponseQCStatus,
    InteractionNature,
    NetworkRelationshipType,
)
from src.features.coresponse.within_study import (
    analyze_within_study_coexpression,
    compute_benjamini_hochberg_fdr,
    validate_coexpression_matrix_eligibility,
)


def test_coexpression_gating_rejects_insufficient_biological_samples() -> None:
    """Test that cohorts with insufficient independent observations (N < 10) are honestly rejected."""
    # GSE232218 style: 8 total samples (4 trt, 4 ctrl)
    expr = pd.DataFrame(
        np.random.lognormal(2.0, 0.5, size=(100, 8)),
        index=[f"ENSG{i:06d}" for i in range(100)],
        columns=[f"Sample_{i}" for i in range(8)],
    )

    net, summary = analyze_within_study_coexpression(
        expression_matrix=expr,
        study_id="study_pilot_small",
        min_independent_samples=10,
    )

    # CRITICAL: Gating prevents reporting correlation on under-replicated data
    assert summary.is_eligible is False
    assert summary.status == CoresponseQCStatus.INSUFFICIENT_SAMPLES
    assert summary.inferred_edges_count == 0
    assert net.edge_count == 0
    assert "Insufficient independent biological observations" in summary.reasons[0]


def test_technical_replicates_cannot_inflate_sample_count() -> None:
    """Test that technical replicates do not inflate biological degrees of freedom."""
    # 12 samples total, but from only 3 biological donors (4 technical replicates each)
    samples = [f"Sample_Donor{d}_Tech{t}" for d in range(1, 4) for t in range(1, 5)]
    expr = pd.DataFrame(
        np.random.lognormal(2.0, 0.5, size=(50, 12)),
        index=[f"ENSG{i:06d}" for i in range(50)],
        columns=samples,
    )

    meta = pd.DataFrame({
        "sample_id": samples,
        "biological_replicate_id": [f"Donor_{d}" for d in range(1, 4) for _ in range(4)],
    }).set_index("sample_id")

    # Even though total columns = 12, independent biological N = 3
    is_eligible, n_bio, reason = validate_coexpression_matrix_eligibility(
        expression_matrix=expr,
        sample_metadata=meta,
        min_independent_samples=8,
    )

    assert is_eligible is False
    assert n_bio == 3  # Correctly collapsed to biological donors
    assert "found N=3" in reason


def test_coexpression_statistical_inference_on_eligible_cohort() -> None:
    """Test statistically rigorous within-study co-expression on an eligible cohort (N=20 samples)."""
    np.random.seed(42)
    n_samples = 20
    n_genes = 30

    # Create synthetic expression where Gene 0 and Gene 1 are strongly correlated (r > 0.85)
    base_signal = np.random.normal(5.0, 1.0, size=n_samples)
    g0 = base_signal + np.random.normal(0.0, 0.2, size=n_samples)
    g1 = base_signal + np.random.normal(0.0, 0.2, size=n_samples)

    # Gene 2 is strongly anti-correlated with Gene 0
    g2 = -base_signal + np.random.normal(0.0, 0.2, size=n_samples)

    # Remaining genes are uncorrelated random noise
    noise = np.random.normal(5.0, 1.0, size=(n_genes - 3, n_samples))

    mat = np.vstack([g0, g1, g2, noise])
    # Shift to positive continuous scale
    mat = mat - np.min(mat) + 1.0

    gene_ids = [f"ENSG{i:06d}" for i in range(n_genes)]
    expr = pd.DataFrame(mat, index=gene_ids, columns=[f"Sample_{j}" for j in range(n_samples)])

    symbols = {gene_ids[0]: "GENE_A", gene_ids[1]: "GENE_B", gene_ids[2]: "GENE_C"}

    net, summary = analyze_within_study_coexpression(
        expression_matrix=expr,
        study_id="study_cohort_eligible",
        gene_symbols=symbols,
        min_independent_samples=15,
        method="spearman",
        min_abs_correlation=0.70,
        max_fdr=0.01,
    )

    assert summary.is_eligible is True
    assert summary.status == CoresponseQCStatus.VALID
    assert net.edge_count >= 1

    # Check that Gene 0 - Gene 1 edge was identified
    g0_edges = net.get_edges_for_gene(gene_ids[0])
    pair_targets = {e.gene_b_canonical_id if e.gene_a_canonical_id == gene_ids[0] else e.gene_a_canonical_id for e in g0_edges}
    assert gene_ids[1] in pair_targets

    edge_0_1 = [e for e in g0_edges if gene_ids[1] in (e.gene_a_canonical_id, e.gene_b_canonical_id)][0]
    assert edge_0_1.relationship_type == NetworkRelationshipType.WITHIN_STUDY_CO_EXPRESSION
    assert edge_0_1.interaction_nature == InteractionNature.STATISTICAL_CO_EXPRESSION
    assert edge_0_1.statistical_estimate is not None
    assert edge_0_1.statistical_estimate > 0.80
    assert edge_0_1.p_value is not None
    assert edge_0_1.p_value < 1e-4
    assert edge_0_1.adjusted_p_value is not None
    assert edge_0_1.adjusted_p_value <= 0.01
    assert edge_0_1.sample_size == 20
    assert edge_0_1.is_chemical_specific is False


def test_benjamini_hochberg_fdr_computation() -> None:
    """Test mathematical accuracy of Benjamini-Hochberg FDR correction."""
    p_vals = [0.001, 0.01, 0.03, 0.04, 0.80]
    fdr = compute_benjamini_hochberg_fdr(p_vals)
    assert len(fdr) == 5
    # Smallest p-val adjusted: 0.001 * 5 / 1 = 0.005
    assert abs(fdr[0] - 0.005) < 1e-6
    # Monotonicity preserved
    assert fdr[0] <= fdr[1] <= fdr[2] <= fdr[3] <= fdr[4]


def test_multi_condition_heterogeneity_rejected_to_prevent_simpsons_paradox() -> None:
    """Test that pooling samples across heterogeneous perturbation conditions is rejected.

    Simpson's Paradox / Bimodal Mean-Shift Confounding:
    When genes X and Y both shift up in response to treatment, pooling controls and treated
    samples creates an artificial positive correlation even when the genes are uncorrelated
    within each condition. Within-study co-expression must reject multi-condition pooling.
    """
    np.random.seed(101)
    ctrl_samples = [f"Ctrl_rep_{i}" for i in range(6)]
    trt_samples = [f"VPA_rep_{i}" for i in range(6)]
    all_samples = ctrl_samples + trt_samples

    # Two genes with zero within-group correlation (pure independent noise)
    g1_ctrl = np.random.normal(5.0, 0.2, size=6)
    g2_ctrl = np.random.normal(5.0, 0.2, size=6)
    # Both shift upward by +3.0 under treatment
    g1_trt = np.random.normal(8.0, 0.2, size=6)
    g2_trt = np.random.normal(8.0, 0.2, size=6)

    g1_all = np.concatenate([g1_ctrl, g1_trt])
    g2_all = np.concatenate([g2_ctrl, g2_trt])

    # In pooled data, correlation is extremely high (> 0.90) due purely to the treatment mean shift
    pooled_r = float(np.corrcoef(g1_all, g2_all)[0, 1])
    assert pooled_r > 0.90

    # Within control and within treatment, correlation is negligible (< 0.50)
    ctrl_r = float(np.corrcoef(g1_ctrl, g2_ctrl)[0, 1])
    trt_r = float(np.corrcoef(g1_trt, g2_trt)[0, 1])
    assert abs(ctrl_r) < 0.50
    assert abs(trt_r) < 0.50

    expr = pd.DataFrame(
        [g1_all, g2_all],
        index=["ENSG00000000001", "ENSG00000000002"],
        columns=all_samples,
    )

    metadata = pd.DataFrame({
        "sample_id": all_samples,
        "group": ["control"] * 6 + ["treatment"] * 6,
        "compound": ["Vehicle"] * 6 + ["Valproic Acid"] * 6,
    }).set_index("sample_id")

    # When metadata indicates heterogeneous conditions, eligibility check must reject pooling
    is_eligible, n_obs, reason = validate_coexpression_matrix_eligibility(
        expression_matrix=expr,
        sample_metadata=metadata,
        min_independent_samples=10,
    )

    assert is_eligible is False
    assert "Heterogeneous experimental conditions detected" in reason
    assert "Simpson's paradox" in reason

    # Within-study analysis should return empty network with INSUFFICIENT_EVIDENCE
    net, summary = analyze_within_study_coexpression(
        expression_matrix=expr,
        study_id="PRJNA802072_pooled_simulation",
        sample_metadata=metadata,
        min_independent_samples=10,
    )

    assert summary.is_eligible is False
    assert summary.status == CoresponseQCStatus.INSUFFICIENT_EVIDENCE
    assert net.edge_count == 0
    assert "Simpson's paradox" in summary.scientific_limitations


def test_prjna802072_under_replicated_arms_honestly_rejected() -> None:
    """Test that PRJNA802072 control and treatment arms (N=6 each) honestly emit INSUFFICIENT_SAMPLES."""
    np.random.seed(42)
    # Control arm: N=6 samples
    ctrl_expr = pd.DataFrame(
        np.random.lognormal(3.0, 0.5, size=(50, 6)),
        index=[f"ENSG{i:06d}" for i in range(50)],
        columns=[f"SRR_ctrl_{i}" for i in range(6)],
    )

    net_ctrl, summary_ctrl = analyze_within_study_coexpression(
        expression_matrix=ctrl_expr,
        study_id="PRJNA802072_control_arm",
        min_independent_samples=10,
    )

    assert summary_ctrl.is_eligible is False
    assert summary_ctrl.status == CoresponseQCStatus.INSUFFICIENT_SAMPLES
    assert summary_ctrl.inferred_edges_count == 0
    assert net_ctrl.edge_count == 0
    assert "found N=6, but minimum required is min_independent_samples=10" in summary_ctrl.reasons[0]

    # Treatment arm: N=6 samples
    trt_expr = pd.DataFrame(
        np.random.lognormal(3.0, 0.5, size=(50, 6)),
        index=[f"ENSG{i:06d}" for i in range(50)],
        columns=[f"SRR_trt_{i}" for i in range(6)],
    )

    net_trt, summary_trt = analyze_within_study_coexpression(
        expression_matrix=trt_expr,
        study_id="PRJNA802072_treatment_arm",
        min_independent_samples=10,
    )

    assert summary_trt.is_eligible is False
    assert summary_trt.status == CoresponseQCStatus.INSUFFICIENT_SAMPLES
    assert summary_trt.inferred_edges_count == 0
    assert net_trt.edge_count == 0

