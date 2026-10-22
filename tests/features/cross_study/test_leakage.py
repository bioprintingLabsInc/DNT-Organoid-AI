"""Tests verifying complete absence of DNT label leakage in Stage 8A.1.

Stage 8A is a molecular-feature stage. It operates strictly without:
- DNT outcome labels (positive/negative/uncertain)
- Regulatory toxicology classifications (EPA, OECD, EFSA, NTP)
- Expected biological answers or target endpoints
- Future training/test split information
"""

import sys
import pytest

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
)
from src.features.cross_study.registry import (
    build_condition_gene_response_registry,
)


def test_zero_import_of_labels_module() -> None:
    """Verify that cross_study module does not import or depend on src.labels."""
    import src.features.cross_study.builder
    import src.features.cross_study.models
    import src.features.cross_study.qc
    import src.features.cross_study.registry
    import src.features.cross_study.serialization

    for mod_name in sys.modules:
        if mod_name.startswith("src.features.cross_study"):
            mod = sys.modules[mod_name]
            # Ensure src.labels is never in module dict
            for attr_name, attr_val in vars(mod).items():
                if hasattr(attr_val, "__name__"):
                    assert "src.labels" not in attr_val.__name__, (
                        f"Leakage detected: {mod_name} imports {attr_val.__name__}"
                    )


def test_records_contain_zero_dnt_labels(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    """Verify that records and registry contain zero DNT outcome label keys or values."""
    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200],
        registry_id="reg_leakage_check",
    )

    forbidden_terms = {
        "dnt_positive",
        "dnt_negative",
        "is_dnt",
        "dnt_label",
        "toxicology_ground_truth",
        "reference_toxicant",
        "positive_control_label",
        "oecd_guideline_endpoint",
        "ntp_dnt_diver_label",
        "model_target",
        "y_true",
    }

    for record in registry.records:
        rec_dict = record.to_dict()
        for k in rec_dict:
            assert k.lower() not in forbidden_terms
            v = rec_dict[k]
            if isinstance(v, str):
                assert v.lower() not in forbidden_terms

    # Check registry provenance
    assert registry.provenance.get("dnt_labels_used") is False
    for k in registry.provenance:
        assert k.lower() not in forbidden_terms


def test_functions_seamlessly_without_any_labels(
    synthetic_contrast_vpa_200: DifferentialExpressionContrastResult,
) -> None:
    """Demonstrate that Stage 8A.1 runs to completion without any Stage 9 label artifacts."""
    registry = build_condition_gene_response_registry(
        contrasts=[synthetic_contrast_vpa_200],
    )
    assert registry.total_records > 0
    assert len(registry.unique_canonical_genes) > 0
