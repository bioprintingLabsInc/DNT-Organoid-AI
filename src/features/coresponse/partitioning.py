"""Stage 8B: AI/ML Partitioning and Data Leakage Safeguards.

Enforces strict separation between:
1. UNRESTRICTED_REFERENCE: Fixed external biological interactomes independent of empirical cohorts.
2. DESCRIPTIVE_EXPLORATION: Observational response summaries not directly used for ML evaluation.
3. TRAINING_ONLY: Empirical data-adaptive networks fitted strictly within training partitions.
4. HELD_OUT_EVALUATION: Validation and test chemicals prohibited from entering network fitting.
"""

from typing import Any

from src.features.coresponse.models import (
    BiologicalNetwork,
    CoresponseQCStatus,
    GeneNetworkEdge,
    NetworkRelationshipType,
    PartitionAssignment,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRecord,
    ConditionGeneResponseRegistry,
)


FORBIDDEN_LEAKAGE_TERMS = {
    "dnt_positive",
    "dnt_negative",
    "regulatory_positive",
    "regulatory_negative",
    "gold_standard_positive",
    "gold_standard_negative",
    "neurotoxic_label",
    "non_neurotoxic_label",
    "training_class",
    "test_class",
}


def filter_registry_by_training_partition(
    registry: ConditionGeneResponseRegistry,
    training_chemical_ids: set[str],
    training_study_family_ids: set[str] | None = None,
) -> ConditionGeneResponseRegistry:
    """Filter registry to include ONLY records from approved training chemicals and study families.

    Guarantees that validation/test chemicals or unapproved study families never participate
    in data-adaptive network construction or module inference.
    """
    clean_chem_ids = {c.strip() for c in training_chemical_ids}
    clean_fam_ids = {f.strip() for f in training_study_family_ids} if training_study_family_ids else None

    filtered_records: list[ConditionGeneResponseRecord] = []
    for r in registry.records:
        if r.compound_id not in clean_chem_ids and r.compound_name not in clean_chem_ids:
            continue
        if clean_fam_ids and r.study_family_id not in clean_fam_ids:
            continue
        filtered_records.append(r)

    return ConditionGeneResponseRegistry(
        registry_id=f"{registry.registry_id}_training_partition",
        records=tuple(filtered_records),
        registry_version=registry.registry_version,
    )


def verify_partition_isolation(
    network: BiologicalNetwork,
    held_out_chemicals: set[str],
) -> tuple[bool, list[str]]:
    """Audit network to verify zero contamination from held-out validation or test chemicals.

    Returns:
        tuple of (is_isolated: bool, violations: list[str])
    """
    violations: list[str] = []
    clean_held_out = {c.strip().upper() for c in held_out_chemicals}

    if network.relationship_type == NetworkRelationshipType.DATA_ADAPTIVE_PREDICTION_FEATURE:
        for edge in network.edges:
            for chem in edge.compound_context:
                if chem.strip().upper() in clean_held_out:
                    violations.append(
                        f"Edge {edge.edge_id}: Held-out evaluation chemical '{chem}' "
                        f"leaked into training-fitted prediction network!"
                    )

    # Check for forbidden supervisory outcome leakage
    for edge in network.edges:
        corpus = f"{edge.scientific_limitations} {edge.provenance_notes or ''}".lower()
        for term in FORBIDDEN_LEAKAGE_TERMS:
            if term in corpus:
                violations.append(
                    f"Edge {edge.edge_id}: Detected forbidden DNT supervisory term '{term}'."
                )

    is_isolated = (len(violations) == 0)
    return is_isolated, violations
