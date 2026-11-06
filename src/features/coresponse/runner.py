"""Stage 8B: Orchestration Runner for Gene Co-Response and Biological Network Features.

Coordinates:
1. Ingestion of established biological reference networks.
2. Gated within-study gene co-expression analysis.
3. Gated cross-condition gene co-response analysis.
4. Topology extraction, induced subnetwork isolation, and module detection.
5. Invariant auditing and canonical serialization.
"""

from pathlib import Path
from typing import Any, Sequence
import pandas as pd

from src.features.coresponse.cross_condition import (
    analyze_cross_condition_coresponse,
    validate_cross_condition_eligibility,
)
from src.features.coresponse.models import (
    BiologicalNetwork,
    CoresponseEvaluationSummary,
    CoresponseQCStatus,
    GeneNetworkTopologyMetrics,
    NetworkRelationshipType,
    PartitionAssignment,
)
from src.features.coresponse.qc import (
    NetworkQCSummary,
    audit_biological_network_qc,
)
from src.features.coresponse.reference_network import (
    filter_network_by_confidence,
    filter_network_by_genes,
    import_reference_network_from_table,
)
from src.features.coresponse.serialization import (
    export_network_to_csv,
    export_topology_metrics_to_csv,
)
from src.features.coresponse.topology import (
    compute_gene_network_topology,
    detect_network_modules,
    extract_chemical_responsive_subnetwork,
)
from src.features.coresponse.within_study import (
    analyze_within_study_coexpression,
    validate_coexpression_matrix_eligibility,
)
from src.features.cross_study.models import (
    ConditionGeneResponseRegistry,
)


def run_stage8b_pipeline(
    reference_network: BiologicalNetwork | None = None,
    registry: ConditionGeneResponseRegistry | None = None,
    expression_matrix: pd.DataFrame | None = None,
    sample_metadata: pd.DataFrame | None = None,
    study_id: str | None = None,
    responsive_genes: Sequence[str] | None = None,
    training_chemicals: set[str] | None = None,
    held_out_chemicals: set[str] | None = None,
    output_edges_csv: str | Path | None = None,
    output_topology_csv: str | Path | None = None,
) -> dict[str, Any]:
    """Execute end-to-end Stage 8B pipeline across available inputs.

    Gated honestly:
    - If registry has insufficient conditions/chemicals, cross-condition status is INSUFFICIENT_EVIDENCE.
    - If expression matrix has insufficient samples, within-study status is INSUFFICIENT_SAMPLES.
    - Reference network topology is evaluated without attributing chemical specificity.
    """
    outcomes: dict[str, Any] = {}

    # 1. Reference Network Processing (if provided)
    if reference_network is not None:
        ref_qc = audit_biological_network_qc(reference_network, held_out_chemicals=held_out_chemicals)
        if not ref_qc.is_valid:
            raise ValueError(f"Reference network failed QC: {ref_qc.findings}")

        # Compute topology across responsive genes or all genes
        topo = compute_gene_network_topology(reference_network, target_genes=responsive_genes)
        outcomes["reference_network"] = reference_network
        outcomes["reference_qc"] = ref_qc
        outcomes["reference_topology"] = topo

        if responsive_genes:
            sub_net = extract_chemical_responsive_subnetwork(reference_network, responsive_genes)
            outcomes["chemical_responsive_subnetwork"] = sub_net
            outcomes["modules"] = detect_network_modules(sub_net)

        if output_topology_csv and topo:
            export_topology_metrics_to_csv(topo, output_topology_csv)

    # 2. Within-Study Co-Expression (if matrix provided)
    if expression_matrix is not None:
        s_id = study_id or "unspecified_study"
        within_net, within_summary = analyze_within_study_coexpression(
            expression_matrix=expression_matrix,
            study_id=s_id,
            sample_metadata=sample_metadata,
        )
        outcomes["within_study_network"] = within_net
        outcomes["within_study_summary"] = within_summary

    # 3. Cross-Condition Co-Response (if registry provided)
    if registry is not None:
        cross_net, cross_summary = analyze_cross_condition_coresponse(
            registry=registry,
            training_chemicals=training_chemicals,
        )
        outcomes["cross_condition_network"] = cross_net
        outcomes["cross_condition_summary"] = cross_summary

        if output_edges_csv and cross_net.edge_count > 0:
            export_network_to_csv(cross_net, output_edges_csv)

    return outcomes
