"""Stage 8B: Gene Co-Response and Biological Network Features Package.

Provides universal, metadata-driven modules distinguishing:
1. Established biological reference networks.
2. Observed within-study gene co-expression.
3. Cross-condition gene co-response.
4. Data-adaptive prediction feature governance and leakage prevention.
"""

from src.features.coresponse.cross_condition import (
    analyze_cross_condition_coresponse,
    validate_cross_condition_eligibility,
)
from src.features.coresponse.models import (
    BiologicalNetwork,
    CoresponseEvaluationSummary,
    CoresponseQCStatus,
    EvidenceCategory,
    GeneNetworkEdge,
    GeneNetworkTopologyMetrics,
    InteractionNature,
    NetworkRelationshipType,
    PartitionAssignment,
    derive_edge_id,
)
from src.features.coresponse.partitioning import (
    filter_registry_by_training_partition,
    verify_partition_isolation,
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
from src.features.coresponse.runner import (
    run_stage8b_pipeline,
)
from src.features.coresponse.serialization import (
    CANONICAL_NETWORK_EDGE_COLUMNS,
    CANONICAL_TOPOLOGY_METRICS_COLUMNS,
    export_network_to_csv,
    export_network_to_dataframe,
    export_topology_metrics_to_csv,
    export_topology_metrics_to_dataframe,
)
from src.features.coresponse.topology import (
    REFERENCE_DEGREE_BIAS_DISCLAIMER,
    build_networkx_graph,
    compute_gene_network_topology,
    detect_network_modules,
    extract_chemical_responsive_subnetwork,
)
from src.features.coresponse.within_study import (
    analyze_within_study_coexpression,
    compute_benjamini_hochberg_fdr,
    validate_coexpression_matrix_eligibility,
)

__all__ = [
    "BiologicalNetwork",
    "CANONICAL_NETWORK_EDGE_COLUMNS",
    "CANONICAL_TOPOLOGY_METRICS_COLUMNS",
    "CoresponseEvaluationSummary",
    "CoresponseQCStatus",
    "EvidenceCategory",
    "GeneNetworkEdge",
    "GeneNetworkTopologyMetrics",
    "InteractionNature",
    "NetworkQCSummary",
    "NetworkRelationshipType",
    "PartitionAssignment",
    "REFERENCE_DEGREE_BIAS_DISCLAIMER",
    "analyze_cross_condition_coresponse",
    "analyze_within_study_coexpression",
    "audit_biological_network_qc",
    "build_networkx_graph",
    "compute_benjamini_hochberg_fdr",
    "compute_gene_network_topology",
    "derive_edge_id",
    "detect_network_modules",
    "export_network_to_csv",
    "export_network_to_dataframe",
    "export_topology_metrics_to_csv",
    "export_topology_metrics_to_dataframe",
    "extract_chemical_responsive_subnetwork",
    "filter_network_by_confidence",
    "filter_network_by_genes",
    "filter_registry_by_training_partition",
    "import_reference_network_from_table",
    "run_stage8b_pipeline",
    "validate_coexpression_matrix_eligibility",
    "validate_cross_condition_eligibility",
    "verify_partition_isolation",
]
