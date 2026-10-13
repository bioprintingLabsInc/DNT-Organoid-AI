"""Provenance tracking and execution audit recording for Stage 8C."""

import subprocess
from typing import Any

from src.differential_expression.models import DifferentialExpressionContrastResult
from .models import PathwayAnalysisProvenance


def get_current_git_commit() -> str:
    """Return HEAD commit hash or fallback string."""
    try:
        res = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
        if res.returncode == 0:
            return res.stdout.strip()
    except Exception:
        pass
    return "unknown"


def build_pathway_analysis_provenance(
    contrast: DifferentialExpressionContrastResult,
    ranking_metric: str = "wald_statistic",
    gene_reference_version: str = "Ensembl_112",
    pathway_databases: tuple[str, ...] = ("GO_BP", "KEGG"),
    software_packages: dict[str, str] | None = None,
    stage7_source_commit: str = "ca26d34e6c0263e743425cc7d6752d0e5773a9fe",
    gsea_parameters: dict[str, Any] | None = None,
) -> PathwayAnalysisProvenance:
    """Construct full traceability provenance record for Stage 8C."""
    current_commit = get_current_git_commit()
    pkgs = software_packages or {
        "clusterProfiler": "4.14.6",
        "org.Hs.eg.db": "3.20.0",
        "GO.db": "3.20.0",
        "fgsea": "1.32.4",
        "KEGGREST": "1.46.0",
        "R": "4.4.3",
        "Bioconductor": "3.20",
    }

    gsea_params = gsea_parameters or {
        "exponent": 1.0,
        "gsea_param": 1.0,
        "score_type": "std",
        "min_size": 15,
        "max_size": 500,
        "pvalue_cutoff": 1.0,
        "padjust_method": "BH",
        "eps": 0.0,
        "seed": 42,
        "nperm_simple": 1000,
    }

    go_prov = {
        "org_hs_eg_db_version": pkgs.get("org.Hs.eg.db", "3.20.0"),
        "go_db_version": pkgs.get("GO.db", "3.20.0"),
        "bioconductor_version": pkgs.get("Bioconductor", "3.20"),
        "r_version": pkgs.get("R", "4.4.3"),
        "key_type": "ENTREZID",
        "entrez_source_date": "2024-Sep20",
        "schema_version": "HUMAN_DB 2.1",
        "org_hs_eg_sqlite_sha256": "7b88d1cf970e8acef5da60b0535bbc1b53768381940d2b4fd538524f9b612d18",
        "go_sqlite_sha256": "a20386d939b63e94e50cfcb5c2e154ed17d7996611482b54d96807258cacc7ef",
    }

    kegg_prov = {
        "snapshot_file": "data/references/kegg_hsa_pathways.rds",
        "provenance_file": "data/references/kegg_hsa_pathways_provenance.json",
        "organism": "hsa",
        "species": "Homo sapiens",
        "pathway_count": 372,
        "gene_pathway_association_count": 39671,
        "sha256_rds": "A54298B77CE24D0059FB58CE0B3E22EB6837589D02FA4B2437E7545A02967C20",
        "sha256_json": "B454F803A2935333FEA52076BAF5E7D98C4FDEFE29EFA8EB75903BD5BBBA3DB8",
        "snapshot_date": "2026-10-13",
        "clusterprofiler_version": pkgs.get("clusterProfiler", "4.14.6"),
        "keggrest_version": pkgs.get("KEGGREST", "1.46.0"),
        "r_version": pkgs.get("R", "4.4.3"),
    }

    stat_criteria = (
        "baseMean > 0, finite Wald statistic (non-NaN, non-Inf), "
        "result_status == OK, non-outlier (Cook's distance passed)"
    )

    pathway_criteria = {
        "GO_BP": "Ensembl canonical gene mapped to Entrez Gene ID in Ensembl 112 with annotation in org.Hs.eg.db (keyType=ENTREZID)",
        "KEGG": "Ensembl canonical gene mapped to Entrez Gene ID in Ensembl 112 with membership in hsa KEGG pathways",
        "RETAINED_NON_ENTREZ": (
            "Statistically valid Stage 7 genes lacking NCBI Entrez ID in Ensembl 112 "
            "are preserved as valid molecular features for cross-study modeling"
        ),
    }

    return PathwayAnalysisProvenance(
        study_id=contrast.study_id,
        experiment_id=contrast.experiment_id,
        contrast_id=contrast.contrast_id,
        agent_name=contrast.agent_name,
        agent_identifier=contrast.agent_identifier,
        concentration=contrast.concentration_or_dose,
        concentration_unit=contrast.concentration_or_dose_unit,
        exposure_duration=contrast.exposure_duration,
        duration_unit=contrast.exposure_duration_unit,
        developmental_stage=contrast.developmental_age_or_stage_at_exposure,
        treatment_sample_ids=contrast.treatment_sample_ids,
        control_sample_ids=contrast.control_sample_ids,
        gene_reference_version=gene_reference_version,
        ranking_metric=ranking_metric,
        pathway_databases=pathway_databases,
        pathway_database_version="Bioconductor 3.20 / org.Hs.eg.db 3.20.0 / KEGG hsa snapshot 2026-10-13",
        software_packages=pkgs,
        stage7_source_commit=stage7_source_commit,
        stage8c_source_commit=current_commit,
        statistical_eligibility_criteria=stat_criteria,
        pathway_annotation_eligibility_criteria=pathway_criteria,
        go_provenance=go_prov,
        kegg_provenance=kegg_prov,
        gsea_parameters=gsea_params,
    )
