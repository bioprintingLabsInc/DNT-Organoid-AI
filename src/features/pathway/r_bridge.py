"""Subprocess execution bridge invoking headless R clusterProfiler for Stage 8C GSEA."""

from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Any

from src.normalization.r_bridge import resolve_r_binary
from .models import PathwayEnrichmentResult, RankedGeneSet


R_GSEA_SCRIPT_PATH = Path("scripts/r/run_preranked_gsea.R")
DEFAULT_KEGG_RDS_PATH = Path("data/references/kegg_hsa_pathways.rds")


@dataclass(frozen=True, slots=True)
class GSEARunnerResult:
    """Raw result container returned by the R GSEA runner bridge."""

    status: str
    r_version: str
    packages: dict[str, str]
    go_results: tuple[PathwayEnrichmentResult, ...]
    kegg_results: tuple[PathwayEnrichmentResult, ...]
    gsea_parameters: dict[str, Any] | None = None
    error_message: str | None = None


def execute_preranked_gsea_r(
    ranked_genes: RankedGeneSet,
    r_binary_path: str = "Rscript",
    r_script_path: Path | str = R_GSEA_SCRIPT_PATH,
    kegg_rds_path: Path | str = DEFAULT_KEGG_RDS_PATH,
    min_size: int = 15,
    max_size: int = 500,
    pvalue_cutoff: float = 1.0,
    seed: int = 42,
    gsea_param: float = 1.0,
    score_type: str = "std",
    nperm_simple: int = 1000,
    eps: float = 0.0,
    padjust_method: str = "BH",
    timeout_seconds: int = 600,
) -> GSEARunnerResult:
    """Execute preranked GSEA in headless R subprocess using clusterProfiler and fgsea."""
    r_bin = resolve_r_binary(r_binary_path)
    if not r_bin:
        return GSEARunnerResult(
            status="ERROR",
            r_version="unknown",
            packages={},
            go_results=(),
            kegg_results=(),
            error_message="Could not resolve R/Rscript binary on host system.",
        )

    script_path = Path(r_script_path)
    if not script_path.exists():
        return GSEARunnerResult(
            status="ERROR",
            r_version="unknown",
            packages={},
            go_results=(),
            kegg_results=(),
            error_message=f"R GSEA script not found at {script_path}",
        )

    # Prepare temporary files for communication
    with tempfile.TemporaryDirectory(prefix="dnt_gsea_") as tmpdir:
        tmp_path = Path(tmpdir)
        ranks_file = tmp_path / "ranks.tsv"
        output_file = tmp_path / "gsea_output.json"

        # Write ranks TSV: entrez_gene_id, ranking_value, canonical_gene_id, approved_symbol
        with open(ranks_file, "w", encoding="utf-8") as f:
            f.write("entrez_gene_id\tranking_value\tcanonical_gene_id\tapproved_symbol\n")
            for g in ranked_genes.ranked_genes:
                symbol = g.approved_symbol or ""
                f.write(f"{g.entrez_gene_id}\t{g.ranking_value}\t{g.canonical_gene_id}\t{symbol}\n")

        cmd = [
            r_bin,
            str(script_path),
            "--ranks",
            str(ranks_file),
            "--output",
            str(output_file),
            "--min_size",
            str(min_size),
            "--max_size",
            str(max_size),
            "--pvalue_cutoff",
            str(pvalue_cutoff),
            "--seed",
            str(seed),
            "--gsea_param",
            str(gsea_param),
            "--score_type",
            str(score_type),
            "--nperm_simple",
            str(nperm_simple),
            "--eps",
            str(eps),
            "--padjust_method",
            str(padjust_method),
        ]

        kegg_path = Path(kegg_rds_path)
        if kegg_path.exists():
            cmd.extend(["--kegg_rds", str(kegg_path)])

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return GSEARunnerResult(
                status="TIMEOUT",
                r_version="unknown",
                packages={},
                go_results=(),
                kegg_results=(),
                error_message=f"R GSEA runner timed out after {timeout_seconds} seconds.",
            )
        except Exception as e:
            return GSEARunnerResult(
                status="ERROR",
                r_version="unknown",
                packages={},
                go_results=(),
                kegg_results=(),
                error_message=f"Subprocess invocation failure: {str(e)}",
            )

        if not output_file.exists():
            return GSEARunnerResult(
                status="ERROR",
                r_version="unknown",
                packages={},
                go_results=(),
                kegg_results=(),
                error_message=f"R script failed to produce output file. Returncode: {proc.returncode}.\nSTDOUT: {proc.stdout}\nSTDERR: {proc.stderr}",
            )

        try:
            with open(output_file, "r", encoding="utf-8") as f:
                data = json.load(f)
        except Exception as e:
            return GSEARunnerResult(
                status="ERROR",
                r_version="unknown",
                packages={},
                go_results=(),
                kegg_results=(),
                error_message=f"Failed to parse R runner output JSON: {str(e)}",
            )

        if data.get("status") != "OK":
            return GSEARunnerResult(
                status="ERROR",
                r_version=data.get("r_version", "unknown"),
                packages=data.get("packages", {}),
                go_results=(),
                kegg_results=(),
                gsea_parameters=data.get("gsea_parameters"),
                error_message=data.get("error_message", "Unknown error from R runner."),
            )

        contrast_id = ranked_genes.contrast_id

        def _to_results(raw_list: list[dict[str, Any]]) -> tuple[PathwayEnrichmentResult, ...]:
            res: list[PathwayEnrichmentResult] = []
            for item in raw_list:
                res.append(
                    PathwayEnrichmentResult(
                        contrast_id=contrast_id,
                        source_database=item["source_database"],
                        pathway_id=item["pathway_id"],
                        pathway_name=item["pathway_name"],
                        enrichment_score=float(item["enrichment_score"]),
                        normalized_enrichment_score=float(item["normalized_enrichment_score"]),
                        p_value=float(item["p_value"]),
                        adjusted_p_value=float(item["adjusted_p_value"]),
                        pathway_size=int(item["pathway_size"]),
                        direction=item["direction"],
                        leading_edge_gene_ids=tuple(str(x) for x in item.get("leading_edge_gene_ids", [])),
                        leading_edge_symbols=tuple(str(x) for x in item.get("leading_edge_symbols", [])),
                        ranking_metric=item.get("ranking_metric", "wald_statistic"),
                        analysis_version=item.get("analysis_version", "clusterProfiler"),
                        is_significant=bool(item.get("is_significant", False)),
                    )
                )
            return tuple(res)

        go_results = _to_results(data.get("go_results", []))
        kegg_results = _to_results(data.get("kegg_results", []))

        return GSEARunnerResult(
            status="OK",
            r_version=data.get("r_version", "unknown"),
            packages=data.get("packages", {}),
            go_results=go_results,
            kegg_results=kegg_results,
            gsea_parameters=data.get("gsea_parameters"),
            error_message=None,
        )
