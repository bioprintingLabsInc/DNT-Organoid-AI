"""Route execution dispatch layer connecting validated routing decisions to execution engines.

Bridges:
- Route.RAW_READ_PIPELINE -> Governed Quantification Bridge Executor
- Route.BULK_RAW_COUNT_PIPELINE -> Frozen Bulk Count Intake Pipeline
Preserves route provenance and explicit execution readiness boundaries.
"""

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from src.ingestion.models import AssetFinding, AssetResult, AssetStatus, RegisteredAsset, Route
from src.routing.router import route_metadata_assets


@dataclass(frozen=True, slots=True)
class RouteExecutionPlan:
    """Execution plan determined by canonical routing decision and modality analysis."""

    route: Route
    status: str
    executor_id: str | None
    target_format: str | None
    downstream_route: Route | None
    registered_assets: tuple[RegisteredAsset, ...]
    routing_results: tuple[AssetResult, ...]
    provenance_chain: dict[str, Any]
    findings: tuple[AssetFinding, ...]
    metadata: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "route": str(self.route),
            "status": self.status,
            "executor_id": self.executor_id,
            "target_format": self.target_format,
            "downstream_route": str(self.downstream_route) if self.downstream_route else None,
            "provenance_chain": self.provenance_chain,
            "findings": [asdict(f) for f in self.findings],
        }


@dataclass(frozen=True, slots=True)
class RouteExecutionResult:
    """Outcome of dispatching and executing a route plan."""

    plan: RouteExecutionPlan
    is_success: bool
    execution_state: str
    emitted_asset_path: Path | None
    emitted_format: str | None
    provenance_record: dict[str, Any]
    error_message: str | None = None
    artifacts: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "is_success": self.is_success,
            "execution_state": self.execution_state,
            "emitted_asset_path": str(self.emitted_asset_path) if self.emitted_asset_path else None,
            "emitted_format": self.emitted_format,
            "provenance_record": self.provenance_record,
            "error_message": self.error_message,
            "artifacts": self.artifacts,
        }


def plan_route_execution(metadata: dict[str, Any]) -> RouteExecutionPlan:
    """Inspect canonical metadata, route registered assets, and construct an execution plan.

    Enforces:
    1. If any asset is INVALID -> plan status INVALID_ASSETS (execution blocked).
    2. If any asset is MANUAL_REVIEW or NEEDS_REVIEW -> plan status MANUAL_REVIEW_REQUIRED.
    3. If all assets route to BULK_RAW_COUNT_PIPELINE -> executor is bulk_count_pipeline_runner.
    4. If all assets route to RAW_READ_PIPELINE:
       - If any asset format is SRA -> status RAW_READ_SRA_EXTRACTION_NOT_IMPLEMENTED.
       - If all asset formats are FASTQ -> executor is quantification_bridge, downstream is BULK_RAW_COUNT_PIPELINE.
    5. Other routes (single cell, processed) -> status EXECUTION_NOT_IMPLEMENTED.
    """
    from src.ingestion.registry import register_assets

    registrations = register_assets(metadata)
    routing_results = route_metadata_assets(metadata)

    findings: list[AssetFinding] = []
    for res in routing_results:
        findings.extend(res.findings)

    # 1. Check for asset invalidity
    if any(r.status == AssetStatus.INVALID for r in routing_results):
        return RouteExecutionPlan(
            route=Route.MANUAL_REVIEW,
            status="INVALID_ASSETS",
            executor_id=None,
            target_format=None,
            downstream_route=None,
            registered_assets=registrations,
            routing_results=routing_results,
            provenance_chain={
                "status": "INVALID_ASSETS",
                "error": "One or more registered input assets failed inspection or metadata validation.",
            },
            findings=tuple(findings),
            metadata=metadata,
        )

    # 2. Check for manual review
    if any(r.route == Route.MANUAL_REVIEW or r.status == AssetStatus.NEEDS_REVIEW for r in routing_results):
        return RouteExecutionPlan(
            route=Route.MANUAL_REVIEW,
            status="MANUAL_REVIEW_REQUIRED",
            executor_id=None,
            target_format=None,
            downstream_route=None,
            registered_assets=registrations,
            routing_results=routing_results,
            provenance_chain={
                "status": "MANUAL_REVIEW_REQUIRED",
                "reason": "Routing ambiguity, declared/detected mismatch, or unsupported modality combination.",
            },
            findings=tuple(findings),
            metadata=metadata,
        )

    # Determine consensus route
    routes = {r.route for r in routing_results}
    if len(routes) != 1:
        # Mixed routes in a single dataset require manual resolution
        findings.append(
            AssetFinding("WARNING", "heterogeneous_routes", "routing", f"Dataset contains mixed routes: {routes}")
        )
        return RouteExecutionPlan(
            route=Route.MANUAL_REVIEW,
            status="MANUAL_REVIEW_REQUIRED",
            executor_id=None,
            target_format=None,
            downstream_route=None,
            registered_assets=registrations,
            routing_results=routing_results,
            provenance_chain={"status": "MANUAL_REVIEW_REQUIRED", "reason": "Mixed routes detected."},
            findings=tuple(findings),
            metadata=metadata,
        )

    selected_route = next(iter(routes))

    # Case A: Bulk Raw Count Pipeline (Direct count matrix)
    if selected_route == Route.BULK_RAW_COUNT_PIPELINE:
        provenance = {
            "input_format": "raw_counts",
            "selected_route": "BULK_RAW_COUNT_PIPELINE",
            "executor": "bulk_count_pipeline_runner",
            "downstream_route": "BULK_RAW_COUNT_PIPELINE",
        }
        return RouteExecutionPlan(
            route=Route.BULK_RAW_COUNT_PIPELINE,
            status="READY_FOR_EXECUTION",
            executor_id="bulk_count_pipeline_runner",
            target_format="normalized_dataset",
            downstream_route=None,
            registered_assets=registrations,
            routing_results=routing_results,
            provenance_chain=provenance,
            findings=tuple(findings),
            metadata=metadata,
        )

    # Case B: Raw Read Pipeline (FASTQ or SRA)
    if selected_route == Route.RAW_READ_PIPELINE:
        detected_formats = {r.detected_format for r in routing_results}
        if "SRA" in detected_formats:
            provenance = {
                "detected_format": "SRA",
                "selected_route": "RAW_READ_PIPELINE",
                "sra_extraction_adapter": "fasterq-dump",
                "raw_read_executor": "quantification_bridge",
                "emitted_format": "raw_counts",
                "downstream_route": "BULK_RAW_COUNT_PIPELINE",
            }
            return RouteExecutionPlan(
                route=Route.RAW_READ_PIPELINE,
                status="READY_FOR_EXECUTION",
                executor_id="quantification_bridge",
                target_format="raw_counts",
                downstream_route=Route.BULK_RAW_COUNT_PIPELINE,
                registered_assets=registrations,
                routing_results=routing_results,
                provenance_chain=provenance,
                findings=tuple(findings),
                metadata=metadata,
            )

        # All FASTQ
        provenance = {
            "detected_format": "FASTQ",
            "selected_route": "RAW_READ_PIPELINE",
            "raw_read_executor": "quantification_bridge",
            "emitted_format": "raw_counts",
            "downstream_route": "BULK_RAW_COUNT_PIPELINE",
        }
        return RouteExecutionPlan(
            route=Route.RAW_READ_PIPELINE,
            status="READY_FOR_EXECUTION",
            executor_id="quantification_bridge",
            target_format="raw_counts",
            downstream_route=Route.BULK_RAW_COUNT_PIPELINE,
            registered_assets=registrations,
            routing_results=routing_results,
            provenance_chain=provenance,
            findings=tuple(findings),
            metadata=metadata,
        )

    # Case C: Other planned routes
    return RouteExecutionPlan(
        route=selected_route,
        status="EXECUTION_NOT_IMPLEMENTED",
        executor_id=None,
        target_format=None,
        downstream_route=None,
        registered_assets=registrations,
        routing_results=routing_results,
        provenance_chain={
            "selected_route": str(selected_route),
            "status": "EXECUTION_NOT_IMPLEMENTED",
        },
        findings=tuple(findings),
        metadata=metadata,
    )


def dispatch_route_execution(
    plan: RouteExecutionPlan,
    executor: Any = None,
    sra_extractor: Any = None,
    **kwargs: Any,
) -> RouteExecutionResult:
    """Dispatch a validated execution plan to its registered executor engine.

    Enforces that:
    - Plans not in READY_FOR_EXECUTION halt immediately with their governance status.
    - Route.RAW_READ_PIPELINE groups paired FASTQs deterministically and dispatches to the quantification bridge.
    - Downstream provenance preserves both the original raw-read route and the derived count matrix transition.
    """
    if plan.status != "READY_FOR_EXECUTION":
        return RouteExecutionResult(
            plan=plan,
            is_success=False,
            execution_state=plan.status,
            emitted_asset_path=None,
            emitted_format=None,
            provenance_record=plan.provenance_chain,
            error_message=f"Execution blocked: route plan status is '{plan.status}'.",
        )

    # Route: RAW_READ_PIPELINE -> quantification_bridge
    if plan.executor_id == "quantification_bridge":
        has_sra = any(r.local_path.suffix.lower() == ".sra" for r in plan.registered_assets)
        if has_sra:
            if sra_extractor is not None:
                extractor_fn = sra_extractor
            else:
                from src.quantification_bridge.sra_runner import extract_sra_to_fastq_pair
                extractor_fn = extract_sra_to_fastq_pair
            pairs_list = []
            extracted_dir = Path("data/interim/sra_extracted")
            temp_scratch_dir = Path("data/interim/sra_temp")
            for asset in plan.registered_assets:
                if asset.local_path.suffix.lower() == ".sra":
                    sample_id = asset.assay_ids[0] if asset.assay_ids else asset.asset_id
                    try:
                        try:
                            p = extractor_fn(
                                sra_path=asset.local_path,
                                output_dir=extracted_dir / asset.local_path.stem,
                                sample_id=sample_id,
                                temp_dir=temp_scratch_dir / asset.local_path.stem,
                            )
                        except TypeError:
                            p = extractor_fn(
                                sra_path=asset.local_path,
                                output_dir=extracted_dir / asset.local_path.stem,
                                sample_id=sample_id,
                            )
                        pairs_list.append(p)
                    except Exception as e:
                        return RouteExecutionResult(
                            plan=plan,
                            is_success=False,
                            execution_state="SRA_EXTRACTION_FAILED",
                            emitted_asset_path=None,
                            emitted_format=None,
                            provenance_record=plan.provenance_chain,
                            error_message=f"SRA extraction failed: {e}",
                        )
            pairs = tuple(pairs_list)
        else:
            from src.quantification_bridge.pairing import group_fastq_pairs
            try:
                pairs = group_fastq_pairs(plan.registered_assets, plan.metadata)
            except Exception as e:
                return RouteExecutionResult(
                    plan=plan,
                    is_success=False,
                    execution_state="FASTQ_PAIRING_FAILED",
                    emitted_asset_path=None,
                    emitted_format=None,
                    provenance_record=plan.provenance_chain,
                    error_message=f"FASTQ pairing validation failed: {e}",
                )

        if executor is not None:
            # Custom / mock executor provided (e.g. for testing)
            res = executor(pairs, plan=plan, **kwargs)
            return res

        # Default: instantiate QuantificationBridgePipeline
        from src.quantification_bridge.pipeline import QuantificationBridgePipeline
        pipeline = QuantificationBridgePipeline()
        study_id = plan.metadata.get("studies", [{}])[0].get("study_id", "generic_study")
        gate_state = pipeline.execute_raw_read_route(pairs=pairs, study_id=study_id, **kwargs)

        return RouteExecutionResult(
            plan=plan,
            is_success=True,
            execution_state=gate_state,
            emitted_asset_path=None,
            emitted_format=None,
            provenance_record={
                **plan.provenance_chain,
                "execution_state": gate_state,
                "paired_libraries": len(pairs),
            },
            artifacts={"pairs": [p.to_dict() for p in pairs]},
        )

    # Route: BULK_RAW_COUNT_PIPELINE -> bulk count pipeline runner
    if plan.executor_id == "bulk_count_pipeline_runner":
        if executor is not None:
            return executor(plan.registered_assets, plan=plan, **kwargs)

        return RouteExecutionResult(
            plan=plan,
            is_success=True,
            execution_state="READY_FOR_BULK_COUNT_PIPELINE",
            emitted_asset_path=plan.registered_assets[0].local_path if plan.registered_assets else None,
            emitted_format="raw_counts",
            provenance_record=plan.provenance_chain,
        )

    return RouteExecutionResult(
        plan=plan,
        is_success=False,
        execution_state="UNKNOWN_EXECUTOR",
        emitted_asset_path=None,
        emitted_format=None,
        provenance_record=plan.provenance_chain,
        error_message=f"No executor registered for executor_id '{plan.executor_id}'.",
    )
