"""Pipeline runner orchestrating end-to-end Bulk RNA-seq differential expression.

Connects:
1. Count Matrix & Sample Info Parsing
2. Step 6A: Bulk Raw-Count QC (`validate_bulk_counts`)
3. Step 6B: Gene Harmonization v1 (`harmonize_bulk_counts`)
4. Step 6C: Bulk Normalization v1 (`normalize_bulk_dataset` with official DESeq2 size factors)
5. Step 7A: Treatment-versus-Matched-Control Contrast Builder v1 (`build_response_contrasts`)
6. Step 7B: Bulk Differential Expression v1 (`run_differential_expression` with official DESeq2 Wald test)
"""

from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
import io
import json
import os
from pathlib import Path
import random
import re
import tempfile
from typing import Any

import pandas as pd

from src.differential_expression.models import (
    DifferentialExpressionContrastResult,
    DifferentialExpressionDataset,
    ExcludedGeneAudit,
    GeneDifferentialExpressionResult,
)
from src.differential_expression.pipeline import run_differential_expression
from src.gene_harmonization.config import load_config as load_harmonization_config
from src.gene_harmonization.harmonizer import harmonize_bulk_counts
from src.gene_harmonization.models import HarmonizedDataset
from src.gene_harmonization.reference import GeneReference
from src.normalization.config import NormalizationConfig
from src.normalization.models import NormalizedDataset
from src.normalization.normalizer import normalize_bulk_dataset
from src.qc.matrix import CountMatrix, read_matrix
from src.qc.result import DatasetResult, Status as QCStatus
from src.qc.validator import validate_bulk_counts
from src.response_builder.builder import build_response_contrasts
from src.response_builder.models import (
    ContrastStatus,
    MolecularResponseContrast,
    ResponseContrastDataset,
)

# Global cached GeneReference to avoid re-reading 5MB TSV on every run
_CACHED_GENE_REFERENCE: GeneReference | None = None


def get_gene_reference() -> GeneReference:
    """Retrieve or load cached GeneReference for Ensembl Release 112."""
    global _CACHED_GENE_REFERENCE
    if _CACHED_GENE_REFERENCE is None:
        cfg = load_harmonization_config()
        _CACHED_GENE_REFERENCE = GeneReference.load_from_config(cfg)
    return _CACHED_GENE_REFERENCE


@dataclass(frozen=True, slots=True)
class ChemicalParams:
    """User-entered exposure, biological context, and developmental timing parameters."""

    agent_name: str = "Rotenone"
    agent_identifier: str = "CID:6758"
    concentration_or_dose: float = 0.5
    concentration_or_dose_unit: str = "uM"
    developmental_age: str = "day_35"
    exposure_duration: float = 24.0
    exposure_duration_unit: str = "h"
    vehicle: str = "0.1% DMSO"
    organoid_context_id: str = "cerebral_organoid_ctx"
    organoid_type: str = "cerebral_organoid"
    brain_region: str = "cerebral_cortex"
    biological_source_id: str = "ipsc_donor_source"
    species: str = "Homo sapiens"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True, slots=True)
class PipelineRunResult:
    """Complete results of the end-to-end differential expression pipeline."""

    is_success: bool
    contrast_result: DifferentialExpressionContrastResult | None
    de_genes_df: pd.DataFrame
    excluded_genes_df: pd.DataFrame
    samples_df: pd.DataFrame
    summary_metrics: dict[str, Any]
    scientific_findings: list[dict[str, str]]
    error_message: str | None = None
    qc_result: DatasetResult | None = None
    harmonized_dataset: HarmonizedDataset | None = None
    normalized_dataset: NormalizedDataset | None = None
    contrast_dataset: ResponseContrastDataset | None = None
    de_dataset: DifferentialExpressionDataset | None = None

    def to_summary_dict(self) -> dict[str, Any]:
        return {
            "is_success": self.is_success,
            "error_message": self.error_message,
            "summary_metrics": self.summary_metrics,
            "scientific_findings": self.scientific_findings,
            "contrast_info": self.contrast_result.to_dict() if self.contrast_result else None,
        }


def parse_count_matrix(file_or_content: Any) -> tuple[pd.DataFrame, list[str]]:
    """Parse uploaded file, string, or DataFrame into a normalized count DataFrame.

    Returns:
        (df, sample_ids) where df has 'gene_id' column followed by sample columns.
    """
    if isinstance(file_or_content, pd.DataFrame):
        df = file_or_content.copy()
    elif isinstance(file_or_content, (str, Path)) and os.path.exists(str(file_or_content)):
        path = str(file_or_content)
        sep = "\t" if path.endswith((".tsv", ".tab", ".txt")) else ","
        df = pd.read_csv(path, sep=sep)
    elif hasattr(file_or_content, "read"):
        # Streamlit UploadedFile or file-like object
        content = file_or_content.read()
        if isinstance(content, bytes):
            content = content.decode("utf-8-sig")
        sep = "\t" if ("\t" in content and "," not in content.split("\n")[0]) else ","
        df = pd.read_csv(io.StringIO(content), sep=sep)
    elif isinstance(file_or_content, str):
        sep = "\t" if ("\t" in file_or_content and "," not in file_or_content.split("\n")[0]) else ","
        df = pd.read_csv(io.StringIO(file_or_content), sep=sep)
    else:
        raise ValueError(f"Unsupported count matrix input type: {type(file_or_content)}")

    # Ensure gene identifier is first column
    if df.empty or len(df.columns) < 2:
        raise ValueError("Count matrix must contain at least a gene column and at least 2 sample columns.")

    # Normalize gene column name
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "gene_id"})
    df["gene_id"] = df["gene_id"].astype(str).str.strip()

    # Verify sample columns
    sample_ids = [str(c).strip() for c in df.columns[1:]]
    df.columns = ["gene_id"] + sample_ids

    # Verify counts are numeric non-negative integers
    for sid in sample_ids:
        df[sid] = pd.to_numeric(df[sid], errors="coerce")
        if df[sid].isna().any():
            raise ValueError(f"Sample column '{sid}' contains missing or non-numeric count values.")
        if (df[sid] < 0).any():
            raise ValueError(f"Sample column '{sid}' contains negative count values.")
        df[sid] = df[sid].round().astype(int)

    return df, sample_ids


def parse_sample_info(
    file_or_content: Any | None,
    matrix_sample_ids: list[str],
) -> pd.DataFrame:
    """Parse or automatically infer sample condition and biological replicate metadata."""
    if file_or_content is not None:
        if isinstance(file_or_content, pd.DataFrame):
            df = file_or_content.copy()
        elif isinstance(file_or_content, (str, Path)) and os.path.exists(str(file_or_content)):
            path = str(file_or_content)
            sep = "\t" if path.endswith((".tsv", ".tab", ".txt")) else ","
            df = pd.read_csv(path, sep=sep)
        elif hasattr(file_or_content, "read"):
            content = file_or_content.read()
            if isinstance(content, bytes):
                content = content.decode("utf-8-sig")
            sep = "\t" if ("\t" in content and "," not in content.split("\n")[0]) else ","
            df = pd.read_csv(io.StringIO(content), sep=sep)
        elif isinstance(file_or_content, str):
            sep = "\t" if ("\t" in file_or_content and "," not in file_or_content.split("\n")[0]) else ","
            df = pd.read_csv(io.StringIO(file_or_content), sep=sep)
        else:
            df = None

        if df is not None and not df.empty:
            # Reconcile columns
            col_map = {str(c).strip().lower(): str(c).strip() for c in df.columns}
            sid_col = col_map.get("sample_id") or col_map.get("sample") or df.columns[0]
            cond_col = (
                col_map.get("condition_type")
                or col_map.get("condition")
                or col_map.get("group")
                or col_map.get("status")
            )
            rep_col = (
                col_map.get("biological_replicate_id")
                or col_map.get("replicate_id")
                or col_map.get("replicate")
                or col_map.get("rep")
            )

            df = df.rename(columns={sid_col: "sample_id"})
            df["sample_id"] = df["sample_id"].astype(str).str.strip()

            # Ensure all matrix samples are in metadata
            meta_sids = set(df["sample_id"])
            if all(s in meta_sids for s in matrix_sample_ids):
                # Filter to matrix samples in exact matrix order
                df = df.set_index("sample_id").loc[matrix_sample_ids].reset_index()

                if cond_col and cond_col in df.columns:
                    # Normalize condition strings to 'control' or 'treatment'
                    def norm_cond(val: Any) -> str:
                        s = str(val).strip().lower()
                        if any(w in s for w in ("ctrl", "control", "veh", "dmso", "untreated", "mock", "baseline")):
                            return "control"
                        return "treatment"

                    df["condition_type"] = df[cond_col].apply(norm_cond)
                else:
                    df["condition_type"] = [_infer_sample_condition(s) for s in df["sample_id"]]

                if rep_col and rep_col in df.columns:
                    df["biological_replicate_id"] = df[rep_col].astype(str).str.strip()
                else:
                    df["biological_replicate_id"] = [f"rep_{i+1}" for i in range(len(df))]

                return df[["sample_id", "condition_type", "biological_replicate_id"]]

    # Auto-infer condition and replicate assignments from sample names
    records = []
    ctrl_count = 0
    trt_count = 0

    for sid in matrix_sample_ids:
        cond = _infer_sample_condition(sid)
        if cond == "control":
            ctrl_count += 1
            rep = f"rep_ctrl_{ctrl_count}"
        else:
            trt_count += 1
            rep = f"rep_trt_{trt_count}"
        records.append({"sample_id": sid, "condition_type": cond, "biological_replicate_id": rep})

    # If all inferred as one class, split evenly (first half control, second half treatment)
    all_ctrl = all(r["condition_type"] == "control" for r in records)
    all_trt = all(r["condition_type"] == "treatment" for r in records)
    if (all_ctrl or all_trt) and len(records) >= 4:
        mid = len(records) // 2
        for i, r in enumerate(records):
            if i < mid:
                r["condition_type"] = "control"
                r["biological_replicate_id"] = f"rep_ctrl_{i+1}"
            else:
                r["condition_type"] = "treatment"
                r["biological_replicate_id"] = f"rep_trt_{i - mid + 1}"

    return pd.DataFrame(records)


def _infer_sample_condition(sample_id: str) -> str:
    """Infer control vs treatment from sample name substrings."""
    s = sample_id.lower()
    if any(w in s for w in ("ctrl", "control", "dmso", "veh", "vehicle", "mock", "baseline", "untreated", "c_")):
        return "control"
    return "treatment"


def build_canonical_metadata(
    samples_df: pd.DataFrame,
    chemical_params: ChemicalParams,
    asset_id: str = "asset_app_01",
    assay_id: str = "assay_app_01",
    gene_identifier_type: str = "ensembl_gene_id",
) -> dict[str, Any]:
    """Construct a full canonical metadata dictionary adhering to schema and backend stage rules."""
    sample_records = []
    for _, row in samples_df.iterrows():
        sid = str(row["sample_id"]).strip()
        ctype = str(row["condition_type"]).strip().lower()
        cid = "cond_treatment" if ctype == "treatment" else "cond_control"
        brep = str(row["biological_replicate_id"]).strip()

        sample_records.append(
            {
                "sample_id": sid,
                "experiment_id": "exp_app_01",
                "condition_id": cid,
                "biological_replicate_id": brep,
                "biological_source_id": chemical_params.biological_source_id,
                "organoid_context_id": chemical_params.organoid_context_id,
                "developmental_age_or_stage_at_collection": chemical_params.developmental_age,
                "developmental_age_or_stage_at_exposure": chemical_params.developmental_age,
            }
        )

    all_sids = [s["sample_id"] for s in sample_records]

    metadata: dict[str, Any] = {
        "studies": [
            {
                "study_id": "study_app_01",
                "source_type": "in_house",
                "metadata_validation_status": "valid",
            }
        ],
        "experiments": [
            {
                "experiment_id": "exp_app_01",
                "study_id": "study_app_01",
                "metadata_validation_status": "valid",
            }
        ],
        "conditions": [
            {
                "condition_id": "cond_treatment",
                "experiment_id": "exp_app_01",
                "treatment_control_status": "treatment",
                "developmental_age_or_stage": chemical_params.developmental_age,
            },
            {
                "condition_id": "cond_control",
                "experiment_id": "exp_app_01",
                "treatment_control_status": "control",
                "developmental_age_or_stage": chemical_params.developmental_age,
            },
        ],
        "samples": sample_records,
        "organoid_contexts": [
            {
                "organoid_context_id": chemical_params.organoid_context_id,
                "organoid_type": chemical_params.organoid_type,
                "brain_region": chemical_params.brain_region,
            }
        ],
        "biological_sources": [
            {
                "biological_source_id": chemical_params.biological_source_id,
                "species": chemical_params.species,
            }
        ],
        "exposures": [
            {
                "condition_id": "cond_treatment",
                "agent_name": chemical_params.agent_name,
                "agent_identifier": chemical_params.agent_identifier,
                "concentration_or_dose": chemical_params.concentration_or_dose,
                "concentration_or_dose_unit": chemical_params.concentration_or_dose_unit,
                "exposure_duration": chemical_params.exposure_duration,
                "exposure_duration_unit": chemical_params.exposure_duration_unit,
                "vehicle": chemical_params.vehicle,
            }
        ],
        "sequencing_assays": [
            {
                "assay_id": assay_id,
                "sample_ids": all_sids,
                "rna_seq_modality": "bulk_rna_seq",
                "library_strategy": "polyA",
                "sequencing_platform": "Illumina NovaSeq 6000",
                "strandedness": "reverse",
                "read_layout": "paired_end",
            }
        ],
        "treatment_control_relationships": [
            {
                "relationship_id": "rel_app_01",
                "treatment_condition_ids": ["cond_treatment"],
                "matched_control_condition_ids": ["cond_control"],
            }
        ],
        "input_data_assets": [
            {
                "asset_id": asset_id,
                "assay_ids": [assay_id],
                "gene_identifier_type": gene_identifier_type,
            }
        ],
    }

    return metadata


def generate_demo_data() -> tuple[pd.DataFrame, pd.DataFrame, ChemicalParams]:
    """Generate realistic demonstration human neural organoid RNA-seq data (Rotenone 0.5 uM).

    Contains 40 authentic human genes covering:
    - Mitochondrial and ER stress markers (upregulated)
    - Neuronal differentiation and synaptic transmission markers (downregulated)
    - Canonical housekeeping genes (unchanged)
    """
    sids = [
        "ctrl_organoid_rep1",
        "ctrl_organoid_rep2",
        "ctrl_organoid_rep3",
        "rotenone_0.5uM_rep1",
        "rotenone_0.5uM_rep2",
        "rotenone_0.5uM_rep3",
    ]

    # Authentic human genes from Ensembl 112
    # 12 Upregulated stress / apoptosis markers
    up_genes = [
        "ENSG00000100292",  # HMOX1 (Heme oxygenase 1)
        "ENSG00000128272",  # ATF4 (Activating transcription factor 4)
        "ENSG00000175197",  # DDIT3 (CHOP)
        "ENSG00000001084",  # GCLC (Glutamate-cysteine ligase)
        "ENSG00000181019",  # NQO1 (NAD(P)H dehydrogenase)
        "ENSG00000044574",  # HSPA5 (BiP / GRP78)
        "ENSG00000116717",  # GADD45A (Growth arrest and DNA damage)
        "ENSG00000087088",  # BAX (Apoptosis regulator BAX)
        "ENSG00000164305",  # CASP3 (Caspase 3)
        "ENSG00000124762",  # CDKN1A (p21)
        "ENSG00000161011",  # SQSTM1 (p62)
        "ENSG00000291237",  # SOD2 (Superoxide dismutase 2)
    ]

    # 14 Downregulated synaptic and mitochondrial markers
    down_genes = [
        "ENSG00000090266",  # NDUFS1 (Complex I subunit)
        "ENSG00000167792",  # NDUFV1 (Complex I subunit)
        "ENSG00000078018",  # MAP2 (Microtubule associated protein 2)
        "ENSG00000258947",  # TUBB3 (Beta-III tubulin)
        "ENSG00000008056",  # SYN1 (Synapsin 1)
        "ENSG00000132639",  # SNAP25 (Synaptosome associated protein)
        "ENSG00000167281",  # RBFOX3 (NeuN)
        "ENSG00000102003",  # SYP (Synaptophysin)
        "ENSG00000106089",  # STX1A (Syntaxin 1A)
        "ENSG00000176884",  # GRIN1 (NMDA receptor 1)
        "ENSG00000155511",  # GRIA1 (AMPA receptor 1)
        "ENSG00000070808",  # CAMK2A (CaMKII alpha)
        "ENSG00000001617",  # SEMA3F (Semaphorin 3F)
        "ENSG00000001631",  # KRIT1 (KRIT1 ankyrin repeat)
    ]

    # 14 Stable housekeeping genes
    housekeeping_genes = [
        "ENSG00000111640",  # GAPDH
        "ENSG00000075624",  # ACTB
        "ENSG00000166710",  # B2M
        "ENSG00000089157",  # RPLP0
        "ENSG00000150991",  # UBC
        "ENSG00000196262",  # PPIA
        "ENSG00000112592",  # TBP
        "ENSG00000168488",  # HPRT1
        "ENSG00000102144",  # PGK1
        "ENSG00000164924",  # YWHAZ
        "ENSG00000073578",  # SDHA
        "ENSG00000144713",  # RPS18
        "ENSG00000000003",  # TSPAN6
        "ENSG00000000419",  # DPM1
    ]

    all_genes = up_genes + down_genes + housekeeping_genes
    mean_levels = [40, 80, 150, 300, 600, 1200, 2500, 5000]

    rng = random.Random(2026)
    rows = []

    for i, gid in enumerate(all_genes):
        base_mu = mean_levels[i % len(mean_levels)]
        if gid in up_genes:
            fc = 3.8
        elif gid in down_genes:
            fc = 0.28
        else:
            fc = 1.0

        row = [gid]
        for sid in sids:
            is_trt = "rotenone" in sid
            target_mu = base_mu * (fc if is_trt else 1.0)
            # Add realistic biological replicate variation
            val = max(10, int(round(rng.gauss(target_mu, max(3.0, target_mu * 0.12)))))
            row.append(val)
        rows.append(row)

    counts_df = pd.DataFrame(rows, columns=["gene_id"] + sids)

    samples_records = []
    for sid in sids:
        ctype = "treatment" if "rotenone" in sid else "control"
        rep_num = sid.split("_rep")[-1]
        rep_id = f"rep_{'trt' if ctype == 'treatment' else 'ctrl'}_{rep_num}"
        samples_records.append(
            {
                "sample_id": sid,
                "condition_type": ctype,
                "biological_replicate_id": rep_id,
            }
        )
    samples_df = pd.DataFrame(samples_records)

    params = ChemicalParams(
        agent_name="Rotenone",
        agent_identifier="CID:6758",
        concentration_or_dose=0.5,
        concentration_or_dose_unit="uM",
        developmental_age="day_35",
        exposure_duration=24.0,
        exposure_duration_unit="h",
        vehicle="0.1% DMSO",
        organoid_context_id="ctx_cerebral_cortex_day35",
        organoid_type="cerebral_organoid",
        brain_region="cerebral_cortex",
        biological_source_id="ipsc_donor_line_1",
        species="Homo sapiens",
    )

    return counts_df, samples_df, params


def run_full_analysis_pipeline(
    counts_input: Any,
    samples_input: Any | None,
    chemical_params: ChemicalParams,
    r_binary_path: str = "Rscript",
    strict_version_check: bool = True,
    r_timeout_seconds: int = 180,
) -> PipelineRunResult:
    """Execute the complete end-to-end differential expression analysis pipeline."""
    temp_matrix_path: str | None = None
    try:
        # 1. Parse and validate count matrix
        counts_df, sample_ids = parse_count_matrix(counts_input)

        # 2. Parse or infer sample table
        samples_df = parse_sample_info(samples_input, sample_ids)

        # Verify replicate counts
        trt_samples = samples_df[samples_df["condition_type"] == "treatment"]
        ctrl_samples = samples_df[samples_df["condition_type"] == "control"]
        if len(trt_samples) < 2 or len(ctrl_samples) < 2:
            return PipelineRunResult(
                is_success=False,
                contrast_result=None,
                de_genes_df=pd.DataFrame(),
                excluded_genes_df=pd.DataFrame(),
                samples_df=samples_df,
                summary_metrics={},
                scientific_findings=[],
                error_message=(
                    f"Differential expression requires at least 2 biological replicates in both treatment and control. "
                    f"Found: {len(trt_samples)} treatment, {len(ctrl_samples)} control."
                ),
            )

        # Detect gene identifier type
        first_few_genes = counts_df["gene_id"].head(20).tolist()
        ensembl_count = sum(1 for g in first_few_genes if str(g).startswith("ENSG"))
        gene_id_type = "ensembl_gene_id" if (ensembl_count / max(1, len(first_few_genes)) > 0.5) else "gene_symbol"

        # 3. Build canonical metadata
        asset_id = "asset_app_upload"
        assay_id = "assay_app_bulk"
        metadata = build_canonical_metadata(
            samples_df=samples_df,
            chemical_params=chemical_params,
            asset_id=asset_id,
            assay_id=assay_id,
            gene_identifier_type=gene_id_type,
        )

        # Write counts to temporary CSV for QC and Harmonization
        with tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, newline="") as tf:
            counts_df.to_csv(tf, index=False)
            temp_matrix_path = tf.name

        # 4. Step 6A: Bulk Raw-Count QC
        qc_result = validate_bulk_counts(temp_matrix_path, metadata, asset_id, assay_id)
        if qc_result.status == QCStatus.FAIL:
            error_msgs = [f.message for f in qc_result.findings if f.severity.value == "ERROR"]
            return PipelineRunResult(
                is_success=False,
                contrast_result=None,
                de_genes_df=pd.DataFrame(),
                excluded_genes_df=pd.DataFrame(),
                samples_df=samples_df,
                summary_metrics={"qc_status": str(qc_result.status)},
                scientific_findings=[{"severity": f.severity.value, "rule": f.rule_id, "message": f.message} for f in qc_result.findings],
                error_message=f"Dataset failed bulk raw-count QC: {'; '.join(error_msgs)}",
                qc_result=qc_result,
            )

        # 5. Step 6B: Gene Harmonization v1
        gene_ref = get_gene_reference()
        harm_result = harmonize_bulk_counts(
            matrix=temp_matrix_path,
            metadata=metadata,
            asset_id=asset_id,
            assay_id=assay_id,
            qc_result=qc_result,
            reference=gene_ref,
        )
        if harm_result.status.value == "FAIL":
            error_msgs = [f.message for f in harm_result.findings if f.severity.value == "ERROR"]
            return PipelineRunResult(
                is_success=False,
                contrast_result=None,
                de_genes_df=pd.DataFrame(),
                excluded_genes_df=pd.DataFrame(),
                samples_df=samples_df,
                summary_metrics={"harmonization_status": str(harm_result.status)},
                scientific_findings=[{"severity": f.severity.value, "rule": f.rule_id, "message": f.message} for f in harm_result.findings],
                error_message=f"Gene Harmonization failed: {'; '.join(error_msgs)}",
                qc_result=qc_result,
                harmonized_dataset=harm_result,
            )

        # 6. Step 6C: Bulk Normalization v1 (Locked DESeq2 Size Factors)
        norm_result = normalize_bulk_dataset(harm_result, metadata)
        if norm_result.status.value == "FAIL":
            error_msgs = [f.message for f in norm_result.findings if f.severity.value == "ERROR"]
            return PipelineRunResult(
                is_success=False,
                contrast_result=None,
                de_genes_df=pd.DataFrame(),
                excluded_genes_df=pd.DataFrame(),
                samples_df=samples_df,
                summary_metrics={"normalization_status": str(norm_result.status)},
                scientific_findings=[{"severity": f.severity.value, "rule": f.rule_id, "message": f.message} for f in norm_result.findings],
                error_message=f"Bulk Normalization failed: {'; '.join(error_msgs)}",
                qc_result=qc_result,
                harmonized_dataset=harm_result,
                normalized_dataset=norm_result,
            )

        # 7. Step 7A: Treatment-versus-Matched-Control Contrast Builder v1
        contrast_dataset = build_response_contrasts(norm_result, metadata)
        if not contrast_dataset.contrasts:
            return PipelineRunResult(
                is_success=False,
                contrast_result=None,
                de_genes_df=pd.DataFrame(),
                excluded_genes_df=pd.DataFrame(),
                samples_df=samples_df,
                summary_metrics={"contrast_status": "NO_CONTRASTS_BUILT"},
                scientific_findings=[{"severity": f.severity.value, "rule": f.rule_id, "message": f.message} for f in contrast_dataset.findings],
                error_message="No treatment-control contrasts could be constructed from metadata.",
                qc_result=qc_result,
                harmonized_dataset=harm_result,
                normalized_dataset=norm_result,
                contrast_dataset=contrast_dataset,
            )

        primary_contrast = contrast_dataset.contrasts[0]
        if primary_contrast.status != ContrastStatus.ELIGIBLE:
            error_msgs = [f.message for f in primary_contrast.findings if f.severity.value in ("ERROR", "REVIEW")]
            return PipelineRunResult(
                is_success=False,
                contrast_result=None,
                de_genes_df=pd.DataFrame(),
                excluded_genes_df=pd.DataFrame(),
                samples_df=samples_df,
                summary_metrics={"contrast_status": str(primary_contrast.status)},
                scientific_findings=[{"severity": f.severity.value, "rule": f.rule_id, "message": f.message} for f in primary_contrast.findings],
                error_message=f"Contrast is not eligible for differential expression: {'; '.join(error_msgs)}",
                qc_result=qc_result,
                harmonized_dataset=harm_result,
                normalized_dataset=norm_result,
                contrast_dataset=contrast_dataset,
            )

        # 8. Step 7B: Bulk Differential Expression v1 (DESeq2 Wald test)
        de_dataset = run_differential_expression(
            contrast_dataset=contrast_dataset,
            normalized_dataset=norm_result,
            harmonized_dataset=harm_result,
            sample_metadata=metadata["samples"],
            r_binary_path=r_binary_path,
            strict_version_check=strict_version_check,
            r_timeout_seconds=r_timeout_seconds,
        )

        cr = de_dataset.contrast_results[0]
        if not cr.is_inferentially_eligible or cr.status != ContrastStatus.ELIGIBLE:
            error_msgs = [f.message for f in cr.findings if f.severity.value == "ERROR"]
            return PipelineRunResult(
                is_success=False,
                contrast_result=cr,
                de_genes_df=pd.DataFrame(),
                excluded_genes_df=pd.DataFrame(),
                samples_df=samples_df,
                summary_metrics={"de_status": str(cr.status)},
                scientific_findings=[{"severity": f.severity.value, "rule": f.rule_id, "message": f.message} for f in cr.findings],
                error_message=f"Differential expression analysis failed or was blocked: {'; '.join(error_msgs)}",
                qc_result=qc_result,
                harmonized_dataset=harm_result,
                normalized_dataset=norm_result,
                contrast_dataset=contrast_dataset,
                de_dataset=de_dataset,
            )

        # 9. Format output DataFrames and Metrics
        de_gene_records = []
        for g in cr.gene_results:
            is_sig = bool(g.adjusted_p_value_bh is not None and g.adjusted_p_value_bh < 0.05)
            direction = "Not Significant"
            if is_sig and g.log2_fold_change is not None:
                direction = "Upregulated" if g.log2_fold_change > 0 else "Downregulated"

            de_gene_records.append(
                {
                    "canonical_gene_id": g.canonical_gene_id,
                    "approved_symbol": g.approved_symbol or g.canonical_gene_id,
                    "base_mean": g.base_mean,
                    "log2_fold_change": g.log2_fold_change,
                    "lfc_standard_error": g.lfc_standard_error,
                    "wald_statistic": g.wald_statistic,
                    "p_value": g.p_value,
                    "adjusted_p_value_bh": g.adjusted_p_value_bh,
                    "significance": "Significant (FDR < 0.05)" if is_sig else "Not Significant",
                    "direction": direction,
                    "result_status": g.result_status,
                }
            )

        de_genes_df = pd.DataFrame(de_gene_records)
        if not de_genes_df.empty:
            de_genes_df = de_genes_df.sort_values(
                by=["adjusted_p_value_bh", "p_value"],
                ascending=[True, True],
                na_position="last",
            ).reset_index(drop=True)

        excluded_records = [e.to_dict() for e in cr.excluded_genes]
        excluded_genes_df = pd.DataFrame(excluded_records)

        # Enriched samples dataframe with size factors and library sizes
        enriched_samples = []
        sf_dict = cr.size_factors_used
        for _, s in samples_df.iterrows():
            sid = s["sample_id"]
            enriched_samples.append(
                {
                    "sample_id": sid,
                    "condition_type": s["condition_type"],
                    "biological_replicate_id": s["biological_replicate_id"],
                    "deseq2_size_factor": sf_dict.get(sid, 1.0),
                    "raw_library_size": int(counts_df[sid].sum()),
                }
            )
        samples_result_df = pd.DataFrame(enriched_samples)

        sig_up_count = sum(1 for g in cr.gene_results if g.adjusted_p_value_bh is not None and g.adjusted_p_value_bh < 0.05 and g.log2_fold_change is not None and g.log2_fold_change > 0)
        sig_down_count = sum(1 for g in cr.gene_results if g.adjusted_p_value_bh is not None and g.adjusted_p_value_bh < 0.05 and g.log2_fold_change is not None and g.log2_fold_change < 0)

        summary_metrics = {
            "contrast_id": cr.contrast_id,
            "status": str(cr.status),
            "treatment_replicates": cr.treatment_replicate_count,
            "control_replicates": cr.control_replicate_count,
            "total_input_genes": cr.total_input_genes_count,
            "evaluated_canonical_genes": cr.eligible_canonical_genes_count,
            "excluded_genes_count": cr.excluded_genes_count,
            "is_gene_universe_conserved": cr.is_gene_universe_conserved,
            "significant_upregulated": sig_up_count,
            "significant_downregulated": sig_down_count,
            "design_formula": cr.design_formula,
            "residual_degrees_of_freedom": cr.residual_degrees_of_freedom,
            "r_version": cr.r_environment_info.get("r_version", "4.4.3"),
            "bioc_version": cr.r_environment_info.get("bioc_version", "3.20"),
            "deseq2_version": cr.r_environment_info.get("deseq2_version", "1.46.0"),
        }

        all_findings = [
            {"severity": f.severity.value, "rule": f.rule_id, "message": f.message}
            for f in cr.findings
        ]

        return PipelineRunResult(
            is_success=True,
            contrast_result=cr,
            de_genes_df=de_genes_df,
            excluded_genes_df=excluded_genes_df,
            samples_df=samples_result_df,
            summary_metrics=summary_metrics,
            scientific_findings=all_findings,
            qc_result=qc_result,
            harmonized_dataset=harm_result,
            normalized_dataset=norm_result,
            contrast_dataset=contrast_dataset,
            de_dataset=de_dataset,
        )

    except Exception as exc:
        return PipelineRunResult(
            is_success=False,
            contrast_result=None,
            de_genes_df=pd.DataFrame(),
            excluded_genes_df=pd.DataFrame(),
            samples_df=pd.DataFrame(),
            summary_metrics={},
            scientific_findings=[],
            error_message=f"Pipeline execution encountered an unexpected error: {str(exc)}",
        )
    finally:
        if temp_matrix_path and os.path.exists(temp_matrix_path):
            try:
                os.remove(temp_matrix_path)
            except OSError:
                pass
