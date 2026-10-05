"""DNT Organoid AI Web Interface: Treatment-Response Explorer.

Statistically rigorous Bulk RNA-seq differential expression analysis web interface.
Architecture:
Browser / Web Interface
        ↓
Existing scientific backend
        ↓
QC → Gene Harmonization → Normalization → Treatment-Control Contrast → DESeq2
        ↓
Molecular Response displayed in the browser (Gene table + visual plot + scientific checks + CSV downloads)

Strictly real-data analysis website. Only one count-matrix upload required.
Sample assignments are entered explicitly by the user in the editable table.
Automatic inference of treatment/control status and biological replicates from sample names is prohibited.
"""

import io
import json
import math
from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from src.app.geo_resolver import fetch_geo_metadata, resolve_geo_samples
from src.app.pipeline_runner import (
    ChemicalParams,
    PipelineRunResult,
    parse_count_matrix,
    parse_sample_info,
    prepare_plot_data,
    run_full_analysis_pipeline,
)

# -----------------------------------------------------------------------------
# Streamlit Page Configuration & Styling
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="DNT Organoid AI Web Interface | Treatment-Response Explorer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    .main-title {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E3A8A;
        margin-bottom: 0.2rem;
    }
    .subtitle {
        font-size: 1.05rem;
        color: #4B5563;
        margin-bottom: 1.5rem;
    }
    .kpi-card {
        background-color: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        padding: 16px;
        text-align: center;
    }
    .kpi-value {
        font-size: 1.8rem;
        font-weight: 700;
        color: #1E293B;
    }
    .kpi-label {
        font-size: 0.85rem;
        color: #64748B;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        margin-top: 4px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# -----------------------------------------------------------------------------
# Header
# -----------------------------------------------------------------------------
st.markdown('<div class="main-title">🧬 DNT Organoid AI Web Interface: Treatment-Response Explorer</div>', unsafe_allow_html=True)
st.markdown(
    '<div class="subtitle">'
    "Statistically rigorous human developmental neurotoxicity (DNT) organoid treatment-versus-matched-control response analysis website. "
    "Orchestrates locked scientific backend: Raw-Count QC → Ensembl 112 Harmonization → DESeq2 Normalization → Step 7A Contrast Gate → Step 7B DESeq2 Wald Test."
    "</div>",
    unsafe_allow_html=True,
)

# -----------------------------------------------------------------------------
# Sidebar: File Uploads & Exposure Inputs
# -----------------------------------------------------------------------------
with st.sidebar:
    st.header("⚙️ Data Ingestion & Setup")

    st.subheader("1. Count Matrix (Required)")
    uploaded_counts = st.file_uploader(
        "Upload raw count matrix (.csv, .tsv, .txt, .gz)",
        type=["csv", "tsv", "txt", "gz"],
        help="Column 1: Gene identifier (Ensembl gene ID or gene symbol). Remaining columns: Sample IDs containing non-negative raw integer counts.",
    )

    counts_df: pd.DataFrame | None = None
    sample_ids: list[str] = []
    excluded_annotations: list[str] = []
    matrix_parse_error: str | None = None

    if uploaded_counts is not None:
        try:
            counts_df, sample_ids, excluded_annotations = parse_count_matrix(uploaded_counts, return_excluded_columns=True)
            if excluded_annotations:
                st.info(f"ℹ️ Recognized and excluded {len(excluded_annotations)} annotation column(s): {', '.join(excluded_annotations)}. Sample count columns preserved.")
        except Exception as e:
            matrix_parse_error = str(e)
            counts_df, sample_ids, excluded_annotations = None, [], []

    st.subheader("2. Dataset Provenance & Resolution (Optional)")
    source_acc = st.text_input(
        "Source Accession / Public ID",
        value=st.session_state.get("source_accession_val", ""),
        placeholder="e.g., GSE232218",
        help="Enter a public GEO Series accession (e.g. GSE232218) to automatically resolve sample records and match conditions.",
        key="source_accession_input",
    )
    if source_acc.strip():
        st.session_state["source_accession_val"] = source_acc.strip()

    resolve_geo_clicked = st.button("🔍 Resolve GEO Metadata", use_container_width=True)

    if resolve_geo_clicked:
        if not sample_ids:
            st.warning("Please upload a count matrix before resolving GEO metadata.")
        elif not source_acc.strip():
            st.warning("Please enter a GEO Series accession (e.g., GSE232218).")
        else:
            try:
                with st.spinner(f"Querying NCBI GEO for {source_acc.strip()}..."):
                    geo_data = fetch_geo_metadata(source_acc.strip())
                    assignments, suggestions = resolve_geo_samples(
                        matrix_sample_ids=sample_ids,
                        geo_data=geo_data,
                        target_agent=st.session_state.get("chem_name_input", "") or "VPA",
                        target_dose=float(st.session_state.get("chem_dose_input", 200.0)),
                        target_dose_unit=st.session_state.get("chem_unit_input", "uM"),
                        vehicle_keyword=st.session_state.get("chem_veh_input", "") or "vehicle",
                    )
                    st.session_state["sample_assignment_init"] = assignments[["sample_id", "condition_type", "biological_replicate_id"]]
                    if "sample_assignment_editor" in st.session_state:
                        del st.session_state["sample_assignment_editor"]
                    if suggestions.get("developmental_age"):
                        st.session_state["dev_age_val"] = suggestions["developmental_age"]
                    if suggestions.get("biological_source"):
                        st.session_state["bio_src_val"] = suggestions["biological_source"]
                    if suggestions.get("exposure_duration") is not None:
                        st.session_state["exp_dur_val"] = float(suggestions["exposure_duration"])
                    if suggestions.get("exposure_duration_unit"):
                        st.session_state["exp_dur_unit_val"] = suggestions["exposure_duration_unit"]
                    if suggestions.get("vehicle"):
                        st.session_state["chem_veh_val"] = suggestions["vehicle"]
                    if suggestions.get("organoid_context"):
                        st.session_state["organoid_ctx_val"] = suggestions["organoid_context"]
                    if suggestions.get("duration_citation"):
                        st.session_state["duration_citation"] = suggestions["duration_citation"]
                    if suggestions.get("duration_status"):
                        st.session_state["duration_status"] = suggestions["duration_status"]
                    trt_n = (assignments["condition_type"] == "treatment").sum()
                    ctrl_n = (assignments["condition_type"] == "control").sum()
                    excl_n = (assignments["condition_type"] == "exclude").sum()
                    st.success(f"✅ Resolved GEO metadata: {trt_n} treatment, {ctrl_n} control, {excl_n} excluded. Please review below.")
            except Exception as exc:
                st.error(f"❌ Failed to resolve GEO metadata: {exc}")

    st.subheader("3. Sample Assignments (Required)")
    if matrix_parse_error:
        st.error(f"❌ Error reading count matrix: {matrix_parse_error}")
        edited_samples = None
    elif sample_ids:
        current_matrix_id = getattr(uploaded_counts, "file_id", getattr(uploaded_counts, "name", str(sample_ids)))
        if st.session_state.get("last_uploaded_matrix_id") != current_matrix_id:
            st.session_state["last_uploaded_matrix_id"] = current_matrix_id
            st.session_state["sample_assignment_init"] = pd.DataFrame(
                {
                    "sample_id": sample_ids,
                    "condition_type": ["" for _ in sample_ids],
                    "biological_replicate_id": ["" for _ in sample_ids],
                }
            )
            if "sample_assignment_editor" in st.session_state:
                del st.session_state["sample_assignment_editor"]

        st.caption("Assign treatment/control status ('control', 'treatment', or 'exclude') and biological replicate ID for each sample:")
        edited_samples = st.data_editor(
            st.session_state["sample_assignment_init"],
            key="sample_assignment_editor",
            column_config={
                "sample_id": st.column_config.TextColumn(
                    "sample_id",
                    disabled=True,
                    help="Sample ID extracted from count matrix column",
                ),
                "condition_type": st.column_config.SelectboxColumn(
                    "condition_type",
                    options=["", "control", "treatment", "exclude"],
                    required=True,
                    help="Assign 'control', 'treatment', or 'exclude' for non-contrast samples",
                ),
                "biological_replicate_id": st.column_config.TextColumn(
                    "biological_replicate_id",
                    required=False,
                    help="Explicit biological replicate identifier (e.g. rep_1, rep_2)",
                ),
            },
            hide_index=True,
            use_container_width=True,
            num_rows="fixed",
        )
    else:
        st.info("Upload a raw count matrix above to configure sample assignments.")
        edited_samples = None

    st.subheader("4. Exposure & Experimental Design")
    with st.expander("🔬 Exposure Parameters", expanded=True):
        chem_name = st.text_input(
            "Chemical / Agent Name *",
            value=st.session_state.get("chem_name_val", "Valproic Acid" if "GSE232218" in source_acc.upper() else ""),
            placeholder="e.g., Valproic Acid",
            help="Name of the chemical compound or stressor tested.",
            key="chem_name_input",
        )
        chem_id = st.text_input(
            "Agent Identifier / CID",
            value="CID:3121" if "GSE232218" in source_acc.upper() else "",
            placeholder="e.g., CID:3121",
            help="PubChem CID, CAS registry number, or standard chemical identifier.",
            key="chem_id_input",
        )

        col_dose1, col_dose2 = st.columns([2, 1])
        with col_dose1:
            chem_dose = st.number_input(
                "Concentration / Dose",
                value=200.0 if "GSE232218" in source_acc.upper() else 10.0,
                step=0.1,
                format="%.2f",
                key="chem_dose_input",
            )
        with col_dose2:
            chem_unit = st.selectbox("Unit", options=["uM", "nM", "mM", "ug/mL", "mg/kg", "ppm"], index=0, key="chem_unit_input")

        col_age1, col_dur = st.columns([1, 1])
        with col_age1:
            dev_age = st.text_input(
                "Developmental Age *",
                value=st.session_state.get("dev_age_val", "Day 6" if "GSE232218" in source_acc.upper() else ""),
                placeholder="e.g., Day 6",
                help="Strictly matches developmental stage between treatment and control (e.g. 'Day 6').",
                key="dev_age_input",
            )
        with col_dur:
            default_dur = st.session_state.get(
                "exp_dur_val",
                48.0 if "GSE232218" in source_acc.upper() else 0.0,
            )
            exp_dur = st.number_input(
                "Duration *",
                value=float(default_dur),
                step=1.0,
                min_value=0.0,
                help="Duration of exposure (e.g. 48.0 h). Must be verified from official records.",
                key="exp_dur_input",
            )

        col_dur_unit, col_veh = st.columns([1, 1])
        with col_dur_unit:
            unit_opts = ["h", "days", "weeks"]
            default_unit = st.session_state.get("exp_dur_unit_val", "h")
            unit_idx = unit_opts.index(default_unit) if default_unit in unit_opts else 0
            exp_dur_unit = st.selectbox("Duration Unit", options=unit_opts, index=unit_idx, key="exp_dur_unit_input")
        with col_veh:
            chem_vehicle = st.text_input(
                "Vehicle",
                value=st.session_state.get("chem_veh_val", "vehicle control" if "GSE232218" in source_acc.upper() else ""),
                placeholder="e.g., vehicle control",
                key="chem_veh_input",
            )

        dur_citation = st.session_state.get("duration_citation")
        dur_status = st.session_state.get("duration_status")
        if dur_citation:
            if dur_status == "verified_publication":
                st.caption(f"📚 **Verified Official Publication Evidence**: {dur_citation}")
            elif dur_status == "verified_geo_sample":
                st.caption(f"📋 **Verified GEO Sample Evidence**: {dur_citation}")
            else:
                st.warning("⚠️ **Unresolved Duration**: Exposure duration was not declared in GEO sample records. Please enter verified experimental duration.")
        elif exp_dur == 0.0:
            st.warning("⚠️ **Duration Required**: Exposure duration cannot be assumed. Please enter verified experimental duration (e.g. 48.0 h).")

        organoid_ctx = st.text_input(
            "Organoid Context",
            value=st.session_state.get("organoid_ctx_val", "Cerebral Cortex Organoid (SOSR-CO)" if "GSE232218" in source_acc.upper() else "Cerebral Cortex Organoid"),
            placeholder="e.g., Cerebral Cortex Organoid",
            key="organoid_ctx_input",
        )
        brain_reg = st.text_input("Brain Region", value="cerebral_cortex", placeholder="e.g., cerebral_cortex", key="brain_reg_input")
        bio_src = st.text_input(
            "Biological Source / Line",
            value=st.session_state.get("bio_src_val", "iPSC (AICS-0023)" if "GSE232218" in source_acc.upper() else ""),
            placeholder="e.g., iPSC donor line",
            key="bio_src_input",
        )

    analyze_button = st.button("🚀 ANALYZE", type="primary", use_container_width=True)

    st.markdown("---")
    st.caption("🔒 **Locked Scientific Architecture**")
    st.caption("• **R Runtime**: R 4.4.3 / Bioc 3.20 / DESeq2 1.46.0")
    st.caption("• **Reference**: Ensembl Human Release 112")
    st.caption("• **Replication**: ≥ 2 Biological Replicates per arm")
    st.caption("• **Inference**: Negative-binomial Wald Test")

# -----------------------------------------------------------------------------
# Pipeline Execution Trigger
# -----------------------------------------------------------------------------
if "pipeline_result" not in st.session_state:
    st.session_state["pipeline_result"] = None

if analyze_button:
    validation_errors = []
    if uploaded_counts is None:
        validation_errors.append("Please upload an RNA-seq count matrix (.csv, .tsv, .txt, .gz) file.")
    elif matrix_parse_error:
        validation_errors.append(f"Invalid count matrix: {matrix_parse_error}")
    elif edited_samples is None or edited_samples.empty:
        validation_errors.append("No samples detected in count matrix. Please verify uploaded file.")
    else:
        # Validate sample assignments
        active_rows = []
        for _, row in edited_samples.iterrows():
            sid = str(row.get("sample_id", "")).strip()
            cond = str(row.get("condition_type", "")).strip().lower() if pd.notna(row.get("condition_type")) else ""
            rep = str(row.get("biological_replicate_id", "")).strip() if pd.notna(row.get("biological_replicate_id")) else ""

            if cond in ("treatment", "control"):
                if not rep or rep.lower() in ("none", "nan"):
                    validation_errors.append(f"Missing biological_replicate_id for active sample '{sid}'. Explicit biological replicate assignment is required.")
                active_rows.append({"sample_id": sid, "condition_type": cond, "biological_replicate_id": rep})
            elif cond not in ("exclude", "excluded", "none", "", "unassigned"):
                validation_errors.append(f"Invalid condition '{cond}' for sample '{sid}'. Please assign 'control', 'treatment', or 'exclude'.")

        trt_count = sum(1 for r in active_rows if r["condition_type"] == "treatment")
        ctrl_count = sum(1 for r in active_rows if r["condition_type"] == "control")
        if trt_count < 2 or ctrl_count < 2:
            validation_errors.append(
                f"Differential expression requires at least 2 biological replicates in both treatment and control. "
                f"Currently assigned active: {trt_count} treatment, {ctrl_count} control."
            )

    if not chem_name.strip():
        validation_errors.append("Please specify Chemical / Agent Name.")
    if not dev_age.strip():
        validation_errors.append("Please specify Developmental Age (e.g., 'Day 6').")

    if validation_errors:
        for err in validation_errors:
            st.error(f"❌ {err}")
    else:
        samples_df = pd.DataFrame(
            {
                "sample_id": [str(s).strip() for s in edited_samples["sample_id"]],
                "condition_type": [str(c).strip().lower() if pd.notna(c) else "exclude" for c in edited_samples["condition_type"]],
                "biological_replicate_id": [str(r).strip() if pd.notna(r) else "" for r in edited_samples["biological_replicate_id"]],
            }
        )

        chem_params = ChemicalParams(
            agent_name=chem_name.strip(),
            agent_identifier=chem_id.strip(),
            concentration_or_dose=chem_dose,
            concentration_or_dose_unit=chem_unit,
            developmental_age=dev_age.strip(),
            exposure_duration=exp_dur,
            exposure_duration_unit=exp_dur_unit,
            vehicle=chem_vehicle.strip(),
            organoid_context_id=f"ctx_{organoid_ctx.lower().replace(' ', '_')}" if organoid_ctx.strip() else "organoid_context",
            organoid_type="cerebral_organoid",
            brain_region=brain_reg.strip() or "cerebral_cortex",
            biological_source_id=bio_src.strip() or "biological_source",
            species="Homo sapiens",
        )

        with st.status("Executing differential expression pipeline...", expanded=True) as status_box:
            def on_progress(stage: str, msg: str):
                status_box.write(f"**{stage}**: {msg}")

            res = run_full_analysis_pipeline(
                counts_input=counts_df if counts_df is not None else uploaded_counts,
                samples_input=samples_df,
                chemical_params=chem_params,
                source_accession=source_acc.strip(),
                progress_callback=on_progress,
                r_timeout_seconds=600,
            )
            if res.is_success:
                status_box.update(label="✅ Differential expression pipeline completed successfully!", state="complete", expanded=False)
            else:
                status_box.update(label=f"❌ Pipeline execution halted: {res.error_message}", state="error", expanded=True)
            st.session_state["pipeline_result"] = res

# -----------------------------------------------------------------------------
# Results Presentation
# -----------------------------------------------------------------------------
res: PipelineRunResult | None = st.session_state.get("pipeline_result")

if res is None:
    # Initial landing screen guidance
    if uploaded_counts is not None and counts_df is not None:
        st.success(f"📁 Count matrix loaded: **{len(counts_df):,}** genes across **{len(sample_ids)}** samples.")
        st.info("👉 Please assign `condition_type` ('control' or 'treatment') and `biological_replicate_id` in the sidebar table, enter exposure details, and click **🚀 ANALYZE**.")
    else:
        st.info("👋 Upload a raw-count matrix in the sidebar, assign sample conditions and biological replicates, enter your exposure parameters, and click **🚀 ANALYZE** to run the pipeline.")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown(
            """
            ### 1. Ingestion & Quality Control
            - Validates non-negative integer count matrix.
            - Reconciles samples with explicit canonical metadata.
            - Evaluates library sizes and zero-count fractions.
            """
        )
    with col2:
        st.markdown(
            """
            ### 2. Harmonization & Normalization
            - Maps identifiers to Ensembl 112 canonical IDs.
            - Prevents silent aggregation of colliding genes.
            - Calculates locked DESeq2 median-of-ratios size factors.
            """
        )
    with col3:
        st.markdown(
            """
            ### 3. Contrast & DESeq2 Modeling
            - Gating on exact developmental age matching.
            - Enforces $\ge 2$ true biological replicates per arm.
            - Official DESeq2 negative-binomial Wald test.
            """
        )

elif not res.is_success:
    st.error(f"❌ Analysis Pipeline Blocked / Failed: {res.error_message}")
    if res.scientific_findings:
        st.subheader("Auditable Findings")
        findings_df = pd.DataFrame(res.scientific_findings)
        st.dataframe(findings_df, use_container_width=True)

else:
    # Analysis Success Banner
    cr = res.contrast_result
    st.success(
        f"✅ **Differential Expression Profile Generated**: `{cr.contrast_id}` | "
        f"Exposure: **{cr.agent_name}** ({cr.concentration_or_dose} {cr.concentration_or_dose_unit}, {cr.exposure_duration} {cr.exposure_duration_unit}) | "
        f"Stage: **{cr.treatment_collection_age_or_stage}** | "
        f"Model: **{cr.design_formula}**"
    )

    # KPI Metrics
    m = res.summary_metrics
    col_m1, col_m2, col_m3, col_m4, col_m5, col_m6 = st.columns(6)
    with col_m1:
        st.markdown(
            f'<div class="kpi-card"><div class="kpi-value">{m.get("evaluated_canonical_genes", 0)}</div>'
            f'<div class="kpi-label">Evaluated Genes</div></div>',
            unsafe_allow_html=True,
        )
    with col_m2:
        st.markdown(
            f'<div class="kpi-card"><div class="kpi-value" style="color: #DC2626;">{m.get("significant_upregulated", 0)}</div>'
            f'<div class="kpi-label">Sig. Upregulated</div></div>',
            unsafe_allow_html=True,
        )
    with col_m3:
        st.markdown(
            f'<div class="kpi-card"><div class="kpi-value" style="color: #2563EB;">{m.get("significant_downregulated", 0)}</div>'
            f'<div class="kpi-label">Sig. Downregulated</div></div>',
            unsafe_allow_html=True,
        )
    with col_m4:
        st.markdown(
            f'<div class="kpi-card"><div class="kpi-value">{m.get("treatment_replicates", 0)} vs {m.get("control_replicates", 0)}</div>'
            f'<div class="kpi-label">Replicates (Trt vs Ctrl)</div></div>',
            unsafe_allow_html=True,
        )
    with col_m5:
        st.markdown(
            f'<div class="kpi-card"><div class="kpi-value">{m.get("residual_degrees_of_freedom", 0)}</div>'
            f'<div class="kpi-label">Residual DF</div></div>',
            unsafe_allow_html=True,
        )
    with col_m6:
        st.markdown(
            f'<div class="kpi-card"><div class="kpi-value" style="font-size: 1.2rem; margin-top: 8px;">{m.get("deseq2_version", "1.46.0")}</div>'
            f'<div class="kpi-label">DESeq2 Release</div></div>',
            unsafe_allow_html=True,
        )

    # Scientific Checks Callout
    with st.expander("🛡️ Verified Scientific Invariants & Quality Gates", expanded=False):
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(
                f"- **Biological Replication**: Verified {cr.treatment_replicate_count} treatment and {cr.control_replicate_count} control biological replicates (minimum 2 required in each arm).\n"
                f"- **Gene-Universe Conservation**: Verified {cr.total_input_genes_count} input rows = {cr.eligible_canonical_genes_count} evaluated canonical genes + {cr.excluded_genes_count} excluded genes.\n"
                f"- **Age Equivalence**: Treatment collection age (`{cr.treatment_collection_age_or_stage}`) and control collection age (`{cr.control_collection_age_or_stage}`) strictly match."
            )
        with c2:
            st.markdown(
                f"- **Size Factors**: Locked DESeq2 ratio size factors applied without alteration.\n"
                f"- **Locked Runtime**: Validated against R {m.get('r_version')} / Bioconductor {m.get('bioc_version')} / DESeq2 {m.get('deseq2_version')}.\n"
                f"- **Statistical Test**: Two-sided negative-binomial Wald test with Benjamini-Hochberg FDR control."
            )

    # Main Tabs
    tab_plots, tab_table, tab_audit, tab_samples, tab_provenance = st.tabs(
        [
            "📊 Visual Plots (Volcano & MA)",
            "🧬 Differential Expression Gene Table",
            "🛡️ Gene-Universe Conservation Audit",
            "📋 Samples & Size Factors",
            "📜 Provenance & Findings Log",
        ]
    )

    # -------------------------------------------------------------------------
    # TAB 1: Visual Plots
    # -------------------------------------------------------------------------
    with tab_plots:
        de_df = res.de_genes_df.copy()
        if not de_df.empty:
            de_df["neg_log10_padj"] = -np.log10(de_df["adjusted_p_value_bh"].replace(0, 1e-300))
            de_df["log10_base_mean"] = np.log10(de_df["base_mean"].replace(0, 1))

            plot_df, is_sampled, total_sig, total_bg = prepare_plot_data(de_df)

            col_p1, col_p2 = st.columns(2)

            with col_p1:
                st.subheader("Volcano Plot")
                color_map = {
                    "Upregulated": "#EF4444",
                    "Downregulated": "#3B82F6",
                    "Not Significant": "#9CA3AF",
                }
                volcano_fig = px.scatter(
                    plot_df,
                    x="log2_fold_change",
                    y="neg_log10_padj",
                    color="direction",
                    color_discrete_map=color_map,
                    hover_name="approved_symbol",
                    hover_data={
                        "canonical_gene_id": True,
                        "log2_fold_change": ":.3f",
                        "adjusted_p_value_bh": ":.3e",
                        "base_mean": ":.1f",
                        "neg_log10_padj": False,
                        "direction": False,
                    },
                    labels={
                        "log2_fold_change": "log2(Fold Change)",
                        "neg_log10_padj": "-log10(Adjusted P-Value)",
                        "direction": "Status",
                    },
                    title=f"Volcano Plot: {cr.agent_name} vs Matched Control",
                )
                volcano_fig.add_hline(
                    y=-math.log10(0.05),
                    line_dash="dash",
                    line_color="#6B7280",
                    annotation_text="FDR = 0.05",
                    annotation_position="bottom right",
                )
                volcano_fig.add_vline(x=1.0, line_dash="dot", line_color="#D1D5DB")
                volcano_fig.add_vline(x=-1.0, line_dash="dot", line_color="#D1D5DB")
                volcano_fig.update_layout(
                    height=520,
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                    template="plotly_white",
                )
                st.plotly_chart(volcano_fig, use_container_width=True)

            with col_p2:
                st.subheader("MA Plot (Mean vs Log2FC)")
                ma_fig = px.scatter(
                    plot_df,
                    x="log10_base_mean",
                    y="log2_fold_change",
                    color="direction",
                    color_discrete_map=color_map,
                    hover_name="approved_symbol",
                    hover_data={
                        "canonical_gene_id": True,
                        "log2_fold_change": ":.3f",
                        "base_mean": ":.1f",
                        "adjusted_p_value_bh": ":.3e",
                        "log10_base_mean": False,
                        "direction": False,
                    },
                    labels={
                        "log10_base_mean": "log10(Base Mean)",
                        "log2_fold_change": "log2(Fold Change)",
                        "direction": "Status",
                    },
                    title=f"MA Plot: {cr.agent_name} Treatment Effect",
                )
                ma_fig.add_hline(y=0.0, line_dash="solid", line_color="#9CA3AF")
                ma_fig.update_layout(
                    height=520,
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                    template="plotly_white",
                )
                st.plotly_chart(ma_fig, use_container_width=True)

            if is_sampled:
                st.caption(
                    f"ℹ️ **Performance-Optimized Display**: Plot displays 100% of statistically significant genes ({total_sig:,}) "
                    f"and a deterministic representative sample of 5,000 of {total_bg:,} background genes (seed 42) "
                    "to maintain responsive, lag-free browser visualization. Complete evaluated canonical genes "
                    f"({len(de_df):,}) are preserved without sampling in the Differential Expression Gene Table, Conservation Audit, and CSV export."
                )
        else:
            st.warning("No differential expression results available to plot.")

    # -------------------------------------------------------------------------
    # TAB 2: Gene Table
    # -------------------------------------------------------------------------
    with tab_table:
        st.subheader("Differential Expression Results")

        filter_col1, filter_col2, filter_col3, filter_col4 = st.columns([2, 2, 2, 3])
        with filter_col1:
            dir_filter = st.selectbox(
                "Filter Direction",
                options=["All Genes", "Significant (FDR < 0.05)", "Upregulated Only", "Downregulated Only"],
            )
        with filter_col2:
            p_cutoff = st.number_input("Max FDR Cutoff", value=0.05, min_value=0.0001, max_value=1.0, step=0.01)
        with filter_col3:
            lfc_cutoff = st.number_input("Min |log2FC|", value=0.0, min_value=0.0, max_value=10.0, step=0.5)
        with filter_col4:
            search_query = st.text_input("🔍 Search Gene Symbol or Ensembl ID", value="")

        display_df = res.de_genes_df.copy()

        if dir_filter == "Significant (FDR < 0.05)":
            display_df = display_df[display_df["adjusted_p_value_bh"] < 0.05]
        elif dir_filter == "Upregulated Only":
            display_df = display_df[(display_df["adjusted_p_value_bh"] < p_cutoff) & (display_df["log2_fold_change"] > 0)]
        elif dir_filter == "Downregulated Only":
            display_df = display_df[(display_df["adjusted_p_value_bh"] < p_cutoff) & (display_df["log2_fold_change"] < 0)]

        if lfc_cutoff > 0:
            display_df = display_df[display_df["log2_fold_change"].abs() >= lfc_cutoff]

        if search_query.strip():
            sq = search_query.strip().lower()
            display_df = display_df[
                display_df["approved_symbol"].str.lower().str.contains(sq)
                | display_df["canonical_gene_id"].str.lower().str.contains(sq)
            ]

        st.caption(f"Showing **{len(display_df)}** of **{len(res.de_genes_df)}** genes evaluated.")

        formatted_table = display_df.copy()
        formatted_table["base_mean"] = formatted_table["base_mean"].map(lambda x: f"{x:.2f}" if pd.notnull(x) else "")
        formatted_table["log2_fold_change"] = formatted_table["log2_fold_change"].map(lambda x: f"{x:.3f}" if pd.notnull(x) else "")
        formatted_table["lfc_standard_error"] = formatted_table["lfc_standard_error"].map(lambda x: f"{x:.3f}" if pd.notnull(x) else "")
        formatted_table["wald_statistic"] = formatted_table["wald_statistic"].map(lambda x: f"{x:.2f}" if pd.notnull(x) else "")
        formatted_table["p_value"] = formatted_table["p_value"].map(lambda x: f"{x:.3e}" if pd.notnull(x) else "")
        formatted_table["adjusted_p_value_bh"] = formatted_table["adjusted_p_value_bh"].map(lambda x: f"{x:.3e}" if pd.notnull(x) else "")

        st.dataframe(
            formatted_table[
                [
                    "canonical_gene_id",
                    "approved_symbol",
                    "base_mean",
                    "log2_fold_change",
                    "lfc_standard_error",
                    "wald_statistic",
                    "p_value",
                    "adjusted_p_value_bh",
                    "direction",
                    "result_status",
                ]
            ],
            use_container_width=True,
            height=450,
        )

    # -------------------------------------------------------------------------
    # TAB 3: Gene Conservation Audit
    # -------------------------------------------------------------------------
    with tab_audit:
        st.subheader("Gene-Universe Conservation Audit")
        st.markdown(
            """
            In accordance with **Gene-Universe Conservation**, every input row in the expression matrix
            is deterministically audited. No genes are silently dropped or imputed.
            """
        )

        c_audit1, c_audit2, c_audit3 = st.columns(3)
        with c_audit1:
            st.metric("Total Input Rows", cr.total_input_genes_count)
        with c_audit2:
            st.metric("DE-Evaluated Canonical Genes", cr.eligible_canonical_genes_count)
        with c_audit3:
            st.metric("Excluded Genes", cr.excluded_genes_count)

        if cr.is_gene_universe_conserved:
            st.success("✅ **Conservation Invariant Holds**: Total Input Rows = DE-Evaluated Genes + Excluded Genes.")
        else:
            st.error("❌ Conservation Invariant Mismatch! Audit failed.")

        if not res.excluded_genes_df.empty:
            st.write("### Excluded Gene Detail")
            st.dataframe(res.excluded_genes_df, use_container_width=True)
        else:
            st.info("Zero genes were excluded from canonical analysis. 100% of input features were uniquely mapped and evaluated.")

    # -------------------------------------------------------------------------
    # TAB 4: Samples & Size Factors
    # -------------------------------------------------------------------------
    with tab_samples:
        st.subheader("Sample Metadata & DESeq2 Size Factors")
        st.markdown(
            "Size factors were estimated during **Step 6C Bulk Normalization v1** using the median-of-ratios method "
            "and locked for differential expression modeling."
        )

        st.dataframe(res.samples_df, use_container_width=True)

        if not res.samples_df.empty:
            sample_fig = px.bar(
                res.samples_df,
                x="sample_id",
                y="deseq2_size_factor",
                color="condition_type",
                color_discrete_map={"control": "#6B7280", "treatment": "#2563EB"},
                labels={"deseq2_size_factor": "Locked DESeq2 Size Factor", "sample_id": "Sample"},
                title="DESeq2 Size Factors by Sample",
            )
            sample_fig.update_layout(height=350, template="plotly_white")
            st.plotly_chart(sample_fig, use_container_width=True)

    # -------------------------------------------------------------------------
    # TAB 5: Provenance & Audit Log
    # -------------------------------------------------------------------------
    with tab_provenance:
        st.subheader("Scientific Provenance & Audit Log")

        m = res.summary_metrics
        if m.get("source_accession"):
            st.info(f"🌐 **Public Repository Provenance**: NCBI GEO `{m['source_accession']}` (Source type: `public`)")
        else:
            st.info("🏢 **In-House Study Provenance**: Source type: `in_house`")

        if m.get("excluded_annotation_columns"):
            st.markdown(f"- **Excluded Annotation Columns**: `{', '.join(m['excluded_annotation_columns'])}`")

        st.write("#### Execution Findings")
        if res.scientific_findings:
            f_df = pd.DataFrame(res.scientific_findings)
            st.dataframe(f_df, use_container_width=True)
        else:
            st.info("No warnings or review findings encountered during analysis.")

        st.write("#### Software Runtime Environment")
        st.json(
            {
                "r_version": m.get("r_version"),
                "bioconductor_version": m.get("bioc_version"),
                "deseq2_version": m.get("deseq2_version"),
                "pipeline_version": "Step 7B v1.0.0",
                "design_formula": m.get("design_formula"),
                "residual_degrees_of_freedom": m.get("residual_degrees_of_freedom"),
                "source_accession": m.get("source_accession", ""),
                "source_type": m.get("source_type", "in_house"),
            }
        )

    # -------------------------------------------------------------------------
    # Downloads Section
    # -------------------------------------------------------------------------
    st.markdown("---")
    st.subheader("📥 Download Analysis Results")

    col_d1, col_d2, col_d3 = st.columns(3)

    # 1. DE Results CSV
    csv_buffer = io.StringIO()
    res.de_genes_df.to_csv(csv_buffer, index=False)
    with col_d1:
        st.download_button(
            label="📄 Download Differential Expression (CSV)",
            data=csv_buffer.getvalue(),
            file_name=f"de_results_{cr.agent_name}_{cr.developmental_age_or_stage_at_exposure or 'contrast'}.csv",
            mime="text/csv",
            use_container_width=True,
        )

    # 2. Excluded Genes CSV
    ex_buffer = io.StringIO()
    res.excluded_genes_df.to_csv(ex_buffer, index=False)
    with col_d2:
        st.download_button(
            label="🛡️ Download Excluded Genes Audit (CSV)",
            data=ex_buffer.getvalue(),
            file_name=f"excluded_genes_audit_{cr.agent_name}.csv",
            mime="text/csv",
            use_container_width=True,
        )

    # 3. Complete JSON Summary
    summary_json = json.dumps(res.to_summary_dict(), indent=2)
    with col_d3:
        st.download_button(
            label="📦 Download Complete Summary (JSON)",
            data=summary_json,
            file_name=f"analysis_summary_{cr.agent_name}.json",
            mime="application/json",
            use_container_width=True,
        )
