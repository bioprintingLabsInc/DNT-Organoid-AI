"""JSON source metadata adapter for structured study manifests and BioStudies JSON."""

import json
from pathlib import Path
from typing import Any

from .base import MappingFinding, MetadataSourceAdapter, RawSampleRecord, StudyContext


class JsonMetadataAdapter(MetadataSourceAdapter):
    """Adapter for JSON study manifests and repository metadata."""

    def can_read(self, source: Any) -> bool:
        if isinstance(source, dict):
            return True
        if isinstance(source, (str, Path)):
            p = Path(source)
            return p.name.lower().endswith(".json") and p.is_file()
        return False

    def parse_study(
        self,
        source: Any,
        study_id: str | None = None,
        experiment_id: str = "exp_01",
        **kwargs: Any,
    ) -> tuple[StudyContext, list[RawSampleRecord], list[MappingFinding]]:
        findings: list[MappingFinding] = []
        source_name = str(source) if isinstance(source, (str, Path)) else "in_memory_json"

        data: dict[str, Any]
        if isinstance(source, dict):
            data = source
        else:
            p = Path(source)
            try:
                data = json.loads(p.read_text(encoding="utf-8-sig", errors="replace"))
            except Exception as e:
                findings.append(
                    MappingFinding(
                        severity="ERROR",
                        rule="json_parse_error",
                        field_path=str(p),
                        message=f"Failed to parse JSON file: {e}",
                    )
                )
                return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

        # Check if this is already a canonical document
        if "samples" in data and isinstance(data["samples"], list):
            sample_list = data["samples"]
        elif "study" in data and "samples" in data["study"]:
            sample_list = data["study"]["samples"]
        else:
            findings.append(
                MappingFinding(
                    severity="ERROR",
                    rule="missing_samples_array",
                    field_path=source_name,
                    message="JSON document must contain a top-level 'samples' array.",
                )
            )
            return StudyContext(study_id=study_id or "study_unknown", findings=tuple(findings)), [], findings

        raw_records: list[RawSampleRecord] = []
        for idx, s_obj in enumerate(sample_list):
            if not isinstance(s_obj, dict):
                continue
            sid = s_obj.get("sample_id") or s_obj.get("sample_alias") or s_obj.get("id") or f"sample_{idx+1}"
            key_prov = {k: k for k in s_obj.keys()}
            raw_records.append(
                RawSampleRecord(
                    sample_id=str(sid),
                    attributes=dict(s_obj),
                    source_file=source_name,
                    source_locator=f"samples[{idx}]",
                    source_key_provenance=key_prov,
                )
            )

        inferred_study_id = (
            study_id
            or data.get("study_id")
            or data.get("bioproject")
            or f"study_{Path(source_name).stem}"
        )

        pub_accessions: list[dict[str, str]] = []
        if data.get("bioproject"):
            pub_accessions.append({"namespace": "BioProject", "identifier": data["bioproject"]})
        if data.get("sra_study"):
            pub_accessions.append({"namespace": "SRA", "identifier": data["sra_study"]})

        context = StudyContext(
            study_id=inferred_study_id,
            study_title=data.get("study_title") or data.get("contrast"),
            source_type="public",
            public_accessions=tuple(pub_accessions),
            experiment_id=experiment_id,
            attributes={k: v for k, v in data.items() if k != "samples"},
            findings=tuple(findings),
        )

        return context, raw_records, findings
