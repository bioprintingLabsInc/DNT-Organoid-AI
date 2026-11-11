"""RSeQC infer_experiment.py strandedness diagnostic wrapper and evidence-based decision engine."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
from typing import Any

from .adapter import sha256_file
from .env_exec import run_wsl, to_wsl_path
from .library_type import LibraryType


@dataclass(frozen=True, slots=True)
class StrandednessDecision:
    """Evidence-based strandedness determination outcome."""

    strandedness_flag: int  # 0: Unstranded, 1: Forward, 2: Reverse, -1: Ambiguous
    strand_mode: str  # "unstranded", "forward", "reverse", "ambiguous"
    confidence: str  # "HIGH", "MEDIUM", "LOW", "AMBIGUOUS", "CONFLICT"
    mean_forward: float
    mean_reverse: float
    mean_unexplained: float
    consensus_indication: str
    rationale: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "strandedness_flag": self.strandedness_flag,
            "strand_mode": self.strand_mode,
            "confidence": self.confidence,
            "mean_forward": round(self.mean_forward, 4),
            "mean_reverse": round(self.mean_reverse, 4),
            "mean_unexplained": round(self.mean_unexplained, 4),
            "consensus_indication": self.consensus_indication,
            "rationale": self.rationale,
        }


def infer_strandedness_decision(
    metrics_or_summary: dict[str, Any],
    library_type: str | LibraryType | None = None,
    configured_strandedness: int | None = None,
) -> StrandednessDecision:
    """Evaluate empirical RSeQC evidence to determine the featureCounts strandedness flag.

    Accepts either:
    1. A strandedness summary dictionary with keys:
       'mean_forward_fraction', 'mean_reverse_fraction', 'mean_unexplained_fraction'
    2. A dictionary mapping sample IDs to infer_experiment result dictionaries.

    If configured_strandedness is provided by investigator/protocol, it is verified
    against empirical RSeQC evidence and MUST NOT silently bypass conflicting observations.

    If library_type orientation is unknown (e.g. LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND),
    verified empirical evidence is required to select counting orientation, halting
    if evidence is ambiguous or non-stranded.

    Scientific Decision Invariants:
    - Reverse-stranded (featureCounts -s 2):
      Dominant antisense reads (mean_rev >= 0.70) with low sense signal (mean_fwd <= 0.20) or diff >= 0.50.
    - Forward-stranded (featureCounts -s 1):
      Dominant sense reads (mean_fwd >= 0.70) with low antisense signal (mean_rev <= 0.20) or diff >= 0.50.
    - Unstranded (featureCounts -s 0):
      Balanced forward/reverse fractions (|mean_fwd - mean_rev| <= 0.15, each >= 0.25).
    - Ambiguous / Conflict:
      Neither criterion met, or empirical evidence conflicts with declared library prep or configured flag;
      halts for investigator review (flag = -1).
    """
    lib_type = LibraryType.from_string(library_type) if library_type else LibraryType.UNKNOWN

    if "mean_forward_fraction" in metrics_or_summary and "mean_reverse_fraction" in metrics_or_summary:
        mean_fwd = float(metrics_or_summary["mean_forward_fraction"])
        mean_rev = float(metrics_or_summary["mean_reverse_fraction"])
        mean_unexp = float(metrics_or_summary.get("mean_unexplained_fraction", 0.0))
    else:
        # Per-sample dictionary
        fwd_list: list[float] = []
        rev_list: list[float] = []
        unexp_list: list[float] = []
        for val in metrics_or_summary.values():
            if isinstance(val, dict):
                f = val.get("forward_fraction")
                r = val.get("reverse_fraction")
                u = val.get("unexplained_fraction", 0.0)
                if f is not None:
                    fwd_list.append(float(f))
                if r is not None:
                    rev_list.append(float(r))
                if u is not None:
                    unexp_list.append(float(u))

        mean_fwd = sum(fwd_list) / len(fwd_list) if fwd_list else 0.0
        mean_rev = sum(rev_list) / len(rev_list) if rev_list else 0.0
        mean_unexp = sum(unexp_list) / len(unexp_list) if unexp_list else 0.0

    delta = abs(mean_fwd - mean_rev)

    # 1. Evaluate empirical signal
    empirical_flag: int
    empirical_mode: str
    if (mean_rev >= 0.70 and mean_fwd <= 0.20) or (mean_rev - mean_fwd >= 0.50):
        empirical_flag = 2
        empirical_mode = "reverse"
    elif (mean_fwd >= 0.70 and mean_rev <= 0.20) or (mean_fwd - mean_rev >= 0.50):
        empirical_flag = 1
        empirical_mode = "forward"
    elif delta <= 0.15 and (mean_fwd >= 0.25 and mean_rev >= 0.25):
        empirical_flag = 0
        empirical_mode = "unstranded"
    else:
        empirical_flag = -1
        empirical_mode = "ambiguous"

    # 2. Investigator-configured strandedness verification
    if configured_strandedness is not None:
        if configured_strandedness == 1:
            is_conflict = (empirical_flag == 2) or (empirical_flag == 0 and mean_fwd >= 0.25 and mean_rev >= 0.25)
            if is_conflict:
                return StrandednessDecision(
                    strandedness_flag=-1,
                    strand_mode="ambiguous",
                    confidence="CONFLICT",
                    mean_forward=mean_fwd,
                    mean_reverse=mean_rev,
                    mean_unexplained=mean_unexp,
                    consensus_indication=f"CONFLICT_CONFIGURED_FORWARD_EMPIRICAL_{empirical_mode.upper()}",
                    rationale=(
                        f"Investigator configured forward strandedness (-s 1), but empirical RSeQC evidence "
                        f"indicates {empirical_mode} (forward={mean_fwd:.1%}, reverse={mean_rev:.1%}). "
                        "Halting downstream quantification because configured strandedness cannot silently bypass "
                        "conflicting empirical observations."
                    ),
                )
            else:
                return StrandednessDecision(
                    strandedness_flag=1,
                    strand_mode="forward",
                    confidence="HIGH",
                    mean_forward=mean_fwd,
                    mean_reverse=mean_rev,
                    mean_unexplained=mean_unexp,
                    consensus_indication="CONSENSUS_VERIFIED_USER_CONFIGURED (featureCounts -s 1)",
                    rationale=(
                        f"Configured forward strandedness (-s 1) accepted and consistent with empirical observations "
                        f"(forward={mean_fwd:.1%}, reverse={mean_rev:.1%})."
                    ),
                )
        elif configured_strandedness == 2:
            is_conflict = (empirical_flag == 1) or (empirical_flag == 0 and mean_fwd >= 0.25 and mean_rev >= 0.25)
            if is_conflict:
                return StrandednessDecision(
                    strandedness_flag=-1,
                    strand_mode="ambiguous",
                    confidence="CONFLICT",
                    mean_forward=mean_fwd,
                    mean_reverse=mean_rev,
                    mean_unexplained=mean_unexp,
                    consensus_indication=f"CONFLICT_CONFIGURED_REVERSE_EMPIRICAL_{empirical_mode.upper()}",
                    rationale=(
                        f"Investigator configured reverse strandedness (-s 2), but empirical RSeQC evidence "
                        f"indicates {empirical_mode} (forward={mean_fwd:.1%}, reverse={mean_rev:.1%}). "
                        "Halting downstream quantification because configured strandedness cannot silently bypass "
                        "conflicting empirical observations."
                    ),
                )
            else:
                return StrandednessDecision(
                    strandedness_flag=2,
                    strand_mode="reverse",
                    confidence="HIGH",
                    mean_forward=mean_fwd,
                    mean_reverse=mean_rev,
                    mean_unexplained=mean_unexp,
                    consensus_indication="CONSENSUS_VERIFIED_USER_CONFIGURED (featureCounts -s 2)",
                    rationale=(
                        f"Configured reverse strandedness (-s 2) accepted and consistent with empirical observations "
                        f"(forward={mean_fwd:.1%}, reverse={mean_rev:.1%})."
                    ),
                )
        elif configured_strandedness == 0:
            if empirical_flag in (1, 2):
                return StrandednessDecision(
                    strandedness_flag=-1,
                    strand_mode="ambiguous",
                    confidence="CONFLICT",
                    mean_forward=mean_fwd,
                    mean_reverse=mean_rev,
                    mean_unexplained=mean_unexp,
                    consensus_indication=f"CONFLICT_CONFIGURED_UNSTRANDED_EMPIRICAL_{empirical_mode.upper()}",
                    rationale=(
                        f"Investigator configured unstranded setting (-s 0), but empirical RSeQC evidence "
                        f"indicates {empirical_mode} (forward={mean_fwd:.1%}, reverse={mean_rev:.1%}). "
                        "Halting downstream quantification because configured strandedness cannot silently bypass "
                        "conflicting empirical observations."
                    ),
                )
            else:
                return StrandednessDecision(
                    strandedness_flag=0,
                    strand_mode="unstranded",
                    confidence="HIGH",
                    mean_forward=mean_fwd,
                    mean_reverse=mean_rev,
                    mean_unexplained=mean_unexp,
                    consensus_indication="CONSENSUS_VERIFIED_USER_CONFIGURED (featureCounts -s 0)",
                    rationale=(
                        f"Configured unstranded setting (-s 0) accepted and consistent with empirical observations "
                        f"(forward={mean_fwd:.1%}, reverse={mean_rev:.1%})."
                    ),
                )

    # 3. Consensus Gating against declared library type
    if lib_type.is_fwd:
        # QuantSeq FWD requires forward sense reads (-s 1)
        if empirical_flag == 1:
            return StrandednessDecision(
                strandedness_flag=1,
                strand_mode="forward",
                confidence="HIGH",
                mean_forward=mean_fwd,
                mean_reverse=mean_rev,
                mean_unexplained=mean_unexp,
                consensus_indication="CONSENSUS_VERIFIED_FORWARD (featureCounts -s 1)",
                rationale=(
                    f"Declared library type {lib_type.value} verified by empirical RSeQC evidence: "
                    f"dominant forward/sense orientation ({mean_fwd:.1%}) with low antisense signal ({mean_rev:.1%})."
                ),
            )
        else:
            return StrandednessDecision(
                strandedness_flag=-1,
                strand_mode="ambiguous",
                confidence="CONFLICT",
                mean_forward=mean_fwd,
                mean_reverse=mean_rev,
                mean_unexplained=mean_unexp,
                consensus_indication=f"CONFLICT_DECLARED_FWD_EMPIRICAL_{empirical_mode.upper()}",
                rationale=(
                    f"Library preparation is declared as {lib_type.value} (which requires forward strand -s 1), "
                    f"but empirical RSeQC fractions indicate {empirical_mode} "
                    f"(forward={mean_fwd:.1%}, reverse={mean_rev:.1%}, unexplained={mean_unexp:.1%}). "
                    "Halting downstream quantification to prevent erroneous strandedness assignment."
                ),
            )

    if lib_type.is_rev:
        # QuantSeq REV requires reverse antisense reads (-s 2)
        if empirical_flag == 2:
            return StrandednessDecision(
                strandedness_flag=2,
                strand_mode="reverse",
                confidence="HIGH",
                mean_forward=mean_fwd,
                mean_reverse=mean_rev,
                mean_unexplained=mean_unexp,
                consensus_indication="CONSENSUS_VERIFIED_REVERSE (featureCounts -s 2)",
                rationale=(
                    f"Declared library type {lib_type.value} verified by empirical RSeQC evidence: "
                    f"dominant reverse/antisense orientation ({mean_rev:.1%}) with low forward signal ({mean_fwd:.1%})."
                ),
            )
        else:
            return StrandednessDecision(
                strandedness_flag=-1,
                strand_mode="ambiguous",
                confidence="CONFLICT",
                mean_forward=mean_fwd,
                mean_reverse=mean_rev,
                mean_unexplained=mean_unexp,
                consensus_indication=f"CONFLICT_DECLARED_REV_EMPIRICAL_{empirical_mode.upper()}",
                rationale=(
                    f"Library preparation is declared as {lib_type.value} (which requires reverse strand -s 2), "
                    f"but empirical RSeQC fractions indicate {empirical_mode} "
                    f"(forward={mean_fwd:.1%}, reverse={mean_rev:.1%}, unexplained={mean_unexp:.1%}). "
                    "Halting downstream quantification to prevent erroneous strandedness assignment."
                ),
            )

    # 4. QuantSeq with unknown orientation: requires verified empirical evidence
    if lib_type.is_quantseq:
        if empirical_flag == 1:
            return StrandednessDecision(
                strandedness_flag=1,
                strand_mode="forward",
                confidence="HIGH",
                mean_forward=mean_fwd,
                mean_reverse=mean_rev,
                mean_unexplained=mean_unexp,
                consensus_indication="EMPIRICAL_QUANTSEQ_FORWARD (featureCounts -s 1)",
                rationale=(
                    f"QuantSeq 3' orientation was not documented in metadata; verified empirical RSeQC evidence "
                    f"resolved orientation to forward/sense ({mean_fwd:.1%} forward vs {mean_rev:.1%} reverse)."
                ),
            )
        elif empirical_flag == 2:
            return StrandednessDecision(
                strandedness_flag=2,
                strand_mode="reverse",
                confidence="HIGH",
                mean_forward=mean_fwd,
                mean_reverse=mean_rev,
                mean_unexplained=mean_unexp,
                consensus_indication="EMPIRICAL_QUANTSEQ_REVERSE (featureCounts -s 2)",
                rationale=(
                    f"QuantSeq 3' orientation was not documented in metadata; verified empirical RSeQC evidence "
                    f"resolved orientation to reverse/antisense ({mean_rev:.1%} reverse vs {mean_fwd:.1%} forward)."
                ),
            )
        else:
            return StrandednessDecision(
                strandedness_flag=-1,
                strand_mode="ambiguous",
                confidence="CONFLICT",
                mean_forward=mean_fwd,
                mean_reverse=mean_rev,
                mean_unexplained=mean_unexp,
                consensus_indication=f"CONFLICT_QUANTSEQ_NON_STRANDED_EMPIRICAL_{empirical_mode.upper()}",
                rationale=(
                    f"QuantSeq 3' orientation was not documented in metadata, and empirical RSeQC fractions "
                    f"indicate {empirical_mode} (forward={mean_fwd:.1%}, reverse={mean_rev:.1%}). QuantSeq libraries "
                    "are strand-specific and cannot be unstranded. Halting quantification to prevent erroneous assignment."
                ),
            )

    # 5. Standard conventional / unspecified handling
    if empirical_flag == 2:
        return StrandednessDecision(
            strandedness_flag=2,
            strand_mode="reverse",
            confidence="HIGH",
            mean_forward=mean_fwd,
            mean_reverse=mean_rev,
            mean_unexplained=mean_unexp,
            consensus_indication="REVERSE_STRANDED (featureCounts -s 2)",
            rationale=(
                f"Dominant reverse strand orientation ({mean_rev:.1%}) with low forward signal "
                f"({mean_fwd:.1%}) indicates reverse-stranded (dUTP / stranded-reverse) RNA-seq."
            ),
        )
    elif empirical_flag == 1:
        return StrandednessDecision(
            strandedness_flag=1,
            strand_mode="forward",
            confidence="HIGH",
            mean_forward=mean_fwd,
            mean_reverse=mean_rev,
            mean_unexplained=mean_unexp,
            consensus_indication="FORWARD_STRANDED (featureCounts -s 1)",
            rationale=(
                f"Dominant forward strand orientation ({mean_fwd:.1%}) with low reverse signal "
                f"({mean_rev:.1%}) indicates forward-stranded (directional) RNA-seq."
            ),
        )
    elif empirical_flag == 0:
        return StrandednessDecision(
            strandedness_flag=0,
            strand_mode="unstranded",
            confidence="HIGH",
            mean_forward=mean_fwd,
            mean_reverse=mean_rev,
            mean_unexplained=mean_unexp,
            consensus_indication="UNSTRANDED (featureCounts -s 0)",
            rationale=(
                f"Balanced forward ({mean_fwd:.1%}) and reverse ({mean_rev:.1%}) orientations "
                f"(diff = {delta:.2%}) confirm standard unstranded RNA-seq library preparation."
            ),
        )
    else:
        return StrandednessDecision(
            strandedness_flag=-1,
            strand_mode="ambiguous",
            confidence="AMBIGUOUS",
            mean_forward=mean_fwd,
            mean_reverse=mean_rev,
            mean_unexplained=mean_unexp,
            consensus_indication="AMBIGUOUS_EVIDENCE_REQUIRES_REVIEW",
            rationale=(
                f"Ambiguous orientation proportions (forward={mean_fwd:.1%}, reverse={mean_rev:.1%}, "
                f"unexplained={mean_unexp:.1%}) do not match unstranded or stranded signatures. "
                "Manual review required."
            ),
        )


def run_infer_experiment(
    bam_path: str | Path,
    bed_path: str | Path,
    sample_size: int = 200000,
    mapq_min: int = 10,
    output_txt_path: str | Path | None = None,
) -> dict[str, Any]:
    """Execute RSeQC infer_experiment.py with governed parameters."""
    b_p = Path(bam_path).resolve()
    bed_p = Path(bed_path).resolve()

    if not b_p.exists():
        raise FileNotFoundError(f"Input BAM not found: {b_p}")
    if not bed_p.exists():
        raise FileNotFoundError(f"Input BED not found: {bed_p}")

    bam_sha = sha256_file(b_p)
    bed_sha = sha256_file(bed_p)

    wsl_bam = to_wsl_path(b_p)
    wsl_bed = to_wsl_path(bed_p)

    cmd = f"infer_experiment.py -r {wsl_bed} -i {wsl_bam} -s {sample_size} -q {mapq_min}"
    res = run_wsl(cmd)
    if res.returncode != 0:
        raise RuntimeError(f"infer_experiment.py failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    raw_output = res.stdout.strip()
    if output_txt_path:
        out_p = Path(output_txt_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        out_p.write_text(raw_output, encoding="utf-8")

    # Parse fractions:
    # Paired-end patterns:
    #   Fraction of reads failed to determine: 0.0412
    #   Fraction of reads explained by "1++,1--,2+-,2-+": 0.0125
    #   Fraction of reads explained by "1+-,1-+,2++,2--": 0.9463
    # Single-end patterns:
    #   Fraction of reads failed to determine: 0.0412
    #   Fraction of reads explained by "++,--": 0.9463
    #   Fraction of reads explained by "+-,-+": 0.0125
    failed_frac = None
    sense_frac = None
    antisense_frac = None

    for line in raw_output.splitlines():
        if "failed to determine" in line:
            m = re.search(r":\s*([0-9.]+)", line)
            if m:
                failed_frac = float(m.group(1))
        # Forward/sense strand:
        elif "1++,1--,2+-,2-+" in line or '"++,--"' in line or '++,--' in line:
            m = re.search(r":\s*([0-9.]+)", line)
            if m:
                sense_frac = float(m.group(1))
        # Reverse/antisense strand:
        elif "1+-,1-+,2++,2--" in line or '"+-,-+"' in line or '+-,-+' in line:
            m = re.search(r":\s*([0-9.]+)", line)
            if m:
                antisense_frac = float(m.group(1))

    return {
        "status": "INFER_EXPERIMENT_SUCCESS",
        "bam_path": str(b_p),
        "bam_sha256": bam_sha,
        "bed_path": str(bed_p),
        "bed_sha256": bed_sha,
        "sample_size": sample_size,
        "mapq_min": mapq_min,
        "raw_output": raw_output,
        "forward_fraction": sense_frac,
        "reverse_fraction": antisense_frac,
        "unexplained_fraction": failed_frac,
    }
