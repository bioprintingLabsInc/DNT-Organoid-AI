"""Governed library preparation typing and validation for bulk RNA-seq pipelines.

Distinguishes:
1. Conventional bulk RNA-seq (unstranded, dUTP reverse-stranded, directional forward).
2. Lexogen QuantSeq 3' mRNA-seq FWD (sense strand, 3' polyA trimmed).
3. Lexogen QuantSeq 3' mRNA-seq REV (antisense strand, 5' polyT trimmed).
4. Lexogen QuantSeq 3' mRNA-seq with UMIs (FWD_UMI, REV_UMI).
5. Unknown / unspecified library types (fail-safe fallback).
"""

from __future__ import annotations

from enum import Enum
import re
from typing import Any


class LibraryType(str, Enum):
    """Recognized RNA-seq library preparation types."""

    CONVENTIONAL_BULK = "CONVENTIONAL_BULK"
    LEXOGEN_QUANTSEQ_3UTR_FWD = "LEXOGEN_QUANTSEQ_3UTR_FWD"
    LEXOGEN_QUANTSEQ_3UTR_REV = "LEXOGEN_QUANTSEQ_3UTR_REV"
    LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND = "LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND"
    LEXOGEN_QUANTSEQ_3UTR_FWD_UMI = "LEXOGEN_QUANTSEQ_3UTR_FWD_UMI"
    LEXOGEN_QUANTSEQ_3UTR_REV_UMI = "LEXOGEN_QUANTSEQ_3UTR_REV_UMI"
    LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND_UMI = "LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND_UMI"
    UNKNOWN = "UNKNOWN"

    @classmethod
    def from_string(cls, val: str | LibraryType | None) -> LibraryType:
        """Parse arbitrary user or repository metadata string into canonical LibraryType.
        
        Guarantees lossless classification without inventing library preparation
        details when absent.
        """
        if val is None:
            return cls.UNKNOWN
        if isinstance(val, cls):
            return val

        raw = str(val).strip()
        if not raw:
            return cls.UNKNOWN

        norm = raw.upper().replace("-", "_").replace(" ", "_").replace("'", "")

        # Exact enum name match
        try:
            return cls(norm)
        except ValueError:
            pass

        # Check for UMI indicators as isolated word/token (avoid false match in 'ILLUMINA')
        has_umi = bool(re.search(r"(?:^|[_ -])UMI(?:$|[_ -])", norm))

        # Check for non-QuantSeq Lexogen kits (e.g. CORALL is full-length bulk RNA-seq)
        if "CORALL" in norm:
            return cls.CONVENTIONAL_BULK

        # Lexogen QuantSeq patterns: MUST explicitly mention QuantSeq or 3' mRNA-Seq.
        # General 'LEXOGEN' alone or unspecified Lexogen kits MUST NOT be assumed to be QuantSeq.
        is_quantseq = (
            "QUANTSEQ" in norm
            or "QUANT_SEQ" in norm
            or "3_MRNA_SEQ" in norm
            or "3PRIME_MRNA" in norm
            or "3_PRIME_MRNA" in norm
            or "3UTR_SEQ" in norm
            or "3_UTR_SEQ" in norm
            or ("LEXOGEN" in norm and ("3UTR" in norm or "3_UTR" in norm or "3PRIME" in norm or "3_PRIME" in norm or "3_MRNA" in norm))
        )

        if is_quantseq:
            # Distinguish REV vs FWD vs UNKNOWN_STRAND using word-boundary tokens
            has_rev = bool(re.search(r"(?:^|[_ -])(?:REV|REVERSE|ANTISENSE)(?:$|[_ -])", norm))
            has_fwd = bool(re.search(r"(?:^|[_ -])(?:FWD|FORWARD|SENSE)(?:$|[_ -])", norm))

            if has_rev and not has_fwd:
                return cls.LEXOGEN_QUANTSEQ_3UTR_REV_UMI if has_umi else cls.LEXOGEN_QUANTSEQ_3UTR_REV
            if has_fwd and not has_rev:
                return cls.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI if has_umi else cls.LEXOGEN_QUANTSEQ_3UTR_FWD

            # If QuantSeq orientation is not documented, keep it unknown!
            # Do NOT assume FWD. Require verified empirical evidence before selecting counting orientation.
            return cls.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND_UMI if has_umi else cls.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND

        # Conventional bulk patterns
        conventional_tokens = {
            "CONVENTIONAL",
            "STANDARD",
            "BULK",
            "BULK_RNA_SEQ",
            "TRUSEQ",
            "NEBNEXT",
            "POLYA",
            "TOTAL_RNA",
            "RIBO_ZERO",
            "SMART_SEQ",
            "ILLUMINA_STRANDED",
            "DUTP",
        }
        for token in conventional_tokens:
            if token in norm:
                return cls.CONVENTIONAL_BULK

        return cls.UNKNOWN

    @property
    def is_quantseq(self) -> bool:
        """Whether this library type is a Lexogen QuantSeq 3' mRNA-seq variant."""
        return self in (
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV_UMI,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND_UMI,
        )

    @property
    def is_fwd(self) -> bool:
        """Whether reads represent the sense (forward) transcript strand."""
        return self in (
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
        )

    @property
    def is_rev(self) -> bool:
        """Whether reads represent the antisense (reverse) transcript strand."""
        return self in (
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV_UMI,
        )

    @property
    def has_umi(self) -> bool:
        """Whether library incorporates unique molecular identifiers (UMIs)."""
        return self in (
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV_UMI,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_UNKNOWN_STRAND_UMI,
        )

    @property
    def expected_strandedness_flag(self) -> int | None:
        """Expected featureCounts strandedness flag (-s): 1=forward, 2=reverse, None=variable/unspecified."""
        if self.is_fwd:
            return 1
        if self.is_rev:
            return 2
        return None

    @property
    def requires_length_normalization(self) -> bool:
        """Whether within-sample abundance comparisons (e.g. TPM) require transcript-length adjustment.
        
        Note: Differential expression engines like DESeq2 ALWAYS consume raw unnormalized integer counts
        for both conventional bulk and QuantSeq. However, for within-sample descriptive abundance metrics:
        - Conventional bulk RNA-seq fragments across the entire transcript, requiring transcript-length
          normalization (TPM/FPKM) to compare relative abundance across genes.
        - QuantSeq 3' mRNA-seq yields exactly 1 read per transcript molecule at the 3' end, so counts
          are inherently length-independent (proportional to CPM without gene length division).
        """
        return not self.is_quantseq
