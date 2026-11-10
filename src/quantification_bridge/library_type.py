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
    LEXOGEN_QUANTSEQ_3UTR_FWD_UMI = "LEXOGEN_QUANTSEQ_3UTR_FWD_UMI"
    LEXOGEN_QUANTSEQ_3UTR_REV_UMI = "LEXOGEN_QUANTSEQ_3UTR_REV_UMI"
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

        # Lexogen QuantSeq patterns
        is_quantseq = (
            "QUANTSEQ" in norm
            or "LEXOGEN" in norm
            or "3_MRNA_SEQ" in norm
            or "3PRIME" in norm
            or "3_PRIME" in norm
            or "3UTR" in norm
        )

        if is_quantseq:
            # Distinguish REV vs FWD
            if "REV" in norm or "ANTISENSE" in norm:
                return cls.LEXOGEN_QUANTSEQ_3UTR_REV_UMI if has_umi else cls.LEXOGEN_QUANTSEQ_3UTR_REV
            if "FWD" in norm or "FORWARD" in norm or "SENSE" in norm:
                return cls.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI if has_umi else cls.LEXOGEN_QUANTSEQ_3UTR_FWD
            # If QuantSeq is mentioned without explicit FWD or REV, Lexogen standard kit is FWD
            # but if UMI is present, assign FWD_UMI
            return cls.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI if has_umi else cls.LEXOGEN_QUANTSEQ_3UTR_FWD

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
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_FWD_UMI,
            LibraryType.LEXOGEN_QUANTSEQ_3UTR_REV_UMI,
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
        """Whether counts must undergo transcript-length normalization (e.g. TPM).
        
        QuantSeq 3' mRNA-seq produces exactly 1 read per transcript molecule at the 3' end,
        meaning counts are inherently independent of transcript length.
        """
        return not self.is_quantseq
