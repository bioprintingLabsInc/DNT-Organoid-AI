"""Governed Raw-Read Quantification Bridge.

Upstream preprocessing and quantification module operating strictly outside
frozen Stages 1-7, Stage 8C, and Stage 8A.1. Registered as the execution engine
for Route.RAW_READ_PIPELINE.
"""

from .pairing import FastqPair, FastqPairingError, determine_read_mate, extract_sample_stem, group_fastq_pairs
from .pipeline import QuantificationBridgePipeline

__version__ = "0.2.0"

__all__ = [
    "FastqPair",
    "FastqPairingError",
    "QuantificationBridgePipeline",
    "determine_read_mate",
    "extract_sample_stem",
    "group_fastq_pairs",
]
