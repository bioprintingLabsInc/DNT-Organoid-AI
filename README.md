# DNT Organoid AI

DNT Organoid AI is a research platform for developing and validating AI models that assess developmental-neurotoxicity-associated molecular responses in human brain organoid experiments. It is intended to support both public human brain-organoid datasets and future in-house experiments while preserving provenance, experimental context, and matched-control relationships.

The intended workflow is:

> Data sources → ingestion → format detection → metadata validation → format-specific QC/preprocessing → gene harmonization → normalization → treatment-versus-matched-control response construction → standardized features → model-ready dataset → AI training → calibration/validation → locked model → prospective prediction

The architecture supports bulk RNA-seq, scRNA-seq, snRNA-seq, raw FASTQ/SRA inputs, and scientifically usable processed expression matrices through modality- and format-appropriate routes. See [System Architecture](docs/architecture/SYSTEM_ARCHITECTURE.md) for the intended boundaries and responsibilities.

## Status and intended use

This repository is under research development. It is not a clinical, diagnostic, or regulatory decision system, and its outputs must not be used as such.
