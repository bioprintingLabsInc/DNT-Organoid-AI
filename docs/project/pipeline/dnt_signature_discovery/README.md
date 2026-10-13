# DNT Signature Discovery

**Historical Identifier**: Stage 11  
**Scientific Name**: Developmental Neurotoxicity Transcriptomic Signature Discovery  
**Status**: **STATUS: PLANNED — NO IMPLEMENTATION YET**

---

## 1. Scientific Objective & Leakage Controls

The purpose of DNT Signature Discovery is to extract sparse, interpretable, cross-study reproducible gene and pathway signatures that reliably distinguish DNT-active compounds from negative reference chemicals across early human organoid neurodevelopment.

> [!CRITICAL]
> **Strict Partitioning Invariant**:
> * Signature discovery occurs **strictly after compound/study-safe dataset partitioning**;
> * Signature discovery uses the **training partition only**;
> * Calibration compounds and untouched test chemicals **cannot influence signature selection, ranking, thresholding, or weight optimization**.

---

## 2. Intended Inputs & Outputs

* **Intended Upstream Input**: Training partition feature vectors (from `pipeline/compound_safe_dataset_splitting/`) + independent training DNT reference labels.
* **Intended Output**: Candidate biomarker gene panels, pathway signature weights, and biological mechanistic perturbation profiles learned without test-set leakage.

---

## 3. Implementation Plan

* Implementation will reside in `src/models/signatures/` and `tests/models/signatures/`.
* No code paths or algorithms are currently implemented.
