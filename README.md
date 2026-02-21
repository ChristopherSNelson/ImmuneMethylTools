# ImmuneMethylTools

A DNA Methylation QC Pipeline for Autoimmune Disease Research.

ImmuneMethylTools is a robust pipeline designed to identify and mitigate common wet-lab and bioinformatics artifacts in DNA methylation data. It provides clinical-grade quality control, surgical data cleaning, and leverages advanced statistical modeling (OLS/Wilcoxon) to isolate true biological signals.

## The 6-Gate Architecture

The pipeline processes data through six rigorous validation gates to prevent error propagation:

1. **Sample QC:** Bisulfite conversion rate & mean depth filters.
2. **Signal Quality:** Contamination detection (Bimodality) & Sex-metadata (XCI) verification.
3. **Sample Integrity:** Technical duplicate detection (Pearson correlation).
4. **Surgical Remediation:** Site-level depth filtering & VDJ-locus beta masking for clonal samples.
5. **Normalization:** Confound detection (Batch/Age/Sex) followed by median-centring.
6. **Inference:** Covariate-adjusted DMR calling (OLS) & Leak-free ML validation.

## Quick Start

```bash
# 1. Setup Environment
bash setup.sh
source venv/bin/activate

# 2. Generate Mock Data
python data/generate_mock_data.py

# 3. Run Pipeline (Generates PDF Report)
python core/pipeline.py --report
```

## The Artifact Map: What we catch

ImmuneMethylTools is built to survive "stumper" artifacts common in immune-cell studies:

| # | Artifact | Key Signal | Remediation |
|---|----------|------------|-------------|
| 1 | **Bisulfite Failure** | non-CpG meth > 2% | Drop Sample |
| 2 | **Contamination** | "Muddy" beta (near 0.5) | Drop Sample |
| 3 | **Sex-Metadata Mixup** | XCI signal mismatch | Drop Sample |
| 4 | **Duplicates** | Pearson r > 0.99 | Drop Replicate |
| 5 | **Clonal VDJ** | High Beta + Long Fragment | **Mask Locus** (Keep Sample) |
| 6 | **Batch Confound** | High Cramér's V (Batch × Disease) | Median-Centring |

---

**For the full Technical Manual, including SOPs on Masking vs. Dropping, Cell-Type Deconvolution, and Configurable Thresholds, see [TECHNICAL_GUIDE.md](TECHNICAL_GUIDE.md).**
