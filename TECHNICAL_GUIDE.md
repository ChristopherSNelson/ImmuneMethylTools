# ImmuneMethylTools: Technical Guide & SOP

This document provides a comprehensive technical guide to ImmuneMethylTools, a DNA Methylation QC pipeline designed for autoimmune disease research. It details the project's features, design choices, setup instructions, and the analytical workflow.

ImmuneMethylTools identifies and documents common wet-lab and bioinformatics artifacts before any differential analysis is attempted, ensuring data quality and robust biological conclusions. An ElasticNet regression is used to find real signals after data cleanup.

---

## Key Features & Design Choices

ImmuneMethylTools is built around several core principles and features to ensure high-quality methylation data analysis:

- Comprehensive Artifact Detection: Proactively identifies and mitigates over 7 classes of common issues that can confound methylation studies, including:
    - Bisulfite conversion failures
    - Low sequencing depth at sample and site levels
    - Sample contamination
    - Technical duplicates
    - Sex-metadata mixups (XCI signal discrepancies)
    - Clonal VDJ artifacts
    - Batch effects confounded with disease labels
- Surgical Data Cleaning: Employs intelligent, precise strategies for data remediation. Rather than outright dropping samples, it utilizes methods like masking clonal VDJ loci (setting `beta_value` to NaN) to preserve statistical power and sample integrity where possible.
- Robust Normalization: Incorporates median-centring normalization to correct for inter-sample variability and performs rigorous checks for batch/disease confounding to avoid introducing bias into results.
- Safe Machine Learning: Validates detected biological signals using ElasticNet Logistic Regression with GroupKFold cross-validation. This explicitly prevents data leakage by ensuring that no patient appears in both training and test folds, leading to more trustworthy AUC metrics.
- Configurable Thresholds: All critical analysis thresholds are exposed in `config.json` at the project root, allowing users to tune sensitivity for their specific dataset without modifying the source code.
- Architectural Mandate (Multi-Stage Pipeline): The analytical workflow follows a strict sequential processing order (QC -> Normalization -> DMR -> ML) to ensure that upstream data integrity is established before downstream analyses run, preventing error propagation.

---

## Setup

To get started with ImmuneMethylTools, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/ChristopherSNelson/ImmuneMethylTools.git
    cd ImmuneMethylTools
    ```
2. Initialize Environment and Generate Mock Data:
    The `setup.sh` script creates a virtual environment, installs all necessary Python packages (including `jupyter` and `ipykernel`), and `generate_mock_data.py` creates a synthetic dataset (`mock_methylation.csv`) containing various simulated artifacts for testing the pipeline.
    ```bash
    bash setup.sh
    source venv/bin/activate
    python data/generate_mock_data.py
    ```
    Ensure your virtual environment is activated before running any scripts or notebooks.

---

## Workflow: Running the Pipeline

The core functionality is orchestrated by `core/pipeline.py`, which runs all artifact detector modules in the correct sequence. Outputs are generated in the `output/` directory.

### Run the full pipeline

```bash
python core/pipeline.py              # run all stages; writes outputs to `output/`
```

Common Pipeline Flags:

- Generate PDF Report:
    ```bash
    python core/pipeline.py --report     # Also generates an 8-section PDF report
    ```
- Skip Figure Generation: (Useful for faster CI/debugging)
    ```bash
    python core/pipeline.py --no-figures
    ```
- Use Custom Configuration: Override default analysis thresholds defined in `config.json` (see "Configurable Analysis Thresholds" section below).
    ```bash
    python core/pipeline.py --config my_study_config.json
    ```

### Pipeline Outputs

All generated files are saved to the `output/` directory:

| File | Contents |
|------|----------|
| `output/clean_methylation.csv` | QC-passed, processed methylation data. |
| `output/audit_log_pipeline_{ts}.csv` | A detailed, timestamped log of all `DETECTED` and `INFO` events from the pipeline run. |
| `output/flagged_samples.csv` | A cumulative registry of all samples flagged across successive pipeline runs. |
| `output/logs/pipeline_{ts}.log` | A full mirror of the pipeline's `stdout` output, timestamped. |
| `output/figures/` | Contains various exploratory data analysis (EDA) figures, including QC metrics, PCA plots, KDEs, and volcano plots. |
| `output/report_{ts}.pdf` | The optional 8-section PDF report, generated when the `--report` flag is used. |

The `data/` directory is input-only and primarily contains `mock_methylation.csv`.

### Run individual modules

While `core/pipeline.py` orchestrates the full workflow, individual modules can be run or imported independently for specific analyses or debugging:

```python
from core.io_utils import data_path, load_methylation
from core.qc_guard import audit_quality, detect_contamination, filter_site_quality
from core.sample_audit import detect_duplicates
from core.dmr_hunter import find_dmrs

df            = load_methylation(data_path("mock_methylation.csv"))
clean_samples = audit_quality(df)
df_clean      = df[df["sample_id"].isin(clean_samples)]

dmrs = find_dmrs(df_clean, clean_samples)
```

---

## Artifact Map: The 7 Artifacts + True Biological Signal

ImmuneMethylTools identifies and addresses a range of artifacts to isolate true biological signals. The mock dataset includes these simulated artifacts and a positive control:

| # | Sample | Artifact | Key Signal | Detector Module |
|---|--------|----------|------------|-----------------|
| 1 | Batch_01 Cases | Confounded batch | +0.1 mean beta shift | `normalizer` (Cramér's V) |
| 2 | S003 / P003 | Clonal VDJ | beta > 0.8, fragment > 3 SD | `repertoire_clonality` |
| 3 | S001, S002 | Bisulfite failure | non_cpg_meth_rate ≈ 0.05 | `qc_guard` |
| 4 | S010 ↔ S_DUP | Sample duplication | Pearson r = 0.9999 | `sample_audit` |
| 5 | S020 | Contamination | muddy beta, low Bimodality Coefficient | `qc_guard` |
| 6 | S030 | Low coverage | mean depth ≈ 5x | `qc_guard` |
| 7 | S035, S036 | Sex metadata mixup | X-linked beta contradicts reported sex | `xci_guard` |
| — | cg00000300–310 | True biological DMR (positive control) | +0.25 beta shift in all Case samples | `dmr_hunter`, `ml_guard` |

Important Notes on Artifacts:
- S001/S002 are designed to fail bisulfite QC early (Stage 1a) to ensure the clonal VDJ artifact (S003/P003) is exercised later in the pipeline.
- S035 (true female XX, reported male XY) and S036 (true male XY, reported female XX) demonstrate sex-metadata discrepancies detected by X-linked methylation patterns.
- The "true" DMR signal (cg00000300–310) is designed to be robust, surviving all artifact filters to serve as a positive control.

---

## Design Choices: Masking vs. Dropping Clonal Samples

A critical design choice in ImmuneMethylTools is to mask clonally expanded samples in VDJ loci rather than outright dropping them. This is a deliberate, principled decision based on statistical rigor:

| Approach | Impact on Non-VDJ Sites | Statistical Impact | Rationale |
|----------|--------------------------|--------------------|-----------|
| Drop Sample | All rows deleted (non-VDJ biology lost) | Reduces N; can silently bias Case/Control ratios and statistical power. | Too aggressive; loses valuable data. |
| Mask VDJ Loci (Default) | Non-VDJ sites retained; VDJ rows set to `NaN` | Minimal N loss; `NaN` values are handled by downstream imputation to cohort mean. | Surgical precision: only artifactual loci are suppressed. Power preservation: avoids unnecessary reduction in statistical power. Ratio integrity: maintains Case/Control ratios. Safe imputation: `NaN`s are imputed to values that minimally distort downstream analyses. Full auditability: masking is logged, allowing review and override. |

Analyst Override: If you wish to exclude clonally expanded VDJ DMR windows from downstream analysis without removing the sample:

```python
clean_dmrs = dmrs[~dmrs["clonal_risk"] & dmrs["significant"]]
```

---

## Multi-Stage Data Integrity Pipeline

The analytical workflow is structured into sequential processing phases, each rigorously validating data quality and mitigating specific classes of artifacts. This layered approach is critical for preventing error propagation and ensuring the reliability of downstream analyses.

| Phase | Stage | Module | Artifact Class Addressed |
|-------|-------|--------|--------------------------|
| 1 | 1a–1c | `qc_guard`, `xci_guard` | Bisulfite failure, low depth, contamination, sex-metadata mismatch |
| 2 | 2 | `sample_audit` | Technical duplicates (Pearson r ≥ 0.99) |
| 2.5 | 2.5 | `qc_guard` | Site-level low depth (< 5 reads per CpG row) |
| 3 | 3–3.5 | `repertoire_clonality` | Clonal VDJ artifacts — surgical `NaN` masking at VDJ loci |
| 4 | 4 | `normalizer` | Batch × disease confound — median-centring after masking |
| 5 | 7 | `ml_guard` | Data leakage — GroupKFold CV where patient (not sample) is the CV unit |

Stage Ordering Rationale:
- Sample-level QC precedes batch correction: Correcting artifact-contaminated samples would absorb the artifact signal into the correction model, obscuring it rather than removing it.
- Deduplication precedes site-level depth filtering: A duplicate pair with low coverage at one replicate should not survive to inflate site-level QC counts.
- VDJ masking precedes normalization: Median-centring inflated clonal betas before masking would shift the cohort median upward; masking must occur first to ensure accurate normalization.
- GroupKFold cross-validation groups by patient: This prevents data leakage where paired samples from the same donor might be split across train/test folds, artificially inflating AUC metrics.

### Detailed Pipeline Stage Order

| Stage | Module | Action |
|-------|--------|--------|
| 1a | `qc_guard` | Bisulfite/depth sample QC → removes failed samples |
| 1b | `qc_guard` | Contamination detection → removes contaminated samples |
| 1c | `xci_guard` | XCI sex-signal check → removes sex-metadata mismatches |
| 2 | `sample_audit` | Duplicate removal → drops sample_b from each flagged pair |
| — | `visuals` | Exclusion accounting figure saved after Stage 2 |
| 2.5 | `qc_guard` | Site-level depth filter → removes rows with depth < 5 from df_clean |
| 3 | `repertoire_clonality` | Clonal VDJ scan → identifies `clonal_samples` |
| 3.5 | `repertoire_clonality` | `mask_clonal_vdj_sites` → sets `beta_value=NaN` at VDJ rows for clonal samples |
| 4 | `normalizer` | Confound check + median-centring → produces `df_norm` |
| 5 | `deconvolution` | Cell fractions + per-sample T/B/Treg table (Case-first) |
| 6 | `dmr_hunter` | Sliding-window DMR caller on `df_norm` (per-sample median per window, then Wilcoxon rank-sum); generates volcano plots |
| 7 | `ml_guard` | ElasticNet GroupKFold CV for robust classification |
| — | `report_gen` | Optional 8-section PDF report (`--report` flag) |

---

## VDJ / Clonal-Risk DMR Windows

The `find_dmrs()` function no longer strictly excludes VDJ-region CpGs from sliding windows. Instead, it provides granular annotation to empower the analyst's decision-making:

| Column | Type | Meaning |
|--------|------|---------|
| `n_vdj_cpgs` | int | Number of CpGs within the window that fall in a VDJ locus |
| `clonal_risk` | bool | `True` if any CpG in the window is located in a VDJ locus |
| `mean_gc` | float | Mean GC content across CpGs in the window |

Significant DMRs with `clonal_risk=True` are flagged as `HIGH CLONALITY` in the console output and the audit log. This is crucial because clonal expansion can drive hypermethylation and long fragments in VDJ loci, meaning a significant DMR in these regions might reflect clonality rather than true disease-associated differential methylation.

Analyst Filter Examples:

- Exclude VDJ-overlapping windows:
    ```python
    clean_dmrs = dmrs[~dmrs["clonal_risk"] & dmrs["significant"]]
    ```
- Explicitly review flagged windows:
    ```python
    flagged = dmrs[dmrs["clonal_risk"] & dmrs["significant"]]
    print(flagged[["window_id", "delta_beta", "p_adj", "n_vdj_cpgs", "mean_gc"]])
    ```

---

## VDJ-Locus Beta Masking (Stage 3.5)

Following the identification of clonal samples by `repertoire_clonality`, the `mask_clonal_vdj_sites()` function surgically sets `beta_value = NaN` for every VDJ-region CpG row belonging to a flagged sample. Crucially, non-VDJ rows in those samples and all rows in clean (non-clonal) samples remain untouched.

Downstream modules are designed to handle these `NaN` values gracefully:

| Module | NaN Handling Mechanism | Net Effect |
|--------|--------------------------|------------|
| `normalizer` | `groupby.median()` skips `NaN`s; `NaN`s propagate into `beta_normalized` | Accurately normalized data; masked values remain `NaN` after normalization. |
| `dmr_hunter` | `pivot.fillna(global_mean)` imputes masked sites to approximately zero after median-centring. | Artificially inflated clonal VDJ betas (~0.97) are suppressed to approximately zero in the feature matrices, preventing them from driving DMR signals. |
| `ml_guard` | `pivot.fillna(per_cpg_mean)` imputes masked sites to approximately zero after median-centring. | Similar to `dmr_hunter`, prevents artifactual signals from influencing ML feature matrices. |

This process ensures that artificially inflated clonal VDJ betas (which can be as high as ~0.97) are effectively suppressed to approximately zero in the analytical feature matrices, thereby preventing them from spuriously driving DMR or ML signals.

---

## Site-Level Depth QC

The `filter_site_quality(df, min_depth=5)` function is employed to remove individual CpG rows whose read depth falls below a specified `min_depth` threshold (defaulting to 5). This is a site-level filter, meaning a single sample can contribute both passing and failing rows. Excluding low-depth sites is crucial because they exhibit inflated binomial variance, thereby reducing noise in the `DMR-hunter` and `ML_guard` feature matrices.

Both `find_dmrs` and `run_safe_model` also accept `min_site_depth=5` as a parameter, providing an additional layer of validation and consistency across the pipeline stages.

---

## Cell-Type Deconvolution (Stage 5)

Methylation data from whole-blood or PBMC samples represents a mixture of cell types. This is a significant consideration because a Differential Methylation Region (DMR) call in a study, for instance, nominally targeting B cells, would be invalid if the Case cohort contained substantially more Regulatory T cells (Tregs) or fewer B cells than the controls. In such scenarios, the observed methylation differences would merely reflect differences in cell-type composition, not true disease biology.

Stage 5 of the pipeline provides two critical checks to address this:

| Function | Purpose |
|----------|---------|
| `estimate_cell_fractions(df)` | Estimates per-sample fractions for T cells, B cells, and Regulatory T cells (Tregs). |
| `detect_lineage_shift(df)` | Flags samples exhibiting unexpected FoxP3 or PAX5 methylation patterns, indicative of lineage shifts. |

### Algorithm

The mock implementation in `deconvolution.py` generates random but biologically plausible cell fractions (e.g., B: 55–80%, T: 10–30%, Treg: 2–8%) and validates against two canonical lineage marker loci:

| Marker | Locus (GRCh38) | Interpretation |
|--------|---------------|----------------|
| FoxP3 | chrX:49,250,136–49,255,701 | Typically hypomethylated in Tregs; hypermethylation suggests loss of Treg identity. |
| PAX5 | chr9:36,896,702–37,175,926 | Typically hypomethylated in B cells; hypermethylation indicates a B-cell lineage shift. |

Mock CpG Proxies: The mock dataset uses `cg00000001`–`cg00000005` as proxies for FoxP3 and `cg00000006`–`cg00000010` for PAX5.

### Production Adaptation for Real Data

For real-world applications, `estimate_cell_fractions()` would typically be replaced or augmented with established reference-based deconvolution methods:

- Houseman 2012 (constrained projection): Utilizes a 6-cell-type reference matrix derived from purified blood fractions (Granulocytes, CD4T, CD8T, B cells, NK, Monocytes). This method is implemented in R via `minfi::estimateCellCounts()`.
- EpiDISH: Offers both reference-free and reference-based deconvolution, capable of handling a broader range of immune subtypes, including Tregs and plasma cells, which might not be covered by other panels.
- MethylResolver: Performs deconvolution using constrained least squares with bootstrap confidence intervals, supporting the use of custom reference matrices.

The expected input schema remains consistent: a long-format DataFrame with `sample_id`, `cpg_id`, and `beta_value` (as used throughout the pipeline). No additional columns are required for these functions to operate.

---

## Configurable Analysis Thresholds

A key design feature of ImmuneMethylTools is the externalization of all critical analysis thresholds into `config.json` at the project root. This allows users and analysts to easily tune the pipeline's sensitivity and behavior for their specific datasets without needing to modify any source code.

Usage Examples:

- Run with the default `config.json`:
    ```bash
    python core/pipeline.py --report
    ```
- Run with a custom configuration file:
    ```bash
    python core/pipeline.py --config my_study_config.json
    ```

Default Thresholds:

The default thresholds are chosen to be biologically conservative, offering a robust starting point. However, real-world datasets may require adjustment for optimal performance:

| Section | Parameter | Default | Notes |
|---------|-----------|---------|-------|
| `qc` | `bisulfite_fail_thresh` | 0.01 | Non-CpG methylation rate; threshold for bisulfite conversion failure. Can be increased to 0.02 for less stringent QC. |
| `qc` | `depth_fail_thresh` | 10 | Mean sample depth cutoff; samples below this are considered low coverage. |
| `qc` | `site_depth_thresh` | 5 | Per-CpG row depth; individual CpG rows below this are removed at Stage 2.5. |
| `qc` | `bc_sigma_thresh` | 2.0 | Bimodality Coefficient (BC) z-score; samples with BC below this many standard deviations from cohort median are flagged as contaminated. |
| `duplicates` | `corr_thresh` | 0.99 | Pearson correlation coefficient threshold; sample pairs with `r` at or above this are flagged as technical duplicates. Can be raised to 0.995 for noisier data. |
| `clonality` | `beta_min` | 0.80 | VDJ hypermethylation cutoff; beta values above this in VDJ loci indicate clonal expansion. |
| `clonality` | `frag_sd_thresh` | 3.0 | Fragment length outlier threshold; flags fragments > `sample_mean + this * sample_std`. This is RRBS-safe and adaptive. |
| `clonality` | `min_locus_hits` | 3 | Minimum number of dual-criteria hits required within a single VDJ locus to flag a sample as clonal. |
| `dmr` | `p_adj_thresh` | 0.05 | Benjamini-Hochberg (BH)-corrected p-value cutoff for calling significant DMRs. |
| `dmr` | `delta_beta_min` | 0.10 | Minimum absolute ΔBeta (Case − Control) required for a biologically meaningful DMR. |
| `dmr` | `chunk_size` | `null` | Number of CpGs to pivot per chunk for DMR analysis (e.g., 50000 for EPIC arrays). `null` (default) means load all CpGs into memory. |
| `ml` | `n_top_cpgs` | 200 | Number of top-variance CpGs selected for the ElasticNet feature matrix. |
| `ml` | `l1_ratio` | 0.5 | ElasticNet mixing parameter (0 = Ridge, 1 = Lasso). |
| `ml` | `c_param` | 1.0 | Inverse of regularization strength. |
| `ml` | `chunk_size` | `null` | Number of CpGs to load per chunk when computing CpG variance for feature selection. `null` (default) means in-memory. |

Missing keys in `config.json` will automatically fall back to the default values defined in `core/config_loader.py`.

---

## Interactive Notebook

The `notebooks/ImmuneMethylTools_Validation.ipynb` is a fully executable and interactive end-to-end demonstration of the pipeline's capabilities. It provides a transparent view into each stage of the analysis, with all underlying logic imported directly from the `core/` modules (ensuring no math or complex logic is reimplemented in the notebook cells).

To launch the notebook:

```bash
jupyter notebook notebooks/ImmuneMethylTools_Validation.ipynb
```

The notebook guides you through six key steps:

| Step | Demonstrates |
|------|-------------|
| Step 0 | Data loading and schema validation using `io_utils.load_methylation()`, ensuring data integrity from the start. |
| Step 1 | Sample-level QC checks, including bisulfite conversion rates, detection of contamination via bimodality, and sex-metadata XCI signal consistency. |
| Step 2 | Technical duplicate detection, illustrated with a Pearson `r` heatmap highlighting the `S010 ↔ S_DUP` simulated pair. |
| Step 3 | Clonal artifact handling, showcasing VDJ beta values before and after masking (`df_pre_site` vs. `df_masked` heatmaps) to visually confirm the surgical precision of the masking operation. |
| Step 4 | Normalization process, including the detection of batch × disease confounding using Cramér's V, and the impact of median-centring on beta distributions. |
| Step 5 | DMR calling and Machine Learning validation, featuring a volcano plot of Differential Methylation Regions and an ElasticNet AUC of `1.0000` on the true biological signal (positive control). |

---

## Tests

The project includes a comprehensive test suite to ensure the reliability and correctness of all modules.

To run the tests:

```bash
pytest tests/ -v
```

The suite comprises 75 tests covering:
- All seven simulated artifacts.
- All pipeline stages and the interactions between them.
- All eight detector modules.
- The `io_utils` safe loader and schema validator.
- The XCI sex-signal filter.

The mock data schema, including its 14 columns, is detailed in `CLAUDE.md` under the "Key Variable Glossary."

---

### Covariate-Adjusted DMR Calling (Stage 6)

The `find_dmrs()` function supports two statistical modes:

1.  **Wilcoxon Rank-Sum (Default):** A non-parametric test robust to non-Gaussian distributions, used when no covariates are provided.
2.  **OLS Linear Modeling (Covariate-Adjusted):** When `covariate_cols` (e.g., `["age", "sex"]`) are provided via `config.json`, the pipeline fits an Ordinary Least Squares model per cluster.
    - **M-Value Scale:** Beta values are logit-transformed to M-values to satisfy OLS normality and homoscedasticity assumptions.
    - **Formula:** `M-value ~ disease_label + age + sex + ...`
    - **Effect Size:** Delta-beta is reported on the original [0,1] scale for biological interpretability, while p-values are derived from the M-value model.

---

## TODO / Future Work

The ImmuneMethylTools pipeline is continuously evolving. Here are key areas identified for future development:

### Adapting to Real EPIC / WGBS Data

The provided mock data (`generate_mock_data.py`) simulates a simplified cohort. Applying ImmuneMethylTools to real EPIC (Illumina MethylationEPIC arrays) or WGBS (Whole-Genome Bisulfite Sequencing) data requires the following adaptations:

- Input Format: Replace `mock_methylation.csv` with a tidy long-format CSV exported from standard tools like `minfi` (for EPIC arrays) or `Bismark`/`Bismark2bedGraph` (for WGBS). The robust schema validator in `io_utils.load_methylation()` will automatically catch missing or mis-typed columns during loading.
- CpG Count Management: Real EPIC arrays contain approximately 850,000 CpGs, while WGBS datasets can exceed 25 million. Both `find_dmrs` and `run_safe_model` are designed to support chunked processing to avoid materializing the entire CpG matrix in memory. This is configured by setting `chunk_size` in `config.json`:
    ```json
    "dmr": { "chunk_size": 50000 },
    "ml":  { "chunk_size": 50000 }
    ```
    At 100 samples, each 50K-CpG chunk typically uses around 40 MB of memory. Without chunking, a full in-memory EPIC matrix would require approximately 680 MB, and a WGBS dataset with 25M CpGs would necessitate a `chunk_size` of up to 100,000 to manage memory effectively.
- VDJ Coordinates: Update the `IS_VDJ_REGION` logic within `repertoire_clonality.py` to utilize accurate IGH/IGK/IGL/TRA/TRB locus coordinates from the GRCh38 human reference genome. The current mock flag in the simulated data will need to be replaced with real genomic intervals.
- Bisulfite Threshold: The default 2% non-CpG methylation rate cutoff is appropriate for WGBS data. For EPIC arrays, bisulfite conversion quality is typically assessed via control probes. This will require replacing the `non_cpg_meth_rate` metric with array control-probe metrics for accurate QC.
- Sex Inference: The `xci_guard` module infers sex based on X-linked CpG beta values. It is essential to ensure that the input data includes a representative set of X-linked probes (which are present on EPIC arrays by default) and to confirm adequate chrX coverage for WGBS datasets.
