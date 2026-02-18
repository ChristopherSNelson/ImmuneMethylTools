# ImmuneMethylTools

Immune-cell DNA methylation QC pipeline for autoimmune disease research.
Detects and documents common wet-lab and bioinformatics artifacts before any
differential analysis is attempted. An elasticNet regression is used to find real signal after data clean up.

---

## Setup

```bash
bash setup.sh          # create venv, install requirements (includes jupyter, ipykernel)
source venv/bin/activate
python data/generate_mock_data.py   # generate mock_methylation.csv
```

---

## Usage

### Run the full pipeline

```bash
python core/pipeline.py              # run all stages; write outputs
python core/pipeline.py --report     # also generate 8-section PDF report
python core/pipeline.py --no-figures # skip figure generation (faster CI)
```

Runs all eight artifact detector modules and writes all outputs to `output/`:

| File | Contents |
|------|----------|
| `output/clean_methylation.csv` | 34 QC-passed samples |
| `output/audit_log_pipeline_{ts}.csv` | Per-run DETECTED/INFO event log |
| `output/flagged_samples.csv` | Cumulative flagged-sample registry |
| `output/logs/pipeline_{ts}.log` | Full stdout mirror |
| `output/figures/` | QC, PCA, KDE, volcano, exclusion plots |
| `output/report_{ts}.pdf` | 8-section PDF report (`--report` flag) |

The `data/` directory is input-only (`mock_methylation.csv`). <!-- contains injected artifact simulations -->

### Run individual modules

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

| # | Sample | Artifact | Key Signal | Detector Module |
|---|--------|----------|------------|-----------------|
| 1 | Batch_01 Cases | Confounded batch | +0.1 mean beta shift | `normalizer` (Cramér's V) |
| 2 | S003 / P003 | Clonal VDJ | beta > 0.8, fragment > 180 bp | `repertoire_clonality` |
| 3 | S001, S002 | Bisulfite failure | non_cpg_meth_rate ≈ 0.05 | `qc_guard` |
| 4 | S010 ↔ S_DUP | Sample duplication | Pearson r = 0.9999 | `sample_audit` |
| 5 | S020 | Contamination | muddy beta, BC = 0.78 | `qc_guard` |
| 6 | S030 | Low coverage | mean depth ≈ 5x | `qc_guard` |
| 7 | S035, S036 | Sex metadata mixup | X-linked beta contradicts reported sex | `xci_guard` |
| — | cg00000300–310 | True biological DMR (positive control) | +0.25 beta shift in all Case samples | `dmr_hunter`, `ml_guard` |

Notes:
- S001/S002 fail bisulfite QC (Stage 1a) before the clonal scan (Stage 3), so the
  clonal VDJ artifact is placed at S003/P003 to ensure detection is exercised.
- S035 (true female XX, reported male XY) and S036 (true male XY, reported female XX):  detected by X-linked methylation in females. The finding contradicts sample metadata, a common simple sample ID check method.
- The "true" DMR signal (cg00000300–310) is autosomal, non-VDJ, non-X-linked;
  it survives every artifact filter and provides a positive control.

---

## SOP: Masking vs. Dropping Clonal Samples

When a sample is flagged as clonally expanded (VDJ locus, high methylation beta, longer fragments due to excluded allele methylation protection),
the pipeline masks rather than drops that sample. This is a deliberate,
principled choice:

| Approach | What happens to non-VDJ sites? | Statistical impact |
|----------|--------------------------------|--------------------|
| Drop sample | All rows deleted — non-VDJ biology lost | Reduces N; can bias case/control ratio silently |
| Mask VDJ loci (default) | Non-VDJ sites retained; VDJ rows set to NaN | Minimal N loss; downstream imputation to cohort mean |

Why masking was chosen:

1. Surgical precision: only the artifact loci (VDJ rows) are suppressed to NaN;
   non-VDJ CpGs, cell-fraction estimates, and sample metadata remain in analysis.
2. Power preservation: dropping removes every observation from a sample,
   reducing statistical power at all CpG sites, not just the contaminated ones.
3. Ratio integrity: dropping a Case sample silently changes the Case:Control
   ratio, which can shift model calibration without any audit trail.
4. Safe imputation: after median-centring, masked NaN values impute to
   approximately zero (cohort mean), causing minimal distortion to downstream
   DMR calling and ML feature matrices.
5. Full auditability: the audit log records exactly which sample × locus
   combinations were masked, so analysts can review and override if needed.

Analyst override: if you want to exclude clonally expanded VDJ DMR windows
from downstream analysis without removing the sample:

```python
clean_dmrs = dmrs[~dmrs["clonal_risk"] & dmrs["significant"]]
```

---

## Multi-Stage Data Integrity Pipeline

The analytical workflow is structured into sequential processing phases, each
validating data quality and mitigating a specific class of artifact. This layered
approach ensures that upstream data integrity is established before downstream
analyses run, preventing error propagation.

| Phase | Stage | Module | Artifact Class Addressed |
|-------|-------|--------|--------------------------|
| 1 | 1a–1c | `qc_guard`, `xci_guard` | Bisulfite failure, low depth, contamination, sex-metadata mismatch |
| 2 | 2 | `sample_audit` | Technical duplicates (Pearson r ≥ 0.995) |
| 2.5 | 2.5 | `qc_guard` | Site-level low depth (< 5 reads per CpG row) |
| 3 | 3–3.5 | `repertoire_clonality` | Clonal VDJ artifacts — surgical NaN masking at VDJ loci |
| 4 | 4 | `normalizer` | Batch × disease confound — median-centring after masking |
| 5 | 7 | `ml_guard` | Data leakage — GroupKFold CV where patient (not sample) is the CV unit |

Stage ordering rationale:

- Sample-level QC precedes batch correction: correcting artifact-contaminated
  samples absorbs the artifact signal into the correction model, obscuring it.
- Deduplication precedes site-level depth filtering: a duplicate pair with low
  coverage at one replicate should not survive to inflate site-level QC counts.
- VDJ masking precedes normalization: median-centring inflated clonal betas before
  masking would shift the cohort median upward; masking must occur first.
- GroupKFold cross-validation groups by patient, not sample: splitting on sample_id
  leaks paired samples from the same donor across train/test folds, inflating AUC.

---

## Pipeline Stage Order

| Stage | Module | Action |
|-------|--------|--------|
| 1a | `qc_guard` | Bisulfite/depth sample QC → removes failed samples |
| 1b | `qc_guard` | Contamination detection → removes contaminated samples |
| 1c | `xci_guard` | XCI sex-signal check → removes sex-metadata mismatches |
| 2 | `sample_audit` | Duplicate removal → drops sample_b from each flagged pair |
| — | `visuals` | Exclusion accounting figure saved after Stage 2 |
| 2.5 | `qc_guard` | Site-level depth filter → removes rows with depth < 5 from df_clean |
| 3 | `repertoire_clonality` | Clonal VDJ scan → `clonal_samples` list |
| 3.5 | `repertoire_clonality` | `mask_clonal_vdj_sites` → `beta_value=NaN` at VDJ rows for clonal samples |
| 4 | `normalizer` | Confound check + median-centring → df_norm |
| 5 | `deconvolution` | Cell fractions + per-sample T/B/Treg table (Case-first) |
| 6 | `dmr_hunter` | Sliding-window DMR caller on df_norm; volcano plots saved |
| 7 | `ml_guard` | ElasticNet GroupKFold CV |
| — | `report_gen` | Optional 8-section PDF report (`--report` flag) |

---

## VDJ / Clonal-Risk DMR Windows

`find_dmrs()` no longer excludes VDJ-region CpGs from sliding windows.
Instead, every window is annotated with two columns:

| Column | Type | Meaning |
|--------|------|---------|
| `n_vdj_cpgs` | int | Number of CpGs in the window that fall in a VDJ locus |
| `clonal_risk` | bool | `True` when any CpG in the window is a VDJ CpG |

Significant DMRs with `clonal_risk=True` are flagged HIGH CLONALITY in
stdout and the audit log. The analyst decides whether to accept or
discard them — clonal expansion drives hypermethylation and long fragments in
VDJ loci, so a significant DMR there may reflect clonality rather than true
disease-associated differential methylation.

```python
# Analyst filter example — exclude VDJ-overlapping windows:
clean_dmrs = dmrs[~dmrs["clonal_risk"] & dmrs["significant"]]

# Or review them explicitly:
flagged = dmrs[dmrs["clonal_risk"] & dmrs["significant"]]
print(flagged[["window_id", "delta_beta", "p_adj", "n_vdj_cpgs"]])
```

---

## VDJ-Locus Beta Masking (Stage 3.5)

After clonal samples are identified, `mask_clonal_vdj_sites()` sets
`beta_value = NaN` at every VDJ-region row belonging to a flagged sample.
Non-VDJ rows in those samples and all rows in clean samples are untouched.

Downstream modules handle the NaN safely:

| Module | NaN handling |
|--------|-------------|
| `normalizer` | `groupby.median()` skips NaN; NaN propagates into `beta_normalized` |
| `dmr_hunter` | `pivot.fillna(global_mean)` imputes masked sites to ~0 after median-centring |
| `ml_guard` | `pivot.fillna(per_cpg_mean)` imputes masked sites to ~0 after median-centring |

Net effect: artificially inflated clonal VDJ betas (~0.97) are suppressed to
approximately zero in the feature matrices rather than driving DMR or ML signal.

---

## Site-Level Depth QC

`filter_site_quality(df, min_depth=5)` removes individual CpG rows below the
depth threshold before downstream analysis. `find_dmrs` and `run_safe_model`
also accept `min_site_depth=5` for an additional validation layer.

---

## Notebook

`notebooks/ImmuneMethylTools_Validation.ipynb` is a fully executable end-to-end
demonstration of the pipeline. All logic is imported from `core/`; no math is
reimplemented in notebook cells.

```bash
jupyter notebook notebooks/ImmuneMethylTools_Validation.ipynb
```

The notebook walks through six steps:

| Step | Demonstrates |
|------|-------------|
| Step 0 | Data loading and schema validation via `io_utils.load_methylation()` |
| Step 1 | Sample QC — bisulfite rates, contamination bimodality, sex-metadata XCI signal |
| Step 2 | Duplicate detection — Pearson r heatmap; S010 ↔ S_DUP pair highlighted |
| Step 3 | Clonal artifact — VDJ beta before/after masking (`df_pre_site` vs. `df_masked` heatmaps) |
| Step 4 | Normalization — batch × disease Cramér's V; beta distributions before and after median-centring |
| Step 5 | DMR calling and ML — volcano plot; ElasticNet AUC = 1.0000 on the true biological signal |

The before/after VDJ heatmap (Step 3) uses `df_pre_site` (raw clonal betas ≈ 0.97)
and `df_masked` (NaN imputed to ~0) to visually confirm the masking operation.

---

## Tests

```bash
pytest tests/ -v
```

75 tests covering all seven simulated artifacts, all pipeline stages,
all eight detector modules, the io_utils safe loader, and the XCI sex-signal filter.

---

## TODO / Future Work

### Age and sex as methylation covariates

Age and sex are currently collected in sample metadata but not modeled anywhere
in the pipeline. Both are known confounders:

- DNA methylation drifts systematically with age (the basis of epigenetic clocks
  such as Horvath 2013), so a Case/Control DMR call without age adjustment risks
  flagging age-associated CpGs as disease signal if case/control cohorts are not 
  age-matched.
- Sex drives a clean PC2 axis in the PCA even in this synthetic dataset — in real
  EPIC data, sex-dimorphic autosomal methylation is strong enough that any CpG
  mildly correlated with sex and also imbalanced across case/control groups will
  appear as a false positive without explicit adjustment.

Planned fix:

- Add Cramér's V / ANOVA checks for age × disease_label and sex × disease_label
  imbalance in `normalizer`, analogous to the existing batch × disease check.
- Replace the nonparametric Wilcoxon in `dmr_hunter` with a linear model
  (`beta ~ disease + age + sex + batch`) so both are proper covariates rather
  than uncontrolled nuisance variables.
- Expose `covariate_cols` as a parameter to `find_dmrs()` so analysts can
  pass additional covariates (age, sex, smoking status) without code changes.

### Adapting to real EPIC / WGBS data

The mock data (`generate_mock_data.py`) simulates a simplified cohort.
Applying ImmuneMethylTools to real data requires the following adaptations:

- Input format: replace `mock_methylation.csv` with a tidy long-format CSV
  exported from `minfi` (EPIC array) or `Bismark`/`Bismark2bedGraph` (WGBS).
  The schema validator in `io_utils.load_methylation()` will catch missing or
  mis-typed columns at load time.
- CpG count: real EPIC arrays have ~850,000 CpGs; WGBS can exceed 25M.
  The `dmr_hunter` sliding window and `ml_guard` pivot are memory-resident —
  profile and consider chunked processing or sparse matrices for WGBS scale.
- VDJ coordinates: update `IS_VDJ_REGION` logic in `repertoire_clonality.py`
  to use real IGH/IGK/IGL/TRA/TRB locus coordinates from GRCh38 (currently a
  mock flag in the simulated data).
- Bisulfite threshold: the 2% non-CpG methylation rate cutoff is appropriate
  for WGBS; for EPIC arrays bisulfite conversion is assessed via control probes
  — replace `non_cpg_meth_rate` with array control-probe metrics.
- Sex inference: `xci_guard` uses X-linked CpG beta values; ensure the input
  includes a representative set of X-linked probes (present on EPIC by default;
  confirm chrX coverage for WGBS).
