# ImmuneMethylTools

Rigorous immune-cell DNA methylation QC pipeline for autoimmune disease research.
Detects and documents common wet-lab and bioinformatics artifacts before any
differential analysis is attempted.

---

## Setup

```bash
bash setup.sh          # create venv, install requirements
source venv/bin/activate
python data/generate_mock_data.py   # generate mock_methylation.csv
```

---

## Usage

### Run the full pipeline

```bash
python core/pipeline.py
```

Runs all eight artifact detector modules in the correct order and writes
all outputs to `output/`:

| File | Contents |
|------|----------|
| `output/clean_methylation.csv` | 34 QC-passed samples |
| `output/audit_log_pipeline_{ts}.csv` | Per-run DETECTED/INFO event log |
| `output/flagged_samples.csv` | Cumulative flagged-sample registry |
| `output/logs/pipeline_{ts}.log` | Full stdout mirror |
| `output/figures/` | QC, PCA, KDE, volcano, exclusion plots |
| `output/report_{ts}.pdf` | 8-section PDF report (`--report` flag) |

The `data/` directory is input-only (`mock_methylation.csv`).

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

### VDJ / Clonal-Risk DMR windows

`find_dmrs()` no longer excludes VDJ-region CpGs from sliding windows.
Instead, every window is annotated with two columns:

| Column | Type | Meaning |
|--------|------|---------|
| `n_vdj_cpgs` | int | Number of CpGs in the window that fall in a VDJ locus |
| `clonal_risk` | bool | `True` when any CpG in the window is a VDJ CpG |

Significant DMRs with `clonal_risk=True` are flagged **HIGH CLONALITY** in
stdout and the audit log. **The Analyst decides** whether to accept or
discard them — clonal expansion drives hypermethylation + long fragments in
VDJ loci, so a significant DMR there may reflect clonality rather than true
disease-associated differential methylation.

```python
# Analyst filter example — exclude VDJ-overlapping windows:
clean_dmrs = dmrs[~dmrs["clonal_risk"] & dmrs["significant"]]

# Or review them explicitly:
flagged = dmrs[dmrs["clonal_risk"] & dmrs["significant"]]
print(flagged[["window_id", "delta_beta", "p_adj", "n_vdj_cpgs"]])
```

### VDJ-locus beta masking (Stage 3.5)

After clonal samples are identified, `mask_clonal_vdj_sites()` sets
`beta_value = NaN` at every VDJ-region row belonging to a flagged sample.
Non-VDJ rows in those samples and all rows in clean samples are untouched.

Downstream modules handle the NaN safely:
- **normalizer** — `groupby.median()` skips NaN; NaN propagates into `beta_normalized`
- **dmr_hunter** — `pivot.fillna(global_mean)` imputes masked sites to ~0 (after median-centring)
- **ml_guard** — `pivot.fillna(per_cpg_mean)` imputes masked sites to ~0

Net effect: artificially inflated clonal VDJ betas (~0.97) are suppressed to
approximately zero in the feature matrices rather than driving DMR or ML signal.

---

### Site-level depth QC

`filter_site_quality(df, min_depth=5)` removes individual CpG rows below the
depth threshold before downstream analysis.  `find_dmrs` and `run_safe_model`
also accept `min_site_depth=5` for defense-in-depth.

---

## Tests

```bash
pytest tests/ -v
```

75 tests covering all seven simulated artifacts, all pipeline stages,
all eight detector modules, the io_utils safe loader, and the XCI guard.

---

## Artifact Map

| # | Sample | Artifact | Key signal |
|---|--------|----------|------------|
| 1 | Batch_01 Cases | Confounded batch | +0.1 mean beta shift |
| 2 | S001 / P001 | Clonal VDJ | beta > 0.8, fragment > 180 bp |
| 3 | S001, S002 | Bisulfite failure | non_cpg_meth_rate ≈ 0.05 |
| 4 | S010 ↔ S_DUP | Sample duplication | Pearson r = 0.9999 |
| 5 | S020 | Contamination | muddy beta, BC = 0.78 |
| 6 | S030 | Low coverage | mean depth ≈ 5x |
| 7 | S035, S036 | Sex metadata mixup | X-linked beta contradicts reported sex |
| — | cg00000300–310 | True biological DMR (positive control) | +0.25 beta shift in all Case samples |
