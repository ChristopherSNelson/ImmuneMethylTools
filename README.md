# ImmuneMethylTools

Rigorous immune-cell DNA methylation QC pipeline for autoimmune disease research.
Detects and documents common wet-lab and bioinformatics artefacts before any
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

Runs all eight artifact detector modules in the correct order, writes
`data/clean_methylation.csv` (36 clean samples), a per-run audit log
`data/audit_log_pipeline_{ts}.csv`, and a run log `logs/pipeline_{ts}.log`.

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

### Site-level depth QC

`filter_site_quality(df, min_depth=5)` removes individual CpG rows below the
depth threshold before downstream analysis.  `find_dmrs` and `run_safe_model`
also accept `min_site_depth=5` for defense-in-depth.

---

## Tests

```bash
pytest tests/ -v
```

64 tests covering all six simulated artifacts, all eight modules, the
io_utils safe loader, and the end-to-end pipeline.

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
