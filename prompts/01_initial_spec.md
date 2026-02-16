# Initial Specification — ImmuneMethylTools

**Date:** 2026-02-15
**Roles:** Senior Python ML Engineer (executor) | Principal Bioinformatics Architect (reviewer)
**Architect:** Christopher S. Nelson <christopher.s.nelson.01@gmail.com>

**Commit authorship:** Every commit must list both co-authors:
  `Co-Authored-By: Christopher S. Nelson <christopher.s.nelson.01@gmail.com>`
  `Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>`

---

## Project Goal

Build a specialized repository `ImmuneMethylTools` demonstrating a rigorous, IC-level approach to analyzing B-cell/T-cell methylation data for autoimmune disease. Must address common pitfalls: batch effects, clonal expansion artifacts (VDJ), confounding, and sample integrity.

---

## Core Directives

- **Bisulfite Intuition:** Check for good bisulfite methylation intuition at every step.
- **Execution Rule:** Follow the plan strictly. Stop after each phase for review.
- **Visuals:** Generate "Before vs. After" visualizations for every data manipulation.
- **Mandatory `__main__` Block:** Every Python script in `core/` MUST include a `if __name__ == "__main__":` block that loads `data/mock_methylation.csv`, runs the module's logic, and prints a timestamped detection log to stdout for every issue found. Format: `[YYYY-MM-DD HH:MM:SS] [MODULE] ✓ DETECTED | <artifact> | <key metric>`

---

## PHASE 0: PROMPT ARCHIVING & SETUP

- Create `prompts/` directory.
- Save this prompt to `prompts/01_initial_spec.md`.
- Create `prompts/LOG.md` to track every instruction (local machine execution).

---

## PHASE 1: INFRASTRUCTURE

1. Initialize repo: `ImmuneMethylTools`.
2. **Config files:**
   - `CLAUDE.md`: Project context, architecture decisions, variable names.
   - `.gitignore` & `LICENSE` (MIT).
   - `requirements.txt`: pandas, numpy, scipy, scikit-learn, matplotlib, seaborn.
   - `setup.sh`: Script to create venv and install requirements.
3. **Directory structure:**
   - `data/` — mock inputs
   - `core/` — detective modules
   - `notebooks/` — final demo
   - `tests/` — sanity checks

---

## PHASE 2: THE "LAB SIMULATOR" (Mock Data)

Write `data/generate_mock_data.py` to produce `mock_methylation.csv`.

### Columns

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `patient_id` | Patient identifier (multiple samples per patient) |
| `batch_id` | Batch label (Batch_01, Batch_02, ...) |
| `age` | Patient age |
| `disease_label` | Case / Control |
| `cpg_id` | CpG site identifier |
| `beta_value` | Methylation beta value [0, 1] |
| `depth` | Sequencing depth at site |
| `fragment_length` | Insert/fragment length (bp) |
| `is_vdj_region` | Boolean — site in VDJ recombination region |
| `non_cpg_meth_rate` | Non-CpG methylation rate (proxy for bisulfite conversion) |

### Required "Stumper" Artifacts

1. **Confounded Batch:** `Batch_01` contains 80% Case samples; shift mean methylation +0.1.
2. **Clonal Artifact:** One Case patient has a VDJ region (`is_vdj_region=True`) with Beta > 0.8 AND Fragment Length > 180 bp.
3. **Bisulfite Failure:** 2 samples with `non_cpg_meth_rate` > 0.02 (2%).
4. **Duplication:** 2 samples with different IDs but inter-sample correlation > 0.99.
5. **Contamination:** 1 sample with a "muddy" beta distribution (peaks shifting toward 0.5).

---

---

## PHASE 3: THE "DETECTIVE" MODULES (Core Logic)

Implement as clean, functional Python scripts in `core/`. Add docstrings explaining biological intent.

### 3.1 `visuals.py` — The EDA Suite
- `plot_qc_metrics(df)`: Histograms for Conversion Rate, Mean Depth, and Global Methylation Levels.
- `plot_beta_distribution(df)`: KDE plot of Beta values per sample (spot "muddy" contamination).
- `plot_pca(df, title, color_by)`: Standard PCA scatter plot.

### 3.2 `qc_guard.py` — The Gatekeeper
- `audit_quality(df)`: Flag samples with `non_cpg_meth_rate > 0.01` or `depth < 10`.
- `detect_contamination(df)`: Cohort-relative bimodality coefficient (BC) check — flag samples > 2σ below cohort BC median AND mean beta in muddy range [0.40, 0.65].
- **Output:** Returns `clean_samples_list`.

### 3.3 `sample_audit.py` — The Integrity Check
- `detect_duplicates(df)`: Pairwise correlation on top-100 high-variance CpGs. Flag pairs > 0.99.
- `snv_concordance_placeholder(df)`: Interface stub for SNV-based identity check.

### 3.4 `normalizer.py` — The Normalizer
- `check_confounding(df, col1, col2)`: Chi-square / Cramér's V check between two categorical columns.
- `robust_normalize(df)`: Median centering — subtract sample median from each beta value. Save Before/After figure.

### 3.5 `repertoire_clonality.py` — The Lineage Guard
- GRCh38 coordinates for B-cell (IGH, IGK, IGL) and T-cell (TRA, TRB, TRG, TRD) loci documented.
- `flag_clonal_artifacts(df)`: Flag VDJ rows with beta > 0.8 AND fragment > 180 bp.
- `get_vdj_summary(df)`: Per-patient VDJ methylation summary.

### 3.6 `deconvolution.py` — The Cell-Type Guard
- `estimate_cell_fractions(df)`: Mock function generating T/B/Treg fractions.
- `detect_lineage_shift(df)`: Check methylation at FoxP3 (Treg) and PAX5 (B-cell) proxy loci.

### 3.7 `dmr_hunter.py` — The Strict Analyst
- **Safety:** Assert `df` contains ONLY `clean_samples`.
- **Filter:** Exclude `is_vdj_region` CpGs.
- **Stats:** Sliding window (size=5, step=1) Wilcoxon rank-sum + BH FDR correction.
- **Criteria:** p_adj < 0.05, |ΔBeta| > 0.10, ≥ 3 CpGs per window.

### 3.8 `ml_guard.py` — The Validator
- `run_safe_model(df)`: `LogisticRegression` (ElasticNet, l1_ratio=0.5) with `GroupKFold` (n=5). Uses top-200 high-variance CpGs. Grouped by `patient_id` to prevent data leakage.

---

## Global Rules (added during Phase 3)

### Authorship — every commit
```
Co-Authored-By: Christopher S. Nelson <christopher.s.nelson.01@gmail.com>
Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

### Style
- **American English spelling** throughout. normalize (not normalise), color (not colour).

---

## Status

- [x] Phase 0 complete — commit e91d292
- [x] Phase 1 complete — commit e91d292
- [x] Phase 2 complete — commit e91d292 (10/10 tests passing)
- [x] Phase 3 complete — commit b7fb85e (8 modules) / 6307753 (36 tests) / 569cb35 (spelling)
- [ ] Phase 4 pending — architect review required
