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
- **Timestamped Audit Log:** Every artifact detector module MUST append its findings to a persistent, timestamped CSV file at `output/audit_log_YYYYMMDD_HHMMSS.csv`. This log serves as the single source of truth for all data quality flags and exclusions for a specific analysis run.

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
   - `data/` — input data only (`mock_methylation.csv`)
   - `output/` — all generated outputs (logs, figures, audit CSVs, reports, clean data)
   - `core/` — artifact detector modules
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

## PHASE 3: THE ARTIFACT DETECTOR MODULES (Core Logic)

Implement as clean, functional Python scripts in `core/`. Add docstrings explaining biological intent.

### 3.1 `visuals.py` — The EDA Suite
- `plot_qc_metrics(df)`: Histograms for Conversion Rate, Mean Depth, and Global Methylation Levels.
- `plot_beta_distribution(df)`: KDE plot of Beta values per sample (spot "muddy" contamination).
- `plot_pca(df, title, color_by)`: Standard PCA scatter plot.

### 3.2 `qc_guard.py` — Sample QC Filter
- `audit_quality(df)`: Flag samples with `non_cpg_meth_rate > 0.01` or `depth < 10`.
- `detect_contamination(df)`: Cohort-relative bimodality coefficient (BC) check — flag samples > 2σ below cohort BC median AND mean beta in muddy range [0.40, 0.65].
- **Output:** Returns `clean_samples_list`.

### 3.3 `sample_audit.py` — Sample Integrity Auditor
- `detect_duplicates(df)`: Pairwise correlation on top-100 high-variance CpGs. Flag pairs > 0.99.
- `snv_concordance_placeholder(df)`: Interface stub for SNV-based identity check.

### 3.4 `normalizer.py` — The Normalizer
- `check_confounding(df, col1, col2)`: Chi-square / Cramér's V check between two categorical columns.
- `robust_normalize(df)`: Median centering — subtract sample median from each beta value. Save Before/After figure.

### 3.5 `repertoire_clonality.py` — Clonal Expansion & VDJ Artifact Detector
- GRCh38 coordinates for B-cell (IGH, IGK, IGL) and T-cell (TRA, TRB, TRG, TRD) loci documented.
- `flag_clonal_artifacts(df)`: Flag VDJ rows with beta > 0.8 AND fragment > 180 bp. Returns `(clonal_rows, flagged_samples)`.
- `mask_clonal_vdj_sites(df, clonal_samples)`: Set `beta_value=NaN` at VDJ-region rows for flagged samples only. Non-VDJ rows and non-clonal samples unchanged. Returns `(df_masked, n_masked)`. Called as Stage 3.5 in pipeline.
- `get_vdj_summary(df)`: Per-patient VDJ methylation summary.

### 3.6 `deconvolution.py` — The Cell-Type Guard
- `estimate_cell_fractions(df)`: Mock function generating T/B/Treg fractions.
- `detect_lineage_shift(df)`: Check methylation at FoxP3 (Treg) and PAX5 (B-cell) proxy loci.

### 3.7 `dmr_hunter.py` — Differential Methylation Region Caller
- **SampleQC:** Assert `df` contains ONLY `clean_samples`.
- **VDJ annotation (not exclusion):** All windows include VDJ CpGs; each window annotated with `n_vdj_cpgs` (int) and `clonal_risk` (bool). Analyst filters via `dmrs[~dmrs["clonal_risk"] & dmrs["significant"]]`.
- **Stats:** Sliding window (size=5, step=1); per-sample mean beta within each window, then Wilcoxon rank-sum on sample-level means + BH FDR correction. Each sample contributes one observation per window (no pseudoreplication).
- **Criteria:** p_adj < 0.05, |ΔBeta| > 0.10, ≥ 3 CpGs per window.
- **NaN handling:** `pivot.fillna(global_mean_beta_normalized)` before Wilcoxon — masked clonal VDJ sites imputed to ~0 (cohort mean after median-centring).

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
- [x] Phase 3 extensions (Sessions 7–11) — see details below
- [ ] Phase 4 pending — architect review required

## Phase 3 Extensions (added post-initial-spec, architect-approved)

### Session 7 — Artifact 6 + Infrastructure
- **Artifact 6 (Low Coverage):** S030 simulated at Poisson λ=5 (mean ≈ 5x depth); `audit_quality()` flags depth < 10 threshold; test added
- **`core/io_utils.py`:** `project_root()`, `data_path()`, `load_methylation()` (6-point schema validator), `Tee` (stdout mirror to log file), `append_flagged_samples()`, `write_audit_log()`
- **`core/pipeline.py`:** End-to-end runner; passes `clean_samples` list automatically stage-to-stage; exports `output/clean_methylation.csv`; audit log `output/audit_log_pipeline_{ts}.csv`
- **`tests/test_io_utils_and_pipeline.py`:** 15 integration tests; module-scoped pipeline fixture (64/64 tests total)

### Session 8 — Site QC, Clonal Annotation, Visualizations, T/B Logging
- **`filter_site_quality()`** in `qc_guard.py`: site-level depth filter (depth < 5 removes rows, not samples); returns `(df_clean, stats_dict)`; wired as Stage 2.5
- **`detect_duplicates` tuple return:** now returns `(pairs_df, ids_to_drop)` matching `flag_clonal_artifacts` pattern
- **VDJ clonal_risk annotation:** `dmr_hunter.find_dmrs()` no longer blanket-excludes VDJ CpGs; adds `n_vdj_cpgs` + `clonal_risk` columns; analyst decides via `dmrs[~dmrs["clonal_risk"]]`
- **Stage 5 T/B/Treg table:** full per-sample table sorted Case-first, B:T ratio, per-group means logged to pipeline stdout
- **`plot_exclusion_accounting()`:** two-panel (pie + waterfall) exclusion figure saved after Stage 2
- **`plot_volcano()`:** standard + clonal-risk coloring variants; saved after Stage 6
- **`README.md`** created at project root (Setup, Usage, VDJ analyst note, artifact map)

### Session 9 — PDF Report Generator
- **`core/report_gen.py`:** 8-section A4 PDF via fpdf2 (no TTF dependency)
  - Sections: Executive Summary, Exclusion Accounting, QC Metrics, Beta KDE, PCA (×2), Volcano (×2), DETECTED Events table, Significant DMR table
  - Git hash + run timestamp in header/footer on every page
  - `_safe()` helper ensures Latin-1 compliance (handles em-dash, arrows, ±, ⚠)
- **`run_pipeline()` gains `save_report=False`** param; `__main__` gains `--report` / `--no-figures` argparse flags
- **Result dict expanded:** adds `n_total`, `audit_csv`, `run_ts`, optional `report_path`
- 64/64 tests passing; 188 KB PDF smoke-tested

### Session 10 — Stage 3.5 VDJ Beta Masking
- **`mask_clonal_vdj_sites(df, clonal_samples)`** added to `core/repertoire_clonality.py`
  - Sets `beta_value=NaN` at VDJ-region rows for flagged clonal samples only
  - Non-VDJ rows and non-clonal samples untouched; returns `(df_masked, n_masked)`
- **Stage 3.5** inserted in `core/pipeline.py` between Stage 3 (clonality scan) and Stage 4 (normalization)
  - Guarded by `if clonal_samples:` — no-op when no clonal samples survive QC
  - Audit entry: `CLONALITY / cohort / INFO / VDJ-locus beta values masked to NaN`
- **NaN propagation chain:**
  - `robust_normalize`: `groupby.median()` skips NaN → NaN carries into `beta_normalized`
  - `dmr_hunter` pivot: `fillna(global_mean)` ≈ 0 after median-centring
  - `ml_guard` pivot: `fillna(per_cpg_mean)` ≈ 0 after median-centring
  - Net: inflated clonal VDJ signal (~0.97) suppressed to ~0 in feature matrices
- **3 new tests** in `test_phase3_modules.py`: sets_nan / preserves_non_vdj / preserves_other_samples
- 67/67 tests passing (commits 30b08e6, d0e1eaa)

### Session 11 — Output Directory Consolidation
- **`output/` directory** introduced as the single destination for all generated files
  - `output/logs/` — run logs (was `logs/`)
  - `output/figures/` — PNGs (was `figures/`)
  - `output/audit_log_*.csv` — per-run audit logs (was `data/audit_log_*.csv`)
  - `output/flagged_samples.csv` — cumulative flagged-sample registry (was `data/`)
  - `output/clean_methylation.csv` — clean data export (was `data/`)
  - `output/report_*.pdf` — PDF reports (was `data/`)
- **`data/` is now input-only** — contains only `mock_methylation.csv`
- **`core/io_utils.py`**: `output_path(filename)` helper added (analogous to `data_path()`); `append_flagged_samples` default updated
- **13 files updated**: FIGURES_DIR, audit log paths, log paths, flagged_samples paths across all `core/` modules and `data/generate_mock_data.py`
- **`.gitignore`**: single `output/` entry replaces `figures/`, `*.png`, `*.pdf`, `*.log`
- **`tests/test_phase3_modules.py`**: `test_audit_log_creation` searches `output/` instead of `data/`
- Old `logs/` and `figures/` directories and stale generated files deleted
- 67/67 tests passing; pipeline verified (commit b8448e7)

### Session 12–13 — XCI Guard, Clonal Artifact Placement Fix, Phase 3 Wrap-Up
- **`core/xci_guard.py`**: `compute_xci_signal()` + `detect_sex_mixups()` (Stage 1c); returns `(flagged_list, report_df)`
- **Clonal artifact moved to P003/S003** — P001/S001 and P002/S002 carry bisulfite failure and are removed at Stage 1a before Stage 3 runs; P003/S003 is the first case sample to survive QC and receive the VDJ inject
- **`detect_duplicates`** and **`flag_clonal_artifacts`** both return `(data_df, id_list)` for consistent pipeline wiring
- **`detect_sex_mixups`** returns `(flagged_list, report_df)` matching XCI guard convention
- **75/75 tests passing** across `test_mock_data.py`, `test_phase3_modules.py`, `test_io_utils_and_pipeline.py`
- Pipeline: 41 input samples → 34 clean (−3 bisulfite/depth, −1 contamination, −2 XCI mixup, −1 duplicate)
- 1 significant DMR (w00305, ΔBeta ≈ +0.112, p_adj ≈ 1.3e-04); AUC = 1.0000

### Phase 4 — Validation Notebook
- **`notebooks/ImmuneMethylTools_Validation.ipynb`**: 6-step end-to-end demo (Steps 0–5)
  - All logic imported from `core/`; no math reimplemented in cells
  - Before/after VDJ heatmaps (`df_pre_site` vs. `df_masked`) confirm Stage 3.5 masking
  - `requirements.txt` updated: added `jupyter>=7.0.0`, `ipykernel>=6.0.0`
- Notebook verified headless via `jupyter nbconvert --to notebook --execute`

### Post-Phase 5 — PCA Covariate Panel
- **`plot_pca_covariates(df)`** added to `core/visuals.py`: 4-panel figure (batch, disease label, sex, age) on normalized data
- Wired into `core/pipeline.py` after Stage 4 PCA figures → saved as `output/figures/pca_covariates.png`
- Notebook: markdown + code cells added after Step 3 normalization (Visual 4)
- Finding confirmed: PC1 = batch × disease confound; PC2 = sex (7M vs 7F Batch_01 Cases, zero overlap, mean PC2 F=+0.79 M=−0.69)
- CLAUDE.md module map and pipeline stage table updated

### Post-Phase 5 — TODO / Future Improvements
- **Age and sex covariates**: both are metadata-only; no model adjustment anywhere
  - Risk (age): age-confounded DMRs in unmatched cohorts
  - Risk (sex): sex drives a clean PC2 axis in PCA; sex-dimorphic autosomal methylation produces false positives if sex is imbalanced across case/control groups
  - Planned: Cramér's V checks for age × disease and sex × disease imbalance in `normalizer`; linear model `beta ~ disease + age + sex + batch` in `dmr_hunter`; `covariate_cols` parameter on `find_dmrs()`
- **Real EPIC / WGBS adaptation**: five blockers documented in README TODO and CLAUDE.md Future Improvements
  - Input format (minfi/Bismark → tidy CSV), CpG scale (chunked/sparse for WGBS), VDJ coordinates (GRCh38 loci), bisulfite QC (control probes for EPIC), sex inference (chrX coverage check)

### Phase 5 — Final Audit
- **`README.md`**: Full clinical-grade rewrite (~280 lines)
  - Artifact map corrected: Clonal VDJ row updated to S003/P003 (was stale S001/P001)
  - New section: SOP — Masking vs. Dropping (rationale for surgical NaN masking)
  - New section: Defense in Depth Strategy (6-gate layered guard architecture with rationale)
  - New section: Pipeline Stage Order table (mirrors CLAUDE.md)
  - New section: Notebook walkthrough (Step 0–5 descriptions)
  - CLI flags documented: `--report`, `--no-figures`
- **`core/report_gen.py`**: Executive Summary now shows "Sex metadata mixups removed: 2"; `n_total` fallback includes `n_sex_flagged`
- **Terminology cleanup**: 14 British "artefact" spellings corrected to "artifact" across `core/` and `data/`
- 75/75 tests passing; all outputs verified

### Session 14 — Package Structure, Configurable Thresholds, Chunked Processing

**Commits: e4d40a7, 37baff1**

#### JSON-Configurable Analysis Thresholds
- **`config.json`** (new, project root): human-editable thresholds for all pipeline stages. Auto-discovered by pipeline; override with `--config PATH` flag. Null = use default.
- **`core/config_loader.py`** (new): `load_config(path=None)` merges user file with hard-coded defaults; missing keys fall back silently.
- All key threshold functions now accept optional kwargs (backward-compatible; defaults = module constants):
  - `audit_quality(bisulfite_thresh, depth_thresh)`
  - `detect_contamination(bc_sigma_thresh, contamination_mean_lo/hi)`
  - `detect_duplicates(corr_thresh)`
  - `flag_clonal_artifacts(beta_min, frag_min)`
  - `find_dmrs(p_adj_thresh, delta_beta_min)`
  - `run_safe_model(n_top_cpgs, l1_ratio, c_param)`
- `pipeline.py`: loads config at startup; logs config source in banner; passes all threshold values as kwargs; `--config` CLI flag added.

#### Chunked CpG Processing
- `find_dmrs(chunk_size=N)`: sliding-window scan in overlapping N-CpG chunks; right border = WINDOW_SIZE−1 CpGs for boundary continuity; BH correction applied globally post-chunks.
- `run_safe_model(chunk_size=N)`: variance computed per chunk without materialising full Sample×CpG matrix; selects global top-N CpGs then builds final small feature matrix.
- Verified: chunk_size=50 on 500-CpG mock → identical 496 windows, correct DMR, AUC=1.0000.
- Recommended settings: `chunk_size=50_000` for EPIC (~40 MB/chunk at 100 samples); `chunk_size=100_000` for WGBS.

#### Package Infrastructure
- **`pyproject.toml`**: package metadata, pinned dependencies, optional `[notebook]` extra, `[tool.pytest.ini_options]` testpaths.
- **`requirements.txt`**: scientific core pinned to exact installed versions (pandas==3.0.0, numpy==2.4.2, scipy==1.17.0, scikit-learn==1.8.0, matplotlib==3.10.8, seaborn==0.13.2, fpdf2==2.8.5); jupyter/ipykernel kept as `>=`.

#### README Documentation
- **Cell-Type Deconvolution section**: algorithm (custom mock, not Houseman/EpiDISH), FoxP3/PAX5 marker CpGs and GRCh38 loci, production adaptation notes (EpiDISH, MethylResolver, Houseman 2012), expected inputs.
- **Configurable Thresholds section**: full parameter table with defaults and tuning notes; `--config` usage examples.
- **Adapting to real data**: CpG scale note updated — chunked processing now implemented (not just planned).

#### Final State
- 75/75 tests passing
- Pipeline: 41 → 34 clean samples; AUC = 1.0000; 1 significant DMR (w00305)
- GitHub: https://github.com/ChristopherSNelson/ImmuneMethylTools (main branch, current)

### Session 15 — Pseudoreplication Fix, Sub-Threshold Negative Controls, gc_content

**Commits: 3093411, 8b89691, 2ff34a6, a5dfb94, c9de513, 622b1d1**

#### DMR Pseudoreplication Fix
- `find_dmrs()`: eliminated pseudoreplication — per-sample mean per window computed before Wilcoxon rank-sum (not per-row, which erroneously gave each CpG-per-sample an independent observation)
- `ml_guard.py`: NaN imputation + feature selection moved inside the CV loop to prevent data leakage

#### Sub-Threshold Negative Controls
- `inject_borderline_signal()`: +0.09 raw shift at cg00001500–cg00001507 (8 CpGs, Case only); after median-centring expected ΔBeta ~0.08, just below `DELTA_BETA_MIN = 0.10` — must not appear as a significant DMR
- `inject_subtle_signal()`: +0.08 raw shift at cg00002000–cg00002005 (6 CpGs, Case only); expected ΔBeta ~0.04, well below threshold
- Both ranges autosomal, non-VDJ, non-X-linked; purpose: confirm the pipeline does not over-call near the detection boundary

#### gc_content Column
- `data/generate_mock_data.py`: `gc_content` column added — per-CpG `Uniform(0.30, 0.70)`, constant across samples, rounded to 2 decimal places
- `core/dmr_hunter.py`: `mean_gc` annotation per DMR window
- `core/io_utils.py`: `gc_content` added to `REQUIRED_COLUMNS` (14 total); validation check 8 added

#### Miscellaneous
- `core/visuals.py`: PCA caption numerics replaced with generic description

**75/75 tests passing**

---

### Session 16 — GRCh38 Coordinates, 10K CpG Scale, Coordinate-Based VDJ Annotation

**Commit: 0cd98c4**

#### Scale-Up
- `N_PATIENTS=100` (50 Case, 50 Control), `N_CPGS=10_000`, `N_X_CPGS=600`
- All test count assertions updated; clean sample count 34 → 94

#### GRCh38 Coordinate Assignment
- `assign_cpg_coordinates()`: distributes 10K CpGs across chr1–22 + chrX using `CHROM_SIZES_GRCH38`; ~282 placed inside real VDJ loci proportional to locus size; reserved signal CpG ranges excluded from VDJ placement
- `_SIGNAL_CLUSTERS` dict: places signal CpGs (true DMR, borderline, subtle) as tight genomic clusters (200–400 bp spacing) so they form clusters under distance-based detection
- `is_x_chromosome` derived from `chrom == "chrX"` (was random assignment)
- Schema expanded to 16 columns: `chrom` and `pos` added; `REQUIRED_COLUMNS` = 16; validators for valid chromosomes and non-negative positions

#### Coordinate-Based VDJ Annotation
- `annotate_vdj_regions(df)` added to `repertoire_clonality.py`: coordinate-based intersection with GRCh38 VDJ loci (6 loci, 2 kb buffer); adds `vdj_locus` column; backward-compatible for legacy data without chrom/pos
- `pipeline.py`: calls `annotate_vdj_regions(df)` after data loading, before QC

#### True Biological DMR Recalibration
- Signal shift set to +0.25; baseline reset to `Normal(0.35, 0.04)` clipped [0.15, 0.55] for all samples to avoid beta ceiling clipping
- Signal CpGs placed as tight genomic cluster on chr6:30000000–30002000

**80/80 tests passing.**
Pipeline results at 10K scale: 101 input → 94 clean; 1 significant DMR cluster (chr6:30000000-30002000, ΔBeta ~+0.18); AUC ~1.0

---

### Session 17 — RRBS-Safe Clonal Detector, Marker-Based Deconvolution, Distance-Based DMR Caller

**Commit: 30d6fe0**

#### RRBS-Safe Clonal Detector
- `flag_clonal_artifacts()` rewritten: SD-based fragment outlier detection (>3 SD from sample mean) replaces hardcoded 180 bp threshold; requires ≥3 hits in a single VDJ locus before flagging
- `annotate_vdj_regions()` now also returns `vdj_locus` (str/None) column
- Constants: `CLONAL_FRAG_MIN=180` → `CLONAL_FRAG_SD_THRESH=3.0` + `CLONAL_MIN_LOCUS_HITS=3`
- `config.json`: `clonality.frag_min` → `clonality.frag_sd_thresh` + `clonality.min_locus_hits`

#### Marker-Based Deconvolution
- `estimate_cell_fractions()` rewritten: T/B/Treg fractions derived from actual FoxP3/PAX5 proxy CpG beta values (cg1–5 = FoxP3, cg6–10 = PAX5); no longer random-seed based
- XCI correction: +0.25 beta offset at X-linked markers for female samples
- Sex-specific FoxP3 lineage shift thresholds: Male=0.15, Female=0.30
- `detect_lineage_shift()` output includes `sex` column
- 3 new tests: `test_clonal_requires_min_locus_hits`, `test_lineage_shift_sex_aware_foxp3`, `test_cell_fractions_marker_based`

#### Distance-Based CpG Cluster DMR Caller
- `find_dmrs()` fully rewritten:
  - `_build_clusters()`: sorts by (chrom, pos), groups consecutive CpGs within `MAX_GAP_BP=1000`, filters `min_cpgs ≥ 3`
  - Per-sample median per cluster before Wilcoxon (eliminates pseudoreplication at cluster level)
  - Output columns: `chrom`, `start_pos`, `end_pos`; cluster IDs format `cl00001` (was `w00001`)
  - Chunked processing path removed (clustering naturally limits pivot size)
  - Annotations: `n_vdj_cpgs`, `clonal_risk`, `mean_gc`, `test_method` per cluster
- `data/generate_mock_data.py`: `inject_true_biological_signal()` rewritten — baseline reset, Case +0.25

**83/83 tests passing.**

---

### Session 18 — Gemini Review Follow-Up, Centralized Logging, Deconvolution Refactor

**Commits: 7b5d4ad, 8690947, 072e5ff, 2790ef7, 532f448**

#### Gemini Review Integration
- `tldr_readme.md`: fixed placeholder clone URL; stripped 67 inline bold instances from prose
- `README.md`: stripped all inline bold from prose per style rule
- `core/visuals.py`: exclusion accounting pie chart — replaced overlapping inline labels with legend; hid percentage labels for slices < 10%
- `data/generate_mock_data.py`: rounded `gc_content` to 2 decimal places

#### Centralized Logging Helpers
- `core/io_utils.py`: `ts()` and `audit_entry()` added as centralized helpers (9 facilities total)
- All 9 core `__main__` blocks refactored: import `ts` and `audit_entry` from `io_utils`; removed 9 duplicated `ts()` definitions and 8 duplicated `ae()` definitions (net −39 lines across 12 files)
- Standalone modules: `ae = lambda sid, st, d, m: audit_entry(MODULE_TAG, sid, st, d, m)`
- Pipeline: `ae = audit_entry` (direct 5-arg alias)

#### Deconvolution Refactor
- `_get_mean_beta(df_rows, default)` helper extracted (replaces 4 inline `if len(df_rows)` patterns)
- 7 magic numbers promoted to named module-level constants (`TREG_FRAC_SCALE`, `B_FRAC_SCALE`, etc.)
- `core/config_loader.py`: docstring example updated to match all current `_DEFAULTS`

**83/83 tests passing.**

---

### Session 19 — Age/Sex Covariates, OLS DMR Model, statsmodels BH Correction

**Commit: 005c3c3**

#### Confound Checks in normalizer.py
- `check_continuous_confound(df, continuous_col, group_col)`: one-way ANOVA; reports F-stat and p-value for continuous covariate vs disease_label
- Pipeline Stage 4: batch (Cramér's V) + sex (Cramér's V) + age (ANOVA) all run before normalization

#### OLS DMR Model
- `find_dmrs(covariate_cols)`: when covariates provided, fits OLS per cluster via `statsmodels`
  - Dependent variable: logit-transformed M-values (addresses beta heteroscedasticity)
  - Model: `M-value ~ disease + covariates`; disease coefficient → p-value and t-statistic
  - Categorical covariates auto-detected and encoded with `C()`
  - `delta_beta` reported on original beta scale for biological interpretability
  - Wilcoxon rank-sum fallback when `covariate_cols=None` (backward compatible)
- BH correction: `statsmodels.stats.multitest.multipletests` (replaced custom implementation)
- Output: `wilcoxon_stat` → `test_stat`; new `test_method` column (`"OLS"` or `"Wilcoxon"`)
- Default `covariate_cols: ["age", "sex"]` in `config.json` and `config_loader.py`

**86/86 tests passing** (3 new: OLS path, Wilcoxon backward compat, age ANOVA confound).

---

### Session 20 — Artifact 8 Lineage Composition Anomaly

**Commit: 4fb01be**

#### Artifact 8 — Lineage Composition
- `inject_artifact8_lineage(df)` added to `data/generate_mock_data.py`:
  - Treg-enriched: S045/S046 (Case) — FoxP3 proxy cg00000001–cg00000005 (chr1) driven to beta ~0.06; below M=0.15 and F=0.30 thresholds → flagged by `detect_lineage_shift()`
  - B-cell depleted: S065/S066 (Control) — PAX5 proxy cg00000006–cg00000010 (chr1) driven to beta ~0.65; above 0.50 PAX5 threshold → flagged by `detect_lineage_shift()`
  - All four samples survive all QC stages (autosomal chr1, non-VDJ, non-X-linked)
- `CLAUDE.md`: Artifact 8 added to Known Stumper Artifacts; FoxP3/PAX5 proxy placement rules documented
- `core/deconvolution.py`: docstring fix for `_get_mean_beta()` helper

**86/86 tests passing.**

---

### Session 21 — plot_before_after Expanded to All 8 Artifacts

**Working changes (no commit yet)**

#### plot_before_after — 3×4 → 5×4 GridSpec
- `data/generate_mock_data.py` — `plot_before_after()`:
  - Docstring updated: "Six-panel" → "Ten-panel figure"; full layout table for all 5 rows
  - `figsize`: `(20, 16)` → `(20, 28)`; `GridSpec`: `(3, 4)` → `(5, 4, hspace=0.60)`
  - Rows 0–2 (Artifacts 1–5): unchanged
  - **Panel F (Artifact 6):** row 3 cols 0-1 — per-sample mean depth histogram; S030 marked with orange line; 10x threshold in red
  - **Panel G (Artifact 7):** row 3 cols 2-3 — jittered strip of mean X-linked beta by sex label; left = true sex (df_before, clean bimodal); right = reported sex (df_after, S035/S036 as black-star outliers in wrong group); local `np.random.default_rng(99)` for reproducible jitter
  - **Panel H (Artifact 8):** row 4 cols 0:2 + 2:4 — FoxP3 × PAX5 proxy scatter per sample; before = uniform cloud; after = S045/S046 (red triangles, low FoxP3) and S065/S066 (purple squares, high PAX5) separated; threshold lines on both panels
- Module docstring: "all 7 artifacts" → "all 8 artifacts"

**86/86 tests passing.**

---

### Session 22 — Report Polish, venv Fix, Test Rename

#### Split artifact simulation figure
- `data/generate_mock_data.py` — `plot_before_after()` split into:
  - `_plot_artifacts_1_to_5()` → `qc_before_after_1.png` (GridSpec 3×4, figsize 20×17)
  - `_plot_artifacts_6_to_8()` → `qc_before_after_2.png` (GridSpec 2×4, figsize 20×12)
  - `plot_before_after()` thin wrapper calling both; prints both paths

#### Report Section 2 — glob-based, future-proof
- `core/report_gen.py`:
  - `_ARTIFACT_FIG_META` dict keyed by filename stem; `glob("qc_before_after_*.png")` drives embedding
  - `idx > 0` page break before each subsequent artifact figure
  - `figure()` helper uses `multi_cell()` for caption wrapping
  - Caption strings: short, `\n`-separated row descriptions

#### Executive Summary — standalone mode fix
- `_parse_audit_summary(audit_df)` reconstructs all 11 pipeline KV fields from audit log CSV
- `__main__` block calls it instead of passing `pipeline_result={}`

#### venv statsmodels fix
- `core/dmr_hunter.py`: `smf` and `smm` moved to top-level imports; lazy import inside `find_dmrs()` removed

#### Test file rename
- `tests/test_phase3_modules.py` → `tests/test_core_integration.py` (git mv)

**86/86 tests passing.**
