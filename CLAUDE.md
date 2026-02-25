# CLAUDE.md — ImmuneMethylTools Project Context

## Project Purpose
Rigorous IC-level analysis of B-cell/T-cell DNA methylation data for autoimmune disease research. Demonstrates detection of common wet-lab and bioinformatics pitfalls.

## Roles
- **Principal Bioinformatics Architect** (user): Reviews each phase, approves design decisions.
- **Senior ML Engineer** (Claude): Implements strictly per spec; stops after each phase for review.

## Architecture Decisions

| Decision | Rationale |
|----------|-----------|
| Beta-value representation | Standard [0,1] bounded methylation measure; logit-transform to M-values for DMR OLS modelling (applied to raw `beta_value`, not `beta_normalized`, since median-centring can push normalized values outside [0,1]).  ML Guard uses `beta_normalized` with `logit_transform=False` + StandardScaler to avoid learning batch identity. |
| Fragment length as clonal proxy | Long fragments in VDJ loci indicate clonal expansion, not true methylation |
| Non-CpG meth rate as bisulfite QC | >2% indicates incomplete bisulfite conversion; discard sample |
| Batch correction AFTER QC | QC filters must run before any batch normalization to avoid correcting artifacts |
| VDJ regions annotated, not excluded | `clonal_risk` + `n_vdj_cpgs` columns added per DMR window; analyst decides via `dmrs[~dmrs["clonal_risk"]]` |
| VDJ-locus beta masking (Stage 3.5) | `mask_clonal_vdj_sites()` sets `beta_value=NaN` at VDJ rows for clonal samples; downstream fillna imputes to ~0 (cohort mean after median-centring) |
| Site-level depth filter (Stage 2.5) | `filter_site_quality()` removes rows with depth < 5 from `df_clean` after sample dedup |
| XCI sex-signal filter at Stage 1c | Sex-metadata mismatches (X-linked beta vs. reported sex) removed before dedup; `detect_sex_mixups()` returns `(flagged_list, report_df)` |
| Detection functions return (data, ids) | `detect_duplicates` and `flag_clonal_artifacts` both return `(data_df, id_list)` for consistent pipeline wiring |
| fpdf2 for PDF reports; no TTF | Helvetica core font; `_safe()` helper replaces non-Latin-1 chars (em-dash, arrows, etc.) |

## Key Variable Glossary

| Variable | Type | Description |
|----------|------|-------------|
| `sample_id` | str | Unique per-sample label (e.g. `S001`) |
| `patient_id` | str | Patient; multiple samples possible per patient |
| `batch_id` | str | Sequencing batch (Batch_01, Batch_02, ...) |
| `age` | int | Patient age at collection |
| `disease_label` | str | `Case` or `Control` |
| `cpg_id` | str | CpG site identifier (e.g. `cg00001245`) |
| `beta_value` | float [0,1] | Methylation level at CpG site |
| `depth` | int | Read depth (coverage) at site |
| `fragment_length` | int (bp) | Insert size; >180 bp in VDJ = clonal flag |
| `is_vdj_region` | bool | True if site falls in VDJ recombination locus |
| `non_cpg_meth_rate` | float [0,1] | Non-CpG methylation rate; >0.02 = bisulfite failure |
| `sex` | str | `M` or `F` — reported sex from sample metadata |
| `is_x_chromosome` | bool | True if CpG falls on X chromosome |
| `gc_content` | float [0,1] | GC dinucleotide fraction around the CpG site; constant per CpG across samples |
| `chrom` | str | GRCh38 chromosome (e.g. `chr1`, `chrX`) |
| `pos` | int | GRCh38 genomic position |

## Known Stumper Artifacts (Simulated)

1. **Confounded Batch** — Batch_01 enriched for Cases (+0.1 mean beta shift)
2. **Clonal Artifact** — VDJ region, beta > 0.8, fragment length outlier (> 3 SD from sample mean)
3. **Bisulfite Failure** — 2 samples (S001, S002), non_cpg_meth_rate ≈ 0.05
4. **Sample Duplication** — S010 ↔ S_DUP, Pearson r = 0.9999
5. **Contamination** — S020, bimodality coefficient 0.78 (muddy beta near 0.5)
6. **Low Coverage** — S030, mean depth ≈ 5x (Poisson λ=5)
7. **Sex Metadata Mixup** — S035 (true F, reported M) + S036 (true M, reported F); X-linked beta contradicts reported sex
8. **Lineage Composition Anomaly** — S045/S046 (FoxP3 proxy beta ≈ 0.06, Treg-enriched); S065/S066 (PAX5 proxy beta ≈ 0.65, B-cell depleted)

## True Biological DMR Signal (Positive Control)

`inject_true_biological_signal()` in `generate_mock_data.py` adds a **batch-independent, disease-associated methylation shift** at CpGs `cg00003000`-`cg00003010` (11 sites, placed as a tight genomic cluster on chr6:30000000-30002000):
- All samples: baseline reset to Normal(0.35, 0.04) clipped to [0.15, 0.55]
- All **Case** samples: +0.25 added to baseline (post-clip range ~0.40-0.80)
- All **Control** samples: unchanged (baseline only)
- After full pipeline: observed ΔBeta ~+0.18 (Case vs. Control, per-sample median per cluster), p_adj < 0.05, 1 significant DMR cluster
- AUC ~1.0 with the larger cohort (94 clean samples)
- These CpGs are autosomal, non-VDJ, non-X-linked — unaffected by any artifact injection
- Purpose: gives `dmr_hunter` and `ml_guard` a genuine positive control; confirms the pipeline can detect a real signal amid all injected artifacts

### Sub-threshold Negative Controls

Two weaker signals are injected as negative controls to confirm the DMR caller does not over-call:

- `inject_borderline_signal()`: +0.09 raw shift at `cg00001500`-`cg00001507` (8 CpGs, Case only). After median-centring, expected ΔBeta ~0.08 — just below the `DELTA_BETA_MIN = 0.10` threshold. Should not appear as a significant DMR.
- `inject_subtle_signal()`: +0.08 raw shift at `cg00002000`-`cg00002005` (6 CpGs, Case only). After median-centring, expected ΔBeta ~0.04 — well below threshold. Should not appear as a significant DMR.

Both ranges are autosomal, non-VDJ, non-X-linked, and outside all other reserved CpG ranges.

## Module Map

| Path | Purpose |
|------|---------|
| `data/generate_mock_data.py` | Simulate all 8 artifacts + 1 true biological DMR signal + 2 sub-threshold negative controls into mock_methylation.csv |
| `core/infrastructure/io_utils.py` | `project_root()`, `data_path()` (inputs), `output_path()` (all outputs), `ts()`, `audit_entry()` (centralized logging helpers), `load_methylation()` (schema validator), `Tee`, `append_flagged_samples()`, `write_audit_log()` |
| `core/infrastructure/visuals.py` | Public helpers: `compute_sample_qc_stats()` (per-sample QC aggregation), `compute_pca_coords()` (pivot → median-impute → scale → PCA; used by notebook and plot functions). Plot functions: QC metrics, beta KDE, PCA, PCA covariate panel (batch/label/sex/age), exclusion accounting (pie + waterfall), volcano plot, coefficient rank chart |
| `core/infrastructure/report_gen.py` | 8-section A4 PDF report via fpdf2 — figures + audit log + DMR table; git hash in footer |
| `core/qc/qc_guard.py` | Bisulfite/depth sample QC, contamination detection, site-level depth filter |
| `core/qc/xci_guard.py` | X-Chromosome Inactivation sex-signal mismatch detector (`compute_xci_signal`, `detect_sex_mixups`) |
| `core/qc/sample_audit.py` | Technical duplicate detection (pairwise Pearson r) |
| `core/qc/repertoire_clonality.py` | RRBS-safe VDJ clonal artifact flagging (SD-based fragment outliers, min locus hits) + VDJ-locus beta masking (`mask_clonal_vdj_sites`) + coordinate-based VDJ annotation with `vdj_locus` column (`annotate_vdj_regions`) |
| `core/analytics/normalizer.py` | Confound checks: batch/sex × disease (Cramér's V) + age × disease (ANOVA), median-centring |
| `core/analytics/deconvolution.py` | Marker-based T/B/Treg cell-fraction estimation (FoxP3/PAX5 proxy CpGs) + sex-aware FoxP3 lineage shift detection |
| `core/analytics/dmr_hunter.py` | Distance-based CpG cluster DMR caller (per-sample median per cluster, OLS with `covariate_cols` or Wilcoxon rank-sum, BH FDR via statsmodels); annotates `n_vdj_cpgs` + `clonal_risk` + `mean_gc` + `test_method` per cluster |
| `core/analytics/ml_guard.py` | ElasticNet + GroupKFold CV classifier; NaN imputation + top-variance feature selection inside Pipeline per fold (leak-free CV) |
| `core/orchestration/pipeline.py` | End-to-end runner; passes clean_samples through all stages; `--report` / `--no-figures` / `--config` CLI flags |
| `core/orchestration/config_loader.py` | `load_config(path=None)` — loads `config.json` from project root; merges with hard-coded defaults |
| `config.json` | Human-editable analysis thresholds (QC, duplicates, clonality, DMR, ML); `null` = use default |
| `pyproject.toml` | Single source of truth for dependencies; optional `[notebook]` and `[dev]` extras, pytest testpaths |
| `requirements.txt` | Autogenerated lockfile via `pip-compile pyproject.toml --output-file requirements.txt --no-annotate` — do NOT edit by hand |
| `tests/` | 87 unit tests across 3 test files |
| `notebooks/` | End-to-end demo notebook (Phase 4) |

## Pipeline Stage Order (MUST NOT change without architect approval)

| Stage | Module | Action |
|-------|--------|--------|
| 1a | qc_guard | Bisulfite/depth sample QC → removes failed samples |
| 1b | qc_guard | Contamination detection → removes contaminated samples |
| 1c | xci_guard | XCI sex-signal check → removes sex-metadata mismatches |
| 2  | sample_audit | Duplicate removal → drops sample_b from each flagged pair |
|    | visuals | Exclusion accounting figure saved after Stage 2 |
| 2.5 | qc_guard | Site-level depth filter → removes rows with depth < 5 from df_clean |
| 3  | repertoire_clonality | Clonal VDJ scan → `clonal_samples` list |
| 3.5 | repertoire_clonality | `mask_clonal_vdj_sites` → `beta_value=NaN` at VDJ rows for clonal samples in df_clean |
| 4  | normalizer | Confound checks (batch/sex/age × disease) + median-centring → df_norm |
|    | visuals | PCA covariate panel (batch/label/sex/age) saved → `pca_covariates.png` |
| 5  | deconvolution | Cell fractions + per-sample T/B/Treg table (Case-first) |
| 6  | dmr_hunter | Distance-based CpG cluster DMR caller on df_norm (per-sample median per cluster, OLS with age/sex covariates or Wilcoxon rank-sum, BH FDR via statsmodels); volcano plots saved |
| 7  | ml_guard | ElasticNet GroupKFold CV |
|    | report_gen | Optional PDF report (--report flag) |

## Report Generation

- Output: `output/report_{YYYYMMDD_HHMMSS}.pdf` (~188 KB with all 7 figures embedded)
- Trigger: `python core/orchestration/pipeline.py --report`
- Standalone: `python core/infrastructure/report_gen.py` (auto-finds most-recent `audit_log_pipeline_*.csv`)
- fpdf2 core font (Helvetica); `_safe()` strips non-Latin-1 chars — no TTF dependency

## Style Rules
- **American English spelling** throughout all code, comments, docstrings, and documentation.
  - normalize / normalization (not normalise / normalisation)
  - color (not colour)
  - Use `sed -i '' -e 's/normalise/normalize/g' -e 's/colour/color/g'` to catch regressions.

## Authorship — MANDATORY for every git commit
Every commit made by Claude MUST include both co-authors in the commit message trailer:

```
Co-Authored-By: Christopher S. Nelson <christopher.s.nelson.01@gmail.com>
Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
```

This applies to all commits regardless of who initiates them.

## Execution Rules
- Stop after each phase; await architect review.
- Generate Before/After plots for every data manipulation.
- Commit incrementally when running autonomously.
- Never skip bisulfite intuition checks.

## Mandatory: core/ Module Main Blocks
Every Python script in `core/` MUST include a `if __name__ == "__main__":` block at the bottom that:
1. Loads `data/mock_methylation.csv`
2. Runs the module's detection/correction logic
3. Prints timestamped results to stdout for each issue detected
4. Format: `[YYYY-MM-DD HH:MM:SS] [MODULE] DETECTED | <issue description> | <key metric>`

## Centralized Audit Logging
In addition to stdout printing, every detection event MUST be logged to a centralized, timestamped CSV file in the `output/` directory (e.g., `output/audit_log_20260215_143005.csv`). All generated outputs (logs, figures, audit CSVs, reports, clean data) go to `output/`. The `data/` directory is input-only (`mock_methylation.csv`).

Each entry in the log must contain:
- `timestamp`: ISO-8601 format of the detection event.
- `module`: The name of the artifact detector module (e.g., `QC_GUARD`, `SAMPLE_AUDIT`).
- `sample_id`: The affected sample ID (or `cohort` for general findings).
- `status`: `DETECTED` (for artifacts) or `INFO` (for general stats).
- `description`: A concise summary of the issue or check performed.
- `metric`: The key numerical value (e.g., Pearson r, Cramér's V, non-CpG rate).

Note: For standalone runs, the filename should be generated at the start of the execution to ensure all detections from that run are grouped together.

## Documentation Structure
The project now employs a two-tiered README structure to balance conciseness and comprehensive detail:
-   **`README.md`**: Provides a concise, one-page overview of the project's purpose, key features, quick setup, and workflow. It serves as the primary entry point for new users or quick reference.
-   **`TECHNICAL_GUIDE.md`**: A comprehensive, detailed guide covering all aspects of the project, including in-depth explanations of features, design choices, artifact map, pipeline stages, configurable thresholds, and future work. The concise `README.md` links to this document for full details.
This structure ensures that essential information is immediately accessible, while extensive documentation is readily available for deeper dives.

## Future Improvements

### Age and sex as methylation covariates — COMPLETED (Session 19)

Implemented in full:
- `check_continuous_confound()` added to `normalizer.py` — one-way ANOVA for age × disease_label imbalance detection.  Sex × disease reuses the existing `check_confounding()` (Cramér's V).
- Pipeline Stage 4 now runs batch, sex, and age confound checks before normalization.
- `find_dmrs()` accepts `covariate_cols: list[str] | None` parameter.  When provided, an OLS linear model (`M-value ~ disease + covariates`) is fitted per cluster via `statsmodels`, with logit-transformed M-values as the dependent variable (beta values are heteroscedastic; M-values satisfy OLS normality assumptions).  The p-value and t-statistic come from the disease coefficient; `delta_beta` is reported on the original beta scale for biological interpretability.
- Categorical covariates (e.g. `sex`) are auto-detected and encoded with `C()`.
- When no covariates are passed, the Wilcoxon rank-sum test is used (backward compatible).
- Default `covariate_cols: ["age", "sex"]` in `config.json` and `config_loader.py`.
- BH correction now uses `statsmodels.stats.multitest.multipletests` (replaced custom implementation).
- Output columns: `wilcoxon_stat` renamed to `test_stat`; new `test_method` column (`"OLS"` or `"Wilcoxon"`).
- 86/86 tests passing (3 new tests: OLS path, Wilcoxon backward compat, age confound).

### Adapting to real EPIC / WGBS data
Five practical blockers to address before running on real cohorts:

| Item | Current (mock) | Real-data requirement |
|------|---------------|----------------------|
| Input format | `mock_methylation.csv` | Tidy long-format CSV from `minfi` (EPIC) or `Bismark` (WGBS); `load_methylation()` schema validator will catch mismatches |
| CpG scale | ~10,000 CpGs | EPIC ~850K; WGBS >25M — set `chunk_size` in `config.json` (50_000 for EPIC, 100_000 for WGBS); chunked path implemented in `find_dmrs` and `run_safe_model` |
| VDJ coordinates | ~~Mock boolean flag~~ DONE | GRCh38 coordinates (`chrom`, `pos`) assigned to all CpGs; `annotate_vdj_regions()` in `repertoire_clonality.py` intersects with 6 real VDJ loci (IGH, IGK, IGL, TRA_TRD, TRB, TRG) with 2 kb buffer |
| Bisulfite QC | non-CpG methylation rate >2% | EPIC uses control-probe metrics from `minfi`; replace `non_cpg_meth_rate` with array control-probe conversion efficiency |
| Sex inference | X-linked CpG beta values | Confirm chrX coverage for WGBS; EPIC includes X-linked probes by default |
