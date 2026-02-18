# CLAUDE.md — ImmuneMethylTools Project Context

## Project Purpose
Rigorous IC-level analysis of B-cell/T-cell DNA methylation data for autoimmune disease research. Demonstrates detection of common wet-lab and bioinformatics pitfalls.

## Roles
- **Principal Bioinformatics Architect** (user): Reviews each phase, approves design decisions.
- **Senior ML Engineer** (Claude): Implements strictly per spec; stops after each phase for review.

## Architecture Decisions

| Decision | Rationale |
|----------|-----------|
| Beta-value representation | Standard [0,1] bounded methylation measure; logit-transform to M-values before linear modelling (applied to raw `beta_value`, not `beta_normalized`, since median-centring can push normalized values outside [0,1]) |
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

## Known Stumper Artifacts (Simulated)

1. **Confounded Batch** — Batch_01 enriched for Cases (+0.1 mean beta shift)
2. **Clonal Artifact** — VDJ region, beta > 0.8, fragment > 180 bp
3. **Bisulfite Failure** — 2 samples (S001, S002), non_cpg_meth_rate ≈ 0.05
4. **Sample Duplication** — S010 ↔ S_DUP, Pearson r = 0.9999
5. **Contamination** — S020, bimodality coefficient 0.78 (muddy beta near 0.5)
6. **Low Coverage** — S030, mean depth ≈ 5x (Poisson λ=5)
7. **Sex Metadata Mixup** — S035 (true F, reported M) + S036 (true M, reported F); X-linked beta contradicts reported sex

## True Biological DMR Signal (Positive Control)

`inject_true_biological_signal()` in `generate_mock_data.py` adds a **batch-independent, disease-associated methylation shift** at CpGs `cg00000300`–`cg00000310` (11 sites):
- All **Case** samples: +0.25 to beta_value (clipped to [0,1])
- All **Control** samples: unchanged
- After full pipeline: observed ΔBeta ≈ +0.185 (Case vs. Control), p_adj ≈ 1.2e-05, 1 significant DMR window (w00305)
- AUC rises to 1.0000 when this signal is present (ElasticNet correctly identifies the sites)
- These CpGs are autosomal, non-VDJ, non-X-linked — unaffected by any artifact injection
- Purpose: gives `dmr_hunter` and `ml_guard` a genuine positive control; confirms the pipeline can detect a real signal amid all injected artifacts

## Module Map

| Path | Purpose |
|------|---------|
| `data/generate_mock_data.py` | Simulate all 7 artifacts + 1 true biological DMR signal into mock_methylation.csv |
| `core/io_utils.py` | `project_root()`, `data_path()` (inputs), `output_path()` (all outputs), `load_methylation()` (schema validator), `Tee`, `append_flagged_samples()`, `write_audit_log()` |
| `core/visuals.py` | QC metrics, beta KDE, PCA, PCA covariate panel (batch/label/sex/age), exclusion accounting (pie + waterfall), volcano plot |
| `core/qc_guard.py` | Bisulfite/depth sample QC, contamination detection, site-level depth filter |
| `core/xci_guard.py` | X-Chromosome Inactivation sex-signal mismatch detector (`compute_xci_signal`, `detect_sex_mixups`) |
| `core/sample_audit.py` | Technical duplicate detection (pairwise Pearson r) |
| `core/normalizer.py` | Batch × disease confound check (Cramér's V), median-centring |
| `core/repertoire_clonality.py` | VDJ clonal artifact flagging (`flag_clonal_artifacts`) + VDJ-locus beta masking (`mask_clonal_vdj_sites`) |
| `core/deconvolution.py` | Mock T/B/Treg cell-fraction estimation + FoxP3/PAX5 lineage shift detection |
| `core/dmr_hunter.py` | Sliding-window Wilcoxon DMR caller; annotates `n_vdj_cpgs` + `clonal_risk` per window |
| `core/ml_guard.py` | ElasticNet + GroupKFold CV classifier (data-leakage validator) |
| `core/pipeline.py` | End-to-end runner; passes clean_samples through all stages; `--report` / `--no-figures` CLI flags |
| `core/report_gen.py` | 8-section A4 PDF report via fpdf2 — figures + audit log + DMR table; git hash in footer |
| `tests/` | 75 unit tests across 3 test files |
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
| 4  | normalizer | Confound check + median-centring → df_norm |
|    | visuals | PCA covariate panel (batch/label/sex/age) saved → `pca_covariates.png` |
| 5  | deconvolution | Cell fractions + per-sample T/B/Treg table (Case-first) |
| 6  | dmr_hunter | Sliding-window DMR caller on df_norm; volcano plots saved |
| 7  | ml_guard | ElasticNet GroupKFold CV |
|    | report_gen | Optional PDF report (--report flag) |

## Report Generation

- Output: `output/report_{YYYYMMDD_HHMMSS}.pdf` (~188 KB with all 7 figures embedded)
- Trigger: `python core/pipeline.py --report`
- Standalone: `python core/report_gen.py` (auto-finds most-recent `audit_log_pipeline_*.csv`)
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

## Future Improvements

### Age and sex as methylation covariates
Both are metadata-only — collected but not modeled anywhere in the pipeline.

- **Risk (age)**: without age adjustment, `dmr_hunter` may flag age-associated CpGs as disease signal if cases and controls are not age-matched.
- **Risk (sex)**: sex drives a clean PC2 separation even in synthetic data; in real EPIC data, sex-dimorphic autosomal methylation is strong enough to produce false positives if sex is imbalanced across case/control groups.
- **Planned fix**:
  - Add Cramér's V / ANOVA checks for age × disease_label and sex × disease_label imbalance in `normalizer`, analogous to the existing batch × disease check.
  - Replace the Wilcoxon test in `dmr_hunter` with a linear model (`beta ~ disease + age + sex + batch`) so both are proper covariates.
  - Expose `covariate_cols: list[str]` parameter to `find_dmrs()` so analysts can pass additional covariates (age, sex, smoking status) without code changes.

### Adapting to real EPIC / WGBS data
Five practical blockers to address before running on real cohorts:

| Item | Current (mock) | Real-data requirement |
|------|---------------|----------------------|
| Input format | `mock_methylation.csv` | Tidy long-format CSV from `minfi` (EPIC) or `Bismark` (WGBS); `load_methylation()` schema validator will catch mismatches |
| CpG scale | ~500 CpGs | EPIC ~850K; WGBS >25M — profile `dmr_hunter` pivot and `ml_guard` pivot; consider chunked / sparse processing |
| VDJ coordinates | Mock boolean flag | Replace with real IGH/IGK/IGL/TRA/TRB locus coordinates from GRCh38 in `repertoire_clonality.py` |
| Bisulfite QC | non-CpG methylation rate >2% | EPIC uses control-probe metrics from `minfi`; replace `non_cpg_meth_rate` with array control-probe conversion efficiency |
| Sex inference | X-linked CpG beta values | Confirm chrX coverage for WGBS; EPIC includes X-linked probes by default |
