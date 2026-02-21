# ImmuneMethylTools — Instruction Log

**Environment:** Local macOS machine (Darwin 25.2.0)
**Executor Model:** claude-sonnet-4-5-20250929
**Architect:** Principal Bioinformatics Architect

---

## Log Entries

### 2026-02-15 — Session 1

**Instruction received:** Initial full project specification (see `prompts/01_initial_spec.md`).

**Actions taken:**

#### Phase 0
- [x] Created `prompts/` directory
- [x] Saved spec to `prompts/01_initial_spec.md`
- [x] Created `prompts/LOG.md` (this file)

#### Phase 1
- [x] Created `CLAUDE.md` with project context, architecture decisions, variable glossary
- [x] Created `.gitignore`
- [x] Created `LICENSE` (MIT)
- [x] Created `requirements.txt`
- [x] Created `setup.sh`
- [x] Created directory skeleton: `data/`, `core/`, `notebooks/`, `tests/`

#### Phase 2
- [x] Wrote `data/generate_mock_data.py`
  - Artifact 1: Confounded Batch_01 (+0.1 mean shift, 80% Cases)
  - Artifact 2: Clonal VDJ artifact (beta > 0.8, fragment > 180 bp)
  - Artifact 3: Bisulfite failure (2 samples, non_cpg_meth_rate > 0.02)
  - Artifact 4: Sample duplication (correlation > 0.99)
  - Artifact 5: Contaminated sample (muddy beta, peaks near 0.5)

---

### 2026-02-15 — Session 1 (continued)

**Instruction received:** Two follow-up directives:
1. Write automated test suite for Phase 2 artifacts.
2. MANDATORY rule: Every `core/` Python script must have a timestamped `__main__` block that prints detected issues to stdout.

**Actions taken:**
- [x] Updated `CLAUDE.md` with mandatory main-block rule
- [x] Wrote `tests/test_mock_data.py` — 5 artifact assertions + pytest + standalone `__main__`

---

_Append new entries below as instructions are received._

---

### 2026-02-15 — Session 2

**Instruction received:** Phase 3 spec — implement eight artifact detector modules in `core/`.

**Actions taken:**

#### Phase 3 — Core modules (commit b7fb85e)
- [x] `core/visuals.py` — EDA suite: `plot_qc_metrics`, `plot_beta_distribution`, `plot_pca`
- [x] `core/qc_guard.py` — Bisulfite/depth filter + cohort-relative BC contamination detection
- [x] `core/sample_audit.py` — Pairwise correlation duplicate detection (top-100 HV CpGs)
- [x] `core/normalizer.py` — Cramér's V confound check + median-centering normalization
- [x] `core/repertoire_clonality.py` — Clonal VDJ artifact flagging with GRCh38 locus reference
- [x] `core/deconvolution.py` — Mock T/B fractions + FoxP3/PAX5 lineage shift detector
- [x] `core/dmr_hunter.py` — Sliding-window Wilcoxon + BH correction; safety assertion on clean_samples
- [x] `core/ml_guard.py` — ElasticNet LogReg with GroupKFold cross-validation

**Key detection results on mock_methylation.csv:**
- qc_guard: S001/S002 (bisulfite failure), S020 (contamination, BC=0.78 vs cohort median 0.92)
- sample_audit: S010 ↔ S_DUP, r=0.9999
- normalizer: Cramér's V=0.62 (batch_id × disease_label, p=3.6e-04)
- repertoire_clonality: 15 rows in S001/P001, mean β=0.97, mean frag=232 bp
- dmr_hunter: 0 significant DMRs (expected — batch confound must be corrected first)
- ml_guard: AUC≈0.43 (random mock data; GroupKFold working correctly)

**sklearn ≥1.8 compatibility fixes applied:**
- Removed `penalty="elasticnet"` (deprecated); use `l1_ratio=` with `solver="saga"` only
- Replaced `make_scorer(needs_proba=True)` with string scorer `"roc_auc"`

---

### 2026-02-15 — Session 2 (continued)

**Instruction received:** Add Phase 3 tests in same style as `tests/test_mock_data.py`.

**Actions taken (commit 6307753):**
- [x] Wrote `tests/test_phase3_modules.py` — 36 tests across all 8 modules
- [x] Edge cases tested: safety assertion (dmr_hunter), sorted output (sample_audit), exact median-centering, dual-criteria clonality, fraction sum-to-one
- [x] 46/46 total tests passing (Phase 2 + Phase 3) in 5.8 s

**Instruction received:** List Christopher S. Nelson as co-author on all commits.

**Actions taken (commit 6307753):**
- [x] Added Authorship section to `CLAUDE.md` (mandatory co-author trailers)
- [x] Added architect name/email + co-author rule to `prompts/01_initial_spec.md`
- [x] Added co-author reminder to `memory/MEMORY.md`

**Instruction received:** Enforce American English spelling — replace normalise/colour variants.

**Actions taken (commit 569cb35):**
- [x] `sed` replaced all British variants in `core/`, `tests/`, `prompts/`, `CLAUDE.md`
  - normalise/normalisation/normalised → normalize/normalization/normalized
  - colour/coloured → color/colored
- [x] Added Style Rules section to `CLAUDE.md` with sed regression check one-liner
- [x] 46/46 tests still passing after spelling changes

---

### 2026-02-15 — Session 3

**Instruction received:** Series of figure and documentation polish tasks.

**Actions taken (commit 8ada017):**
- [x] `core/visuals.py:139` — clarified KDE subtitle: "muddy/contaminated samples show peaks shifting toward 0.5"
- [x] Regenerated `figures/normalization_before_after.png` (stale after Session 2 spelling fix; code was already correct)

**Actions taken (commit 0d8539c):**
- [x] `data/generate_mock_data.py` — removed letter prefix (A/B/C/D/E) from `qc_before_after.png` panel titles; now "Before: ..." / "After: ..."

**Actions taken (commit b6db52b):**
- [x] `prompts/01_initial_spec.md` — replaced character-metaphor module names with scientific descriptions:
  - "The Gatekeeper" → Sample QC Filter
  - "The Integrity Check" → Sample Integrity Auditor
  - "The Lineage Guard" → Clonal Expansion & VDJ Artifact Detector
  - "The Strict Analyst" → Differential Methylation Region Caller
  - Phase 3 section header and `core/` directory note updated to "artifact detector modules"
- [x] `prompts/LOG.md` — updated Session 2 entry to match ("detective" → "artifact detector")

**Actions taken (commit be626b3):**
- [x] `data/generate_mock_data.py` — further `qc_before_after.png` improvements:
  - "Before"/"After" → "Before Artifact Injection"/"After Artifact Injection" across all 5 panels
  - Non-CpG Meth Rate histogram: added `range=(0, 0.07)` so both Before/After bins start at 0
  - Panel D subtitle now conditional: "Pearson r, 12-sample subset" (before) vs "Pearson r, subset incl. duplicate pair" (after)

**Actions taken (commit b0f74b6):**
- [x] `data/generate_mock_data.py` — annotate per-(batch, disease) sample counts inside Panel A boxplot axes
  - Labels placed at 2% below top edge (`va='top'`) using post-draw ylim
  - ylim locked after annotation to prevent autoscale drift

---

### 2026-02-15 — Session 4

**Instruction received:** Persist flagged sample IDs to `data/flagged_samples.csv` and capture each module's stdout run to a timestamped log file under `logs/`.

**Actions taken (commit 2068ec6):**
- [x] Created `core/io_utils.py` — shared I/O utilities (no pandas dependency):
  - `append_flagged_samples(rows, csv_path)` — appends dicts to CSV via `csv.DictWriter`; creates file with header on first call, appends thereafter
  - `Tee` class — context manager that mirrors `sys.stdout` to a log file simultaneously
- [x] Updated `__main__` block in `core/qc_guard.py`:
  - Wraps run in `Tee(logs/qc_guard_{ts}.log)`
  - Writes bisulfite-failure rows (`flag_type=bisulfite_failure`, `detail=non_cpg_meth_rate=X`) and contamination rows (`flag_type=contamination`, `detail=BC=X mean_beta=Y`) to CSV
- [x] Updated `__main__` block in `core/sample_audit.py`:
  - Wraps run in `Tee(logs/sample_audit_{ts}.log)`
  - Writes both sides of each duplicate pair (`flag_type=duplicate`, `detail=r=X paired with Y`) to CSV
- [x] Updated `__main__` block in `core/repertoire_clonality.py`:
  - Wraps run in `Tee(logs/repertoire_clonality_{ts}.log)`
  - Writes clonal samples (`flag_type=clonal_vdj`, `detail=mean_beta=X mean_frag=Ybp`) to CSV
- [x] Fixed `core/visuals.py` Panel 2 depth histogram: `range=(0, max*1.1)` + `set_xlim(left=0)` so x-axis is always anchored at zero

**Verification results:**
- `data/flagged_samples.csv` — 6 rows from 3 modules (S001/S002 bisulfite, S020 contamination, S010/S_DUP duplicate, S001 clonal_vdj)
- `logs/` — 3 timestamped `.log` files, one per module run
- 46/46 tests passing

---

### 2026-02-15 — Session 5

**Instructions received:**
1. Update `CLAUDE.md` with new "Centralized Audit Logging" directive.
2. Update `prompts/01_initial_spec.md` with matching "Timestamped Audit Log" core directive.
3. Implement audit logging across all 8 `core/` modules.

**Schema** (`data/audit_log_{YYYYMMDD_HHMMSS}.csv`):
```
timestamp, module, sample_id, status, description, metric
```
- `status`: `DETECTED` (artifact) or `INFO` (general stat)
- `sample_id`: affected sample or `cohort` for cohort-wide findings
- Per-run file — filename generated at execution start; history never overwritten

**Actions taken (commit 0c9e956):**
- [x] `CLAUDE.md` — added "Centralized Audit Logging" section with full schema and filename convention
- [x] `prompts/01_initial_spec.md` — added "Timestamped Audit Log" bullet to Core Directives
- [x] `core/io_utils.py` — added `write_audit_log(entries, csv_path)`: writes per-run CSV via `csv.DictWriter`, fresh file each call
- [x] All 8 `__main__` blocks updated:
  - Local `ae(sample_id, status, description, metric)` helper builds each entry with live ISO-8601 timestamp
  - Entries collected throughout the run; `write_audit_log()` called at end; path printed to stdout
  - `MODULE_TAG` uppercase constant used as the `module` field (e.g. `QC_GUARD`, `SAMPLE_AUDIT`, `CLONALITY`)

**Sample audit log output — `qc_guard` standalone run:**
```
timestamp,module,sample_id,status,description,metric
2026-02-15T20:54:38,QC_GUARD,cohort,INFO,Samples evaluated,n=41
2026-02-15T20:54:38,QC_GUARD,S001,DETECTED,Bisulfite conversion failure,non_cpg_meth_rate=0.050
2026-02-15T20:54:38,QC_GUARD,S002,DETECTED,Bisulfite conversion failure,non_cpg_meth_rate=0.050
2026-02-15T20:54:38,QC_GUARD,cohort,INFO,Clean samples retained,n=39
2026-02-15T20:54:38,QC_GUARD,S020,DETECTED,Contamination — low bimodality coefficient,BC=0.7797
```

**Existing outputs preserved:**
- `data/flagged_samples.csv` — cumulative cross-run sample registry unchanged
- `logs/*.log` — Tee stdout capture for the 3 flagging modules unchanged
- 46/46 tests passing

---

### 2026-02-15 — Session 6

**Instruction received:** Add Artifact 6 (low-coverage failure) to `data/generate_mock_data.py` — force S030 depth to `Poisson(λ=5)`.

**Actions taken (commit 35da468):**
- [x] `data/generate_mock_data.py`:
  - Module docstring updated to list Artifact 6 — Low Coverage (S030, mean ~5x)
  - New `inject_artifact6_low_depth(df)` function: applies `RNG.poisson(lam=5)` to all S030 rows; prints mean depth vs 10x threshold
  - Called in `main()` after Artifact 5; summary stats line added for S030
- [x] `tests/test_phase3_modules.py`:
  - `test_audit_quality_clean_count`: expected clean count updated 39 → 38 (S001 + S002 bisulfite + S030 depth = 3 failures)
  - `test_audit_quality_excludes_bisulfite_failures`: added `assert "S030" not in clean`
- 46/46 tests passing

---

### 2026-02-15 — Session 7

**Instructions received:** Sync tests and error messages with current codebase state (four action items).

**Actions taken (commit eb79397):**
- [x] `core/dmr_hunter.py` — `"SAFETY VIOLATION"` → `"DATA VALIDATION ERROR"` in assert; `# Safety assertion` → `# Input validation`
- [x] `tests/test_phase3_modules.py`:
  - `test_dmr_hunter_safety_assertion_triggers`: `match=` updated to `"DATA VALIDATION ERROR"`
  - `test_audit_log_creation` (new, 47 tests total): checks `data/audit_log_*.csv` exist; searches all logs for S001/S002 as DETECTED entries (not just most-recent, which may be VISUALS cohort-only)
  - `log()` / `run()` standalone helpers: `module` param → `artifact_detector`
  - Module docstring: `"detective modules"` → `"artifact detector modules"`
- S030 assertion kept — Artifact 6 confirmed active (S030 mean depth = 4.96x in current mock data)
- 47/47 tests passing

---

### 2026-02-16 — Session 8

**Instructions received (across two messages — context compaction mid-session):**

1. Restore `penalty="elasticnet"` in `core/ml_guard.py` LogisticRegression for older sklearn compat (committed 36cb66b).
2. Create `core/pipeline.py` — end-to-end runner importing all 8 detectors, passing `clean_samples` automatically (committed c9dacb2).
3. Add `df_clean.to_csv("data/clean_methylation.csv")` as final pipeline step; add `clean_csv` to return dict (committed f25d1c8).
4. Add centralized "Safe Loading" to `core/io_utils.py`: `project_root()`, `data_path()`, `load_methylation()` with 6 validation checks; refactor all 9 `__main__` blocks (committed e9e9d51).
5. Label module audit log files with module tag (e.g. `audit_log_QC_GUARD_{ts}.csv`); pipeline gets `audit_log_pipeline_{ts}.csv` (committed bc231bc).
6. Add high-priority tests: `tests/test_io_utils_and_pipeline.py` (15 tests, module-scoped pipeline fixture) + `test_ml_guard_group_kfold_patient_exclusivity` (committed 41f99f9).
7. Implement `filter_site_quality()` in `core/qc_guard.py` — site-level depth QC returning `(df_clean, stats_dict)`.
8. Wire `min_site_depth=5` param into `dmr_hunter.find_dmrs()` and `ml_guard.run_safe_model()`.
9. Add Stage 2.5 (site-level depth QC) to `pipeline.py`.
10. Change `detect_duplicates` to return `tuple[pd.DataFrame, list[str]]` — `ids_to_drop` is `sample_b` from flagged pairs (mirrors `flag_clonal_artifacts` pattern).

**Actions taken (commit e90b864):**
- [x] `core/qc_guard.py`:
  - Added `SITE_DEPTH_THRESH=5`, `SITE_LOW_DEPTH_SAMPLE_WARN=20.0` constants
  - Added `filter_site_quality(df, min_depth=5) → (df_clean, stats_dict)`: removes rows with `depth < min_depth`; stats include n_total, n_low, pct_low, per_sample_pct series, min_depth
  - `__main__` block: logs cohort INFO + per-sample DETECTED when pct_s > 20%
- [x] `core/dmr_hunter.py`: added `min_site_depth=5` to `find_dmrs()`; depth filter applied after input assertion, before VDJ mask
- [x] `core/ml_guard.py`: added `min_site_depth=5` to `run_safe_model()`; depth filter applied before `pivot_table`
- [x] `core/pipeline.py`:
  - Imports: `SITE_DEPTH_THRESH`, `SITE_LOW_DEPTH_SAMPLE_WARN`, `filter_site_quality` from qc_guard
  - Stage 2.5 inserted between Stage 2 (dedup) and Stage 3 (clonality); applies `filter_site_quality(df_clean)` and writes audit entries
  - Stage 2 now uses `ids_to_drop` list directly from `detect_duplicates` to update `clean_samples`
- [x] `core/sample_audit.py`: `detect_duplicates` now returns `(pairs_df, ids_to_drop)` tuple; `__main__` prints `ids_to_drop`; docstring updated
- [x] `tests/test_phase3_modules.py`: 4 existing `detect_duplicates` callers updated to unpack tuple; added `test_detect_duplicates_returns_ids_to_drop`
- **64/64 tests passing**

**Key design decisions:**
- `detect_duplicates` now matches `flag_clonal_artifacts` pattern: both return `(data, id_list)`
- `filter_site_quality` is site-level (not sample-level): individual rows removed, samples kept
- `min_site_depth` on `find_dmrs` / `run_safe_model` provides defense-in-depth even if pipeline Stage 2.5 is bypassed

---

### 2026-02-16 — Session 8 (continued)

**Instruction received:** Replace blanket VDJ exclusion in `dmr_hunter` with a `clonal_risk` flag; let the Analyst decide. Write to stdout + audit log. Add Usage note in README.

**Actions taken (commit b86725a):**
- [x] `core/dmr_hunter.py`:
  - Removed `df[~df["is_vdj_region"]]` blanket exclusion
  - Built `vdj_cpgs` set; annotates each window with `n_vdj_cpgs` (int) and `clonal_risk` (bool)
  - Significant windows with `clonal_risk=True` print `⚠ HIGH CLONALITY` tag; DETECTED audit entries for sig clonal windows
  - Empty-DataFrame fallback and docstring updated
- [x] `core/pipeline.py` Stage 6: `⚠ HIGH CLONALITY` tag on sig DMR lines; DETECTED audit entries
- [x] `tests/test_phase3_modules.py`: `test_dmr_hunter_excludes_vdj` → `test_dmr_hunter_flags_vdj_cpgs_as_clonal_risk`; output columns test updated
- [x] `README.md`: created — Setup, Usage, VDJ clonal_risk analyst note + filter code examples, artifact map
- 64/64 tests passing

**Instruction received:** Log per-sample T/B/Treg balance table sorted Case/Control to pipeline log.

**Actions taken (commit 5892e98):**
- [x] `core/pipeline.py` Stage 5: joined disease_label; added `b_t_ratio`; sorted Case-first; per-group mean summary lines; full 36-row per-sample table printed to stdout (→ log via Tee)
- 64/64 tests passing

**Instruction received:** Exclusion accounting pie+waterfall figure and volcano plot (with and without clonal-risk coloring).

**Actions taken (commit 3c63007):**
- [x] `core/visuals.py`: added `plot_exclusion_accounting()` (two-panel: pie + floating-bar waterfall) and `plot_volcano()` (standard + clonal-risk coloring variant)
- [x] `core/pipeline.py`: calls both functions at appropriate stages; volcano called twice (standard + clonal risk)
- 64/64 tests passing

### 2026-02-16 — Session 9

**Instruction received:** PDF report generator — bundle audit log + figures into a structured PDF. `core/report_gen.py`. Add `--report` flag to pipeline `__main__`. Include git hash + run timestamp in footer on every page.

**Actions taken (commit 98d0306):**
- [x] `requirements.txt`: added `fpdf2>=2.6.0`
- [x] `core/report_gen.py` (new file, 8-section A4 PDF):
  - `_safe()` — Latin-1 sanitizer replaces em-dash, arrows etc. (no TTF required)
  - `_git_info()` — returns (short_hash, commit_date) via subprocess
  - `_Report(FPDF)` — subclass with `header()`, `footer()`, `section()`, `body()`, `kv()`, `figure()` helpers
  - `_detected_table()` — bordered table of DETECTED events, alternating row fill
  - `_sig_dmr_table()` — bordered table of significant DMRs, clonal_risk rows amber
  - `generate_report(pipeline_result, audit_csv, output_path, run_ts)` — 8-section report
  - `__main__` — standalone mode: auto-finds most-recent audit_log_pipeline_*.csv
- [x] `core/pipeline.py`:
  - `run_pipeline()` gains `save_report=False` param
  - Result dict built inside `with Tee` block; adds `n_total`, `audit_csv`, `run_ts` keys
  - `generate_report()` called conditionally; `report_path` added to result dict
  - `__main__` replaced with argparse: `--report` and `--no-figures` flags
- 64/64 tests passing
- Smoke test: `python core/pipeline.py --report` → 188 KB PDF generated cleanly

---

### 2026-02-17 — Session 10

**Instruction received:** Implement Stage 3.5 — surgically mask VDJ-locus `beta_value` to NaN in clonally-expanded samples before normalization and ML, so clonal expansion signal cannot drive the classifier.

**Actions taken (commit 30b08e6):**
- [x] `core/repertoire_clonality.py`:
  - Added `mask_clonal_vdj_sites(df, clonal_samples) → (df_masked, n_masked)`
  - Sets `beta_value=NaN` only at rows where `sample_id in clonal_samples AND is_vdj_region=True`
  - Non-VDJ rows and non-clonal samples untouched; returns a copy + count of masked sites
- [x] `core/pipeline.py`:
  - Import: added `mask_clonal_vdj_sites` alongside `flag_clonal_artifacts`
  - Stage 3 audit entry description updated: "VDJ CpGs masked (beta=NaN) before normalization and ML"
  - Stage 3.5 block inserted between Stage 3 and Stage 4 (guarded by `if clonal_samples:`)
  - Summary line: "Clonal VDJ rows (masked)" (was "excluded")
  - Module docstring updated to reflect Stage 3.5
- [x] `tests/test_phase3_modules.py`:
  - Import updated to include `mask_clonal_vdj_sites`
  - 3 new tests added: `test_mask_clonal_vdj_sites_sets_nan`, `_preserves_non_vdj`, `_preserves_other_samples`
- [x] `CLAUDE.md`: Architecture Decision row added for Stage 3.5 masking; module map + pipeline stage table updated; test count 64 → 67
- [x] `README.md`: test count updated; Stage 3.5 section added explaining NaN propagation and downstream fillna behavior
- **67/67 tests passing**

**Design note:** Stage 3.5 is silent in a normal pipeline run because S001 (the only clonal sample) is expelled in Stage 1a for bisulfite failure — correct by design. The masking fires when a clonally-expanded sample survives QC. Unit tests exercise the function directly on raw df, bypassing QC stages.

**NaN propagation through downstream modules:**
- `robust_normalize`: `groupby.median()` skips NaN → NaN carries into `beta_normalized`
- `dmr_hunter pivot`: `fillna(global_mean_beta_normalized)` ≈ 0 after median-centring
- `ml_guard pivot`: `fillna(per_cpg_mean)` ≈ 0 after median-centring
- Inflated clonal VDJ signal (~0.97) suppressed to ~0 in feature matrices

---

### 2026-02-17 — Session 11

**Instruction received:** Consolidate all generated outputs into a single `output/` directory; `data/` becomes input-only. Clean up stale logs, figures, and CSVs.

**Actions taken (commit b8448e7):**
- [x] `core/io_utils.py`:
  - Added `output_path(filename)` helper — analogous to `data_path()`, resolves to `output/filename`
  - `append_flagged_samples` default updated from `"data/flagged_samples.csv"` to `None` (computes `output_path("flagged_samples.csv")` inside)
  - Module docstring updated to list 7 facilities; `write_audit_log` docstring references `output/`
  - `Tee` docstring example path updated
- [x] `FIGURES_DIR` updated in 4 files: `core/visuals.py`, `core/normalizer.py`, `core/report_gen.py`, `data/generate_mock_data.py` → `output/figures/`
- [x] All 8 `core/` module `__main__` blocks: audit log paths and run log paths → `output/`; flagged_samples paths → `output/`; `makedirs` targets updated
- [x] `core/pipeline.py`: `_log`, `_flag_csv`, `_audit_csv`, `clean_csv`, `_report_path` all → `output/`; `makedirs` updated; summary print strings updated
- [x] `core/report_gen.py`: glob pattern for auto-discovery → `output/`; `report_gen.py` standalone output → `output_path()`; renamed local var `output_path` → `report_file` to avoid shadowing imported function
- [x] `tests/test_phase3_modules.py`: `test_audit_log_creation` searches `output/` instead of `data/`
- [x] `.gitignore`: single `output/` entry replaces `figures/`, `*.png`, `*.pdf`, `*.log`; added `logs/` and `figures/` as legacy patterns
- [x] Deleted old `logs/`, `figures/` directories and all stale generated files from `data/`
- [x] Created `output/logs/` and `output/figures/` skeleton
- 67/67 tests passing; pipeline verified — all outputs land in `output/`

**New directory contract:**
```
data/       ← INPUT ONLY (mock_methylation.csv)
output/
  logs/                        pipeline and module run logs
  figures/                     all PNGs
  audit_log_{MODULE}_{ts}.csv  standalone module runs
  audit_log_pipeline_{ts}.csv  full pipeline runs
  flagged_samples.csv          cumulative cross-run registry
  clean_methylation.csv        QC-passed export
  report_{ts}.pdf              PDF report
```

---

### 2026-02-17 — Session 12

**Instruction received:** Implement Stage 1c — X-Chromosome Inactivation (XCI) Guard. Detect samples whose X-linked methylation signal contradicts reported sex metadata (sample swap / label error).

**Biological rationale:** XX females exhibit Lyonization: one X silenced → mean X-linked beta ~0.50 in bulk tissue. XY males have a single active X → unmethylated baseline ~0.25. A mismatch between `sex` metadata and observed X-linked beta signals a metadata error or case/control label mixup.

**Actions taken (commit fd5cbc6):**

- [x] `data/generate_mock_data.py`:
  - Schema extended to 13 columns: added `sex` ("M"/"F") and `is_x_chromosome` (bool, last 30 of 500 CpGs; indices 471–500)
  - Sex assigned by patient parity: odd patient number → Female, even → Male
  - `inject_xci_signal(df)` — new function; called **after** artifacts 1–6 so batch shifts, bisulfite inflation, and contamination smearing do not corrupt X-linked ground truth; Female clip: `Normal(0.50, 0.05).clip(0.35, 0.65)`; Male clip: `Normal(0.25, 0.04).clip(0.10, 0.33)`
  - `inject_artifact7_sex_mixup(df)` — new function; swaps `sex` metadata only (not beta values): S035 (true F → reported M), S036 (true M → reported F)
  - Verified: S035 X-mean β = 0.499; S036 X-mean β = 0.248; S010↔S_DUP r = 0.9992
  - Module docstring updated to list Artifact 7; summary stats section extended

- [x] `core/io_utils.py`:
  - `REQUIRED_COLUMNS` extended from 11 to 13 columns (added `sex`, `is_x_chromosome`)
  - `load_methylation()`: added check 7 — sex vocabulary validation (`{"M","F"}`)
  - `_minimal_valid_df()` in tests updated to include new columns

- [x] `core/xci_guard.py` (new module):
  - `MODULE = "XCI_GUARD"`; thresholds: `XCI_FEMALE_LO=0.40`, `XCI_FEMALE_HI=0.65`, `XCI_MALE_HI=0.35`, `N_X_CPGS_MIN=5`
  - `compute_xci_signal(df) → summary_df` — groups on `(sample_id, sex)` over X-linked rows; returns `mean_x_beta`, `n_x_cpgs`
  - `detect_sex_mixups(df) → (flagged_list, report_df)` — flags samples with insufficient `n_x_cpgs` as non-mismatch; adds `xci_mismatch` bool column to report
  - `__main__` block: per-sample signal table, DETECTED entries, audit log, flagged_samples log

- [x] `core/pipeline.py`:
  - Import: `from xci_guard import detect_sex_mixups`
  - Stage 1c inserted between Stage 1b (contamination) and Stage 2 (duplicate removal)
  - `n_sex_flagged` counter + DETECTED audit entries + flagged_rows entries
  - Exclusion accounting figure updated: XCI stage added as 3rd exclusion category
  - Summary print + return dict both include `n_sex_flagged`
  - Module docstring updated to list Stage 1c

- [x] `tests/test_phase3_modules.py`:
  - Import: `compute_xci_signal`, `detect_sex_mixups`, `XCI_FEMALE_LO/HI`, `XCI_MALE_HI`
  - 5 new tests: `test_compute_xci_signal_returns_per_sample`, `test_detect_sex_mixups_flags_s035_s036`, `test_detect_sex_mixups_does_not_flag_clean_samples`, `test_female_xci_signal_in_range`, `test_male_xci_signal_below_threshold`

- [x] `tests/test_io_utils_and_pipeline.py`:
  - `_minimal_valid_df()` updated to include `sex` and `is_x_chromosome`
  - `test_pipeline_returns_expected_keys`: added `n_sex_flagged`
  - `test_pipeline_clean_samples_count`: 36 → 34
  - `test_pipeline_exports_clean_csv`: 36 → 34 samples, "13 CLAUDE.md schema columns"
  - New test: `test_pipeline_n_sex_flagged` (asserts == 2)

- [x] `tests/test_mock_data.py`:
  - New test: `test_low_coverage_sample` (Artifact 6 — S030 mean depth < 10x)
  - New test: `test_sex_mixup_artifact` (Artifact 7 — S035 X-beta ≥ 0.40 but sex="M"; S036 X-beta < 0.35 but sex="F")

- [x] `CLAUDE.md`:
  - Architecture Decisions: added XCI Guard row
  - Key Variable Glossary: added `sex` and `is_x_chromosome` rows
  - Known Stumper Artifacts: added Artifact 7
  - Module Map: added `core/xci_guard.py`; updated `generate_mock_data.py` entry; test count 67 → 75
  - Pipeline Stage Order: Stage 1c inserted between 1b and 2

- [x] `.gitignore`: added `*.Rhistory` (commit bc3abec)

**Sample count ledger (updated):**
```
41 input → 38 (Stage 1a: 3 QC) → 37 (Stage 1b: contamination) → 35 (Stage 1c: 2 XCI) → 34 (Stage 2: 1 dup)
```

**Key design decision — XCI signal injection order:**
The XCI `inject_xci_signal()` step must run AFTER artifacts 1–6. Artifact 1 adds +0.10 to all Batch_01 Case betas (including male X-linked CpGs), pushing male X-beta to ~0.35 and triggering false-positive mismatches. Moving XCI injection last (before the metadata swap) gives every sample clean ground-truth X-linked betas, so detection thresholds are reliably met.

**75/75 tests passing.**

---

### 2026-02-17 — Session 12 (continued)

**Instruction received:** Spike in a "True" Biological Signal — inject a genuine Case/Control DMR (not a batch artifact) so `dmr_hunter` has a real positive-control window to find.

**Actions taken (commit e129902):**
- [x] `data/generate_mock_data.py`:
  - Added `inject_true_biological_signal(df)`: shifts `cg00000300`–`cg00000310` (11 CpGs) by +0.25 for all Case samples, regardless of batch
  - Chosen window is autosomal, far from VDJ/FoxP3/PAX5 proxy markers, and not X-linked — unaffected by `inject_xci_signal()`
  - Called in `main()` after Artifact 6, before `inject_xci_signal()`
  - Summary stat added: Case mean β vs Control mean β at DMR window
  - Module docstring updated to note the true biological signal alongside the seven stumper artifacts
- Verified: Case mean β = 0.897, Control mean β = 0.712, Δ = +0.185 (well above dmr_hunter's 0.10 threshold)

**Pipeline smoke-test results (post-signal-spike):**
- Stage 6 DMR Hunter: **1 significant DMR** found — `w00305`, ΔBeta = +0.1123, p_adj = 1.24e-05 ✓
- Stage 7 ML Guard: AUC = 1.0000 (the true DMR is a perfect discriminator on this mock cohort) ✓
- 75/75 tests passing ✓

**Instruction received:** Update LOG.md and wrap up; Phase 4 notebook deferred to next session.

---

### 2026-02-17 — Session 12.5

**Instruction received:** Phase 4 — create end-to-end validation notebook.

**Actions taken:**
- [x] `notebooks/ImmuneMethylTools_Validation.ipynb` — 6-step notebook importing all core modules, walking through QC → dedup → clonal masking → normalization → DMR calling → ML classification
- Phase 4 complete; 75/75 tests passing

---

### 2026-02-18 — Session 13

**Instructions received:** PDF report layout fixes, missing figure generation, caption edits, README updates, and session cleanup.

**Actions taken (commit 4da6ea2):**
- [x] `core/report_gen.py`:
  - Removed horizontal rule from `header()` that was overlapping the run/commit subtitle
  - Increased post-section-header spacing from 4 mm to 10 mm (~2 line heights) — prevents tables/figures from overlapping the blue section band
  - Added `pdf.add_page()` before Section 3 and Section 8 — every section now starts on its own page
- [x] `core/pipeline.py`:
  - Added `plot_qc_metrics`, `plot_beta_distribution`, `plot_pca` to imports from `visuals`
  - Added `plot_qc_metrics(df)` + `plot_beta_distribution(df)` calls on raw df immediately after data load (guarded by `save_figures`)
  - Added `plot_pca(df_norm, color_by="disease_label")` + `plot_pca(df_norm, color_by="batch_id")` calls after Stage 4 normalization
  - All four previously-missing figures now generated on every pipeline run
- [x] `README.md`: `artefacts` → `artifacts` (American English style rule)

**Actions taken (commit ceae60c):**
- [x] `core/report_gen.py` Section 3 caption: replaced inaccurate "bimodality coefficient / flagged samples marked" with accurate description of the three histogram panels and Stage 1a threshold lines
- [x] `core/report_gen.py` Section 4 caption: "mode near 0.5" → "distributions closer to 0.5"

**Actions taken (commit 5a3d7ce):**
- [x] `README.md`:
  - Clean sample count: 36 → 34
  - Test count: 67 → 75
  - Artifact Map: added Artifact 7 (sex mixup) and true biological DMR positive control row

**Session close:**
- [x] `MEMORY.md` condensed — removed all content duplicated in CLAUDE.md; retained only commit authorship reminder, phase summary, current state, env quirks, and terminology
- [x] `prompts/LOG.md` updated with this entry
- **75/75 tests passing**

---

### 2026-02-18 — Session 13 (continued)

**Instructions received:** Clinical README rewrite, PDF report sex-mixup fix, artefact→artifact spelling cleanup.

**Actions taken:**
- [x] `README.md` — full clinical rewrite: Multi-Stage Data Integrity Pipeline section, SOP for masking vs. dropping, deconvolution docs, configurable thresholds table, notebook walkthrough, TODO/Future Work
- [x] `core/report_gen.py` — fixed sex-mixup sample IDs in report text
- [x] Global artefact→artifact spelling pass
- Phase 5 complete

---

### 2026-02-18 — Session 14

**Instructions received:** JSON-configurable thresholds, pyproject.toml, pinned requirements, deconvolution README section, chunked CpG processing, `--config` CLI flag.

**Actions taken:**
- [x] `config.json` — human-editable thresholds (QC, duplicates, clonality, DMR, ML); `null` = use default
- [x] `core/config_loader.py` — `load_config(path=None)` merges JSON with hard-coded defaults
- [x] `pyproject.toml` — package metadata, pinned deps, optional `[notebook]` extra
- [x] `requirements.txt` — exact version pins
- [x] `core/pipeline.py` — `--config PATH` CLI flag; loads config at startup
- [x] `core/dmr_hunter.py` / `core/ml_guard.py` — chunked processing (`chunk_size` param)
- [x] `README.md` — deconvolution section, configurable thresholds docs

---

### 2026-02-18 — Session 15

**Instructions received:** Pseudoreplication fix, PCA caption fix, sub-threshold negative controls, gc_content covariate.

**Actions taken (commits 3093411, 8b89691, 2ff34a6, a5dfb94, c9de513, 622b1d1):**
- [x] `core/dmr_hunter.py` — eliminated pseudoreplication: per-sample mean per window before Wilcoxon rank-sum (not per-row)
- [x] `core/ml_guard.py` — moved feature selection + NaN imputation inside CV loop to prevent data leakage
- [x] `core/visuals.py` — replaced hardcoded PCA caption numerics with generic description
- [x] `data/generate_mock_data.py` — added `inject_borderline_signal()` (+0.09 at cg150-157) and `inject_subtle_signal()` (+0.08 at cg200-205) as DMR negative controls
- [x] `data/generate_mock_data.py` — added `gc_content` column: per-CpG Uniform(0.30, 0.70), constant across samples
- [x] `core/dmr_hunter.py` — added `mean_gc` annotation per DMR window
- [x] `core/io_utils.py` — added `gc_content` to REQUIRED_COLUMNS (14 total), validation check 8
- [x] `CLAUDE.md` — updated Key Variable Glossary, module map, sub-threshold negative controls docs
- [x] Tests updated for new columns; **75/75 tests passing**

---

### 2026-02-19 — Session 16

**Instructions received:** Wire VDJ detector to GRCh38 coordinates at 10K CpG scale. Scale from 40 patients / 500 CpGs to 100 patients / 10,000 CpGs. Add real chromosome coordinates. Coordinate-based VDJ annotation.

**Actions taken (commit 0cd98c4):**
- [x] `data/generate_mock_data.py` — major rewrite:
  - `N_PATIENTS=100` (50 Case, 50 Control), `N_CPGS=10_000`, `N_X_CPGS=600`
  - Added `VDJ_LOCI_GRCH38` (6 loci with 2 kb buffer), `CHROM_SIZES_GRCH38` (chr1-22)
  - New `assign_cpg_coordinates()`: distributes CpGs across chr1-22 + chrX; ~282 placed inside real VDJ loci proportional to locus size; reserved signal ranges excluded from VDJ placement
  - `build_manifest()` outputs `chrom` and `pos` columns; `is_x_chromosome` derived from `chrom == "chrX"`; batch splits parameterized (no more hardcoded 20/16)
  - `add_baseline_methylation()` derives `is_vdj_region` from coordinate map instead of random boolean
  - Signal ranges scaled 10x: true DMR cg3000-3010, borderline cg1500-1507, subtle cg2000-2005
  - True bio DMR shift calibrated to +0.15 (post-centering delta ~0.11, above 0.10 threshold; +0.10 was too small)
  - Artifact 2 restricted to chr14 (IGH locus) via `df["chrom"] == "chr14"` filter
  - Duplication moved after XCI injection to preserve correlation on sex-driven high-variance CpGs
- [x] `core/repertoire_clonality.py`:
  - VDJ_LOCI_GRCH38 updated to 6 entries (merged TRA_TRD, 2 kb buffer)
  - New `annotate_vdj_regions(df)`: coordinate-based VDJ annotation with backward compatibility for legacy data without chrom/pos
- [x] `core/io_utils.py`: REQUIRED_COLUMNS expanded from 14 to 16 (added `chrom`, `pos`); validators for valid chromosomes and non-negative positions
- [x] `core/pipeline.py`: imports and calls `annotate_vdj_regions(df)` after data loading
- [x] `tests/test_mock_data.py`: 5 new tests (`test_chrom_and_pos_columns_exist`, `test_total_cpg_count`, `test_total_sample_count`, `test_vdj_cpgs_in_real_coordinates`, `test_x_linked_cpgs_on_chrx`)
- [x] `tests/test_phase3_modules.py`: QC clean count 38 -> 98
- [x] `tests/test_io_utils_and_pipeline.py`: `_minimal_valid_df` updated to 16 columns; pipeline count assertions updated (34 -> 94 clean samples)
- [x] `CLAUDE.md`: schema updated (16 cols), signal ranges, VDJ coordinates marked done in Future Improvements

**Sample count ledger (updated):**
```
101 input -> 98 (Stage 1a: 3 QC) -> 97 (Stage 1b: contamination) -> 95 (Stage 1c: 2 XCI) -> 94 (Stage 2: 1 dup)
```

**Pipeline results at 10K scale:**
- 3-4 overlapping significant DMRs at true signal (w03003-w03006, p_adj < 1e-13, delta ~+0.11)
- Sub-threshold negative controls correctly non-significant
- Clonal artifact: 122 VDJ rows on chr14 (IGH locus) in P003/S003
- AUC ~0.65 (weaker with diluted signal at 10K scale — expected)

**Known item for next session:** PCA point labels illegible at 100 samples — remove sample_id annotations from `plot_pca` and `plot_pca_covariates`

**80/80 tests passing.**

---

### 2026-02-20 — Session 17

**Instructions received:** Two biological critiques + DMR caller rewrite.

#### Part 1: RRBS-safe clonal detector + marker-based deconvolution (commit 30d6fe0)

**Instruction received:** Make clonal detector RRBS-safe (SD-based fragment outlier detection, min locus hits) and replace fixed random-seed deconvolution with marker-based cell fractions (FoxP3/PAX5 proxies, sex-aware FoxP3 thresholds).

**Actions taken:**
- [x] `core/repertoire_clonality.py`:
  - `annotate_vdj_regions()` now returns both `is_vdj_region` (bool) and `vdj_locus` (str or None) columns
  - `flag_clonal_artifacts()` rewritten: SD-based fragment outlier detection (>3 SD from sample mean) replaces hardcoded 180 bp threshold; requires ≥3 hits in a single VDJ locus before flagging
  - `get_vdj_summary()` updated to accept `frag_sd_thresh` parameter
  - Constants: `CLONAL_FRAG_MIN=180` → `CLONAL_FRAG_SD_THRESH=3.0` + `CLONAL_MIN_LOCUS_HITS=3`
- [x] `core/deconvolution.py`:
  - `estimate_cell_fractions()` rewritten: computes from actual FoxP3/PAX5 proxy CpG beta values instead of random seed
  - XCI correction for females: +0.25 beta offset at X-linked markers
  - Sex-specific FoxP3 lineage shift thresholds: Male=0.15, Female=0.30
  - `detect_lineage_shift()` output includes `sex` column
- [x] `config.json`: `clonality.frag_min` → `clonality.frag_sd_thresh` + `clonality.min_locus_hits`
- [x] `core/config_loader.py`: defaults updated to match
- [x] `core/pipeline.py`: Stage 3 passes `frag_sd_thresh` and `min_locus_hits`; Stage 5 lineage shift output includes sex tag
- [x] `tests/test_phase3_modules.py`: 3 new tests (`test_clonal_requires_min_locus_hits`, `test_lineage_shift_sex_aware_foxp3`, `test_cell_fractions_marker_based`)
- **83/83 tests passing**

#### Part 2: Distance-based CpG cluster DMR caller (uncommitted)

**Instruction received:** Replace fixed-count sliding window with distance-based clustering: sort by genomic position, cluster within 1000 bp gap, filter ≥3 CpGs, per-sample median aggregation, Wilcoxon rank-sum, BH FDR correction.

**Actions taken:**
- [x] `core/dmr_hunter.py` — full rewrite:
  - Replaced `WINDOW_SIZE=5, STEP_SIZE=1` with `MAX_GAP_BP=1000`
  - New `_build_clusters()` helper: sorts CpGs by (chrom, pos), groups consecutive CpGs within max_gap bp, filters by min_cpgs (≥3)
  - `find_dmrs()` rewritten: builds clusters, pivots only clustered CpGs, computes per-sample MEDIAN per cluster, Wilcoxon rank-sum, BH correction
  - Output includes `chrom`, `start_pos`, `end_pos` columns; cluster IDs format `cl00001` (was `w00001`)
  - Removed chunked processing path (clustering naturally reduces pivot size)
  - Annotations preserved: `n_vdj_cpgs`, `clonal_risk`, `mean_gc` per cluster
- [x] `data/generate_mock_data.py`:
  - Added `_SIGNAL_CLUSTERS` to `assign_cpg_coordinates()`: places signal CpGs as tight genomic clusters (200-400 bp spacing) so they form clusters under distance-based detection
  - `inject_true_biological_signal()` rewritten to avoid beta ceiling clipping: baseline reset to Normal(0.35, 0.04) for all samples, then +0.25 for Case only (raw delta ~0.25, post-normalization delta ~0.18)
- [x] `core/pipeline.py`: Stage 6 messages say "clusters" instead of "windows"
- [x] `core/report_gen.py`: DMR table columns include `chrom`; "Window" → "Cluster" in header
- [x] `tests/test_phase3_modules.py`: `test_dmr_hunter_output_columns` updated with `chrom`, `start_pos`, `end_pos`
- [x] `CLAUDE.md`: True Biological DMR Signal section, Module Map, Stage 6 description, test count all updated

**Pipeline results with distance-based clustering:**
- 1 significant DMR cluster (true biological signal, chr6:30000000-30002000, ΔBeta ~+0.18, p_adj < 0.05)
- 2 sub-threshold negative controls correctly non-significant
- AUC ~1.0 (strong signal with reset baseline)
- **83/83 tests passing**

**Key technical challenges solved:**
1. Signal CpGs (cg3000-3010) were randomly scattered across chromosomes → wouldn't form clusters. Fixed by placing them as tight genomic clusters via `_SIGNAL_CLUSTERS`.
2. Beta ceiling clipping: baseline ~0.85 + 0.25 shift clipped to 1.0, killing the delta. Fixed by resetting baseline to 0.35 before injection.

---

### 2026-02-20 — Session 18

**Executor Model:** claude-opus-4-6

**Instructions received:** Gemini review follow-up, style fixes, refactoring across all core modules.

#### Part 1: Gemini commit review + docs style fixes

Reviewed commit `8d746b1` (by gemini-2.5-pro): two-tiered README restructure (`README.md` concise + `tldr_readme.md` comprehensive), `GEMINI.md` mandate file, `CLAUDE.md` docs-structure section.

**Actions taken (commits 7b5d4ad, 8690947, 072e5ff):**
- [x] `tldr_readme.md` — fixed placeholder clone URL (`your-repo` → `ChristopherSNelson/ImmuneMethylTools`)
- [x] `README.md` — stripped all inline bold (`**...**`) from prose per style rule
- [x] `tldr_readme.md` — stripped all inline bold from prose (67 replacements)

#### Part 2: Pie chart fix + gc_content rounding

**Actions taken (commit 2790ef7):**
- [x] `core/visuals.py` — exclusion accounting pie chart: replaced overlapping inline labels with a legend; hid percentage labels for slices < 10%
- [x] `data/generate_mock_data.py` — rounded `gc_content` values to 2 decimal places

#### Part 3: Centralize logging helpers + deconvolution refactor + config docstring

**Actions taken (commit 532f448):**
- [x] `core/io_utils.py` — added `ts()` and `audit_entry()` as centralized helpers (9 facilities now); updated module docstring
- [x] Refactored all 9 core module `__main__` blocks to import `ts` and `audit_entry` from `io_utils`; removed 9 duplicated `ts()` definitions and 8 duplicated `ae()` definitions (net -39 lines across 12 files)
  - Standalone modules use lambda alias: `ae = lambda sid, st, d, m: audit_entry(MODULE_TAG, sid, st, d, m)`
  - Pipeline uses direct alias: `ae = audit_entry`
  - xci_guard: replaced `_ts()` and inline dicts with centralized functions
- [x] `core/deconvolution.py`:
  - Extracted `_get_mean_beta(df_rows, default)` helper (replaces 4 inline `if len(df_rows)` patterns)
  - Promoted 7 magic numbers to module-level constants: `TREG_FRAC_SCALE`, `TREG_FRAC_FLOOR`, `B_FRAC_SCALE`, `B_FRAC_FLOOR`, `T_FRAC_FLOOR`, `OTHER_FRAC_FLOOR`, `OTHER_FRAC_COMPLEMENT`
- [x] `core/config_loader.py` — updated docstring example `config.json` to match current `_DEFAULTS` (added `frag_sd_thresh`, `min_locus_hits`, `site_low_depth_sample_warn`, `contamination_mean_lo/hi`, `chunk_size` for dmr/ml)
- [x] `core/visuals.py` — added sample_ID coloring note to KDE legend-omitted caption

**83/83 tests passing.**
