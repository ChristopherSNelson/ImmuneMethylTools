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
