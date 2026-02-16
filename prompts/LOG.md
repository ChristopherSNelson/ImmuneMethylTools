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
