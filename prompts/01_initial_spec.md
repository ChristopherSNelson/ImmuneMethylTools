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

## Status

- [x] Phase 0 complete
- [x] Phase 1 complete
- [x] Phase 2 complete
- [ ] Phase 3+ pending
