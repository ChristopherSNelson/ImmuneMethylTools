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
