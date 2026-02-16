# CLAUDE.md — ImmuneMethylTools Project Context

## Project Purpose
Rigorous IC-level analysis of B-cell/T-cell DNA methylation data for autoimmune disease research. Demonstrates detection of common wet-lab and bioinformatics pitfalls.

## Roles
- **Principal Bioinformatics Architect** (user): Reviews each phase, approves design decisions.
- **Senior ML Engineer** (Claude): Implements strictly per spec; stops after each phase for review.

## Architecture Decisions

| Decision | Rationale |
|----------|-----------|
| Beta-value representation | Standard [0,1] bounded methylation measure; logit-transform before linear modelling |
| Fragment length as clonal proxy | Long fragments in VDJ loci indicate clonal expansion, not true methylation |
| Non-CpG meth rate as bisulfite QC | >2% indicates incomplete bisulfite conversion; discard sample |
| Batch correction AFTER QC | QC filters must run before any batch normalization to avoid correcting artifacts |

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

## Known Stumper Artifacts (Simulated)

1. **Confounded Batch** — Batch_01 enriched for Cases (+0.1 mean beta shift)
2. **Clonal Artifact** — VDJ region, beta > 0.8, fragment > 180 bp
3. **Bisulfite Failure** — 2 samples, non_cpg_meth_rate > 0.02
4. **Sample Duplication** — 2 samples with Pearson r > 0.99
5. **Contamination** — 1 sample with beta distribution peaked near 0.5

## Module Map

| Path | Purpose |
|------|---------|
| `data/generate_mock_data.py` | Simulate all artifacts into mock_methylation.csv |
| `core/` | QC, batch correction, differential methylation modules |
| `notebooks/` | End-to-end demo notebook |
| `tests/` | Unit tests for each core module |

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
