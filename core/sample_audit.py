"""
core/sample_audit.py — ImmuneMethylTools Integrity Check
=========================================================
Detects technical duplicates via beta-profile correlation; documents the
interface for SNV-based concordance (gold standard, not yet implemented).

Biological intent
-----------------
Sample swap and duplication are among the most damaging errors in
population-level methylation studies — they inflate effective sample size,
introduce spurious correlations that survive batch correction, and produce
false-positive DMRs.

The gold standard for identity verification is SNV fingerprinting
(identity-by-descent / IBD), but high-correlation beta profiles at
high-variance CpGs serve as a fast, sequencing-only proxy:

  1. Retain the top-N high-variance CpGs (most informative for identity).
     Low-variance CpGs carry no sample-discriminating information and inflate
     apparent similarity between unrelated samples.
  2. Compute pairwise Pearson r across the retained CpGs.
  3. Flag pairs with r > 0.99 — profiles this similar are indistinguishable
     from technical duplicates at typical WGBS noise levels.

This approach mirrors the strategy used in MethylAid and minfi's outlier
detection, where informative site selection precedes any correlation step.
"""

import os
from itertools import combinations

import pandas as pd
from scipy.stats import pearsonr

# ── Thresholds ────────────────────────────────────────────────────────────────
DUP_CORR_THRESH = 0.99   # Pearson r ≥ this → flagged as technical duplicate
N_HIGH_VAR_CPGS = 100    # top-variance CpGs used for fingerprinting
MIN_CPGS_FOR_AUDIT = 10     # minimum CpGs required to attempt correlation


# =============================================================================
# Public API
# =============================================================================


def detect_duplicates(
    df: pd.DataFrame,
    corr_thresh: float = DUP_CORR_THRESH,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Identify sample pairs with suspiciously high beta-value correlation.

    Uses top-100 high-variance CpGs (maximum discriminating power) following
    the MethylAid / minfi convention.

    Parameters
    ----------
    df          : long-format methylation DataFrame
    corr_thresh : Pearson r at or above this value → technical duplicate (default 0.99)

    Returns
    -------
    pairs_df   : pd.DataFrame — all sample pairs with columns
                 [sample_a, sample_b, pearson_r, duplicate_flag]
                 sorted by pearson_r descending.
                 Empty DataFrame if fewer than MIN_CPGS_FOR_AUDIT sites available.
    ids_to_drop : list[str] — sample_b IDs from flagged pairs (the samples to
                  remove; sample_a is retained as the canonical representative)
    """
    # Pivot to CpG × Sample
    pivot = df.pivot_table(index="cpg_id", columns="sample_id", values="beta_value")

    # Select high-variance CpGs
    cpg_var = pivot.var(axis=1).sort_values(ascending=False)
    top_cpgs = cpg_var.head(N_HIGH_VAR_CPGS).index

    if len(top_cpgs) < MIN_CPGS_FOR_AUDIT:
        return (
            pd.DataFrame(columns=["sample_a", "sample_b", "pearson_r", "duplicate_flag"]),
            [],
        )

    X = pivot.loc[top_cpgs].fillna(pivot.mean())
    samples = X.columns.tolist()

    records = []
    for s_a, s_b in combinations(samples, 2):
        r, _ = pearsonr(X[s_a].values, X[s_b].values)
        records.append({"sample_a": s_a, "sample_b": s_b, "pearson_r": round(r, 6)})

    result = pd.DataFrame(records)
    result["duplicate_flag"] = result["pearson_r"] >= corr_thresh
    result = result.sort_values("pearson_r", ascending=False).reset_index(drop=True)

    ids_to_drop = result.loc[result["duplicate_flag"], "sample_b"].tolist()
    return result, ids_to_drop


def snv_concordance_placeholder(df: pd.DataFrame) -> str:
    """
    Placeholder for SNV-based identity check (gold standard).

    In a production pipeline, germline SNVs called from WGBS reads would be
    compared between samples using identity-by-descent (IBD) metrics.  Any
    pair with IBD > 0.9 and different patient_id labels represents a swap or
    contamination event.

    Required inputs (not available in this mock dataset):
      - Per-sample BAM files aligned to GRCh38
      - SNV calling: samtools mpileup → bcftools call
      - IBD estimation: plink --genome or king --ibdseg

    Interface stub:
      snv_concordance(bam_paths: list[str], ref_fasta: str) -> pd.DataFrame
      Returns DataFrame with columns [sample_a, sample_b, ibd_fraction, swap_flag]
    """
    return (
        "SNV concordance not yet implemented.  "
        "Required: BAM files → variant calls → pairwise IBD matrix.  "
        "See docstring for interface specification."
    )


# =============================================================================
# __main__ — run sample audit on mock data
# =============================================================================

if __name__ == "__main__":
    import sys
    from datetime import datetime

    # Ensure sibling core/ modules are importable regardless of working directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from io_utils import (  # noqa: E402
        Tee, append_flagged_samples, data_path, load_methylation, project_root, write_audit_log,
    )

    MODULE = "sample_audit"
    MODULE_TAG = "SAMPLE_AUDIT"
    _now = datetime.now()
    run_ts = _now.strftime("%Y-%m-%dT%H:%M:%S")
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base = project_root()
    _log = os.path.join(_base, "output", "logs", f"{MODULE}_{ts_tag}.log")
    _csv = os.path.join(_base, "output", "flagged_samples.csv")
    _audit_csv = os.path.join(_base, "output", f"audit_log_{MODULE_TAG}_{ts_tag}.csv")

    os.makedirs(os.path.join(_base, "output", "logs"), exist_ok=True)

    audit_entries = []

    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def ae(sample_id, status, description, metric):
        return {
            "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            "module": MODULE_TAG,
            "sample_id": sample_id,
            "status": status,
            "description": description,
            "metric": metric,
        }

    with Tee(_log):
        df = load_methylation(data_path("mock_methylation.csv"))
        n_samples = df["sample_id"].nunique()
        audit_entries.append(ae("cohort", "INFO", "Samples evaluated", f"n={n_samples}"))

        result, ids_to_drop = detect_duplicates(df)
        n_pairs = len(result)
        audit_entries.append(ae(
            "cohort", "INFO", "Pairwise correlation evaluated",
            f"n_pairs={n_pairs}",
        ))
        flagged = result[result["duplicate_flag"]]

        if len(flagged):
            for _, row in flagged.iterrows():
                print(
                    f"[{ts()}] [SAMPLE_AUDIT] DETECTED | Duplicate pair | "
                    f"{row.sample_a} ↔ {row.sample_b}  Pearson r={row.pearson_r:.6f}"
                )
                for sid, other in [
                    (row.sample_a, row.sample_b),
                    (row.sample_b, row.sample_a),
                ]:
                    audit_entries.append(ae(
                        sid, "DETECTED",
                        f"Technical duplicate — paired with {other}",
                        f"r={row.pearson_r:.4f}",
                    ))
        else:
            print(
                f"[{ts()}] [SAMPLE_AUDIT]           | No duplicates detected "
                f"(threshold r ≥ {DUP_CORR_THRESH})"
            )
            audit_entries.append(ae(
                "cohort", "INFO", "Duplicate check — none detected",
                f"threshold={DUP_CORR_THRESH}",
            ))

        print(f"\n[{ts()}] [SAMPLE_AUDIT] IDs to drop (sample_b of flagged pairs): {ids_to_drop}")
        print(f"\n[{ts()}] [SAMPLE_AUDIT] Top-5 correlated pairs:")
        print(result.head().to_string(index=False))

        print(f"\n[{ts()}] [SAMPLE_AUDIT] SNV check: {snv_concordance_placeholder(df)}")
        audit_entries.append(ae(
            "cohort", "INFO", "SNV concordance check — stub only",
            "requires_BAM_files",
        ))

        # ── Persist flagged_samples.csv ────────────────────────────────────────────
        flagged_rows = []
        for _, row in flagged.iterrows():
            for sid, other in [
                (row.sample_a, row.sample_b),
                (row.sample_b, row.sample_a),
            ]:
                flagged_rows.append({
                    "run_timestamp": run_ts,
                    "module": MODULE,
                    "sample_id": sid,
                    "flag_type": "duplicate",
                    "detail": f"r={row.pearson_r:.4f} paired with {other}",
                })

        append_flagged_samples(flagged_rows, _csv)
        n_unique = len({r["sample_id"] for r in flagged_rows})
        print(
            f"[{ts()}] [SAMPLE_AUDIT] {n_unique} unique flagged sample(s) written "
            f"→ output/flagged_samples.csv"
        )

        # ── Audit log ──────────────────────────────────────────────────────────────
        write_audit_log(audit_entries, _audit_csv)
        print(f"[{ts()}] [SAMPLE_AUDIT] Audit log → {_audit_csv}")
