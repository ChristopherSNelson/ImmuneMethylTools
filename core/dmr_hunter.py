"""
core/dmr_hunter.py — ImmuneMethylTools Strict Analyst
======================================================
Differentially Methylated Region (DMR) detection using a sliding-window
non-parametric test on clean, normalized data.

Biological intent
-----------------
A DMR is a contiguous stretch of CpGs where the methylation difference between
disease groups is:
  (a) Statistically significant after multiple-testing correction (p-adj < 0.05)
  (b) Biologically meaningful  (|ΔBeta| > 0.10 — below this, differences have
      negligible impact on transcription factor binding affinity)
  (c) Regionally consistent    (≥ 3 CpGs per window — single-site hits are
      likely noise; a region requires multiple concordant signals)

Safety guarantees
-----------------
  1. ASSERTION: the input DataFrame must contain ONLY samples from
     clean_samples_list — preventing artefact-contaminated data from biasing
     DMR calls.  Passing uncleaned data raises AssertionError immediately.
  2. VDJ-region CpGs are EXCLUDED so clonal expansion artefacts (high beta +
     long fragment) are not called as case-associated DMRs.

Statistical approach
--------------------
Sliding window over CpGs sorted by numeric index (proxy for chromosomal order):
  - Window size: 5 CpGs (step = 1)
  - Test:        Wilcoxon rank-sum (non-parametric; robust to non-Gaussian betas)
  - Correction:  Benjamini-Hochberg (BH) FDR
  - Filter:      p_adj < 0.05, |ΔBeta| > 0.10, ≥ 3 CpGs
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import ranksums

# ── Parameters ────────────────────────────────────────────────────────────────
WINDOW_SIZE    = 5      # sliding window width (CpGs)
STEP_SIZE      = 1      # step between windows
P_ADJ_THRESH   = 0.05   # BH-corrected p-value threshold
DELTA_BETA_MIN = 0.10   # minimum |ΔBeta| to qualify as a DMR
MIN_CPGS       = 3      # minimum CpGs per window


# =============================================================================
# Internal helpers
# =============================================================================


def _bh_correction(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction (implemented without statsmodels).

    Returns adjusted p-values clipped to [0, 1], preserving the rank order
    of significance.
    """
    n = len(pvalues)
    if n == 0:
        return np.array([])

    order    = np.argsort(pvalues)
    ranked_p = pvalues[order]
    adjusted = ranked_p * n / (np.arange(1, n + 1))

    # Enforce monotonicity via right-to-left cumulative minimum
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    result         = np.empty(n)
    result[order]  = np.minimum(adjusted, 1.0)
    return result


# =============================================================================
# Public API
# =============================================================================


def find_dmrs(
    df: pd.DataFrame,
    clean_samples: list[str],
    normalized_col: str = "beta_value",
) -> pd.DataFrame:
    """
    Sliding-window Wilcoxon DMR caller on a clean, pre-normalized DataFrame.

    Parameters
    ----------
    df             : long-format methylation DataFrame (QC-filtered & normalized)
    clean_samples  : list of sample IDs that passed QC.  The function ASSERTS
                     that df contains ONLY these samples.
    normalized_col : column to use as methylation signal
                     ('beta_value' or 'beta_normalized' after normalizer.robust_normalize)

    Returns
    -------
    pd.DataFrame with columns:
        window_id, cpgs, n_cpgs, case_mean, ctrl_mean, delta_beta,
        wilcoxon_stat, p_value, p_adj, significant
    Sorted by p_adj ascending.
    """
    # ── Safety assertion ───────────────────────────────────────────────────────
    present = set(df["sample_id"].unique())
    allowed = set(clean_samples)
    contam  = present - allowed
    assert not contam, (
        f"find_dmrs SAFETY VIOLATION: DataFrame contains samples NOT in "
        f"clean_samples_list: {sorted(contam)}.  "
        f"Filter df to clean_samples before calling."
    )

    # ── Exclude VDJ / clonal-flagged CpGs ────────────────────────────────────
    df_clean = df[~df["is_vdj_region"].astype(bool)].copy()

    # ── Pivot: CpG × Sample ───────────────────────────────────────────────────
    pivot = df_clean.pivot_table(
        index="cpg_id", columns="sample_id", values=normalized_col
    ).fillna(df_clean[normalized_col].mean())

    # Sort CpGs by numeric index (proxy for chromosomal order in mock data)
    cpg_order = sorted(pivot.index.tolist(), key=lambda c: int(c.lstrip("cg")))
    pivot = pivot.loc[cpg_order]

    # ── Case / Control sample lists ────────────────────────────────────────────
    meta = df_clean[["sample_id", "disease_label"]].drop_duplicates("sample_id")
    case_sids = meta.loc[meta["disease_label"] == "Case",    "sample_id"].tolist()
    ctrl_sids = meta.loc[meta["disease_label"] == "Control", "sample_id"].tolist()
    case_sids = [s for s in case_sids if s in pivot.columns]
    ctrl_sids = [s for s in ctrl_sids if s in pivot.columns]

    # ── Sliding-window Wilcoxon ────────────────────────────────────────────────
    cpgs     = pivot.index.tolist()
    n_total  = len(cpgs)
    records  = []

    for start in range(0, n_total - WINDOW_SIZE + 1, STEP_SIZE):
        window_cpgs = cpgs[start: start + WINDOW_SIZE]
        if len(window_cpgs) < MIN_CPGS:
            continue

        window_data = pivot.loc[window_cpgs]
        case_vals   = window_data[case_sids].values.flatten()
        ctrl_vals   = window_data[ctrl_sids].values.flatten()

        case_vals = case_vals[~np.isnan(case_vals)]
        ctrl_vals = ctrl_vals[~np.isnan(ctrl_vals)]

        if len(case_vals) < 2 or len(ctrl_vals) < 2:
            continue

        stat, p = ranksums(case_vals, ctrl_vals)
        delta   = float(np.mean(case_vals) - np.mean(ctrl_vals))

        records.append(
            {
                "window_id":      f"w{start:05d}",
                "cpgs":           ",".join(window_cpgs),
                "n_cpgs":         len(window_cpgs),
                "case_mean":      round(float(np.mean(case_vals)), 5),
                "ctrl_mean":      round(float(np.mean(ctrl_vals)), 5),
                "delta_beta":     round(delta, 5),
                "wilcoxon_stat":  round(float(stat), 4),
                "p_value":        float(p),
            }
        )

    if not records:
        return pd.DataFrame(
            columns=[
                "window_id", "cpgs", "n_cpgs", "case_mean", "ctrl_mean",
                "delta_beta", "wilcoxon_stat", "p_value", "p_adj", "significant",
            ]
        )

    result = pd.DataFrame(records)

    # ── BH correction ─────────────────────────────────────────────────────────
    result["p_adj"] = _bh_correction(result["p_value"].values)

    # ── Apply DMR filter criteria ──────────────────────────────────────────────
    result["significant"] = (
        (result["p_adj"]           < P_ADJ_THRESH)
        & (result["delta_beta"].abs() > DELTA_BETA_MIN)
        & (result["n_cpgs"]          >= MIN_CPGS)
    )

    return result.sort_values("p_adj").reset_index(drop=True)


# =============================================================================
# __main__ — run DMR calling on mock data (clean samples only)
# =============================================================================

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))          # for io_utils
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # for core.*
    from core.qc_guard import audit_quality
    from io_utils import write_audit_log  # noqa: E402

    MODULE = "DMR_HUNTER"
    _now   = datetime.now()
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base  = os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    )
    _audit_csv = os.path.join(_base, "data", f"audit_log_{ts_tag}.csv")

    audit_entries = []

    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def ae(sample_id, status, description, metric):
        return {
            "timestamp":   datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            "module":      MODULE,
            "sample_id":   sample_id,
            "status":      status,
            "description": description,
            "metric":      metric,
        }

    CSV = os.path.join(_base, "data", "mock_methylation.csv")
    print(f"[{ts()}] [DMR_HUNTER] Loading {CSV}")
    df = pd.read_csv(CSV)

    # ── QC gate ────────────────────────────────────────────────────────────────
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)].copy()
    print(f"[{ts()}] [DMR_HUNTER]           | Clean samples: n={len(clean_samples)}")
    audit_entries.append(ae(
        "cohort", "INFO", "Clean samples loaded for DMR analysis",
        f"n={len(clean_samples)}",
    ))

    # ── Call DMRs ──────────────────────────────────────────────────────────────
    dmrs = find_dmrs(df_clean, clean_samples)
    sig  = dmrs[dmrs["significant"]]

    print(
        f"[{ts()}] [DMR_HUNTER] DETECTED | Significant DMRs | "
        f"n={len(sig)} of {len(dmrs)} windows tested "
        f"(p_adj<{P_ADJ_THRESH}, |ΔBeta|>{DELTA_BETA_MIN}, n_cpgs≥{MIN_CPGS})"
    )
    audit_entries.append(ae(
        "cohort", "INFO", "Sliding window DMR scan complete",
        f"n_windows={len(dmrs)}",
    ))

    if len(sig):
        for _, row in sig.head(10).iterrows():
            print(
                f"[{ts()}] [DMR_HUNTER]           | {row.window_id} | "
                f"ΔBeta={row.delta_beta:+.4f}  "
                f"p_adj={row.p_adj:.3e}  "
                f"n_cpgs={row.n_cpgs}"
            )
            audit_entries.append(ae(
                "cohort", "DETECTED",
                f"Significant DMR — {row.window_id}",
                f"delta_beta={row.delta_beta:+.4f} p_adj={row.p_adj:.3e}",
            ))
    else:
        print(
            f"[{ts()}] [DMR_HUNTER]           | No significant DMRs — "
            f"check clean_samples, normalization, and thresholds"
        )
        audit_entries.append(ae(
            "cohort", "INFO", "Significant DMRs — none detected",
            f"n_sig=0 of {len(dmrs)} windows",
        ))

    write_audit_log(audit_entries, _audit_csv)
    print(f"[{ts()}] [DMR_HUNTER] Audit log → {_audit_csv}")
