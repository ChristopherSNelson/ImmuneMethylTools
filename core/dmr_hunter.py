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
     clean_samples_list — preventing artifact-contaminated data from biasing
     DMR calls.  Passing uncleaned data raises AssertionError immediately.
  2. VDJ-region CpGs are INCLUDED but every window is annotated with
     `n_vdj_cpgs` (count of VDJ CpGs in the window) and a boolean
     `clonal_risk` flag.  Significant DMRs with clonal_risk=True are logged
     as DETECTED so the Analyst can decide whether to accept or exclude them.
     This replaces the previous blanket exclusion — retaining the data while
     making the risk explicit.

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
WINDOW_SIZE = 5      # sliding window width (CpGs)
STEP_SIZE = 1      # step between windows
P_ADJ_THRESH = 0.05   # BH-corrected p-value threshold
DELTA_BETA_MIN = 0.10   # minimum |ΔBeta| to qualify as a DMR
MIN_CPGS = 3      # minimum CpGs per window


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

    order = np.argsort(pvalues)
    ranked_p = pvalues[order]
    adjusted = ranked_p * n / (np.arange(1, n + 1))

    # Enforce monotonicity via right-to-left cumulative minimum
    for i in range(n - 2, -1, -1):
        adjusted[i] = min(adjusted[i], adjusted[i + 1])

    result = np.empty(n)
    result[order] = np.minimum(adjusted, 1.0)
    return result


# =============================================================================
# Public API
# =============================================================================


def find_dmrs(
    df: pd.DataFrame,
    clean_samples: list[str],
    normalized_col: str = "beta_value",
    min_site_depth: int = 5,
    p_adj_thresh: float = P_ADJ_THRESH,
    delta_beta_min: float = DELTA_BETA_MIN,
    chunk_size: int | None = None,
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
    min_site_depth : per-row minimum read depth; rows below are excluded before
                     the pivot (default: 5 — matches SITE_DEPTH_THRESH in qc_guard)
    p_adj_thresh   : BH-corrected p-value cutoff (default 0.05)
    delta_beta_min : minimum |ΔBeta| to qualify as a DMR (default 0.10)
    chunk_size     : number of CpGs to pivot at once; None = load all CpGs into
                     memory (default, fine for mock data and small EPIC arrays).
                     Set to e.g. 50_000 for EPIC arrays or 100_000 for WGBS to
                     keep per-chunk memory usage below ~500 MB for 100 samples.

    Returns
    -------
    pd.DataFrame with columns:
        window_id, cpgs, n_cpgs, case_mean, ctrl_mean, delta_beta,
        wilcoxon_stat, p_value, p_adj, significant, n_vdj_cpgs, clonal_risk
    Sorted by p_adj ascending.
    `clonal_risk` is True when any CpG in the window overlaps a VDJ locus.
    The Analyst should review these windows before reporting them.
    """
    # ── Input validation ───────────────────────────────────────────────────────
    present = set(df["sample_id"].unique())
    allowed = set(clean_samples)
    contam = present - allowed
    assert not contam, (
        f"find_dmrs DATA VALIDATION ERROR: DataFrame contains samples NOT in "
        f"clean_samples_list: {sorted(contam)}.  "
        f"Filter df to clean_samples before calling."
    )

    # ── Exclude low-depth sites ────────────────────────────────────────────────
    if min_site_depth > 0:
        df = df[df["depth"] >= min_site_depth].copy()

    # ── Identify VDJ CpGs (annotated per-window; not excluded) ────────────────
    df_clean = df.copy()
    vdj_cpgs = set(
        df_clean.loc[df_clean["is_vdj_region"].astype(bool), "cpg_id"].unique()
    )

    # ── Sort CpGs by numeric index (proxy for chromosomal order in mock data) ─
    global_mean = float(df_clean[normalized_col].mean())
    cpg_order = sorted(df_clean["cpg_id"].unique().tolist(), key=lambda c: int(c.lstrip("cg")))

    # ── Case / Control sample lists ────────────────────────────────────────────
    meta = df_clean[["sample_id", "disease_label"]].drop_duplicates("sample_id")
    all_samples_in_df = set(df_clean["sample_id"].unique())
    case_sids = [
        s for s in meta.loc[meta["disease_label"] == "Case", "sample_id"].tolist()
        if s in all_samples_in_df
    ]
    ctrl_sids = [
        s for s in meta.loc[meta["disease_label"] == "Control", "sample_id"].tolist()
        if s in all_samples_in_df
    ]

    # ── Helper: run Wilcoxon over a sorted CpG list given a pre-built pivot ───
    def _window_records(cpg_list, pivot_df, global_start_offset):
        recs = []
        n = len(cpg_list)
        for local_start in range(0, n - WINDOW_SIZE + 1, STEP_SIZE):
            w_cpgs = cpg_list[local_start: local_start + WINDOW_SIZE]
            if len(w_cpgs) < MIN_CPGS:
                continue
            w_data = pivot_df.loc[w_cpgs]
            c_vals = w_data[case_sids].values.flatten()
            t_vals = w_data[ctrl_sids].values.flatten()
            c_vals = c_vals[~np.isnan(c_vals)]
            t_vals = t_vals[~np.isnan(t_vals)]
            if len(c_vals) < 2 or len(t_vals) < 2:
                continue
            stat, p = ranksums(c_vals, t_vals)
            delta = float(np.mean(c_vals) - np.mean(t_vals))
            n_vdj = sum(1 for c in w_cpgs if c in vdj_cpgs)
            recs.append({
                "window_id": f"w{global_start_offset + local_start:05d}",
                "cpgs": ",".join(w_cpgs),
                "n_cpgs": len(w_cpgs),
                "case_mean": round(float(np.mean(c_vals)), 5),
                "ctrl_mean": round(float(np.mean(t_vals)), 5),
                "delta_beta": round(delta, 5),
                "wilcoxon_stat": round(float(stat), 4),
                "p_value": float(p),
                "n_vdj_cpgs": n_vdj,
            })
        return recs

    # ── Sliding-window Wilcoxon — in-memory or chunked ─────────────────────────
    n_total = len(cpg_order)
    records = []

    use_chunks = chunk_size is not None and chunk_size < n_total
    if not use_chunks:
        # In-memory path: build full CpG × Sample pivot once
        pivot = df_clean.pivot_table(
            index="cpg_id", columns="sample_id", values=normalized_col
        ).fillna(global_mean).reindex(cpg_order, fill_value=global_mean)
        case_sids = [s for s in case_sids if s in pivot.columns]
        ctrl_sids = [s for s in ctrl_sids if s in pivot.columns]
        records = _window_records(cpg_order, pivot, global_start_offset=0)
    else:
        # Chunked path: process chunk_size CpGs at a time.
        # Each chunk includes a right border of (WINDOW_SIZE - 1) CpGs so that
        # windows spanning the chunk boundary are processed by the chunk that
        # contains their START position.
        border = WINDOW_SIZE - 1
        for chunk_start in range(0, n_total, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_total)
            ext_end = min(chunk_end + border, n_total)
            ext_cpgs = cpg_order[chunk_start:ext_end]

            chunk_df = df_clean[df_clean["cpg_id"].isin(ext_cpgs)]
            chunk_pivot = chunk_df.pivot_table(
                index="cpg_id", columns="sample_id", values=normalized_col
            ).fillna(global_mean).reindex(ext_cpgs, fill_value=global_mean)

            # Only keep windows whose global start falls within [chunk_start, chunk_end).
            # Windows starting in the border region are skipped here and will be
            # emitted by the next chunk (where they start within its main range).
            chunk_recs = _window_records(ext_cpgs, chunk_pivot, global_start_offset=chunk_start)
            records.extend(r for r in chunk_recs if int(r["window_id"][1:]) < chunk_end)

    if not records:
        return pd.DataFrame(
            columns=[
                "window_id", "cpgs", "n_cpgs", "case_mean", "ctrl_mean",
                "delta_beta", "wilcoxon_stat", "p_value", "p_adj",
                "significant", "n_vdj_cpgs", "clonal_risk",
            ]
        )

    result = pd.DataFrame(records)

    # ── BH correction ─────────────────────────────────────────────────────────
    result["p_adj"] = _bh_correction(result["p_value"].values)

    # ── Apply DMR filter criteria ──────────────────────────────────────────────
    result["significant"] = (
        (result["p_adj"] < p_adj_thresh)
        & (result["delta_beta"].abs() > delta_beta_min)
        & (result["n_cpgs"] >= MIN_CPGS)
    )

    # ── Flag windows that overlap VDJ loci ────────────────────────────────────
    result["clonal_risk"] = result["n_vdj_cpgs"] > 0

    return result.sort_values("p_adj").reset_index(drop=True)


# =============================================================================
# __main__ — run DMR calling on mock data (clean samples only)
# =============================================================================

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))          # for io_utils
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # for core.*
    from core.qc_guard import audit_quality
    from io_utils import data_path, load_methylation, project_root, write_audit_log  # noqa: E402

    MODULE = "DMR_HUNTER"
    _now = datetime.now()
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base = project_root()
    _audit_csv = os.path.join(_base, "output", f"audit_log_{MODULE}_{ts_tag}.csv")

    audit_entries = []

    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def ae(sample_id, status, description, metric):
        return {
            "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            "module": MODULE,
            "sample_id": sample_id,
            "status": status,
            "description": description,
            "metric": metric,
        }

    df = load_methylation(data_path("mock_methylation.csv"))

    # ── QC filter ──────────────────────────────────────────────────────────────
    clean_samples = audit_quality(df)
    df_clean = df[df["sample_id"].isin(clean_samples)].copy()
    print(f"[{ts()}] [DMR_HUNTER]           | Clean samples: n={len(clean_samples)}")
    audit_entries.append(ae(
        "cohort", "INFO", "Clean samples loaded for DMR analysis",
        f"n={len(clean_samples)}",
    ))

    # ── Call DMRs ──────────────────────────────────────────────────────────────
    dmrs = find_dmrs(df_clean, clean_samples)
    sig = dmrs[dmrs["significant"]]

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
            risk_tag = " ⚠ HIGH CLONALITY" if row.clonal_risk else ""
            print(
                f"[{ts()}] [DMR_HUNTER]           | {row.window_id}{risk_tag} | "
                f"ΔBeta={row.delta_beta:+.4f}  "
                f"p_adj={row.p_adj:.3e}  "
                f"n_cpgs={row.n_cpgs}  "
                f"n_vdj_cpgs={row.n_vdj_cpgs}"
            )
            status = "DETECTED" if row.clonal_risk else "DETECTED"
            audit_entries.append(ae(
                "cohort", status,
                f"Significant DMR — {row.window_id}"
                + (" — HIGH CLONALITY (VDJ overlap)" if row.clonal_risk else ""),
                f"delta_beta={row.delta_beta:+.4f} p_adj={row.p_adj:.3e} n_vdj={row.n_vdj_cpgs}",
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

    # ── Clonal risk summary ────────────────────────────────────────────────────
    clonal_windows = dmrs[dmrs["clonal_risk"]]
    sig_clonal = sig[sig["clonal_risk"]] if len(sig) else pd.DataFrame()
    n_cr = len(clonal_windows)
    n_sc = len(sig_clonal)
    print(
        f"[{ts()}] [DMR_HUNTER] {'DETECTED' if n_sc else 'INFO    '} | "
        f"VDJ clonal_risk windows | "
        f"{n_cr} total ({n_sc} significant) — analyst review required for flagged windows"
    )
    audit_entries.append(ae(
        "cohort",
        "DETECTED" if n_sc else "INFO",
        "VDJ clonal_risk window summary",
        f"n_clonal_risk={n_cr} n_sig_clonal={n_sc}",
    ))
    for _, row in sig_clonal.iterrows():
        print(
            f"[{ts()}] [DMR_HUNTER] DETECTED | HIGH CLONALITY — significant DMR in VDJ locus | "
            f"{row.window_id}  ΔBeta={row.delta_beta:+.4f}  n_vdj_cpgs={row.n_vdj_cpgs}"
        )
        audit_entries.append(ae(
            "cohort", "DETECTED",
            f"HIGH CLONALITY — significant DMR overlaps VDJ locus: {row.window_id}",
            f"delta_beta={row.delta_beta:+.4f} n_vdj_cpgs={row.n_vdj_cpgs}",
        ))

    write_audit_log(audit_entries, _audit_csv)
    print(f"[{ts()}] [DMR_HUNTER] Audit log → {_audit_csv}")
