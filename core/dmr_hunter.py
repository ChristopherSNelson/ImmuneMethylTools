"""
core/dmr_hunter.py — ImmuneMethylTools Strict Analyst
======================================================
Differentially Methylated Region (DMR) detection using distance-based CpG
clustering and a non-parametric test on clean, normalized data.

Biological intent
-----------------
A DMR is a contiguous stretch of CpGs where the methylation difference between
disease groups is:
  (a) Statistically significant after multiple-testing correction (p-adj < 0.05)
  (b) Biologically meaningful  (|ΔBeta| > 0.10 — below this, differences have
      negligible impact on transcription factor binding affinity)
  (c) Regionally consistent    (>= 3 CpGs per cluster — single-site hits are
      likely noise; a region requires multiple concordant signals)

Safety guarantees
-----------------
  1. ASSERTION: the input DataFrame must contain ONLY samples from
     clean_samples_list — preventing artifact-contaminated data from biasing
     DMR calls.  Passing uncleaned data raises AssertionError immediately.
  2. VDJ-region CpGs are INCLUDED but every cluster is annotated with
     `n_vdj_cpgs` (count of VDJ CpGs in the cluster) and a boolean
     `clonal_risk` flag.  Significant DMRs with clonal_risk=True are logged
     as DETECTED so the Analyst can decide whether to accept or exclude them.

Statistical approach
--------------------
Distance-based CpG clustering using genomic coordinates (chrom, pos):
  - Sort CpGs by (chrom, pos)
  - Cluster: group consecutive CpGs where each is within MAX_GAP_BP (1000 bp)
    of the next; a new chromosome always starts a new cluster
  - Filter:  only analyze clusters with >= MIN_CPGS (3) CpGs
  - Aggregation:  Per-sample median beta within the cluster (avoids
                  pseudoreplication; each sample contributes one observation)
  - Test:        Wilcoxon rank-sum (non-parametric; robust to non-Gaussian betas)
  - Correction:  Benjamini-Hochberg (BH) FDR on cluster-level p-values
  - Filter:      p_adj < 0.05, |ΔBeta| > 0.10, >= 3 CpGs

This approach reduces the number of statistical tests compared to a
sliding window, improving FDR power and producing biologically coherent
regions that correspond to actual genomic neighborhoods.
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import ranksums

# ── Parameters ────────────────────────────────────────────────────────────────
MAX_GAP_BP = 1000     # maximum inter-CpG distance (bp) within a cluster
P_ADJ_THRESH = 0.05   # BH-corrected p-value threshold
DELTA_BETA_MIN = 0.10   # minimum |ΔBeta| to qualify as a DMR
MIN_CPGS = 3      # minimum CpGs per cluster


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


def _build_clusters(
    df: pd.DataFrame,
    max_gap: int = MAX_GAP_BP,
    min_cpgs: int = MIN_CPGS,
) -> list[list[str]]:
    """
    Group CpGs into distance-based clusters using genomic coordinates.

    CpGs are sorted by (chrom, pos).  A new cluster starts whenever the gap
    to the next CpG exceeds max_gap bp or the chromosome changes.  Only
    clusters with >= min_cpgs CpGs are returned.

    Parameters
    ----------
    df       : long-format DataFrame with cpg_id, chrom, pos columns
    max_gap  : maximum inter-CpG distance in bp (default 1000)
    min_cpgs : minimum CpGs per cluster (default 3)

    Returns
    -------
    list of lists, each containing cpg_id strings for one cluster
    """
    # Deduplicate to one row per CpG
    cpg_coords = (
        df[["cpg_id", "chrom", "pos"]]
        .drop_duplicates("cpg_id")
        .sort_values(["chrom", "pos"])
        .reset_index(drop=True)
    )

    clusters: list[list[str]] = []
    current_cluster: list[str] = []
    prev_chrom: str | None = None
    prev_pos: int = 0

    for _, row in cpg_coords.iterrows():
        if (
            prev_chrom is not None
            and row["chrom"] == prev_chrom
            and (row["pos"] - prev_pos) <= max_gap
        ):
            current_cluster.append(row["cpg_id"])
        else:
            if len(current_cluster) >= min_cpgs:
                clusters.append(current_cluster)
            current_cluster = [row["cpg_id"]]
        prev_chrom = row["chrom"]
        prev_pos = row["pos"]

    # Don't forget the last cluster
    if len(current_cluster) >= min_cpgs:
        clusters.append(current_cluster)

    return clusters


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
    Distance-based cluster Wilcoxon DMR caller on a clean, pre-normalized
    DataFrame.

    CpGs are clustered by genomic proximity (within 1000 bp), then each
    cluster with >= 3 CpGs is tested for differential methylation between
    Case and Control groups using per-sample median aggregation and the
    Wilcoxon rank-sum test.

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
    chunk_size     : reserved for EPIC/WGBS scaling (not used in clustering mode)

    Returns
    -------
    pd.DataFrame with columns:
        window_id, chrom, start_pos, end_pos, cpgs, n_cpgs, case_mean,
        ctrl_mean, delta_beta, wilcoxon_stat, p_value, p_adj, significant,
        n_vdj_cpgs, clonal_risk, mean_gc
    Sorted by p_adj ascending.
    `clonal_risk` is True when any CpG in the cluster overlaps a VDJ locus.
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

    # ── Identify VDJ CpGs (annotated per-cluster; not excluded) ────────────────
    df_clean = df.copy()
    vdj_cpgs = set(
        df_clean.loc[df_clean["is_vdj_region"].astype(bool), "cpg_id"].unique()
    )

    # ── Build CpG → gc_content lookup (one value per CpG) ────────────────────
    gc_map = {}
    if "gc_content" in df_clean.columns:
        gc_map = (
            df_clean[["cpg_id", "gc_content"]]
            .drop_duplicates("cpg_id")
            .set_index("cpg_id")["gc_content"]
            .to_dict()
        )

    # ── Build CpG coordinate lookup ───────────────────────────────────────────
    cpg_coord_df = (
        df_clean[["cpg_id", "chrom", "pos"]]
        .drop_duplicates("cpg_id")
        .set_index("cpg_id")
    )

    # ── Build distance-based clusters ──────────────────────────────────────────
    clusters = _build_clusters(df_clean)

    # ── Case / Control sample lists ────────────────────────────────────────────
    global_mean = float(df_clean[normalized_col].mean())
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

    # ── Build pivot table (CpG x Sample) ───────────────────────────────────────
    # Only pivot CpGs that belong to clusters (reduces memory for large datasets)
    clustered_cpgs = [cpg for cluster in clusters for cpg in cluster]
    pivot_df = df_clean[df_clean["cpg_id"].isin(clustered_cpgs)]
    pivot = pivot_df.pivot_table(
        index="cpg_id", columns="sample_id", values=normalized_col
    ).fillna(global_mean)
    case_sids = [s for s in case_sids if s in pivot.columns]
    ctrl_sids = [s for s in ctrl_sids if s in pivot.columns]

    # ── Test each cluster ──────────────────────────────────────────────────────
    records = []
    for ci, cluster_cpgs in enumerate(clusters):
        # Extract cluster data from pivot
        present_cpgs = [c for c in cluster_cpgs if c in pivot.index]
        if len(present_cpgs) < MIN_CPGS:
            continue

        c_data = pivot.loc[present_cpgs]

        # Per-sample MEDIAN across CpGs in the cluster
        c_medians = c_data[case_sids].median(axis=0).values
        t_medians = c_data[ctrl_sids].median(axis=0).values
        c_medians = c_medians[~np.isnan(c_medians)]
        t_medians = t_medians[~np.isnan(t_medians)]

        if len(c_medians) < 2 or len(t_medians) < 2:
            continue

        stat, p = ranksums(c_medians, t_medians)
        delta = float(np.mean(c_medians) - np.mean(t_medians))
        n_vdj = sum(1 for c in present_cpgs if c in vdj_cpgs)

        # Mean GC content across CpGs in the cluster
        if gc_map:
            gc_vals = [gc_map[c] for c in present_cpgs if c in gc_map]
            mean_gc = round(float(np.mean(gc_vals)), 4) if gc_vals else None
        else:
            mean_gc = None

        # Cluster genomic span
        coords = cpg_coord_df.loc[present_cpgs]
        chrom = coords["chrom"].iloc[0]
        start_pos = int(coords["pos"].min())
        end_pos = int(coords["pos"].max())

        records.append({
            "window_id": f"cl{ci:05d}",
            "chrom": chrom,
            "start_pos": start_pos,
            "end_pos": end_pos,
            "cpgs": ",".join(present_cpgs),
            "n_cpgs": len(present_cpgs),
            "case_mean": round(float(np.mean(c_medians)), 5),
            "ctrl_mean": round(float(np.mean(t_medians)), 5),
            "delta_beta": round(delta, 5),
            "wilcoxon_stat": round(float(stat), 4),
            "p_value": float(p),
            "n_vdj_cpgs": n_vdj,
            "mean_gc": mean_gc,
        })

    _empty_cols = [
        "window_id", "chrom", "start_pos", "end_pos", "cpgs", "n_cpgs",
        "case_mean", "ctrl_mean", "delta_beta", "wilcoxon_stat", "p_value",
        "p_adj", "significant", "n_vdj_cpgs", "clonal_risk", "mean_gc",
    ]
    if not records:
        return pd.DataFrame(columns=_empty_cols)

    result = pd.DataFrame(records)

    # ── BH correction ─────────────────────────────────────────────────────────
    result["p_adj"] = _bh_correction(result["p_value"].values)

    # ── Apply DMR filter criteria ──────────────────────────────────────────────
    result["significant"] = (
        (result["p_adj"] < p_adj_thresh)
        & (result["delta_beta"].abs() > delta_beta_min)
        & (result["n_cpgs"] >= MIN_CPGS)
    )

    # ── Flag clusters that overlap VDJ loci ────────────────────────────────────
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
        f"n={len(sig)} of {len(dmrs)} clusters tested "
        f"(p_adj<{P_ADJ_THRESH}, |ΔBeta|>{DELTA_BETA_MIN}, n_cpgs>={MIN_CPGS})"
    )
    audit_entries.append(ae(
        "cohort", "INFO", "Distance-based cluster DMR scan complete",
        f"n_clusters={len(dmrs)}",
    ))

    if len(sig):
        for _, row in sig.head(10).iterrows():
            risk_tag = " !! HIGH CLONALITY" if row.clonal_risk else ""
            print(
                f"[{ts()}] [DMR_HUNTER]           | {row.window_id}{risk_tag} | "
                f"{row.chrom}:{row.start_pos:,}-{row.end_pos:,}  "
                f"ΔBeta={row.delta_beta:+.4f}  "
                f"p_adj={row.p_adj:.3e}  "
                f"n_cpgs={row.n_cpgs}  "
                f"n_vdj_cpgs={row.n_vdj_cpgs}"
            )
            audit_entries.append(ae(
                "cohort", "DETECTED",
                f"Significant DMR — {row.window_id} ({row.chrom}:{row.start_pos:,}-{row.end_pos:,})"
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
            f"n_sig=0 of {len(dmrs)} clusters",
        ))

    # ── Clonal risk summary ────────────────────────────────────────────────────
    clonal_clusters = dmrs[dmrs["clonal_risk"]]
    sig_clonal = sig[sig["clonal_risk"]] if len(sig) else pd.DataFrame()
    n_cr = len(clonal_clusters)
    n_sc = len(sig_clonal)
    print(
        f"[{ts()}] [DMR_HUNTER] {'DETECTED' if n_sc else 'INFO    '} | "
        f"VDJ clonal_risk clusters | "
        f"{n_cr} total ({n_sc} significant) — analyst review required for flagged clusters"
    )
    audit_entries.append(ae(
        "cohort",
        "DETECTED" if n_sc else "INFO",
        "VDJ clonal_risk cluster summary",
        f"n_clonal_risk={n_cr} n_sig_clonal={n_sc}",
    ))
    for _, row in sig_clonal.iterrows():
        print(
            f"[{ts()}] [DMR_HUNTER] DETECTED | HIGH CLONALITY — significant DMR in VDJ locus | "
            f"{row.window_id}  {row.chrom}:{row.start_pos:,}-{row.end_pos:,}  "
            f"ΔBeta={row.delta_beta:+.4f}  n_vdj_cpgs={row.n_vdj_cpgs}"
        )
        audit_entries.append(ae(
            "cohort", "DETECTED",
            f"HIGH CLONALITY — significant DMR overlaps VDJ locus: {row.window_id}",
            f"delta_beta={row.delta_beta:+.4f} n_vdj_cpgs={row.n_vdj_cpgs}",
        ))

    write_audit_log(audit_entries, _audit_csv)
    print(f"[{ts()}] [DMR_HUNTER] Audit log → {_audit_csv}")
