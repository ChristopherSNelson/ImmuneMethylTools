"""
core/qc_guard.py — ImmuneMethylTools QC Gatekeeper
====================================================
First-pass sample filter; MUST run before normalization or DMR analysis.

Biological intent
-----------------
Three laboratory failures leave quantitative fingerprints that QC catches:

  1. Incomplete bisulfite conversion — cytosines that should convert to uracil
     remain as cytosine.  The diagnostic: non-CpG methylation rate > 1 %.
     In a healthy conversion, non-CpG cytosines are fully converted (rate ≈ 0.4 %).
     We use 1 % (rather than the widely cited 2 %) for conservative early exclusion.

  2. Low coverage — sites with mean depth < 10 reads have inflated binomial
     sampling variance; beta values at such sites are statistically unreliable.

  3. Contamination / 'muddy' samples — cross-contamination from a foreign cell
     type collapses the characteristic bimodal beta distribution into a unimodal
     hump near 0.5.  Sarle's Bimodality Coefficient (BC) quantifies this:
       BC > 0.555 → bimodal (clean)
       BC < 0.555 → unimodal/muddy (contaminated)
     We require BOTH low BC AND a mean beta in the muddy range [0.40, 0.65]
     to reduce false positives from genuinely unimodal loci.

Architecture note
-----------------
audit_quality() returns a clean_samples_list.  downstream modules
(normalizer, dmr_hunter, ml_guard) MUST consume this list and filter
before proceeding — batch correction of artifact-contaminated samples
will absorb the artifact signal into the correction model.
"""

import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew

# ── Thresholds ────────────────────────────────────────────────────────────────
BISULFITE_FAIL_THRESH = 0.01    # non-CpG meth rate; > 1 % → bisulfite failure
DEPTH_FAIL_THRESH = 10      # reads; mean depth below this → unreliable betas
# Contamination: flag samples whose BC falls > BC_SIGMA_THRESH standard
# deviations below the cohort median AND whose mean beta sits in the muddy
# range [0.40, 0.65].  A relative threshold is more robust than Sarle's fixed
# 0.555 because clean WGBS bimodal distributions cluster tightly near BC ≈ 0.90;
# even modest contamination (50 % mix) reduces BC to ~0.78, well outside 2σ.
BC_SIGMA_THRESH = 2.0     # z-score below cohort BC median to flag
CONTAMINATION_MEAN_LO = 0.40    # contaminated sample mean beta floor
CONTAMINATION_MEAN_HI = 0.65    # contaminated sample mean beta ceiling
# Site-level depth filter (applied before pivot in DMR/ML stages)
SITE_DEPTH_THRESH = 5    # per-site minimum depth; rows below → excluded
SITE_LOW_DEPTH_SAMPLE_WARN = 20.0  # % of sites below threshold → DETECTED audit entry


# =============================================================================
# Internal helpers
# =============================================================================


def _sarle_bimodality_coefficient(values: np.ndarray) -> float:
    """
    Sarle's Bimodality Coefficient (BC).

    BC = (skewness² + 1) / (excess_kurtosis + 3·(n−1)² / ((n−2)(n−3)))

    BC > 0.555  ⟹ bimodal distribution (clean methylation profile)
    BC < 0.555  ⟹ unimodal distribution (contaminated / muddy)

    Returns np.nan for samples with fewer than 4 values or degenerate
    denominators (e.g., all-constant distributions).
    """
    n = len(values)
    if n < 4:
        return np.nan

    m3 = skew(values)
    m4 = kurtosis(values, fisher=True)  # excess kurtosis (Gaussian = 0)
    correction = 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    denom = m4 + correction
    if denom <= 0:
        return np.nan
    return float((m3 ** 2 + 1) / denom)


# =============================================================================
# Public API
# =============================================================================


def audit_quality(df: pd.DataFrame) -> list[str]:
    """
    Flag samples failing bisulfite conversion or coverage thresholds.

    Parameters
    ----------
    df : long-format methylation DataFrame (mock_methylation schema)

    Returns
    -------
    clean_samples_list : list of sample IDs that PASS all QC filters
    """
    sample_stats = df.groupby("sample_id").agg(
        mean_ncpg=("non_cpg_meth_rate", "mean"),
        mean_depth=("depth", "mean"),
    )

    bisulfite_fail = sample_stats["mean_ncpg"] > BISULFITE_FAIL_THRESH
    depth_fail = sample_stats["mean_depth"] < DEPTH_FAIL_THRESH
    any_fail = bisulfite_fail | depth_fail

    return sample_stats[~any_fail].index.tolist()


def detect_contamination(df: pd.DataFrame) -> tuple[list[str], pd.DataFrame]:
    """
    Distribution-based contamination check using Sarle's Bimodality Coefficient.

    A clean WGBS/RRBS sample has a bimodal beta distribution (unmethylated ~0.1,
    methylated ~0.85).  Contamination collapses this into a unimodal hump near
    0.5 — the 'muddy' fingerprint.

    Parameters
    ----------
    df : long-format methylation DataFrame

    Returns
    -------
    flagged_samples : list of sample IDs flagged as contaminated
    report_df       : DataFrame with per-sample bimodality statistics
    """
    # Compute per-sample BC and mean beta
    records = []
    for sid, grp in df.groupby("sample_id"):
        betas = grp["beta_value"].dropna().values
        bc = _sarle_bimodality_coefficient(betas)
        mean_b = float(betas.mean())
        records.append({"sample_id": sid, "bimodality_coeff": bc, "mean_beta": mean_b})

    report = pd.DataFrame(records).set_index("sample_id")

    # Cohort-relative threshold: flag samples > BC_SIGMA_THRESH σ below median
    valid_bc = report["bimodality_coeff"].dropna()
    bc_median = float(valid_bc.median())
    bc_std = float(valid_bc.std())
    bc_threshold = bc_median - BC_SIGMA_THRESH * bc_std

    report["low_bimodality"] = report["bimodality_coeff"] < bc_threshold
    report["muddy_mean"] = report["mean_beta"].between(
        CONTAMINATION_MEAN_LO, CONTAMINATION_MEAN_HI
    )
    report["contamination_flag"] = report["low_bimodality"] & report["muddy_mean"]
    report["bimodality_coeff"] = report["bimodality_coeff"].round(4)
    report["mean_beta"] = report["mean_beta"].round(4)

    flagged_samples = report[report["contamination_flag"]].index.tolist()
    return flagged_samples, report


def filter_site_quality(
    df: pd.DataFrame,
    min_depth: int = SITE_DEPTH_THRESH,
) -> tuple[pd.DataFrame, dict]:
    """
    Remove individual CpG rows whose read depth falls below *min_depth*.

    This is a site-level filter (not sample-level): a single sample can
    contribute both passing and failing rows.  Low-depth sites have inflated
    binomial variance; excluding them before pivot reduces noise in the
    DMR-hunter and ML feature matrices.

    Parameters
    ----------
    df        : long-format methylation DataFrame
    min_depth : per-site minimum read depth (default: SITE_DEPTH_THRESH = 5)

    Returns
    -------
    df_clean : DataFrame with low-depth rows removed (copy, no `low_depth_flag` column)
    stats    : dict with keys:
                 n_total          — rows before filtering
                 n_low            — rows below threshold
                 pct_low          — percentage of rows below threshold
                 per_sample_pct   — pd.Series: % low-depth rows per sample_id
                 min_depth        — threshold used
    """
    df = df.copy()
    df["_low_depth"] = df["depth"] < min_depth
    n_total = len(df)
    n_low = int(df["_low_depth"].sum())
    pct_low = n_low / n_total * 100 if n_total else 0.0
    per_sample_pct = (
        df.groupby("sample_id")["_low_depth"].mean() * 100
    ).rename("pct_low_depth")
    df_clean = df[~df["_low_depth"]].drop(columns=["_low_depth"]).copy()
    stats = {
        "n_total": n_total,
        "n_low": n_low,
        "pct_low": pct_low,
        "per_sample_pct": per_sample_pct,
        "min_depth": min_depth,
    }
    return df_clean, stats


# =============================================================================
# __main__ — run QC on mock data
# =============================================================================

if __name__ == "__main__":
    import os
    import sys
    from datetime import datetime

    # Ensure sibling core/ modules are importable regardless of working directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from io_utils import (  # noqa: E402
        Tee, append_flagged_samples, data_path, load_methylation, project_root, write_audit_log,
    )

    MODULE = "qc_guard"
    MODULE_TAG = "QC_GUARD"
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
        n_total = df["sample_id"].nunique()
        audit_entries.append(ae("cohort", "INFO", "Samples evaluated", f"n={n_total}"))

        # ── Audit quality ──────────────────────────────────────────────────────────
        clean_samples = audit_quality(df)
        n_clean = len(clean_samples)
        n_failed = n_total - n_clean
        flagged_qc = [s for s in df["sample_id"].unique() if s not in clean_samples]
        print(
            f"[{ts()}] [QC_GUARD] DETECTED | Bisulfite/depth QC | "
            f"{n_failed}/{n_total} samples failed → {flagged_qc}"
        )
        print(f"[{ts()}] [QC_GUARD]           | Clean samples retained | n={n_clean}")

        ncpg_per_sample = df.groupby("sample_id")["non_cpg_meth_rate"].mean()
        for sid in flagged_qc:
            audit_entries.append(ae(
                sid, "DETECTED", "Bisulfite conversion failure",
                f"non_cpg_meth_rate={ncpg_per_sample[sid]:.3f}",
            ))
        audit_entries.append(ae("cohort", "INFO", "Clean samples retained", f"n={n_clean}"))

        # ── Detect contamination ───────────────────────────────────────────────────
        contaminated, report = detect_contamination(df)
        if contaminated:
            for sid in contaminated:
                row = report.loc[sid]
                print(
                    f"[{ts()}] [QC_GUARD] DETECTED | Contamination | "
                    f"sample={sid}  BC={row.bimodality_coeff:.4f}  "
                    f"mean_beta={row.mean_beta:.4f}"
                )
                audit_entries.append(ae(
                    sid, "DETECTED", "Contamination — low bimodality coefficient",
                    f"BC={row.bimodality_coeff:.4f}",
                ))
        else:
            print(
                f"[{ts()}] [QC_GUARD]           | No contamination detected "
                f"(cohort-relative BC threshold, mean∈"
                f"[{CONTAMINATION_MEAN_LO},{CONTAMINATION_MEAN_HI}])"
            )
            audit_entries.append(ae(
                "cohort", "INFO", "Contamination check — none detected",
                f"threshold=cohort_median-{BC_SIGMA_THRESH}sigma",
            ))

        print(f"\n[{ts()}] [QC_GUARD] Final clean_samples_list: {clean_samples}")

        # ── Site-level depth filter ────────────────────────────────────────────────
        _, site_stats = filter_site_quality(df)
        pct = site_stats["pct_low"]
        print(
            f"[{ts()}] [QC_GUARD]           | Site-level depth filtering | "
            f"{site_stats['n_low']}/{site_stats['n_total']} rows below "
            f"{site_stats['min_depth']}x ({pct:.1f}%)"
        )
        audit_entries.append(ae(
            "cohort", "INFO", "Site-level depth filtering complete",
            f"{pct:.1f}% sites < {site_stats['min_depth']}x dropped",
        ))
        for sid, pct_s in site_stats["per_sample_pct"].items():
            if pct_s > SITE_LOW_DEPTH_SAMPLE_WARN:
                print(
                    f"[{ts()}] [QC_GUARD] DETECTED | High proportion of low-depth sites | "
                    f"sample={sid}  {pct_s:.1f}% sites < {site_stats['min_depth']}x"
                )
                audit_entries.append(ae(
                    sid, "DETECTED", "High proportion of low-depth sites",
                    f"{pct_s:.1f}% sites < {site_stats['min_depth']}x",
                ))

        # ── Persist flagged_samples.csv ────────────────────────────────────────────
        flagged_rows = []

        for sid in flagged_qc:
            flagged_rows.append({
                "run_timestamp": run_ts,
                "module": MODULE,
                "sample_id": sid,
                "flag_type": "bisulfite_failure",
                "detail": f"non_cpg_meth_rate={ncpg_per_sample[sid]:.3f}",
            })

        for sid in contaminated:
            r = report.loc[sid]
            flagged_rows.append({
                "run_timestamp": run_ts,
                "module": MODULE,
                "sample_id": sid,
                "flag_type": "contamination",
                "detail": f"BC={r.bimodality_coeff:.2f} mean_beta={r.mean_beta:.3f}",
            })

        append_flagged_samples(flagged_rows, _csv)
        n_unique = len({r["sample_id"] for r in flagged_rows})
        print(
            f"[{ts()}] [QC_GUARD] {n_unique} unique flagged sample(s) written "
            f"→ output/flagged_samples.csv"
        )

        # ── Audit log ──────────────────────────────────────────────────────────────
        write_audit_log(audit_entries, _audit_csv)
        print(f"[{ts()}] [QC_GUARD] Audit log → {_audit_csv}")
