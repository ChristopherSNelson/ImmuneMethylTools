"""
core/normalizer.py — ImmuneMethylTools Normalizer
===================================================
Confounding check (Cramér's V) and median-centring normalization.

Biological intent
-----------------
Normalization MUST happen AFTER QC (CLAUDE.md architecture rule) to prevent
correcting artifacts into the model rather than removing them.

Confounding check
-----------------
If a batch is enriched for one disease label (Cramér's V > 0.5), naïve batch
correction will absorb the biological signal of interest.  This check raises
a warning so the analyst can use a model-based approach (e.g., ComBat-Seq
with disease as a covariate) instead of blind batch removal.

Cramér's V interpretation:
  < 0.10  — negligible association
  0.10–0.30 — weak
  0.30–0.50 — moderate
  > 0.50   — strong (confounded; downstream correction is biased)

Median centring
---------------
Preferred over quantile normalization for methylation data because:
  1. Preserves relative CpG-to-CpG differences within a sample.
  2. Robust to outlier sites (e.g., clonal VDJ loci with inflated beta).
  3. Does not distort the bimodal beta distribution shape, which is a
     biologically meaningful feature used by contamination detectors.

Each sample's median beta is subtracted from its own beta values, centering
all samples at zero without assuming any particular distribution shape.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import chi2_contingency  # noqa: E402

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── Thresholds ────────────────────────────────────────────────────────────────
CRAMERS_V_WARN = 0.50   # V > 0.5 → strong association; flag as confounded


# =============================================================================
# Internal helpers
# =============================================================================


def _cramers_v(contingency: np.ndarray) -> float:
    """Cramér's V effect size from a contingency table."""
    chi2, _, _, _ = chi2_contingency(contingency)
    n = contingency.sum()
    r, c = contingency.shape
    return float(np.sqrt(chi2 / (n * (min(r, c) - 1))))


# =============================================================================
# Public API
# =============================================================================


def check_confounding(
    df: pd.DataFrame,
    col1: str = "batch_id",
    col2: str = "disease_label",
) -> dict:
    """
    Chi-square test + Cramér's V for association between two categorical
    sample-level columns.

    Parameters
    ----------
    df   : long-format methylation DataFrame
    col1 : first categorical column  (default: 'batch_id')
    col2 : second categorical column (default: 'disease_label')

    Returns
    -------
    dict with keys:
        chi2, p_value, cramers_v, confounded (bool), contingency_table
    """
    meta = df[["sample_id", col1, col2]].drop_duplicates("sample_id")
    ct = pd.crosstab(meta[col1], meta[col2])

    chi2_stat, p_val, dof, _ = chi2_contingency(ct.values)
    v = _cramers_v(ct.values)

    return {
        "chi2": round(chi2_stat, 4),
        "p_value": round(p_val, 6),
        "dof": int(dof),
        "cramers_v": round(v, 4),
        "confounded": v > CRAMERS_V_WARN,
        "contingency_table": ct,
    }


def robust_normalize(df: pd.DataFrame, save_figure: bool = True) -> pd.DataFrame:
    """
    Median-centring normalization: subtract each sample's median beta value.

    A Before/After figure is saved to figures/ when save_figure=True.

    Parameters
    ----------
    df          : long-format methylation DataFrame
    save_figure : if True, save Before/After bar chart to figures/

    Returns
    -------
    pd.DataFrame : copy of df with additional column 'beta_normalized'
                   (original 'beta_value' column is preserved unchanged)
    """
    df = df.copy()
    sample_medians = df.groupby("sample_id")["beta_value"].median()
    df["_sample_median"] = df["sample_id"].map(sample_medians)
    df["beta_normalized"] = (df["beta_value"] - df["_sample_median"]).clip(-1.0, 1.0)
    df = df.drop(columns=["_sample_median"])

    if save_figure:
        _plot_normalization(df)

    return df


def _plot_normalization(df: pd.DataFrame) -> str:
    """Save a Before/After bar chart of per-sample median beta."""
    samples = sorted(df["sample_id"].unique())
    raw_med = df.groupby("sample_id")["beta_value"].median().reindex(samples)
    norm_med = df.groupby("sample_id")["beta_normalized"].median().reindex(samples)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        "Median-Centring Normalization — Before / After",
        fontsize=12,
        fontweight="bold",
    )

    x = range(len(samples))

    for ax, values, ylabel, stage in [
        (axes[0], raw_med.values, "Median Beta", "Before"),
        (axes[1], norm_med.values, "Normalized Median Beta", "After"),
    ]:
        ax.bar(x, values, color="#3498DB", alpha=0.75)
        ax.axhline(0, color="red", linestyle="--", linewidth=0.9)
        ax.set_title(f"{stage} Normalization\nMedian Beta per Sample")
        ax.set_xlabel("Sample Index")
        ax.set_ylabel(ylabel)
        ax.tick_params(axis="x", labelsize=5, rotation=90)
        ax.set_xticks(list(x))
        ax.set_xticklabels(samples, fontsize=4)

    plt.tight_layout()
    out = os.path.join(FIGURES_DIR, "normalization_before_after.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# __main__ — confounding check + normalize mock data
# =============================================================================

if __name__ == "__main__":
    import sys
    from datetime import datetime

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from io_utils import audit_entry, data_path, load_methylation, project_root, ts, write_audit_log  # noqa: E402

    MODULE = "NORMALIZER"
    _now = datetime.now()
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base = project_root()
    _audit_csv = os.path.join(_base, "output", f"audit_log_{MODULE}_{ts_tag}.csv")

    audit_entries = []
    ae = lambda sid, st, d, m: audit_entry(MODULE, sid, st, d, m)

    df = load_methylation(data_path("mock_methylation.csv"))
    n_samples = df["sample_id"].nunique()
    audit_entries.append(ae("cohort", "INFO", "Samples loaded", f"n={n_samples}"))

    # ── Confounding check ──────────────────────────────────────────────────────
    result = check_confounding(df, "batch_id", "disease_label")
    evt_status = "DETECTED" if result["confounded"] else "OK      "
    print(
        f"[{ts()}] [NORMALIZER] {evt_status} | batch_id × disease_label confounding | "
        f"Cramér's V={result['cramers_v']:.4f}  "
        f"p={result['p_value']:.4e}  chi2={result['chi2']:.2f}"
    )
    print("\n  Contingency table:")
    print(result["contingency_table"].to_string())

    audit_entries.append(ae(
        "cohort",
        "DETECTED" if result["confounded"] else "INFO",
        "Batch × disease confound check"
        + (" — strong confound detected" if result["confounded"] else " — OK"),
        f"V={result['cramers_v']:.4f} p={result['p_value']:.4e}",
    ))

    # ── Normalize ──────────────────────────────────────────────────────────────
    df_norm = robust_normalize(df, save_figure=True)
    post_std = df_norm.groupby("sample_id")["beta_normalized"].std()
    print(
        f"\n[{ts()}] [NORMALIZER] DETECTED | Median-centring applied | "
        f"post-norm std range [{post_std.min():.4f}, {post_std.max():.4f}]"
    )
    print(
        f"[{ts()}] [NORMALIZER]           | "
        f"Figure → output/figures/normalization_before_after.png"
    )
    audit_entries.append(ae(
        "cohort", "INFO", "Median-centring normalization applied",
        f"std_range=[{post_std.min():.4f},{post_std.max():.4f}]",
    ))

    write_audit_log(audit_entries, _audit_csv)
    print(f"[{ts()}] [NORMALIZER] Audit log → {_audit_csv}")
