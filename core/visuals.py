"""
core/visuals.py — ImmuneMethylTools EDA Suite
==============================================
Exploratory data analysis visualisations for methylation QC pipelines.
All functions write to figures/ and return the saved path.

Biological intent
-----------------
Before any statistical analysis, a QC-aware bioinformatician must visually
inspect three distributions that betray the most common wet-lab failures:

  1. Conversion rate  — incomplete bisulfite conversion inflates non-CpG
                        methylation; samples with rate < 99 % are suspect.
  2. Sequencing depth — low-coverage sites have inflated binomial variance;
                        sites with depth < 10 reads are statistically unreliable.
  3. Beta distribution — contamination, batch effects, and clonal artefacts each
                         leave a distinct fingerprint in the per-sample KDE.
  4. PCA               — batch effects and outlier samples separate cleanly in
                         PC space, making this the fastest sanity check before
                         any corrective step.
"""

import os

import matplotlib

matplotlib.use("Agg")  # headless — no display needed
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)


# =============================================================================
# 1. QC Metrics Dashboard
# =============================================================================


def plot_qc_metrics(df: pd.DataFrame, save_path: str | None = None) -> str:
    """
    Three-panel QC dashboard: Conversion Rate, Mean Depth, Global Beta.

    Parameters
    ----------
    df        : long-format methylation DataFrame (mock_methylation schema)
    save_path : optional override for output path

    Returns
    -------
    str : path to saved figure
    """
    sample_stats = (
        df.groupby("sample_id")
        .agg(
            conversion_rate=("non_cpg_meth_rate", lambda x: 1.0 - x.mean()),
            mean_depth=("depth", "mean"),
            mean_beta=("beta_value", "mean"),
        )
        .reset_index()
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle(
        "QC Metrics — Per-Sample Distributions", fontsize=13, fontweight="bold"
    )

    # Panel 1 — Bisulfite Conversion Rate
    axes[0].hist(
        sample_stats["conversion_rate"], bins=20, color="#2ECC71", edgecolor="white"
    )
    axes[0].axvline(
        0.99, color="red", linestyle="--", linewidth=1.2, label="99 % threshold"
    )
    axes[0].set_title("Bisulfite Conversion Rate\n(1 – non-CpG meth rate)")
    axes[0].set_xlabel("Conversion Rate")
    axes[0].set_ylabel("# Samples")
    axes[0].legend(fontsize=8)

    # Panel 2 — Mean Sequencing Depth
    axes[1].hist(
        sample_stats["mean_depth"], bins=20,
        range=(0, max(sample_stats["mean_depth"]) * 1.1),
        color="#3498DB", edgecolor="white"
    )
    axes[1].axvline(
        10, color="red", linestyle="--", linewidth=1.2, label="Depth ≥ 10"
    )
    axes[1].set_xlim(left=0)
    axes[1].set_title("Mean Sequencing Depth\nper Sample")
    axes[1].set_xlabel("Mean Read Depth")
    axes[1].set_ylabel("# Samples")
    axes[1].legend(fontsize=8)

    # Panel 3 — Global Mean Beta
    axes[2].hist(
        sample_stats["mean_beta"], bins=20, color="#9B59B6", edgecolor="white"
    )
    axes[2].set_title("Global Mean Beta\nper Sample")
    axes[2].set_xlabel("Mean Beta Value")
    axes[2].set_ylabel("# Samples")

    plt.tight_layout()
    out = save_path or os.path.join(FIGURES_DIR, "qc_metrics.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# 2. Beta Distribution KDE
# =============================================================================


def plot_beta_distribution(df: pd.DataFrame, save_path: str | None = None) -> str:
    """
    Kernel Density Estimate (KDE) of beta values, one curve per sample.

    Contaminated samples ('muddy') appear as a broad hump centred near 0.5
    rather than the clean bimodal peaks at ~0.10 and ~0.85.

    Parameters
    ----------
    df        : long-format methylation DataFrame
    save_path : optional override for output path

    Returns
    -------
    str : path to saved figure
    """
    samples = sorted(df["sample_id"].unique())
    n = len(samples)
    palette = sns.color_palette("husl", n)

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_title(
        "Beta Value Distribution — KDE per Sample\n"
        "(muddy/contaminated samples show peaks shifting toward 0.5)",
        fontsize=11,
    )

    for sid, color in zip(samples, palette):
        betas = df.loc[df["sample_id"] == sid, "beta_value"].dropna()
        betas.plot.kde(ax=ax, color=color, alpha=0.40, linewidth=0.9, label=sid)

    ax.set_xlabel("Beta Value")
    ax.set_ylabel("Density")
    ax.set_xlim(-0.05, 1.05)

    if n > 15:
        if ax.get_legend():
            ax.get_legend().remove()
        ax.text(
            0.99,
            0.97,
            f"n={n} samples (legend omitted for clarity)",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            color="grey",
        )
    else:
        ax.legend(fontsize=6, ncol=3)

    plt.tight_layout()
    out = save_path or os.path.join(FIGURES_DIR, "beta_distribution_kde.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# 3. PCA Scatter
# =============================================================================


def plot_pca(
    df: pd.DataFrame,
    title: str = "PCA — Sample Clustering",
    color_by: str = "disease_label",
    save_path: str | None = None,
) -> str:
    """
    PCA scatter plot on per-sample beta profiles.

    Batch effects and contamination separate cleanly in PC space, making this
    the fastest sanity check before any corrective step.

    Parameters
    ----------
    df        : long-format methylation DataFrame
    title     : plot title
    color_by  : sample-level metadata column to color points by
                (e.g. 'disease_label', 'batch_id')
    save_path : optional override for output path

    Returns
    -------
    str : path to saved figure
    """
    # Build sample × CpG matrix
    pivot = df.pivot_table(index="sample_id", columns="cpg_id", values="beta_value")
    pivot = pivot.fillna(pivot.mean())

    # Align metadata
    meta = (
        df[["sample_id", color_by]]
        .drop_duplicates("sample_id")
        .set_index("sample_id")
    )
    meta = meta.loc[pivot.index]

    # PCA
    X_scaled = StandardScaler().fit_transform(pivot.values)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_scaled)
    var_exp = pca.explained_variance_ratio_

    # Plot
    categories = sorted(meta[color_by].unique())
    palette = dict(zip(categories, sns.color_palette("Set1", len(categories))))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_title(f"{title}\nColored by: {color_by}", fontsize=11)

    for cat in categories:
        mask = meta[color_by].values == cat
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            label=cat,
            color=palette[cat],
            s=60,
            alpha=0.75,
            edgecolors="white",
            linewidths=0.4,
        )

    for i, sid in enumerate(pivot.index):
        ax.annotate(
            sid,
            (coords[i, 0], coords[i, 1]),
            fontsize=4.5,
            alpha=0.6,
            ha="center",
            va="bottom",
        )

    ax.set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)")
    ax.legend(title=color_by, fontsize=8, title_fontsize=8)
    ax.axhline(0, color="grey", linewidth=0.4, linestyle="--")
    ax.axvline(0, color="grey", linewidth=0.4, linestyle="--")

    plt.tight_layout()
    safe_name = color_by.replace(" ", "_")
    out = save_path or os.path.join(FIGURES_DIR, f"pca_{safe_name}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# __main__ — generate all EDA plots on mock data
# =============================================================================

if __name__ == "__main__":
    import sys
    from datetime import datetime

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from io_utils import data_path, load_methylation, project_root, write_audit_log  # noqa: E402

    MODULE = "VISUALS"
    _now   = datetime.now()
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base  = project_root()
    _audit_csv = os.path.join(_base, "data", f"audit_log_{MODULE}_{ts_tag}.csv")

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

    df = load_methylation(data_path("mock_methylation.csv"))
    n_samples = df["sample_id"].nunique()
    print(f"[{ts()}] [VISUALS] Shape: {df.shape} | Samples: {n_samples}")
    audit_entries.append(ae("cohort", "INFO", "Dataset loaded", f"n_samples={n_samples}"))

    p1 = plot_qc_metrics(df)
    print(f"[{ts()}] [VISUALS] DETECTED | QC metrics dashboard saved | {p1}")
    audit_entries.append(ae("cohort", "INFO", "QC metrics dashboard saved", p1))

    p2 = plot_beta_distribution(df)
    print(f"[{ts()}] [VISUALS] DETECTED | Beta KDE saved | {p2}")
    audit_entries.append(ae("cohort", "INFO", "Beta distribution KDE saved", p2))

    p3 = plot_pca(df, title="PCA — Disease Label", color_by="disease_label")
    print(f"[{ts()}] [VISUALS] DETECTED | PCA (disease_label) saved | {p3}")
    audit_entries.append(ae("cohort", "INFO", "PCA colored by disease_label saved", p3))

    p4 = plot_pca(df, title="PCA — Batch", color_by="batch_id")
    print(f"[{ts()}] [VISUALS] DETECTED | PCA (batch_id) saved | {p4}")
    audit_entries.append(ae("cohort", "INFO", "PCA colored by batch_id saved", p4))

    write_audit_log(audit_entries, _audit_csv)
    print(f"[{ts()}] [VISUALS] All EDA figures written to figures/")
    print(f"[{ts()}] [VISUALS] Audit log → {_audit_csv}")
