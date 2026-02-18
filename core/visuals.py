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
  3. Beta distribution — contamination, batch effects, and clonal artifacts each
                         leave a distinct fingerprint in the per-sample KDE.
  4. PCA               — batch effects and outlier samples separate cleanly in
                         PC space, making this the fastest sanity check before
                         any corrective step.
"""

import os

import matplotlib

matplotlib.use("Agg")  # headless — no display needed
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "figures")
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
# 3b. PCA — Covariate Panel (batch / label / sex / age)
# =============================================================================


def plot_pca_covariates(
    df: pd.DataFrame,
    save_path: str | None = None,
) -> str:
    """
    4-panel PCA scatter colored by batch, disease label, sex, and age.

    Sex drives PC2 in this cohort (confirmed by post-hoc analysis). Showing
    all four covariates in a single figure makes the confound structure
    immediately visible without requiring separate plots.

    Parameters
    ----------
    df        : long-format methylation DataFrame; must contain columns
                batch_id, disease_label, sex, age
    save_path : optional override for output path

    Returns
    -------
    str : path to saved figure
    """
    pivot = df.pivot_table(index="sample_id", columns="cpg_id", values="beta_value")
    pivot = pivot.fillna(pivot.mean())

    meta = (
        df[["sample_id", "batch_id", "disease_label", "sex", "age"]]
        .drop_duplicates("sample_id")
        .set_index("sample_id")
        .loc[pivot.index]
    )

    X_scaled = StandardScaler().fit_transform(pivot.values)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_scaled)
    var1, var2 = pca.explained_variance_ratio_ * 100

    pca_df = pd.DataFrame(
        {"PC1": coords[:, 0], "PC2": coords[:, 1]}, index=pivot.index
    ).join(meta)

    batch_colors = dict(
        zip(sorted(meta["batch_id"].unique()), sns.color_palette("Set1", meta["batch_id"].nunique()))
    )
    label_colors = {"Case": sns.color_palette("Set1")[0], "Control": sns.color_palette("Set1")[2]}
    sex_colors   = {"M": sns.color_palette("Set1")[1], "F": sns.color_palette("Set1")[4]}

    panels = [
        ("batch_id",      batch_colors, "Batch"),
        ("disease_label", label_colors, "Disease label"),
        ("sex",           sex_colors,   "Sex  (PC2 driver)"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle(
        f"PCA — cohort covariate structure (n={len(pca_df)} samples)",
        fontsize=12,
    )

    for ax, (col, cmap, label) in zip(axes.flat, panels):
        for grp, gdf in pca_df.groupby(col):
            ax.scatter(
                gdf.PC1, gdf.PC2,
                label=grp, color=cmap.get(grp, "#aaaaaa"),
                s=55, edgecolors="k", linewidths=0.4,
            )
            for _, row in gdf.iterrows():
                ax.annotate(
                    row.name, (row.PC1, row.PC2),
                    fontsize=5, ha="left", va="bottom", color="#444",
                )
        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.axhline(0, color="grey", lw=0.5, ls="--")
        ax.set_xlabel(f"PC1 ({var1:.1f}% var)")
        ax.set_ylabel(f"PC2 ({var2:.1f}% var)")
        ax.set_title(label)
        ax.legend(fontsize=8, framealpha=0.7)

    ax = axes.flat[3]
    sc = ax.scatter(
        pca_df.PC1, pca_df.PC2,
        c=pca_df.age, cmap="plasma",
        s=55, edgecolors="k", linewidths=0.4,
    )
    for _, row in pca_df.iterrows():
        ax.annotate(
            row.name, (row.PC1, row.PC2),
            fontsize=5, ha="left", va="bottom", color="#444",
        )
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel(f"PC1 ({var1:.1f}% var)")
    ax.set_ylabel(f"PC2 ({var2:.1f}% var)")
    ax.set_title("Age")
    plt.colorbar(sc, ax=ax, label="age (years)")

    plt.tight_layout()
    out = save_path or os.path.join(FIGURES_DIR, "pca_covariates.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# 4. Exclusion Accounting (Pie + Waterfall)
# =============================================================================


def plot_exclusion_accounting(
    n_total: int,
    steps: list[tuple[str, int]],
    save_path: str | None = None,
) -> str:
    """
    Two-panel figure showing sample attrition through the QC pipeline.

    Left panel  — Pie chart: one slice per exclusion reason + final clean count.
    Right panel — Waterfall chart: step-down bars showing cumulative loss at
                  each filtering stage, with connector lines and drop labels.

    Parameters
    ----------
    n_total : int
        Total samples before any filtering.
    steps : list of (label, n_dropped)
        Ordered list of (stage_label, n_samples_dropped) pairs, e.g.
        [("Bisulfite/Depth QC", 3), ("Contamination", 1), ("Duplicate", 1)].
    save_path : optional override for output path

    Returns
    -------
    str : path to saved figure
    """
    n_clean = n_total - sum(n for _, n in steps)
    drop_colors = ["#E74C3C", "#E67E22", "#9B59B6", "#1ABC9C"]   # up to 4 steps
    clean_color = "#2ECC71"

    # Cumulative remaining after each step (includes "Input" and "Clean")
    remaining = [n_total]
    for _, n in steps:
        remaining.append(remaining[-1] - n)

    fig, (ax_pie, ax_wf) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        f"Sample Exclusion Accounting  —  {n_total} input → {n_clean} clean",
        fontsize=13, fontweight="bold",
    )

    # ── Left: Pie ─────────────────────────────────────────────────────────────
    pie_labels = [f"Clean\n(n={n_clean})"] + [
        f"{lbl}\n(n={n})" for lbl, n in steps
    ]
    pie_sizes = [n_clean] + [n for _, n in steps]
    pie_colors = [clean_color] + drop_colors[: len(steps)]
    explode = [0.05] + [0.0] * len(steps)   # pull out the clean slice

    wedges, texts, autotexts = ax_pie.pie(
        pie_sizes,
        labels=pie_labels,
        colors=pie_colors,
        explode=explode,
        autopct="%1.1f%%",
        startangle=90,
        pctdistance=0.75,
        textprops={"fontsize": 9},
    )
    for at in autotexts:
        at.set_fontsize(8)
    ax_pie.set_title("Cohort Composition\n(% of input samples)", fontsize=10)

    # ── Right: Waterfall ──────────────────────────────────────────────────────
    # Columns: Input (solid) | step drops (floating) | Clean (solid)
    n_cols = len(steps) + 2      # Input + n_steps + Clean
    x_pos = list(range(n_cols))
    x_labels = ["Input"] + [lbl.replace(" ", "\n") for lbl, _ in steps] + ["Clean"]

    # Input bar (solid)
    ax_wf.bar(0, n_total, color="#95A5A6", edgecolor="white", linewidth=0.8, zorder=3)
    ax_wf.text(0, n_total + 0.4, str(n_total), ha="center", va="bottom",
               fontsize=10, fontweight="bold")

    # Floating drop bars for each exclusion step
    for i, (lbl, n_drop) in enumerate(steps, 1):
        n_after = remaining[i]
        color = drop_colors[i - 1]
        # Invisible spacer to float the bar
        ax_wf.bar(i, n_after, color="none", zorder=2)
        ax_wf.bar(i, n_drop, bottom=n_after, color=color,
                  edgecolor="white", linewidth=0.8, zorder=3)
        # Annotate drop
        ax_wf.text(i, n_after + n_drop / 2, f"−{n_drop}",
                   ha="center", va="center", fontsize=9,
                   fontweight="bold", color="white", zorder=4)
        # Remaining count above bar bottom
        ax_wf.text(i, n_after - 0.4, str(n_after), ha="center", va="top",
                   fontsize=9, color="#2C3E50")

        # Connector line from right edge of previous remaining level
        ax_wf.plot([i - 0.45, i - 0.05], [n_after, n_after],
                   color="#BDC3C7", linewidth=1.0, linestyle="--", zorder=1)

    # Clean bar (solid)
    last_i = len(steps) + 1
    ax_wf.bar(last_i, n_clean, color=clean_color,
              edgecolor="white", linewidth=0.8, zorder=3)
    ax_wf.text(last_i, n_clean + 0.4, str(n_clean),
               ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Connector from last drop bar base to clean bar top
    ax_wf.plot([last_i - 0.45, last_i - 0.05], [n_clean, n_clean],
               color="#BDC3C7", linewidth=1.0, linestyle="--", zorder=1)

    ax_wf.set_xticks(x_pos)
    ax_wf.set_xticklabels(x_labels, fontsize=9)
    ax_wf.set_ylabel("Number of Samples")
    ax_wf.set_ylim(0, n_total * 1.12)
    ax_wf.set_title("QC Waterfall — Samples Remaining at Each Stage", fontsize=10)
    ax_wf.spines[["top", "right"]].set_visible(False)
    ax_wf.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)
    ax_wf.set_axisbelow(True)

    plt.tight_layout()
    out = save_path or os.path.join(FIGURES_DIR, "exclusion_accounting.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# =============================================================================
# 5. Volcano Plot (DMR Hunter output)
# =============================================================================


def plot_volcano(
    dmrs: pd.DataFrame,
    color_clonal_risk: bool = False,
    save_path: str | None = None,
) -> str:
    """
    Volcano plot of DMR-hunter sliding-window results.

    X-axis : ΔBeta (Case − Control mean methylation difference)
    Y-axis : −log₁₀(BH-adjusted p-value)

    Reference lines mark the significance threshold (p_adj < 0.05) and the
    minimum biological effect size (|ΔBeta| > 0.10).

    Parameters
    ----------
    dmrs             : DataFrame returned by find_dmrs() — must contain
                       delta_beta, p_adj, significant, clonal_risk columns.
    color_clonal_risk: If True, significant windows overlapping VDJ loci are
                       colored orange and labeled "High Clonality" so the
                       Analyst can distinguish them from clean DMRs.
    save_path        : optional override for output path

    Returns
    -------
    str : path to saved figure
    """
    if dmrs.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "No DMR windows to plot",
                transform=ax.transAxes, ha="center", va="center", fontsize=12)
        out = save_path or os.path.join(
            FIGURES_DIR,
            "volcano_clonal_risk.png" if color_clonal_risk else "volcano.png",
        )
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out

    df = dmrs.copy()
    df["neg_log10_padj"] = -np.log10(np.clip(df["p_adj"], 1e-300, 1.0))

    has_clonal = "clonal_risk" in df.columns

    # ── Classify windows ──────────────────────────────────────────────────────
    non_sig = df[~df["significant"]]
    if color_clonal_risk and has_clonal:
        sig_clean = df[df["significant"] & ~df["clonal_risk"]]
        sig_clonal = df[df["significant"] & df["clonal_risk"]]
    else:
        sig_clean = df[df["significant"]]
        sig_clonal = df.iloc[0:0]   # empty

    n_sig = int(df["significant"].sum())
    n_sig_clonal = len(sig_clonal)

    fig, ax = plt.subplots(figsize=(10, 6))

    suffix = " (colored by VDJ clonal risk)" if color_clonal_risk else ""
    title = (
        f"Volcano Plot — DMR Hunter{suffix}\n"
        f"n_windows={len(df)}  n_significant={n_sig}"
        + (f"  n_high_clonality={n_sig_clonal}" if color_clonal_risk else "")
    )
    ax.set_title(title, fontsize=11)

    # Non-significant (grey)
    ax.scatter(
        non_sig["delta_beta"], non_sig["neg_log10_padj"],
        color="#BDC3C7", alpha=0.5, s=18, linewidths=0, label="Not significant",
    )
    # Significant, clean (red)
    if len(sig_clean):
        ax.scatter(
            sig_clean["delta_beta"], sig_clean["neg_log10_padj"],
            color="#E74C3C", alpha=0.85, s=40, linewidths=0.5,
            edgecolors="white", label=f"Significant (n={len(sig_clean)})",
            zorder=3,
        )
    # Significant, clonal risk (orange) — only when color_clonal_risk=True
    if len(sig_clonal):
        ax.scatter(
            sig_clonal["delta_beta"], sig_clonal["neg_log10_padj"],
            color="#E67E22", alpha=0.9, s=55, linewidths=0.6,
            edgecolors="white", marker="D",
            label=f"High Clonality (n={n_sig_clonal})",
            zorder=4,
        )

    # Reference lines
    p_thresh = 0.05
    db_thresh = 0.10
    y_thresh = -np.log10(p_thresh)
    ax.axhline(y_thresh, color="#E74C3C", linestyle="--", linewidth=0.9,
               alpha=0.7, label=f"p_adj = {p_thresh}")
    ax.axvline(db_thresh, color="#2C3E50", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.axvline(-db_thresh, color="#2C3E50", linestyle=":", linewidth=0.8, alpha=0.6,
               label=f"|ΔBeta| = {db_thresh}")
    ax.axvline(0, color="#95A5A6", linewidth=0.5, alpha=0.5)

    ax.set_xlabel("ΔBeta (Case − Control)", fontsize=11)
    ax.set_ylabel("−log₁₀(p_adj)", fontsize=11)
    ax.legend(fontsize=8, framealpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    fname = "volcano_clonal_risk.png" if color_clonal_risk else "volcano.png"
    out = save_path or os.path.join(FIGURES_DIR, fname)
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
    print(f"[{ts()}] [VISUALS] All EDA figures written to output/figures/")
    print(f"[{ts()}] [VISUALS] Audit log → {_audit_csv}")
