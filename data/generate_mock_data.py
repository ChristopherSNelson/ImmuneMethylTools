"""
generate_mock_data.py — ImmuneMethylTools Lab Simulator
========================================================
Generates mock_methylation.csv with seven embedded "stumper" artifacts that
mimic real-world pitfalls in immune-cell WGBS / RRBS analysis:

  Artifact 1 — Confounded Batch:   Batch_01 enriched for Cases (+0.1 beta shift)
  Artifact 2 — Clonal Artifact:    VDJ locus, beta > 0.8, fragment > 180 bp
  Artifact 3 — Bisulfite Failure:  2 samples with non_cpg_meth_rate > 0.02
  Artifact 4 — Sample Duplication: 2 samples with Pearson r > 0.99
  Artifact 5 — Contamination:      1 sample with muddy beta (peak near 0.5)
  Artifact 6 — Low Coverage:       S030 depth forced to Poisson(λ=5), mean ~5x
  Artifact 7 — Sex Metadata Mixup: S035 (true F) reported M; S036 (true M) reported F

In addition to the stumper artifacts, one genuine biological signal is injected:
  True DMR      — cg00000300–cg00000310: all Case samples +0.25 beta (autosomal; batch-independent)
  This gives dmr_hunter a real, non-artifactual DMR to detect after batch correction.

Outputs
-------
  data/mock_methylation.csv
  output/figures/qc_before_after.png   — Before/After visualisation of each artifact
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless — no display needed
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy.stats import pearsonr  # noqa: E402

# ── Reproducibility ──────────────────────────────────────────────────────────
RNG = np.random.default_rng(seed=42)

# ── Study parameters ─────────────────────────────────────────────────────────
N_PATIENTS = 40          # 20 Case, 20 Control
N_CPGS = 500         # CpG sites per sample
N_X_CPGS = 30        # last N_X_CPGS of N_CPGS are X-linked (cg00000471–cg00000500)
N_BATCHES = 3
CASE_LABEL = "Case"
CTRL_LABEL = "Control"

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "figures")
OUT_CSV = os.path.join(os.path.dirname(__file__), "mock_methylation.csv")
os.makedirs(FIGURES_DIR, exist_ok=True)


# =============================================================================
# 1.  PATIENT / SAMPLE MANIFEST
# =============================================================================

def build_manifest() -> pd.DataFrame:
    """
    Create one row per (sample, CpG).  Assigns batch with confounded Case/Batch_01
    enrichment (Artifact 1 setup).
    """
    patient_ids = [f"P{i:03d}" for i in range(1, N_PATIENTS + 1)]
    disease_lbls = ([CASE_LABEL] * 20) + ([CTRL_LABEL] * 20)
    ages = RNG.integers(25, 70, size=N_PATIENTS).tolist()

    # Artifact 1 — Batch_01 gets 80 % of Case patients
    case_patients = [p for p, d in zip(patient_ids, disease_lbls) if d == CASE_LABEL]
    ctrl_patients = [p for p, d in zip(patient_ids, disease_lbls) if d == CTRL_LABEL]

    n_case_b1 = int(0.80 * len(case_patients))      # 16 / 20
    n_ctrl_b1 = int(0.20 * len(ctrl_patients))      # 4  / 20

    # Remaining 4 Case patients split 2/2 across Batch_02 and Batch_03
    # so every batch has at least some Cases — required for batch correction models.
    batch_map = {}
    for p in case_patients[:n_case_b1]:
        batch_map[p] = "Batch_01"
    for p in case_patients[n_case_b1:n_case_b1 + 2]:
        batch_map[p] = "Batch_02"
    for p in case_patients[n_case_b1 + 2:]:
        batch_map[p] = "Batch_03"
    for p in ctrl_patients[:n_ctrl_b1]:
        batch_map[p] = "Batch_01"
    for p in ctrl_patients[n_ctrl_b1:16]:
        batch_map[p] = "Batch_02"
    for p in ctrl_patients[16:]:
        batch_map[p] = "Batch_03"

    age_map = dict(zip(patient_ids, ages))
    disease_map = dict(zip(patient_ids, disease_lbls))

    # Odd patient number → Female; even → Male
    sex_map = {
        f"P{i:03d}": ("F" if i % 2 == 1 else "M")
        for i in range(1, N_PATIENTS + 1)
    }

    rows = []
    sample_counter = 1
    for pid in patient_ids:
        sid = f"S{sample_counter:03d}"
        sample_counter += 1
        for cg_idx in range(1, N_CPGS + 1):
            is_x = cg_idx > N_CPGS - N_X_CPGS   # True for cg_idx 471–500
            rows.append({
                "sample_id": sid,
                "patient_id": pid,
                "batch_id": batch_map[pid],
                "age": age_map[pid],
                "disease_label": disease_map[pid],
                "cpg_id": f"cg{cg_idx:08d}",
                "sex": sex_map[pid],
                "is_x_chromosome": is_x,
            })

    return pd.DataFrame(rows)


# =============================================================================
# 2.  BASELINE METHYLATION VALUES
# =============================================================================

def add_baseline_methylation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Draw beta values from a bimodal distribution typical of CpG methylation:
      ~60 % fully methylated (beta ~ 0.85),  ~40 % fully unmethylated (beta ~ 0.10).
    Small Gaussian noise added per CpG site.

    Bisulfite intuition notes:
    - Fragment length baseline: Normal(150, 12).  P(>180 bp) ≈ 0.6 %,  making
      the VDJ clonal artifact (>200 bp) a genuine outlier, not background noise.
    - VDJ-region baseline beta: set to low/unmethylated (~0.10) after the
      genome-wide draw.  Biologically, active VDJ loci in B/T cells must remain
      accessible (unmethylated) for recombination; hypermethylation there only
      arises AFTER a dominant clone silences the locus — i.e., the artifact.
    """
    n = len(df)
    # Site-level random component (same shape across samples for a given CpG)
    cpg_means = RNG.choice([0.10, 0.85], size=N_CPGS, p=[0.40, 0.60])
    cpg_mean_series = np.tile(cpg_means, N_PATIENTS)  # repeated for each sample

    noise = RNG.normal(0, 0.06, size=n)
    beta = np.clip(cpg_mean_series + noise, 0.0, 1.0)
    df["beta_value"] = beta

    # Sequencing depth: Negative Binomial ~ mean 30, overdispersion
    df["depth"] = RNG.negative_binomial(n=5, p=5 / (5 + 25), size=n).clip(1, None)

    # Fragment length: tight Normal(150, 12) — >180 bp is ~0.6 % of baseline.
    # This ensures the clonal VDJ signal (injected at 200–260 bp) is a clear outlier.
    df["fragment_length"] = RNG.normal(150, 12, size=n).astype(int).clip(80, 220)

    # VDJ region: ~3 % of CpGs flagged as VDJ loci
    vdj_cpgs = set(RNG.choice(range(1, N_CPGS + 1), size=int(0.03 * N_CPGS), replace=False))
    df["is_vdj_region"] = df["cpg_id"].apply(
        lambda c: int(c.lstrip("cg")) in vdj_cpgs
    )

    # Bisulfite intuition: VDJ baseline must be unmethylated.
    # Overwrite the bimodal draw for VDJ CpGs with a low-methylation distribution.
    vdj_mask = df["is_vdj_region"]
    df.loc[vdj_mask, "beta_value"] = RNG.normal(0.10, 0.05, size=vdj_mask.sum()).clip(0.0, 0.35)

    # Non-CpG methylation rate (bisulfite conversion proxy): ~ N(0.004, 0.001)
    df["non_cpg_meth_rate"] = RNG.normal(0.004, 0.001, size=n).clip(0, 1)

    return df


# =============================================================================
# ARTIFACTS
# =============================================================================

def inject_artifact1_confounded_batch(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 1 — Confounded Batch
    Batch_01 Case samples receive a systematic +0.10 beta shift.
    This mimics a plate / reagent batch effect that is collinear with disease.
    """
    mask = (df["batch_id"] == "Batch_01") & (df["disease_label"] == CASE_LABEL)
    df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + 0.10).clip(0, 1)
    print(f"  [Artifact 1] +0.1 shift applied to {mask.sum()} rows "
          f"({df.loc[mask, 'sample_id'].nunique()} samples in Batch_01/Case)")
    return df


def inject_artifact2_clonal_vdj(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 2 — Clonal Expansion Artifact in VDJ Locus
    Pick one Case patient; for all their VDJ-region CpGs set beta > 0.8 and
    fragment_length > 180 bp.  This mimics a dominant clone in which VDJ loci
    are hypermethylated and the original fragment is long (compact chromatin).
    """
    case_patients = df.loc[df["disease_label"] == CASE_LABEL, "patient_id"].unique()
    # Skip index 0 (P001 → S001) and index 1 (P002 → S002): both carry the
    # bisulfite-failure artifact (inject_artifact3_bisulfite_failure targets
    # all_samples[:2]).  Selecting index 2 ensures the clonal patient passes
    # sample-level QC so the VDJ masking stage can be demonstrated.
    clonal_patient = case_patients[2]  # deterministic; avoids bisulfite-failure overlap

    mask = (df["patient_id"] == clonal_patient) & (df["is_vdj_region"])
    n_affected = mask.sum()

    df.loc[mask, "beta_value"] = RNG.uniform(0.82, 0.97, size=n_affected)
    df.loc[mask, "fragment_length"] = RNG.integers(200, 260, size=n_affected)

    print(f"  [Artifact 2] Clonal VDJ artifact injected into patient {clonal_patient}: "
          f"{n_affected} CpG rows (beta > 0.8, fragment > 180 bp)")
    return df


def inject_artifact3_bisulfite_failure(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 3 — Incomplete Bisulfite Conversion
    Two samples receive non_cpg_meth_rate drawn from N(0.05, 0.01) — well above
    the 2 % failure threshold.  High non-CpG methylation = cytosines not converted
    = inflated beta values genome-wide.
    """
    all_samples = df["sample_id"].unique()
    # Choose 2 samples that are NOT already the clonal patient's sample
    bad_samples = all_samples[:2]  # deterministic; S001, S002 (both Case/Batch_01)

    for sid in bad_samples:
        mask = df["sample_id"] == sid
        df.loc[mask, "non_cpg_meth_rate"] = RNG.normal(0.05, 0.008, size=mask.sum()).clip(0.02, 1)
        # Bisulfite failure also artifactually inflates beta
        df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + RNG.normal(0.08, 0.02, size=mask.sum())).clip(0, 1)

    print(f"  [Artifact 3] Bisulfite failure injected into samples: {list(bad_samples)}")
    return df


def inject_artifact4_sample_duplication(df: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 4 — Technical Duplicate
    Clone sample S010's beta values into a new sample S_DUP with tiny noise so
    Pearson r > 0.99.  S_DUP gets a unique patient_id (P_DUP) so it does not
    inflate the CpG count attributed to S010's patient.
    """
    source_sid = "S010"
    dup_sid = "S_DUP"
    dup_pid = "P_DUP"

    source_rows = df[df["sample_id"] == source_sid].copy()
    source_rows["sample_id"] = dup_sid
    source_rows["patient_id"] = dup_pid
    # Add tiny noise (std = 0.002) to keep r > 0.99
    source_rows["beta_value"] = (
        source_rows["beta_value"] + RNG.normal(0, 0.002, size=len(source_rows))
    ).clip(0, 1)

    df = pd.concat([df, source_rows], ignore_index=True)

    # Verify
    orig = df.loc[df["sample_id"] == source_sid, "beta_value"].values
    clone = df.loc[df["sample_id"] == dup_sid, "beta_value"].values
    r, _ = pearsonr(orig, clone)
    print(f"  [Artifact 4] Duplicate {dup_sid} (patient {dup_pid}) cloned from {source_sid}: r = {r:.4f}")
    return df


def inject_artifact5_contamination(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 5 — Sample Contamination
    Pick one sample and collapse its bimodal beta distribution toward 0.5.
    Real-world: cross-contamination from another sample type smears the signal.
    We use a mixture: 50 % original + 50 % uniform(0.3, 0.7) noise.
    """
    contaminated_sid = "S020"
    mask = df["sample_id"] == contaminated_sid
    n = mask.sum()

    original_beta = df.loc[mask, "beta_value"].values
    contaminant = RNG.uniform(0.35, 0.65, size=n)      # peaks near 0.5
    mixed_beta = 0.50 * original_beta + 0.50 * contaminant
    df.loc[mask, "beta_value"] = mixed_beta.clip(0, 1)

    print(f"  [Artifact 5] Contamination injected into {contaminated_sid}: "
          f"mean beta shifted to {mixed_beta.mean():.3f}")
    return df


def inject_artifact6_low_depth(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 6 — Low Coverage (depth failure)
    Force sample S030's read depth to Poisson(λ=5), yielding a mean of ~5x.
    Sites with mean depth < 10 reads have inflated binomial sampling variance;
    beta values from such sites are statistically unreliable and must be
    excluded by the depth filter in qc_guard.audit_quality().
    """
    mask = df["sample_id"] == "S030"
    df.loc[mask, "depth"] = RNG.poisson(lam=5, size=mask.sum())
    mean_depth = df.loc[mask, "depth"].mean()
    print(f"  [Artifact 6] Low-coverage failure injected into S030: "
          f"mean depth = {mean_depth:.1f}x (threshold: 10x)")
    return df


def inject_true_biological_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Inject a legitimate DMR that is NOT a batch effect or artifact.
    We pick a 10-CpG window and shift ALL Cases (regardless of batch)
    by +0.25. This creates a tight, significant cluster.

    Chosen region (cg00000300–cg00000310) is far from the proxy markers
    (FoxP3/PAX5/VDJ) and is autosomal (not X-linked), so it is unaffected
    by the XCI signal re-injection.

    Purpose: gives dmr_hunter at least one genuinely significant, non-artifactual
    DMR to find after batch correction, demonstrating the full detection pipeline.
    """
    target_cpgs = [f"cg{i:08d}" for i in range(300, 311)]
    mask = (df["disease_label"] == "Case") & (df["cpg_id"].isin(target_cpgs))
    df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + 0.25).clip(0, 1)
    print(f"  [Signal Spike] Injected true biological signal into {len(target_cpgs)} CpGs "
          f"({mask.sum()} Case rows; +0.25 shift).")
    return df


def inject_borderline_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Inject a borderline sub-threshold biological signal at cg00000150-cg00000157
    (8 CpGs). Raw shift +0.09 to Case beta values; after median-centring the
    observed ΔBeta should be ~0.08, just below the DELTA_BETA_MIN = 0.10 cutoff.

    Purpose: negative control confirming the pipeline does not over-call near
    the detection boundary.
    """
    target_cpgs = [f"cg{i:08d}" for i in range(150, 158)]
    mask = (df["disease_label"] == "Case") & (df["cpg_id"].isin(target_cpgs))
    df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + 0.09).clip(0, 1)
    print(f"  [Borderline Signal] Injected borderline signal into {len(target_cpgs)} CpGs "
          f"({mask.sum()} Case rows; +0.09 shift, expected ΔBeta ~0.08).")
    return df


def inject_subtle_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Inject a clearly sub-threshold biological signal at cg00000200-cg00000205
    (6 CpGs). Raw shift +0.08 to Case beta values; after median-centring the
    observed ΔBeta should be ~0.04, well below the DELTA_BETA_MIN = 0.10 cutoff.

    Purpose: negative control confirming the pipeline ignores weak signals.
    """
    target_cpgs = [f"cg{i:08d}" for i in range(200, 206)]
    mask = (df["disease_label"] == "Case") & (df["cpg_id"].isin(target_cpgs))
    df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + 0.08).clip(0, 1)
    print(f"  [Subtle Signal] Injected subtle signal into {len(target_cpgs)} CpGs "
          f"({mask.sum()} Case rows; +0.08 shift, expected ΔBeta ~0.04).")
    return df


def inject_xci_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Re-inject XCI-appropriate methylation signal for all X-linked CpGs.

    Called AFTER artifacts 1-6 so that batch shifts, bisulfite failures, and
    contamination do not corrupt the ground-truth X-linked beta used by the
    XCI guard.  Artifact 7 then swaps sex metadata for S035/S036 without
    touching beta values.

    Female (XX): XCI → mean ~0.50 at X-linked sites (one X active, one silenced)
    Male   (XY): single active X → unmethylated baseline, mean ~0.25
    """
    female_x_mask = df["is_x_chromosome"].astype(bool) & (df["sex"] == "F")
    df.loc[female_x_mask, "beta_value"] = (
        RNG.normal(0.50, 0.05, size=int(female_x_mask.sum())).clip(0.35, 0.65)
    )
    male_x_mask = df["is_x_chromosome"].astype(bool) & (df["sex"] == "M")
    df.loc[male_x_mask, "beta_value"] = (
        RNG.normal(0.25, 0.04, size=int(male_x_mask.sum())).clip(0.10, 0.33)
    )
    n_f = int(female_x_mask.sum())
    n_m = int(male_x_mask.sum())
    print(f"  [XCI Signal] X-linked betas re-injected: "
          f"{n_f} female rows (~0.50) | {n_m} male rows (~0.25)")
    return df


def inject_artifact7_sex_mixup(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 7 — Sex Metadata Mixup

    S035 (truly female, P035): X-linked beta was injected as female (~0.50,
      XCI signal), but the 'sex' metadata column is overwritten to "M".
      Reported sex contradicts the observed X-linked methylation signal.

    S036 (truly male, P036): X-linked beta was injected as male (~0.25),
      but the 'sex' metadata column is overwritten to "F".
      Reported sex contradicts the observed X-linked methylation signal.

    detect_sex_mixups() (core/xci_guard.py) will flag both samples for exclusion.
    """
    df.loc[df["sample_id"] == "S035", "sex"] = "M"   # true F, reported M
    df.loc[df["sample_id"] == "S036", "sex"] = "F"   # true M, reported F
    n_s035 = (df["sample_id"] == "S035").sum()
    n_s036 = (df["sample_id"] == "S036").sum()
    print(f"  [Artifact 7] Sex metadata mixup injected: "
          f"S035 ({n_s035} rows, true F → reported M), "
          f"S036 ({n_s036} rows, true M → reported F)")
    return df


# =============================================================================
# 3.  BEFORE / AFTER VISUALISATION
# =============================================================================

def plot_before_after(df_before: pd.DataFrame, df_after: pd.DataFrame) -> None:
    """
    Six-panel figure showing the fingerprint of each artifact before and after
    injection so developers can visually confirm the simulation worked.
    """
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle(
        "ImmuneMethylTools — Mock Data: Before vs. After Artifact Injection",
        fontsize=14, fontweight="bold", y=0.98
    )
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.4)

    palette = {"Case": "#E74C3C", "Control": "#3498DB"}

    # ── Panel A: Batch × Disease (Artifact 1) ──────────────────────────────
    ax_a1 = fig.add_subplot(gs[0, 0])
    ax_a2 = fig.add_subplot(gs[0, 1])
    panels_a = [(ax_a1, df_before, "Before Artifact Injection"), (ax_a2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_a:
        sample_mean = (df_.groupby(["sample_id", "batch_id", "disease_label"])
                       ["beta_value"].mean().reset_index())
        sns.boxplot(data=sample_mean, x="batch_id", y="beta_value",
                    hue="disease_label", palette=palette, ax=ax,
                    linewidth=0.8, fliersize=2)
        # Annotate n per (batch, disease) group inside the top of each box column
        counts = sample_mean.groupby(["batch_id", "disease_label"]).size()
        batches = sorted(sample_mean["batch_id"].unique())
        hue_ord = sorted(sample_mean["disease_label"].unique())
        n_hue = len(hue_ord)
        bw = 0.8 / n_hue
        offsets = [(i - (n_hue - 1) / 2) * bw for i in range(n_hue)]
        ylo, yhi = ax.get_ylim()
        y_top = yhi - (yhi - ylo) * 0.02   # 2% below the top edge
        for bi, batch in enumerate(batches):
            for di, disease in enumerate(hue_ord):
                n = counts.get((batch, disease), 0)
                ax.text(bi + offsets[di], y_top, f"n={n}",
                        ha="center", va="top", fontsize=6)
        ax.set_ylim(ylo, yhi)   # lock ylim so text doesn't trigger autoscale
        ax.set_title(f"{title}: Batch × Disease\nmean beta", fontsize=9)
        ax.set_xlabel("Batch", fontsize=8)
        ax.set_ylabel("Mean beta", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, title_fontsize=6)

    # ── Panel B: VDJ Fragment Length (Artifact 2) ──────────────────────────
    ax_b1 = fig.add_subplot(gs[0, 2])
    ax_b2 = fig.add_subplot(gs[0, 3])
    panels_b = [(ax_b1, df_before, "Before Artifact Injection"), (ax_b2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_b:
        vdj_data = df_[df_["is_vdj_region"]]
        sns.scatterplot(data=vdj_data, x="fragment_length", y="beta_value",
                        hue="disease_label", palette=palette, ax=ax,
                        s=8, alpha=0.5, linewidth=0)
        ax.axvline(180, color="k", linestyle="--", linewidth=0.8, label="180 bp")
        ax.axhline(0.8, color="gray", linestyle="--", linewidth=0.8, label="β=0.8")
        ax.set_title(f"{title}: VDJ Loci\nFragment vs Beta", fontsize=9)
        ax.set_xlabel("Fragment length (bp)", fontsize=8)
        ax.set_ylabel("Beta value", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=5)

    # ── Panel C: Non-CpG Meth Rate (Artifact 3) ────────────────────────────
    ax_c1 = fig.add_subplot(gs[1, 0])
    ax_c2 = fig.add_subplot(gs[1, 1])
    panels_c = [(ax_c1, df_before, "Before Artifact Injection"), (ax_c2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_c:
        sample_ncpg = df_.groupby("sample_id")["non_cpg_meth_rate"].mean().reset_index()
        ax.hist(sample_ncpg["non_cpg_meth_rate"], bins=30, range=(0, 0.07), color="#2ECC71", edgecolor="white")
        ax.axvline(0.02, color="red", linestyle="--", linewidth=1.2, label="2% threshold")
        ax.set_title(f"{title}: Non-CpG Meth Rate\n(bisulfite QC)", fontsize=9)
        ax.set_xlabel("Non-CpG meth rate", fontsize=8)
        ax.set_ylabel("# Samples", fontsize=8)
        ax.set_xlim(0, 0.07)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)

    # ── Panel D: Sample Correlation Heatmap (Artifact 4) ───────────────────
    ax_d1 = fig.add_subplot(gs[1, 2])
    ax_d2 = fig.add_subplot(gs[1, 3])
    panels_d = [(ax_d1, df_before, "Before Artifact Injection"), (ax_d2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_d:
        pivot = df_.pivot_table(index="cpg_id", columns="sample_id", values="beta_value")
        # Subset to a manageable number of samples for visibility
        cols = sorted(pivot.columns)[:12]
        if "S_DUP" in pivot.columns:
            cols = cols[:11] + ["S_DUP"]
        corr = pivot[cols].corr()
        sns.heatmap(corr, ax=ax, cmap="RdYlBu_r", vmin=0.7, vmax=1.0,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"shrink": 0.7})
        subtitle = (
            "Pearson r, subset incl. duplicate pair" if "S_DUP" in pivot.columns
            else "Pearson r, 12-sample subset"
        )
        ax.set_title(f"{title}: Sample Correlation\n({subtitle})", fontsize=9)
        ax.tick_params(labelsize=5, rotation=45)

    # ── Panel E: Beta Distribution (Artifact 5) ────────────────────────────
    ax_e1 = fig.add_subplot(gs[2, 0:2])
    ax_e2 = fig.add_subplot(gs[2, 2:4])
    target_sid = "S020"
    ref_sid = "S005"
    panels_e = [(ax_e1, df_before, "Before Artifact Injection"), (ax_e2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_e:
        for sid, color, label in [
            (target_sid, "#E74C3C", f"{target_sid} (contaminated)"),
            (ref_sid, "#3498DB", f"{ref_sid} (reference)"),
        ]:
            sub = df_[df_["sample_id"] == sid]["beta_value"]
            if len(sub):
                ax.hist(sub, bins=50, alpha=0.55, color=color, label=label,
                        edgecolor="none", density=True)
        ax.set_title(f"{title}: Beta Distribution\n(contamination check)", fontsize=9)
        ax.set_xlabel("Beta value", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)

    out_path = os.path.join(FIGURES_DIR, "qc_before_after.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  [Visualisation] Saved: {out_path}")


# =============================================================================
# 4.  MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("ImmuneMethylTools — Mock Data Generator")
    print("=" * 60)

    # ── Build manifest & baseline data ──────────────────────────────────────
    print("\n[1] Building patient manifest...")
    manifest = build_manifest()
    print(f"    Rows: {len(manifest):,}  |  Samples: {manifest['sample_id'].nunique()}"
          f"  |  CpGs per sample: {N_CPGS}")

    print("[2] Sampling baseline methylation...")
    df = add_baseline_methylation(manifest.copy())

    # ── Snapshot BEFORE artifacts ────────────────────────────────────────────
    df_before = df.copy()

    # ── Inject artifacts ─────────────────────────────────────────────────────
    print("\n[3] Injecting stumper artifacts...")
    df = inject_artifact1_confounded_batch(df)
    df = inject_artifact2_clonal_vdj(df)
    df = inject_artifact3_bisulfite_failure(df)
    df = inject_artifact4_sample_duplication(df, manifest)
    df = inject_artifact5_contamination(df)
    df = inject_artifact6_low_depth(df)
    df = inject_true_biological_signal(df)
    df = inject_borderline_signal(df)
    df = inject_subtle_signal(df)
    df = inject_xci_signal(df)
    df = inject_artifact7_sex_mixup(df)

    # ── Clip & round ─────────────────────────────────────────────────────────
    df["beta_value"] = df["beta_value"].clip(0.0, 1.0).round(4)
    df["non_cpg_meth_rate"] = df["non_cpg_meth_rate"].clip(0.0, 1.0).round(6)

    # ── Save CSV ─────────────────────────────────────────────────────────────
    print(f"\n[4] Saving CSV → {OUT_CSV}")
    df.to_csv(OUT_CSV, index=False)
    print(f"    Final shape: {df.shape}  ({df['sample_id'].nunique()} samples)")

    # ── Summary statistics ────────────────────────────────────────────────────
    print("\n[5] Artifact Summary:")
    print(f"    Batch_01 Case mean beta   : "
          f"{df[(df.batch_id=='Batch_01') & (df.disease_label==CASE_LABEL)]['beta_value'].mean():.3f}")
    print(f"    Batch_01 Control mean beta: "
          f"{df[(df.batch_id=='Batch_01') & (df.disease_label==CTRL_LABEL)]['beta_value'].mean():.3f}")
    bis_fail = df.groupby("sample_id")["non_cpg_meth_rate"].mean()
    print(f"    Samples with non_cpg > 2% : "
          f"{(bis_fail > 0.02).sum()} → {list(bis_fail[bis_fail > 0.02].index)}")
    if "S_DUP" in df["sample_id"].values:
        orig = df[df.sample_id == "S010"]["beta_value"].values
        clone = df[df.sample_id == "S_DUP"]["beta_value"].values
        r, _ = pearsonr(orig, clone)
        print(f"    S010 vs S_DUP Pearson r   : {r:.4f}")
    print(f"    S020 (contaminated) mean β : "
          f"{df[df.sample_id == 'S020']['beta_value'].mean():.3f}")
    print(f"    S030 (low coverage) mean depth: "
          f"{df[df.sample_id == 'S030']['depth'].mean():.1f}x")
    s35_x_beta = df[(df.sample_id == "S035") & df["is_x_chromosome"].astype(bool)]["beta_value"].mean()
    s36_x_beta = df[(df.sample_id == "S036") & df["is_x_chromosome"].astype(bool)]["beta_value"].mean()
    print(f"    S035 reported sex='M', true X mean β: {s35_x_beta:.3f} (expect ~0.50 — female XCI)")
    print(f"    S036 reported sex='F', true X mean β: {s36_x_beta:.3f} (expect ~0.25 — male)")
    print(f"    is_x_chromosome: {df['is_x_chromosome'].sum()} X-linked rows  "
          f"({df['is_x_chromosome'].astype(bool).sum() // df['sample_id'].nunique()} per sample)")
    true_dmr_cpgs = [f"cg{i:08d}" for i in range(300, 311)]
    case_dmr_beta = df[(df["disease_label"] == "Case") & (df["cpg_id"].isin(true_dmr_cpgs))]["beta_value"].mean()
    ctrl_dmr_beta = df[(df["disease_label"] == "Control") & (df["cpg_id"].isin(true_dmr_cpgs))]["beta_value"].mean()
    print(f"    True DMR (cg300–310)  Case mean β: {case_dmr_beta:.3f}  "
          f"Control mean β: {ctrl_dmr_beta:.3f}  Δ={case_dmr_beta - ctrl_dmr_beta:+.3f}")
    border_cpgs = [f"cg{i:08d}" for i in range(150, 158)]
    case_border = df[(df["disease_label"] == "Case") & (df["cpg_id"].isin(border_cpgs))]["beta_value"].mean()
    ctrl_border = df[(df["disease_label"] == "Control") & (df["cpg_id"].isin(border_cpgs))]["beta_value"].mean()
    print(f"    Borderline (cg150–157) Case mean β: {case_border:.3f}  "
          f"Control mean β: {ctrl_border:.3f}  Δ={case_border - ctrl_border:+.3f}")
    subtle_cpgs = [f"cg{i:08d}" for i in range(200, 206)]
    case_subtle = df[(df["disease_label"] == "Case") & (df["cpg_id"].isin(subtle_cpgs))]["beta_value"].mean()
    ctrl_subtle = df[(df["disease_label"] == "Control") & (df["cpg_id"].isin(subtle_cpgs))]["beta_value"].mean()
    print(f"    Subtle (cg200–205)     Case mean β: {case_subtle:.3f}  "
          f"Control mean β: {ctrl_subtle:.3f}  Δ={case_subtle - ctrl_subtle:+.3f}")

    # ── Before/After visualisation ────────────────────────────────────────────
    print("\n[6] Generating Before/After visualisation...")
    # Align df_before to same columns for heatmap panel (no S_DUP row)
    plot_before_after(df_before, df)

    print("\n[Done] Mock data generation complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()
