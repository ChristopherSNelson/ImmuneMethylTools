"""
core/deconvolution.py — ImmuneMethylTools Cell-Type Guard
==========================================================
Estimates immune cell-type fractions and flags unexpected lineage shifts.

Biological intent
-----------------
Whole-blood and PBMC methylation data is a mixture of cell types.  If a
'B-cell' study contains samples with very low B-cell fractions (or high
Treg fractions), the methylation signal is driven by cell-type composition,
not disease-associated methylation — a confound indistinguishable from true
DMRs without deconvolution.

Two canonical lineage markers are checked:

  FoxP3 (FOXP3, Xp11.23) — Regulatory T-cell (Treg) master transcription factor.
    In Tregs:       FoxP3 promoter CpGs are unmethylated (low beta ≈ 0.1–0.2).
    In non-Tregs:   FoxP3 is silenced by methylation (high beta ≈ 0.8–0.9).
    Detection:      Unexpected FoxP3 hypomethylation in non-Treg samples
                    → elevated Treg contamination.

  Pax5 (PAX5, 9p13) — B-cell master regulator.
    In B cells:     PAX5 locus is hypomethylated (low beta ≈ 0.1–0.2).
    In T cells:     PAX5 is methylated (high beta ≈ 0.7–0.9).
    Detection:      PAX5 hypermethylation → loss of B-cell epigenetic identity
                    → lineage shift or mislabelled sample.

GRCh38 locus coordinates
-------------------------
  FOXP3:  chrX:49,250,136–49,255,701
  PAX5:   chr9:36,896,702–37,175,926

Sex-aware FoxP3 thresholds
---------------------------
FoxP3 sits on chrX.  In females, X-Chromosome Inactivation (XCI) elevates
the observed beta by ~0.25 relative to males (single active X).  To avoid
false-positive Treg flags in female samples, separate thresholds are used:
  Male threshold:   0.15 (single X, no XCI offset)
  Female threshold: 0.30 (accounts for XCI-driven beta elevation)

Mock implementation note
------------------------
The mock dataset uses numeric cpg_ids without genomic coordinates.  We
designate specific CpG indices as proxies for these loci:
  FoxP3 marker CpGs: cg00000001–cg00000005 (X-linked in mock data)
  PAX5  marker CpGs: cg00000006–cg00000010 (X-linked in mock data)

Both marker panels fall on chrX in the mock CpG layout, so XCI correction
is applied to both when computing cell fractions for female samples.

In production, these would be replaced with published reference CpG panels
(e.g., Houseman 2012, Reinius 2012, or EpiDISH reference matrices).
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd

# ── Locus reference (GRCh38) ─────────────────────────────────────────────────
FOXP3_LOCUS = {
    "chrom": "chrX", "start": 49_250_136, "end": 49_255_701,
    "lineage": "Treg", "expected_in_Treg": "low beta (<0.30)",
}
PAX5_LOCUS = {
    "chrom": "chr9", "start": 36_896_702, "end": 37_175_926,
    "lineage": "B-cell", "expected_in_Bcell": "low beta (<0.30)",
}

# ── Mock marker CpG proxy identifiers ────────────────────────────────────────
FOXP3_CPG_PROXY: list[str] = [f"cg{i:08d}" for i in range(1, 6)]   # cg00000001–05
PAX5_CPG_PROXY: list[str] = [f"cg{i:08d}" for i in range(6, 11)]  # cg00000006–10

# ── Thresholds ────────────────────────────────────────────────────────────────
FOXP3_TREG_FLAG_THRESH_M = 0.15   # FoxP3 beta below this (males) → elevated Treg signal
FOXP3_TREG_FLAG_THRESH_F = 0.30   # FoxP3 beta below this (females) → elevated Treg signal
                                   # (higher than male threshold because XCI inflates X-linked beta;
                                   #  female X-linked baseline is ~0.50 with SD=0.05, so 0.30
                                   #  is well below normal variation)
PAX5_BCELL_SHIFT_THRESH = 0.50    # PAX5 beta above this → loss of B-cell identity
TREG_HIGH_FRAC_THRESH = 0.12     # estimated Treg fraction above this → flag

# XCI baseline offset for X-linked markers in females
# Female X-linked beta ≈ 0.50 vs male ≈ 0.25 due to X-Chromosome Inactivation
XCI_BETA_OFFSET = 0.25


# =============================================================================
# Public API
# =============================================================================


def estimate_cell_fractions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Marker-based cell-fraction estimator using FoxP3 and PAX5 proxy CpGs.

    Uses the 5-CpG proxy panels to derive per-sample fractions:
      - FoxP3 (cg00000001-05, chrX): low methylation → Treg promoter active → Treg fraction
      - PAX5  (cg00000006-10, chrX in mock data): low methylation → B-cell identity → B fraction

    Sex correction: FoxP3 and PAX5 proxy CpGs fall on chrX in mock data, so
    female samples show elevated beta (~0.50) due to X-Chromosome Inactivation.
    An XCI offset is subtracted for females before converting to fractions.

    In production this would be replaced with a reference-based deconvolution
    method (e.g., EpiDISH, MethylResolver, or CIBERSORT).

    Parameters
    ----------
    df : long-format methylation DataFrame (must include ``sex`` column)

    Returns
    -------
    pd.DataFrame with columns:
        [sample_id, b_fraction, t_fraction, treg_fraction, other_fraction]
    All fractions sum to 1.0 per sample.
    """
    records = []
    for sid, grp in df.groupby("sample_id"):
        sex = grp["sex"].iloc[0] if "sex" in grp.columns else "M"

        foxp3_rows = grp[grp["cpg_id"].isin(FOXP3_CPG_PROXY)]
        pax5_rows = grp[grp["cpg_id"].isin(PAX5_CPG_PROXY)]

        foxp3_beta = float(foxp3_rows["beta_value"].mean()) if len(foxp3_rows) else 0.5
        pax5_beta = float(pax5_rows["beta_value"].mean()) if len(pax5_rows) else 0.5

        # Sex-adjust X-linked markers: subtract XCI offset for females
        foxp3_adj = foxp3_beta
        pax5_adj = pax5_beta
        foxp3_x = (len(foxp3_rows) > 0 and foxp3_rows["is_x_chromosome"].astype(bool).any())
        pax5_x = (len(pax5_rows) > 0 and pax5_rows["is_x_chromosome"].astype(bool).any())
        if sex == "F":
            if foxp3_x:
                foxp3_adj = max(0.0, foxp3_beta - XCI_BETA_OFFSET)
            if pax5_x:
                pax5_adj = max(0.0, pax5_beta - XCI_BETA_OFFSET)

        # FoxP3: low methylation → Treg promoter active → higher Treg fraction
        # Scale [0, 1] adjusted beta to ~[0.01, 0.12] Treg fraction range
        treg_frac = max(0.01, (1.0 - foxp3_adj) * 0.12)

        # PAX5: low methylation → B-cell identity maintained → higher B fraction
        # Scale [0, 1] adjusted beta to ~[0.10, 0.75] B fraction range
        b_frac = max(0.10, (1.0 - pax5_adj) * 0.75)

        # T-cell: complement (with floor)
        t_frac = max(0.05, 1.0 - b_frac - treg_frac - 0.05)

        # Other (monocytes, NK, etc.)
        other_frac = max(0.01, 1.0 - b_frac - t_frac - treg_frac)

        # Normalize to sum to exactly 1.0
        total = b_frac + t_frac + treg_frac + other_frac
        records.append({
            "sample_id": sid,
            "b_fraction": round(b_frac / total, 4),
            "t_fraction": round(t_frac / total, 4),
            "treg_fraction": round(treg_frac / total, 4),
            "other_fraction": round(other_frac / total, 4),
        })

    return pd.DataFrame(records)


def detect_lineage_shift(df: pd.DataFrame) -> pd.DataFrame:
    """
    Check per-sample methylation at FoxP3 and PAX5 proxy loci for lineage shifts.

    FoxP3 hypomethylation → unexpected Treg signature in B-cell samples.
    PAX5  hypermethylation → loss of B-cell epigenetic identity.

    Sex-aware FoxP3 thresholds: FoxP3 sits on chrX, so X-Chromosome
    Inactivation inflates female beta values by ~0.25.  Males and females
    use different thresholds to avoid false-positive Treg flags in females.

    Parameters
    ----------
    df : long-format methylation DataFrame (must include ``sex`` column)

    Returns
    -------
    pd.DataFrame with per-sample lineage flags:
        [sample_id, sex, foxp3_mean_beta, pax5_mean_beta,
         treg_flag, bcell_shift_flag, any_lineage_flag]
    """
    records = []
    for sid, grp in df.groupby("sample_id"):
        sex = grp["sex"].iloc[0] if "sex" in grp.columns else None

        foxp3_rows = grp[grp["cpg_id"].isin(FOXP3_CPG_PROXY)]
        pax5_rows = grp[grp["cpg_id"].isin(PAX5_CPG_PROXY)]

        foxp3_mean = float(foxp3_rows["beta_value"].mean()) if len(foxp3_rows) else np.nan
        pax5_mean = float(pax5_rows["beta_value"].mean()) if len(pax5_rows) else np.nan

        # Sex-specific FoxP3 threshold (chrX → XCI inflates female beta)
        foxp3_thresh = FOXP3_TREG_FLAG_THRESH_F if sex == "F" else FOXP3_TREG_FLAG_THRESH_M
        treg_flag = (not np.isnan(foxp3_mean)) and (foxp3_mean < foxp3_thresh)
        bcell_shift = (not np.isnan(pax5_mean)) and (pax5_mean > PAX5_BCELL_SHIFT_THRESH)

        records.append(
            {
                "sample_id": sid,
                "sex": sex,
                "foxp3_mean_beta": round(foxp3_mean, 4) if not np.isnan(foxp3_mean) else np.nan,
                "pax5_mean_beta": round(pax5_mean, 4) if not np.isnan(pax5_mean) else np.nan,
                "treg_flag": treg_flag,
                "bcell_shift_flag": bcell_shift,
                "any_lineage_flag": treg_flag or bcell_shift,
            }
        )

    return pd.DataFrame(records)


# =============================================================================
# __main__ — run deconvolution checks on mock data
# =============================================================================

if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from io_utils import data_path, load_methylation, project_root, write_audit_log  # noqa: E402

    MODULE = "DECONVOLVE"
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
    audit_entries.append(ae(
        "cohort", "INFO", "Samples loaded",
        f"n={df['sample_id'].nunique()}",
    ))

    # ── Cell fractions ─────────────────────────────────────────────────────────
    fracs = estimate_cell_fractions(df)
    treg_high = fracs[fracs["treg_fraction"] > TREG_HIGH_FRAC_THRESH]
    mean_b = fracs["b_fraction"].mean()
    mean_t = fracs["t_fraction"].mean()
    mean_treg = fracs["treg_fraction"].mean()
    print(
        f"[{ts()}] [DECONVOLVE] DETECTED | Estimated cell fractions | "
        f"mean B={mean_b:.3f}  "
        f"mean T={mean_t:.3f}  "
        f"mean Treg={mean_treg:.3f}"
    )
    audit_entries.append(ae(
        "cohort", "INFO", "Cell fractions estimated",
        f"mean_B={mean_b:.3f} mean_T={mean_t:.3f} mean_Treg={mean_treg:.3f}",
    ))
    if len(treg_high):
        print(
            f"[{ts()}] [DECONVOLVE]           | "
            f"Elevated Treg (>{TREG_HIGH_FRAC_THRESH}) samples: "
            f"{treg_high['sample_id'].tolist()}"
        )
        for sid in treg_high["sample_id"]:
            frac = float(treg_high.loc[treg_high["sample_id"] == sid, "treg_fraction"].iloc[0])
            audit_entries.append(ae(
                sid, "DETECTED", "Elevated Treg fraction",
                f"treg_fraction={frac:.4f}",
            ))

    # ── Lineage shift ──────────────────────────────────────────────────────────
    shifts = detect_lineage_shift(df)
    flagged = shifts[shifts["any_lineage_flag"]]
    if len(flagged):
        for _, row in flagged.iterrows():
            parts = []
            if row.treg_flag:
                parts.append(f"FoxP3 β={row.foxp3_mean_beta:.3f}")
            if row.bcell_shift_flag:
                parts.append(f"PAX5 β={row.pax5_mean_beta:.3f}")
            sex_tag = f" sex={row.sex}" if hasattr(row, "sex") and row.sex else ""
            print(
                f"[{ts()}] [DECONVOLVE] DETECTED | Lineage shift | "
                f"sample={row.sample_id}{sex_tag}  {' | '.join(parts)}"
            )
            audit_entries.append(ae(
                row.sample_id, "DETECTED", "Lineage shift detected",
                " ".join(parts),
            ))
    else:
        print(
            f"[{ts()}] [DECONVOLVE]           | No lineage shifts detected "
            f"at current thresholds"
        )
        audit_entries.append(ae(
            "cohort", "INFO", "Lineage shift check — none detected",
            f"foxp3_thresh_M={FOXP3_TREG_FLAG_THRESH_M} foxp3_thresh_F={FOXP3_TREG_FLAG_THRESH_F} pax5_thresh={PAX5_BCELL_SHIFT_THRESH}",
        ))

    print(
        f"\n[{ts()}] [DECONVOLVE] Locus reference: "
        f"FOXP3={FOXP3_LOCUS['chrom']}:{FOXP3_LOCUS['start']:,}-{FOXP3_LOCUS['end']:,}  "
        f"PAX5={PAX5_LOCUS['chrom']}:{PAX5_LOCUS['start']:,}-{PAX5_LOCUS['end']:,}"
    )

    write_audit_log(audit_entries, _audit_csv)
    print(f"[{ts()}] [DECONVOLVE] Audit log → {_audit_csv}")
