"""
core/repertoire_clonality.py — ImmuneMethylTools Lineage Guard
===============================================================
Detects clonal expansion artifacts in immune receptor gene loci.

Biological intent
-----------------
In B and T cells the VDJ recombination loci MUST remain UNMETHYLATED during
active recombination (open chromatin = low beta).  After a dominant clone is
established, the silenced allele becomes HYPERMETHYLATED and the chromatin
compacts, producing long fragment lengths.

This dual signature (beta > 0.80 AND fragment length outlier in VDJ loci) is
therefore NOT a biological methylation signal; it is a clonal expansion
artifact that must be excluded from differential methylation analysis.

RRBS-safe detection
-------------------
Rather than a fixed fragment-length cutoff (which varies by library prep),
fragment outliers are defined as > 3 SD from each sample's own mean.
A sample is only flagged if at least 3 dual-criteria sites fall within a
single VDJ locus, reducing false positives from sporadic noise.

GRCh38 locus coordinates (with 2 kb buffer)
---------------------------------------------
Locus    | Chrom | Start         | End           | Lineage
---------|-------|---------------|---------------|--------
IGH      | chr14 | 105,584,437   | 106,881,844   | B-cell (heavy chain)
IGK      | chr2  |  88,855,361   |  90,237,368   | B-cell (kappa chain)
IGL      | chr22 |  22,024,076   |  22,924,913   | B-cell (lambda chain)
TRA_TRD  | chr14 |  22,088,057   |  23,023,075   | T-cell (alpha/delta)
TRB      | chr7  | 142,297,011   | 142,815,287   | T-cell (beta chain)
TRG      | chr7  |  38,238,024   |  38,370,055   | T-cell (gamma chain)

The `annotate_vdj_regions()` function intersects CpG coordinates (chrom, pos)
with these intervals to set the `is_vdj_region` flag.  For legacy data without
coordinates, the existing boolean flag is used as-is.
"""

import os
from datetime import datetime

import pandas as pd

# ── Thresholds ────────────────────────────────────────────────────────────────
CLONAL_BETA_MIN = 0.80       # hypermethylation threshold within VDJ loci
CLONAL_FRAG_SD_THRESH = 3.0  # flag fragment length outliers > this many SD from sample mean
CLONAL_MIN_LOCUS_HITS = 3    # require ≥ this many dual-criteria hits in a single VDJ locus

# GRCh38 reference locus table (chrom, start, end, lineage) — with 2 kb buffer
VDJ_LOCI_GRCH38: dict[str, tuple[str, int, int, str]] = {
    "IGH":     ("chr14", 105_584_437, 106_881_844, "B-cell"),
    "IGK":     ("chr2",   88_855_361,  90_237_368, "B-cell"),
    "IGL":     ("chr22",  22_024_076,  22_924_913, "B-cell"),
    "TRA_TRD": ("chr14",  22_088_057,  23_023_075, "T-cell"),
    "TRB":     ("chr7",  142_297_011, 142_815_287, "T-cell"),
    "TRG":     ("chr7",   38_238_024,  38_370_055, "T-cell"),
}


# =============================================================================
# Public API
# =============================================================================


def annotate_vdj_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Set ``is_vdj_region`` and ``vdj_locus`` from GRCh38 coordinate intervals.

    If the DataFrame contains ``chrom`` and ``pos`` columns, the function
    deduplicates by ``cpg_id``, checks each CpG against the VDJ locus
    intervals, and maps the result back to all rows.

    If ``chrom`` or ``pos`` is absent (legacy data), the DataFrame is
    returned unchanged — backward compatible.

    Parameters
    ----------
    df : long-format methylation DataFrame

    Returns
    -------
    pd.DataFrame with ``is_vdj_region`` and ``vdj_locus`` set from coordinates
    (or unchanged for legacy data)
    """
    if "chrom" not in df.columns or "pos" not in df.columns:
        return df

    # Deduplicate: one check per CpG, not per row
    cpg_coords = df[["cpg_id", "chrom", "pos"]].drop_duplicates("cpg_id")

    # Build interval lookup grouped by chromosome: (start, end, locus_name)
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for name, (chrom, start, end, _lineage) in VDJ_LOCI_GRCH38.items():
        by_chrom.setdefault(chrom, []).append((start, end, name))

    def _check(row):
        for start, end, name in by_chrom.get(row.chrom, []):
            if start <= row.pos <= end:
                return name
        return None

    cpg_locus = cpg_coords.apply(_check, axis=1)
    locus_map = dict(zip(cpg_coords["cpg_id"], cpg_locus))
    df = df.copy()
    df["vdj_locus"] = df["cpg_id"].map(locus_map)
    df["is_vdj_region"] = df["vdj_locus"].notna()
    return df


def flag_clonal_artifacts(
    df: pd.DataFrame,
    beta_min: float = CLONAL_BETA_MIN,
    frag_sd_thresh: float = CLONAL_FRAG_SD_THRESH,
    min_locus_hits: int = CLONAL_MIN_LOCUS_HITS,
) -> tuple[pd.DataFrame, list[str]]:
    """
    Identify VDJ-region CpGs bearing the clonal expansion signature:
      High Methylation (beta > beta_min) AND Fragment Length Outlier (> sample
      mean + frag_sd_thresh * sample SD).

    RRBS-safe: instead of a fixed fragment-length cutoff (which varies by
    library prep), fragment outliers are defined relative to each sample's
    own distribution (>3 SD from the sample mean by default).

    A sample is only flagged if at least ``min_locus_hits`` dual-criteria
    sites fall within a single VDJ locus, reducing false positives from
    sporadic noise.

    Parameters
    ----------
    df             : long-format methylation DataFrame with ``is_vdj_region``
    beta_min       : VDJ beta above this → clonal hypermethylation (default 0.80)
    frag_sd_thresh : flag fragment lengths > sample_mean + this * sample_std (default 3.0)
    min_locus_hits : require ≥ this many hits in one locus to flag (default 3)

    Returns
    -------
    clonal_rows     : subset of df rows meeting all criteria
    flagged_samples : list of sample IDs carrying the artifact
    """
    # Per-sample fragment length statistics (across ALL sites)
    frag_stats = df.groupby("sample_id")["fragment_length"].agg(["mean", "std"])

    # Compute per-row fragment threshold: sample_mean + frag_sd_thresh * sample_std
    sample_means = df["sample_id"].map(frag_stats["mean"])
    sample_stds = df["sample_id"].map(frag_stats["std"]).fillna(0)
    frag_threshold = sample_means + frag_sd_thresh * sample_stds

    vdj_mask = df["is_vdj_region"].astype(bool)
    high_beta = df["beta_value"] > beta_min
    long_frag = df["fragment_length"] > frag_threshold
    candidate_mask = vdj_mask & high_beta & long_frag

    candidates = df[candidate_mask].copy()

    if candidates.empty:
        return candidates, []

    # Require ≥ min_locus_hits in a single VDJ locus per sample
    if "vdj_locus" in candidates.columns:
        locus_counts = (
            candidates.groupby(["sample_id", "vdj_locus"])
            .size()
            .reset_index(name="n_hits")
        )
        qualifying = locus_counts[locus_counts["n_hits"] >= min_locus_hits]
        qualifying_pairs = set(
            zip(qualifying["sample_id"], qualifying["vdj_locus"])
        )
        if qualifying_pairs:
            keep_mask = candidates.apply(
                lambda r: (r["sample_id"], r["vdj_locus"]) in qualifying_pairs,
                axis=1,
            )
            clonal_rows = candidates[keep_mask].copy()
        else:
            clonal_rows = candidates.iloc[:0].copy()
    else:
        # Legacy path (no vdj_locus column): require min_locus_hits total
        per_sample = candidates.groupby("sample_id").size()
        qualifying_samples = per_sample[per_sample >= min_locus_hits].index.tolist()
        clonal_rows = candidates[candidates["sample_id"].isin(qualifying_samples)].copy()

    flagged_samples = clonal_rows["sample_id"].unique().tolist()
    return clonal_rows, flagged_samples


def mask_clonal_vdj_sites(
    df: pd.DataFrame,
    clonal_samples: list,
) -> tuple:
    """
    Set beta_value to NaN for VDJ-region rows in clonally-expanded samples.
    Non-VDJ rows and non-clonal samples are unchanged.
    Downstream stages (normalizer, dmr_hunter, ml_guard) all handle NaN safely.

    Parameters
    ----------
    df              : long-format methylation DataFrame
    clonal_samples  : list of sample IDs identified as clonally expanded

    Returns
    -------
    df_masked : copy of df with VDJ-region beta_value set to NaN for clonal samples
    n_masked  : count of sites set to NaN
    """
    df = df.copy()
    mask = df["sample_id"].isin(clonal_samples) & df["is_vdj_region"].astype(bool)
    df.loc[mask, "beta_value"] = float("nan")
    return df, int(mask.sum())


def get_vdj_summary(
    df: pd.DataFrame,
    beta_min: float = CLONAL_BETA_MIN,
    frag_sd_thresh: float = CLONAL_FRAG_SD_THRESH,
) -> pd.DataFrame:
    """
    Per-patient summary of VDJ-locus methylation and fragment statistics.

    Useful for spotting which patients have elevated VDJ methylation before
    applying the strict dual-criterion clonal filter.

    Parameters
    ----------
    df             : long-format methylation DataFrame
    beta_min       : beta threshold for clonal hit counting (default 0.80)
    frag_sd_thresh : SD threshold for fragment outlier detection (default 3.0)

    Returns
    -------
    pd.DataFrame sorted by clonal_hits descending
    """
    # Per-sample fragment length stats (across ALL sites, not just VDJ)
    frag_stats = df.groupby("sample_id")["fragment_length"].agg(["mean", "std"])

    vdj = df[df["is_vdj_region"].astype(bool)].copy()

    # Per-row fragment threshold: sample_mean + frag_sd_thresh * sample_std
    s_means = vdj["sample_id"].map(frag_stats["mean"])
    s_stds = vdj["sample_id"].map(frag_stats["std"]).fillna(0)
    frag_threshold = s_means + frag_sd_thresh * s_stds

    clonal_flag = (vdj["beta_value"] > beta_min) & (
        vdj["fragment_length"] > frag_threshold
    )
    vdj = vdj.assign(is_clonal=clonal_flag)

    summary = (
        vdj.groupby(["patient_id", "disease_label"])
        .agg(
            n_cpgs=("cpg_id", "count"),
            mean_beta=("beta_value", "mean"),
            max_beta=("beta_value", "max"),
            mean_frag=("fragment_length", "mean"),
            max_frag=("fragment_length", "max"),
            clonal_hits=("is_clonal", "sum"),
        )
        .reset_index()
        .sort_values("clonal_hits", ascending=False)
    )
    for col in ["mean_beta", "max_beta", "mean_frag"]:
        summary[col] = summary[col].round(3)
    return summary


# =============================================================================
# __main__ — flag clonal artifacts in mock data
# =============================================================================

if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    from core.infrastructure.io_utils import (  # noqa: E402
        Tee, append_flagged_samples, audit_entry, data_path, load_methylation,
        project_root, ts, write_audit_log,
    )

    MODULE = "repertoire_clonality"
    MODULE_TAG = "CLONALITY"
    _now = datetime.now()
    run_ts = _now.strftime("%Y-%m-%dT%H:%M:%S")
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base = project_root()
    _log = os.path.join(_base, "output", "logs", f"{MODULE}_{ts_tag}.log")
    _csv = os.path.join(_base, "output", "flagged_samples.csv")
    _audit_csv = os.path.join(_base, "output", f"audit_log_{MODULE_TAG}_{ts_tag}.csv")

    os.makedirs(os.path.join(_base, "output", "logs"), exist_ok=True)

    audit_entries = []
    ae = lambda sid, st, d, m: audit_entry(MODULE_TAG, sid, st, d, m)

    with Tee(_log):
        df = load_methylation(data_path("mock_methylation.csv"))
        n_vdj = int(df["is_vdj_region"].astype(bool).sum())
        audit_entries.append(ae(
            "cohort", "INFO", "VDJ clonality scan initiated",
            f"n_vdj_rows={n_vdj}",
        ))

        clonal_rows, flagged_samples = flag_clonal_artifacts(df)

        if flagged_samples:
            print(
                f"[{ts()}] [CLONALITY] DETECTED | Clonal VDJ artifact | "
                f"{len(clonal_rows)} CpG rows flagged | samples={flagged_samples}"
            )
            for sid in flagged_samples:
                sub = clonal_rows[clonal_rows["sample_id"] == sid]
                mean_b = sub["beta_value"].mean()
                mean_f = sub["fragment_length"].mean()
                print(
                    f"[{ts()}] [CLONALITY]           | Sample {sid} | "
                    f"n_rows={len(sub)}  "
                    f"mean_beta={mean_b:.3f}  "
                    f"mean_frag={mean_f:.0f} bp"
                )
                audit_entries.append(ae(
                    sid, "DETECTED", "Clonal VDJ artifact — hypermethylation + long fragment outlier",
                    f"mean_beta={mean_b:.3f} mean_frag={mean_f:.0f}bp",
                ))
        else:
            print(
                f"[{ts()}] [CLONALITY]           | No clonal artifacts detected "
                f"(beta>{CLONAL_BETA_MIN} AND frag>{CLONAL_FRAG_SD_THRESH} SD in VDJ loci)"
            )
            audit_entries.append(ae(
                "cohort", "INFO", "Clonal VDJ check — none detected",
                f"beta_min={CLONAL_BETA_MIN} frag_sd={CLONAL_FRAG_SD_THRESH}",
            ))

        print(f"\n[{ts()}] [CLONALITY] VDJ per-patient summary:")
        summary = get_vdj_summary(df)
        print(summary.to_string(index=False))

        print(f"\n[{ts()}] [CLONALITY] GRCh38 Locus Reference:")
        for name, (chrom, start, end, lineage) in VDJ_LOCI_GRCH38.items():
            print(f"  {name:4s}  {chrom:5s}  {start:>12,} – {end:>12,}  ({lineage})")

        # ── Persist flagged_samples.csv ────────────────────────────────────────────
        flagged_rows = []
        for sid in flagged_samples:
            sub = clonal_rows[clonal_rows["sample_id"] == sid]
            flagged_rows.append({
                "run_timestamp": run_ts,
                "module": MODULE,
                "sample_id": sid,
                "flag_type": "clonal_vdj",
                "detail": (
                    f"mean_beta={sub['beta_value'].mean():.2f} "
                    f"mean_frag={sub['fragment_length'].mean():.0f}bp"
                ),
            })

        append_flagged_samples(flagged_rows, _csv)
        n_unique = len({r["sample_id"] for r in flagged_rows})
        print(
            f"[{ts()}] [CLONALITY] {n_unique} unique flagged sample(s) written "
            f"→ output/flagged_samples.csv"
        )

        # ── Audit log ──────────────────────────────────────────────────────────────
        write_audit_log(audit_entries, _audit_csv)
        print(f"[{ts()}] [CLONALITY] Audit log → {_audit_csv}")
