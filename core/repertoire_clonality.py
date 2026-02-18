"""
core/repertoire_clonality.py — ImmuneMethylTools Lineage Guard
===============================================================
Detects clonal expansion artefacts in immune receptor gene loci.

Biological intent
-----------------
In B and T cells the VDJ recombination loci MUST remain UNMETHYLATED during
active recombination (open chromatin = low beta).  After a dominant clone is
established, the silenced allele becomes HYPERMETHYLATED and the chromatin
compacts — producing long fragment lengths.

This dual signature (beta > 0.80 AND fragment > 180 bp in VDJ loci) is
therefore NOT a biological methylation signal — it is a clonal expansion
ARTEFACT that must be EXCLUDED from differential methylation analysis.

GRCh38 locus coordinates (reference — for BAM-level annotation)
---------------------------------------------------------------
Locus  | Chrom | Start         | End           | Lineage
-------|-------|---------------|---------------|--------
IGH    | chr14 | 105,586,437   | 106,879,844   | B-cell (heavy chain)
IGK    | chr2  |  89,156,874   |  90,274,235   | B-cell (kappa chain)
IGL    | chr22 |  22,026,076   |  22,922,913   | B-cell (lambda chain)
TRA    | chr14 |  21,621,904   |  22,552,132   | T-cell (alpha chain)
TRB    | chr7  | 142,299,011   | 142,813,287   | T-cell (beta chain)
TRG    | chr7  |  38,279,362   |  38,407,656   | T-cell (gamma chain)
TRD    | chr14 |  22,090,057   |  22,551,535   | T-cell (delta chain)

Note: In this mock dataset the `is_vdj_region` flag serves as a proxy for
all of the above loci.  Real-world annotation would intersect CpG coordinates
with the BED intervals defined in VDJ_LOCI_GRCH38 below.
"""

import os
from datetime import datetime

import pandas as pd

# ── Thresholds ────────────────────────────────────────────────────────────────
CLONAL_BETA_MIN = 0.80   # hypermethylation threshold within VDJ loci
CLONAL_FRAG_MIN = 180    # compact-chromatin fragment length (bp)

# GRCh38 reference locus table (chrom, start, end, lineage)
VDJ_LOCI_GRCH38: dict[str, tuple[str, int, int, str]] = {
    "IGH": ("chr14", 105_586_437, 106_879_844, "B-cell"),
    "IGK": ("chr2", 89_156_874, 90_274_235, "B-cell"),
    "IGL": ("chr22", 22_026_076, 22_922_913, "B-cell"),
    "TRA": ("chr14", 21_621_904, 22_552_132, "T-cell"),
    "TRB": ("chr7", 142_299_011, 142_813_287, "T-cell"),
    "TRG": ("chr7", 38_279_362, 38_407_656, "T-cell"),
    "TRD": ("chr14", 22_090_057, 22_551_535, "T-cell"),
}


# =============================================================================
# Public API
# =============================================================================


def flag_clonal_artifacts(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Identify VDJ-region CpGs bearing the clonal expansion signature:
      High Methylation (beta > 0.80) AND Long Fragment (> 180 bp).

    Both conditions must co-occur: high beta alone can reflect true silencing;
    long fragment alone is baseline noise.  Together they are diagnostic.

    Parameters
    ----------
    df : long-format methylation DataFrame with `is_vdj_region` column

    Returns
    -------
    clonal_rows     : subset of df rows meeting both criteria
    flagged_samples : list of sample IDs carrying the artefact
    """
    vdj_mask = df["is_vdj_region"].astype(bool)
    high_beta = df["beta_value"] > CLONAL_BETA_MIN
    long_frag = df["fragment_length"] > CLONAL_FRAG_MIN
    clonal_mask = vdj_mask & high_beta & long_frag

    clonal_rows = df[clonal_mask].copy()
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


def get_vdj_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Per-patient summary of VDJ-locus methylation and fragment statistics.

    Useful for spotting which patients have elevated VDJ methylation before
    applying the strict dual-criterion clonal filter.

    Parameters
    ----------
    df : long-format methylation DataFrame

    Returns
    -------
    pd.DataFrame sorted by clonal_hits descending
    """
    vdj = df[df["is_vdj_region"].astype(bool)].copy()
    clonal_flag = (vdj["beta_value"] > CLONAL_BETA_MIN) & (
        vdj["fragment_length"] > CLONAL_FRAG_MIN
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
# __main__ — flag clonal artefacts in mock data
# =============================================================================

if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from io_utils import (  # noqa: E402
        Tee, append_flagged_samples, data_path, load_methylation, project_root, write_audit_log,
    )

    MODULE = "repertoire_clonality"
    MODULE_TAG = "CLONALITY"
    _now = datetime.now()
    run_ts = _now.strftime("%Y-%m-%dT%H:%M:%S")
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base = project_root()
    _log = os.path.join(_base, "logs", f"{MODULE}_{ts_tag}.log")
    _csv = os.path.join(_base, "data", "flagged_samples.csv")
    _audit_csv = os.path.join(_base, "data", f"audit_log_{MODULE_TAG}_{ts_tag}.csv")

    os.makedirs(os.path.join(_base, "logs"), exist_ok=True)

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
        n_vdj = int(df["is_vdj_region"].astype(bool).sum())
        audit_entries.append(ae(
            "cohort", "INFO", "VDJ clonality scan initiated",
            f"n_vdj_rows={n_vdj}",
        ))

        clonal_rows, flagged_samples = flag_clonal_artifacts(df)

        if flagged_samples:
            print(
                f"[{ts()}] [CLONALITY] DETECTED | Clonal VDJ artefact | "
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
                    sid, "DETECTED", "Clonal VDJ artifact — hypermethylation + long fragment",
                    f"mean_beta={mean_b:.3f} mean_frag={mean_f:.0f}bp",
                ))
        else:
            print(
                f"[{ts()}] [CLONALITY]           | No clonal artefacts detected "
                f"(beta>{CLONAL_BETA_MIN} AND frag>{CLONAL_FRAG_MIN} bp in VDJ loci)"
            )
            audit_entries.append(ae(
                "cohort", "INFO", "Clonal VDJ check — none detected",
                f"beta_min={CLONAL_BETA_MIN} frag_min={CLONAL_FRAG_MIN}bp",
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
            f"→ data/flagged_samples.csv"
        )

        # ── Audit log ──────────────────────────────────────────────────────────────
        write_audit_log(audit_entries, _audit_csv)
        print(f"[{ts()}] [CLONALITY] Audit log → {_audit_csv}")
