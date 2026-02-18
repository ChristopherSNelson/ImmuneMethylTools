"""
xci_guard.py — X-Chromosome Inactivation (XCI) Signal Detector

In XX (female) individuals, one X chromosome is silenced (Lyonization).
This produces a characteristic mean methylation of ~0.5 at X-linked CpG sites
in bulk tissue.  XY (male) samples show low methylation at X-linked loci (~0.25).

A mismatch between reported sex and observed X-linked methylation signals a
sample swap or metadata error — a potential case/control label mixup.

Detection logic
---------------
  Female (sex="F"): mean X-linked beta must be in [XCI_FEMALE_LO, XCI_FEMALE_HI]
  Male   (sex="M"): mean X-linked beta must be < XCI_MALE_HI
  Samples with fewer than N_X_CPGS_MIN X-linked CpGs are skipped (insufficient data).

Outputs
-------
  output/audit_log_XCI_GUARD_{ts}.csv   — per-run DETECTED/INFO entries
  output/flagged_samples.csv            — cumulative flagged-sample log
"""

import os
import sys
from datetime import datetime

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from io_utils import (  # noqa: E402
    Tee,
    append_flagged_samples,
    data_path,
    load_methylation,
    output_path,
    write_audit_log,
)

MODULE = "XCI_GUARD"

# ── Detection thresholds ──────────────────────────────────────────────────────
XCI_FEMALE_LO = 0.40  # female mean X-beta must be >= 0.40
XCI_FEMALE_HI = 0.65  # female mean X-beta must be <= 0.65
XCI_MALE_HI = 0.35    # male mean X-beta must be < 0.35
N_X_CPGS_MIN = 5      # minimum X-linked CpGs required to assess signal


# =============================================================================
# Public API
# =============================================================================


def compute_xci_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute per-sample mean beta at X-chromosome CpGs.

    Parameters
    ----------
    df : pd.DataFrame
        Methylation DataFrame containing at minimum the columns:
        sample_id, sex, is_x_chromosome, beta_value.

    Returns
    -------
    pd.DataFrame with columns:
        sample_id, sex, mean_x_beta, n_x_cpgs
    Only rows where is_x_chromosome is True are used.
    """
    x_df = df[df["is_x_chromosome"].astype(bool)]
    summary = (
        x_df.groupby(["sample_id", "sex"])
        .agg(mean_x_beta=("beta_value", "mean"), n_x_cpgs=("beta_value", "count"))
        .reset_index()
    )
    return summary


def detect_sex_mixups(df: pd.DataFrame) -> tuple[list[str], pd.DataFrame]:
    """
    Flag samples whose X-linked methylation signal contradicts their reported sex.

    Female (sex="F"): mean X-linked beta must be in [XCI_FEMALE_LO, XCI_FEMALE_HI].
    Male   (sex="M"): mean X-linked beta must be < XCI_MALE_HI.
    Samples with fewer than N_X_CPGS_MIN X-linked CpGs are skipped.

    Parameters
    ----------
    df : pd.DataFrame
        Methylation DataFrame with columns: sample_id, sex, is_x_chromosome, beta_value.

    Returns
    -------
    flagged_samples : list[str]
        Sample IDs with sex-signal mismatch (candidates for exclusion).
    report_df : pd.DataFrame
        Per-sample summary with columns:
        sample_id, sex, mean_x_beta, n_x_cpgs, xci_mismatch (bool).
    """
    summary = compute_xci_signal(df)
    if summary.empty or "sex" not in summary.columns:
        return [], summary

    f_mask = summary["sex"] == "F"
    m_mask = summary["sex"] == "M"
    sufficient = summary["n_x_cpgs"] >= N_X_CPGS_MIN

    f_ok = f_mask & sufficient & summary["mean_x_beta"].between(XCI_FEMALE_LO, XCI_FEMALE_HI)
    m_ok = m_mask & sufficient & (summary["mean_x_beta"] < XCI_MALE_HI)

    # A sample is mismatched if it has sufficient data AND does not satisfy its
    # sex-appropriate beta criterion.  Samples with insufficient X-linked CpGs
    # are left as non-mismatch to avoid false positives.
    summary["xci_mismatch"] = ~(f_ok | m_ok | ~sufficient)
    flagged = summary.loc[summary["xci_mismatch"], "sample_id"].tolist()
    return flagged, summary


# =============================================================================
# __main__ — standalone detection run
# =============================================================================

if __name__ == "__main__":
    _now = datetime.now()
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _log_path = output_path(f"logs/xci_guard_{ts_tag}.log")
    _audit_path = output_path(f"audit_log_XCI_GUARD_{ts_tag}.csv")
    _flag_path = output_path("flagged_samples.csv")
    os.makedirs(output_path("logs"), exist_ok=True)
    os.makedirs(output_path(""), exist_ok=True)

    with Tee(_log_path):
        def _ts() -> str:
            return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        print(f"[{_ts()}] [{MODULE}] Loading mock data...")
        df = load_methylation(data_path("mock_methylation.csv"), verbose=True)
        print(f"[{_ts()}] [{MODULE}] Samples: {df['sample_id'].nunique()}")

        flagged, report = detect_sex_mixups(df)

        audit_entries = []
        flagged_rows = []

        print(f"\n[{_ts()}] [{MODULE}] XCI Signal Summary (X-linked CpGs only):")
        for _, row in report.iterrows():
            status = "MISMATCH" if row["xci_mismatch"] else "OK      "
            print(
                f"[{_ts()}] [{MODULE}]   {row['sample_id']:8s}  "
                f"sex={row['sex']}  mean_x_beta={row['mean_x_beta']:.3f}  "
                f"n_x_cpgs={row['n_x_cpgs']:3d}  [{status}]"
            )

        print()
        if flagged:
            for sid in flagged:
                row = report.set_index("sample_id").loc[sid]
                metric = f"mean_x_beta={row['mean_x_beta']:.3f} sex={row['sex']}"
                print(
                    f"[{_ts()}] [{MODULE}] DETECTED | XCI sex-signal mismatch | "
                    f"{sid} | {metric}"
                )
                audit_entries.append({
                    "timestamp": _now.strftime("%Y-%m-%dT%H:%M:%S"),
                    "module": MODULE,
                    "sample_id": sid,
                    "status": "DETECTED",
                    "description": "X-linked methylation signal contradicts reported sex",
                    "metric": metric,
                })
                flagged_rows.append({
                    "run_timestamp": ts_tag,
                    "module": MODULE,
                    "sample_id": sid,
                    "flag_type": "sex_mismatch",
                    "detail": metric,
                })
        else:
            print(f"[{_ts()}] [{MODULE}] No sex-signal mismatches detected.")

        audit_entries.append({
            "timestamp": _now.strftime("%Y-%m-%dT%H:%M:%S"),
            "module": MODULE,
            "sample_id": "cohort",
            "status": "INFO",
            "description": "XCI sex-signal check complete",
            "metric": f"n_flagged={len(flagged)} of {report['sample_id'].nunique()}",
        })

        write_audit_log(audit_entries, _audit_path)
        append_flagged_samples(flagged_rows, _flag_path)

        print(f"\n[{_ts()}] [{MODULE}] Audit log  → {_audit_path}")
        print(f"[{_ts()}] [{MODULE}] Flag log   → {_flag_path}")
        print(f"[{_ts()}] [{MODULE}] Run log    → {_log_path}")
