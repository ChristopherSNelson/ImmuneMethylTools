"""
core/io_utils.py — ImmuneMethylTools Shared I/O Utilities
==========================================================
Provides seven facilities used by core modules and notebooks:

  project_root() -> str
      Absolute path to the project root (parent of core/).
      Resolves relative to this file, so it works from any working directory.

  data_path(filename) -> str
      Absolute path to an INPUT file in the project's data/ directory.
      Use for mock_methylation.csv and other input data only.

  output_path(filename) -> str
      Absolute path to a file in the project's output/ directory.
      Use for all generated files: logs, figures, audit CSVs, reports.

  load_methylation(csv_path, *, verbose=True) -> pd.DataFrame
      Safely load and schema-validate a methylation CSV.
      Raises FileNotFoundError / ValueError on any violation.

  append_flagged_samples(rows, csv_path)
      Persist flagged-sample records to a cumulative CSV log.
      Rows are appended; history is never overwritten.

  write_audit_log(entries, csv_path)
      Write a per-run audit log (DETECTED + INFO events) to a new
      timestamped CSV file under output/.  Schema:
        timestamp, module, sample_id, status, description, metric

  class Tee
      Context manager that mirrors sys.stdout to a log file simultaneously.
      Used via context manager in each __main__ block.
"""

import csv
import os
import sys

import pandas as pd

# ── Schema constants ──────────────────────────────────────────────────────────

# All 16 columns defined in CLAUDE.md Key Variable Glossary.
REQUIRED_COLUMNS: list[str] = [
    "sample_id",
    "patient_id",
    "batch_id",
    "age",
    "disease_label",
    "cpg_id",
    "beta_value",
    "depth",
    "fragment_length",
    "is_vdj_region",
    "non_cpg_meth_rate",
    "sex",              # "M" or "F" — reported sex from sample metadata
    "is_x_chromosome",  # bool — True if CpG falls on X chromosome
    "gc_content",       # float [0,1] — GC fraction around the CpG site
    "chrom",            # str — GRCh38 chromosome (e.g. "chr1", "chrX")
    "pos",              # int — GRCh38 genomic position
]

VALID_DISEASE_LABELS: frozenset[str] = frozenset({"Case", "Control"})

# Column order for the flagged-samples CSV
_FIELDNAMES = ["run_timestamp", "module", "sample_id", "flag_type", "detail"]

# Column order for the per-run audit log CSV
_AUDIT_FIELDNAMES = [
    "timestamp", "module", "sample_id", "status", "description", "metric"
]


# =============================================================================
# project_root / data_path
# =============================================================================


def project_root() -> str:
    """
    Absolute path to the project root (parent of the core/ directory).

    Resolves relative to this source file, so it returns the correct path
    regardless of the caller's working directory — including notebooks and
    scripts outside core/.
    """
    return os.path.normpath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    )


def data_path(filename: str) -> str:
    """
    Absolute path to an input file in the project's data/ directory.

    Parameters
    ----------
    filename : str
        Filename only (e.g. "mock_methylation.csv").

    Returns
    -------
    str — full absolute path, e.g. /path/to/ImmuneMethylTools/data/filename
    """
    return os.path.join(project_root(), "data", filename)


def output_path(filename: str) -> str:
    """
    Absolute path to a file in the project's output/ directory.

    Use this for all generated outputs: audit logs, figures, clean data
    exports, PDF reports, and run logs.  The output/ directory is created
    automatically by each module that writes to it.

    Parameters
    ----------
    filename : str
        Filename or sub-path (e.g. "clean_methylation.csv",
        "logs/pipeline_20260217.log", "figures/pca.png").

    Returns
    -------
    str — full absolute path, e.g. /path/to/ImmuneMethylTools/output/filename
    """
    return os.path.join(project_root(), "output", filename)


# =============================================================================
# load_methylation — safe loader + schema validator
# =============================================================================


def load_methylation(csv_path: str, *, verbose: bool = True) -> pd.DataFrame:
    """
    Safely load and validate a methylation CSV file.

    Checks (in order):
      1. File exists at csv_path.
      2. All 16 required columns present (CLAUDE.md schema).
      3. beta_value        ∈ [0, 1]
      4. non_cpg_meth_rate ∈ [0, 1]
      5. depth             ≥ 0
      6. disease_label     ∈ {"Case", "Control"}
      7. sex               ∈ {"M", "F"}
      8. gc_content        ∈ [0, 1]
      9. chrom             valid chromosome strings
     10. pos               ≥ 0

    Parameters
    ----------
    csv_path : str
        Absolute or relative path to the methylation CSV.
    verbose : bool
        If True (default), print a one-line loading summary to stdout.
        Set to False in pipeline / notebook contexts that manage their own
        logging output.

    Returns
    -------
    pd.DataFrame — validated methylation data, schema unchanged.

    Raises
    ------
    FileNotFoundError
        File not found at csv_path.
    ValueError
        Schema or value-range violation; message names the failing check.
    """
    # 1. File existence
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(
            f"load_methylation: file not found: {csv_path}\n"
            "  Run `python data/generate_mock_data.py` to create mock data."
        )

    df = pd.read_csv(csv_path)

    # 2. Required columns
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"load_methylation: missing required column(s): {missing}\n"
            f"  Found: {df.columns.tolist()}"
        )

    # 3. beta_value range
    bv = df["beta_value"].dropna()
    n_bad = int(((bv < 0) | (bv > 1)).sum())
    if n_bad:
        raise ValueError(
            f"load_methylation: {n_bad} beta_value row(s) outside [0, 1]."
        )

    # 4. non_cpg_meth_rate range
    nr = df["non_cpg_meth_rate"].dropna()
    n_bad = int(((nr < 0) | (nr > 1)).sum())
    if n_bad:
        raise ValueError(
            f"load_methylation: {n_bad} non_cpg_meth_rate row(s) outside [0, 1]."
        )

    # 5. depth non-negative
    n_bad = int((df["depth"].dropna() < 0).sum())
    if n_bad:
        raise ValueError(
            f"load_methylation: {n_bad} depth value(s) < 0."
        )

    # 6. disease_label vocabulary
    bad_labels = set(df["disease_label"].dropna().unique()) - VALID_DISEASE_LABELS
    if bad_labels:
        raise ValueError(
            f"load_methylation: unexpected disease_label value(s): {bad_labels}. "
            f"Expected {set(VALID_DISEASE_LABELS)}."
        )

    # 7. sex vocabulary
    valid_sex = frozenset({"M", "F"})
    bad_sex = set(df["sex"].dropna().unique()) - valid_sex
    if bad_sex:
        raise ValueError(
            f"load_methylation: sex column contains invalid values: {bad_sex}. "
            f"Expected 'M' or 'F'."
        )

    # 8. gc_content range
    gc = df["gc_content"].dropna()
    n_bad = int(((gc < 0) | (gc > 1)).sum())
    if n_bad:
        raise ValueError(
            f"load_methylation: {n_bad} gc_content row(s) outside [0, 1]."
        )

    # 9. chrom — valid chromosome strings
    valid_chroms = frozenset(
        [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    )
    bad_chroms = set(df["chrom"].dropna().unique()) - valid_chroms
    if bad_chroms:
        raise ValueError(
            f"load_methylation: chrom column contains invalid values: {bad_chroms}. "
            f"Expected chr1-chr22, chrX, chrY."
        )

    # 10. pos non-negative
    n_bad = int((df["pos"].dropna() < 0).sum())
    if n_bad:
        raise ValueError(
            f"load_methylation: {n_bad} pos value(s) < 0."
        )

    if verbose:
        n_samples = df["sample_id"].nunique()
        n_cpgs = df["cpg_id"].nunique()
        label_counts = (
            df.drop_duplicates("sample_id")["disease_label"]
            .value_counts()
            .to_dict()
        )
        batches = sorted(df["batch_id"].unique())
        print(
            f"load_methylation: {os.path.basename(csv_path)}"
            f"  shape={df.shape}"
            f"  samples={n_samples} {label_counts}"
            f"  cpgs={n_cpgs}"
            f"  batches={batches}"
        )

    return df


# =============================================================================
# append_flagged_samples
# =============================================================================


def append_flagged_samples(rows, csv_path=None):
    """
    Append flagged-sample records to a persistent CSV log.

    Parameters
    ----------
    rows : list of dict
        Each dict must contain the keys:
        run_timestamp, module, sample_id, flag_type, detail
    csv_path : str
        Path to the CSV file.  Created with a header row if it does not yet
        exist; subsequent calls append without touching the header.

    Notes
    -----
    Uses csv.DictWriter — no pandas dependency.
    Parent directory is created automatically if absent.
    """
    if not rows:
        return

    if csv_path is None:
        csv_path = output_path("flagged_samples.csv")

    parent = os.path.dirname(os.path.abspath(csv_path))
    os.makedirs(parent, exist_ok=True)

    file_exists = os.path.isfile(csv_path)

    with open(csv_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_FIELDNAMES)
        if not file_exists:
            writer.writeheader()
        writer.writerows(rows)


# =============================================================================
# write_audit_log
# =============================================================================


def write_audit_log(entries, csv_path):
    """
    Write a per-run audit log to a new timestamped CSV file.

    Parameters
    ----------
    entries : list of dict
        Each dict must contain the keys:
        timestamp, module, sample_id, status, description, metric
        status must be 'DETECTED' or 'INFO'.
    csv_path : str
        Full path to the output CSV (e.g. data/audit_log_20260215_143005.csv).
        File is created fresh each call — use a timestamp in the filename so
        successive runs do not overwrite each other.

    Notes
    -----
    Uses csv.DictWriter — no pandas dependency.
    Parent directory is created automatically if absent.
    """
    if not entries:
        return

    parent = os.path.dirname(os.path.abspath(csv_path))
    os.makedirs(parent, exist_ok=True)

    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_AUDIT_FIELDNAMES)
        writer.writeheader()
        writer.writerows(entries)


# =============================================================================
# Tee — stdout mirror
# =============================================================================


class Tee:
    """
    Context manager that duplicates sys.stdout to a log file.

    Both the terminal and the log file receive every write; neither is
    buffered independently.  The original sys.stdout is restored on exit.

    Usage::

        with Tee("output/logs/qc_guard_20260215_181404.log"):
            print("Appears on terminal AND in the log file.")

    Parameters
    ----------
    log_path : str
        Path to the log file.  Parent directory is created if absent.
    """

    def __init__(self, log_path: str):
        self._log_path = log_path
        self._log_fh = None
        self._orig = None

    # ── Context manager protocol ───────────────────────────────────────────────

    def __enter__(self) -> "Tee":
        parent = os.path.dirname(os.path.abspath(self._log_path))
        os.makedirs(parent, exist_ok=True)
        self._log_fh = open(self._log_path, "w")
        self._orig = sys.stdout
        sys.stdout = self
        return self

    def __exit__(self, *_) -> None:
        sys.stdout = self._orig
        self._log_fh.close()

    # ── File-like interface ────────────────────────────────────────────────────

    def write(self, data: str) -> None:
        self._orig.write(data)
        self._log_fh.write(data)

    def flush(self) -> None:
        self._orig.flush()
        self._log_fh.flush()
