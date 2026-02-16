"""
core/io_utils.py — ImmuneMethylTools Shared I/O Utilities
==========================================================
Provides two facilities used by sample-flagging modules:

  append_flagged_samples(rows, csv_path)
      Persist flagged-sample records to a cumulative CSV log.
      Rows are appended; history is never overwritten.

  class Tee
      Context manager that mirrors sys.stdout to a timestamped log file so
      every __main__ run is captured for auditing.
"""

import csv
import os
import sys

# Column order for the flagged-samples CSV
_FIELDNAMES = ["run_timestamp", "module", "sample_id", "flag_type", "detail"]


# =============================================================================
# append_flagged_samples
# =============================================================================


def append_flagged_samples(rows, csv_path="data/flagged_samples.csv"):
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

    parent = os.path.dirname(os.path.abspath(csv_path))
    os.makedirs(parent, exist_ok=True)

    file_exists = os.path.isfile(csv_path)

    with open(csv_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_FIELDNAMES)
        if not file_exists:
            writer.writeheader()
        writer.writerows(rows)


# =============================================================================
# Tee — stdout mirror
# =============================================================================


class Tee:
    """
    Context manager that duplicates sys.stdout to a log file.

    Both the terminal and the log file receive every write; neither is
    buffered independently.  The original sys.stdout is restored on exit.

    Usage::

        with Tee("logs/qc_guard_20260215_181404.log"):
            print("Appears on terminal AND in the log file.")

    Parameters
    ----------
    log_path : str
        Path to the log file.  Parent directory is created if absent.
    """

    def __init__(self, log_path: str):
        self._log_path = log_path
        self._log_fh   = None
        self._orig     = None

    # ── Context manager protocol ───────────────────────────────────────────────

    def __enter__(self) -> "Tee":
        parent = os.path.dirname(os.path.abspath(self._log_path))
        os.makedirs(parent, exist_ok=True)
        self._log_fh = open(self._log_path, "w")
        self._orig   = sys.stdout
        sys.stdout   = self
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
