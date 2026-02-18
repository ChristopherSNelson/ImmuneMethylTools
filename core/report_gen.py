"""
core/report_gen.py -- ImmuneMethylTools PDF Report Generator
=============================================================
Bundles pipeline audit log + figures into a structured, single-file PDF report.

Usage (triggered by pipeline):
    python core/pipeline.py --report

Usage (standalone, re-uses most-recent audit log):
    python core/report_gen.py
    python core/report_gen.py --audit data/audit_log_pipeline_20260216_094537.csv

Report sections
---------------
1. Executive Summary        key-value table (counts, confound, AUC, git hash)
2. Sample Exclusion         exclusion_accounting.png
3. Global QC Metrics        qc_metrics.png
4. Beta Distribution        beta_distribution_kde.png
5. PCA                      pca_disease_label.png + pca_batch_id.png (stacked)
6. Volcano Plots            volcano.png + volcano_clonal_risk.png (stacked)
7. DETECTED Events Log      full bordered table from audit_log CSV
8. Significant DMR Table    significant windows, amber rows for clonal_risk=True
"""

import glob
import os
import subprocess
from datetime import datetime

import pandas as pd
from fpdf import FPDF

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from io_utils import data_path, project_root  # noqa: E402

FIGURES_DIR = os.path.join(project_root(), "figures")

# ---------------------------------------------------------------------------
# Latin-1 sanitizer  (Helvetica core font only supports Latin-1)
# ---------------------------------------------------------------------------

_UNICODE_MAP = str.maketrans({
    "\u2014": "--",    # em-dash
    "\u2013": "-",     # en-dash
    "\u2192": "->",    # rightwards arrow
    "\u2190": "<-",    # leftwards arrow
    "\u2194": "<->",   # left-right arrow
    "\u2265": ">=",    # greater-than or equal
    "\u2264": "<=",    # less-than or equal
    "\u2260": "!=",    # not equal
    "\u00b1": "+/-",   # plus-minus
    "\u2248": "~=",    # almost equal
    "\u00d7": "x",     # multiplication sign
    "\u25c6": "[*]",   # black diamond (used in volcano legend)
    "\u26a0": "[!]",   # warning sign
    "\u00a0": " ",     # non-breaking space
})


def _safe(text: str) -> str:
    """Replace non-Latin-1 characters so Helvetica core font can render them."""
    return str(text).translate(_UNICODE_MAP).encode("latin-1", errors="replace").decode("latin-1")


# ---------------------------------------------------------------------------
# Git provenance
# ---------------------------------------------------------------------------

def _git_info() -> tuple:
    """Return (short_hash, commit_date) from git, or fallback strings."""
    base = project_root()
    try:
        short_hash = subprocess.run(
            ["git", "-C", base, "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        commit_date = subprocess.run(
            ["git", "-C", base, "log", "-1", "--format=%ci"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()[:19]      # "YYYY-MM-DD HH:MM:SS"
    except Exception:
        short_hash = "unknown"
        commit_date = "unknown"
    return short_hash, commit_date


# ---------------------------------------------------------------------------
# FPDF subclass
# ---------------------------------------------------------------------------

class _Report(FPDF):
    """FPDF subclass with ImmuneMethylTools report helpers."""

    def __init__(self, run_ts: str, git_hash: str, commit_date: str):
        super().__init__(orientation="P", unit="mm", format="A4")
        self._run_ts = run_ts
        self._git_hash = git_hash
        self._commit_date = commit_date
        self.set_auto_page_break(auto=True, margin=20)
        self.set_margins(15, 15, 15)

    # -- Header / footer -------------------------------------------------

    def header(self) -> None:
        self.set_font("Helvetica", "B", 11)
        self.cell(0, 8, "ImmuneMethylTools -- QC Pipeline Report", align="C")
        self.ln(5)
        self.set_font("Helvetica", "", 8)
        self.set_text_color(100, 100, 100)
        self.cell(
            0, 5,
            _safe(f"Run: {self._run_ts}   |   git: {self._git_hash}   |   commit: {self._commit_date}"),
            align="C",
        )
        self.set_text_color(0, 0, 0)
        self.ln(3)
        self.set_draw_color(200, 200, 200)
        self.set_line_width(0.3)
        self.line(self.l_margin, self.get_y(), self.w - self.r_margin, self.get_y())
        self.ln(4)

    def footer(self) -> None:
        self.set_y(-15)
        self.set_font("Helvetica", "I", 7)
        self.set_text_color(120, 120, 120)
        self.cell(
            0, 5,
            _safe(f"git:{self._git_hash} | {self._run_ts} | Page {self.page_no()}"),
            align="C",
        )
        self.set_text_color(0, 0, 0)

    # -- Content helpers -------------------------------------------------

    def section(self, title: str) -> None:
        """Bold section heading with a thin colored band."""
        if self.get_y() > 240:
            self.add_page()
        self.set_font("Helvetica", "B", 13)
        self.set_fill_color(240, 244, 255)
        self.cell(0, 9, f"  {_safe(title)}", border=0, fill=True)
        self.ln(4)

    def body(self, text: str) -> None:
        """Normal paragraph text."""
        self.set_font("Helvetica", "", 9)
        self.multi_cell(0, 5, _safe(text))
        self.ln(2)

    def kv(self, key: str, value: str, flagged: bool = False) -> None:
        """Key-value table row, optionally highlighted in amber."""
        col_w = 75
        self.set_font("Helvetica", "", 9)
        if flagged:
            self.set_fill_color(255, 235, 180)
        else:
            self.set_fill_color(248, 248, 248)
        self.cell(col_w, 6, f"  {_safe(key)}", border=1, fill=True)
        if flagged:
            self.set_fill_color(255, 245, 210)
        else:
            self.set_fill_color(255, 255, 255)
        self.cell(0, 6, f"  {_safe(value)}", border=1, fill=True)
        self.ln()

    def figure(self, img_path: str, caption: str = "", w: float = 170) -> None:
        """Insert a centered figure. Graceful on missing file."""
        if not os.path.isfile(img_path):
            self.set_text_color(180, 0, 0)
            self.body(f"[Figure not found: {img_path}]")
            self.set_text_color(0, 0, 0)
            return
        x_offset = (self.w - w) / 2
        self.image(img_path, x=x_offset, w=w)
        if caption:
            self.set_font("Helvetica", "I", 8)
            self.set_text_color(100, 100, 100)
            self.cell(0, 5, _safe(caption), align="C")
            self.set_text_color(0, 0, 0)
        self.ln(4)


# ---------------------------------------------------------------------------
# Table helpers
# ---------------------------------------------------------------------------

_DET_COLS = ["timestamp", "module", "sample_id", "description", "metric"]
_DET_HEADS = ["Time", "Module", "Sample", "Description", "Metric"]
_DET_WIDTHS = [33, 22, 20, 65, 40]

_DMR_COLS = ["window_id", "n_cpgs", "delta_beta", "p_adj", "n_vdj_cpgs", "clonal_risk"]
_DMR_HEADS = ["Window", "CpGs", "dBeta", "p_adj", "VDJ CpGs", "Clonal Risk"]
_DMR_WIDTHS = [40, 18, 22, 28, 25, 47]


def _table_header(pdf: _Report, heads: list, widths: list) -> None:
    pdf.set_font("Helvetica", "B", 8)
    pdf.set_fill_color(220, 228, 255)
    for head, w in zip(heads, widths):
        pdf.cell(w, 7, f" {head}", border=1, fill=True)
    pdf.ln()


def _detected_table(pdf: _Report, detected: pd.DataFrame) -> None:
    """Render the DETECTED events as a bordered table."""
    if detected.empty:
        pdf.body("No DETECTED events in this run.")
        return

    _table_header(pdf, _DET_HEADS, _DET_WIDTHS)

    for i, (_, row) in enumerate(detected.iterrows()):
        if pdf.get_y() > 265:
            pdf.add_page()
            _table_header(pdf, _DET_HEADS, _DET_WIDTHS)

        pdf.set_font("Helvetica", "", 7)
        if i % 2 == 0:
            pdf.set_fill_color(255, 255, 255)
        else:
            pdf.set_fill_color(248, 248, 252)

        vals = [
            _safe(str(row.get("timestamp", ""))[:19]),
            _safe(str(row.get("module", ""))),
            _safe(str(row.get("sample_id", ""))),
            _safe(str(row.get("description", ""))[:60]),
            _safe(str(row.get("metric", ""))[:35]),
        ]
        for val, w in zip(vals, _DET_WIDTHS):
            pdf.cell(w, 6, f" {val}", border=1, fill=True)
        pdf.ln()


def _sig_dmr_table(pdf: _Report, sig_dmrs: pd.DataFrame) -> None:
    """Render the significant DMR table, amber rows for clonal_risk=True."""
    _table_header(pdf, _DMR_HEADS, _DMR_WIDTHS)

    for i, (_, row) in enumerate(sig_dmrs.iterrows()):
        if pdf.get_y() > 265:
            pdf.add_page()
            _table_header(pdf, _DMR_HEADS, _DMR_WIDTHS)

        pdf.set_font("Helvetica", "", 8)
        clonal = bool(row.get("clonal_risk", False))
        if clonal:
            pdf.set_fill_color(255, 235, 180)
        elif i % 2 == 0:
            pdf.set_fill_color(255, 255, 255)
        else:
            pdf.set_fill_color(248, 248, 252)

        vals = [
            _safe(str(row.get("window_id", ""))),
            str(int(row.get("n_cpgs", 0))),
            f"{float(row.get('delta_beta', 0)):+.4f}",
            f"{float(row.get('p_adj', 1)):.3e}",
            str(int(row.get("n_vdj_cpgs", 0))),
            "YES [!]" if clonal else "no",
        ]
        for val, w in zip(vals, _DMR_WIDTHS):
            pdf.cell(w, 6, f" {val}", border=1, fill=True)
        pdf.ln()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_report(
    pipeline_result: dict,
    audit_csv: str,
    output_path: str,
    run_ts: str,
) -> str:
    """
    Generate a PDF report from the pipeline result dict and audit log CSV.

    Parameters
    ----------
    pipeline_result : dict returned by run_pipeline()
    audit_csv       : path to the pipeline audit_log CSV for this run
    output_path     : destination PDF path (created if parent dir missing)
    run_ts          : ISO-format run timestamp string (e.g. "2026-02-16T09:45:37")

    Returns
    -------
    str  path to the saved PDF
    """
    git_hash, commit_date = _git_info()
    pdf = _Report(run_ts=run_ts, git_hash=git_hash, commit_date=commit_date)

    # Load audit log
    try:
        audit_df = pd.read_csv(audit_csv)
    except Exception:
        audit_df = pd.DataFrame(
            columns=["timestamp", "module", "sample_id", "status", "description", "metric"]
        )
    detected = audit_df[audit_df["status"] == "DETECTED"].reset_index(drop=True)

    pr = pipeline_result  # alias

    n_clean = len(pr.get("clean_samples", []))
    n_total = pr.get(
        "n_total",
        n_clean
        + pr.get("n_qc_failed", 0)
        + pr.get("n_contaminated", 0)
        + pr.get("n_deduped", 0),
    )

    pdf.add_page()

    # ── Section 1: Executive Summary ────────────────────────────────────
    pdf.section("1. Executive Summary")

    n_qc = pr.get("n_qc_failed", 0)
    n_cont = pr.get("n_contaminated", 0)
    n_dup = pr.get("n_deduped", 0)
    n_sig = pr.get("n_sig_dmrs", 0)
    confound = pr.get("confounded", False)
    cv = pr.get("cramers_v", 0.0)
    auc = pr.get("mean_auc", 0.0)
    std_auc = pr.get("std_auc", 0.0)

    confound_str = (
        f"YES  (Cramer's V = {cv:.4f})" if confound else "No"
    )
    pdf.kv("Input samples", str(n_total))
    pdf.kv("Bisulfite/depth failures", str(n_qc), flagged=n_qc > 0)
    pdf.kv("Contaminated samples removed", str(n_cont), flagged=n_cont > 0)
    pdf.kv("Technical duplicates removed", str(n_dup), flagged=n_dup > 0)
    pdf.kv("Final clean samples", str(n_clean))
    pdf.kv("Batch confound detected", confound_str, flagged=confound)
    pdf.kv("Significant DMRs", str(n_sig), flagged=n_sig > 0)
    pdf.kv("Classification AUC (+/- SD)", f"{auc:.4f} +/- {std_auc:.4f}")
    pdf.kv("Git hash", git_hash)
    pdf.kv("Commit date", commit_date)
    pdf.kv("Run timestamp", run_ts)
    pdf.ln(3)

    # ── Section 2: Sample Exclusion Accounting ───────────────────────────
    pdf.add_page()
    pdf.section("2. Sample Exclusion Accounting")
    pdf.figure(
        os.path.join(FIGURES_DIR, "exclusion_accounting.png"),
        caption=(
            "Sample exclusion waterfall and pie breakdown. "
            "Left: proportional exclusion pie. Right: stepwise sample counts."
        ),
    )

    # ── Section 3: Global QC Metrics ────────────────────────────────────
    pdf.section("3. Global QC Metrics")
    pdf.figure(
        os.path.join(FIGURES_DIR, "qc_metrics.png"),
        caption=(
            "Per-sample non-CpG methylation rate, mean coverage, "
            "and bimodality coefficient. Flagged samples marked."
        ),
    )

    # ── Section 4: Beta Distribution ────────────────────────────────────
    pdf.add_page()
    pdf.section("4. Sample Purity -- Beta Distribution KDE")
    pdf.figure(
        os.path.join(FIGURES_DIR, "beta_distribution_kde.png"),
        caption=(
            "Kernel density estimate of beta values per sample. "
            "Contaminated samples show a mode near 0.5."
        ),
    )

    # ── Section 5: PCA ──────────────────────────────────────────────────
    pdf.add_page()
    pdf.section("5. Principal Component Analysis")
    pdf.figure(
        os.path.join(FIGURES_DIR, "pca_disease_label.png"),
        w=150,
        caption="PCA colored by disease label (Case / Control).",
    )
    pdf.figure(
        os.path.join(FIGURES_DIR, "pca_batch_id.png"),
        w=150,
        caption=(
            "PCA colored by batch ID. "
            "Overlap with disease-label coloring indicates potential confound."
        ),
    )

    # ── Section 6: Volcano Plots ─────────────────────────────────────────
    pdf.add_page()
    pdf.section("6. Volcano Plots -- Differential Methylation Regions")
    pdf.figure(
        os.path.join(FIGURES_DIR, "volcano.png"),
        w=148,
        caption=(
            "Volcano plot: x = delta_beta, y = -log10(p_adj). "
            "Dashed lines: |dBeta| = 0.10, p_adj = 0.05."
        ),
    )
    pdf.figure(
        os.path.join(FIGURES_DIR, "volcano_clonal_risk.png"),
        w=148,
        caption=(
            "Volcano with VDJ clonal-risk windows highlighted in orange. "
            "Analyst should review these windows before reporting."
        ),
    )

    # ── Section 7: DETECTED Events Log ──────────────────────────────────
    pdf.add_page()
    pdf.section("7. DETECTED Events Log")
    pdf.body(
        f"Total DETECTED entries: {len(detected)} "
        f"(source: {os.path.basename(audit_csv)})"
    )
    _detected_table(pdf, detected)

    # ── Section 8: Significant DMR Table ─────────────────────────────────
    pdf.section("8. Significant DMR Table")
    dmrs = pr.get("dmrs", pd.DataFrame())
    if "significant" in dmrs.columns:
        sig_dmrs = dmrs[dmrs["significant"]].reset_index(drop=True)
    else:
        sig_dmrs = pd.DataFrame()

    if len(sig_dmrs) == 0:
        pdf.body(
            "No significant DMRs detected in this run (expected for un-corrected "
            "mock data -- correct batch confound before re-running DMR Hunter)."
        )
    else:
        pdf.body(
            f"{len(sig_dmrs)} significant window(s). "
            "Amber rows indicate VDJ clonal-risk overlap -- review before reporting."
        )
        _sig_dmr_table(pdf, sig_dmrs)

    # ── Save ─────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    pdf.output(output_path)
    return output_path


# =============================================================================
# __main__
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="ImmuneMethylTools -- standalone PDF report generator."
    )
    parser.add_argument(
        "--audit",
        default=None,
        help="Path to audit_log CSV. Defaults to most-recent audit_log_pipeline_*.csv.",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PDF path. Defaults to data/report_<timestamp>.pdf.",
    )
    args = parser.parse_args()

    # Resolve audit CSV
    if args.audit:
        audit_csv = args.audit
    else:
        pattern = os.path.join(project_root(), "data", "audit_log_pipeline_*.csv")
        candidates = sorted(glob.glob(pattern))
        if not candidates:
            raise FileNotFoundError("No audit_log_pipeline_*.csv found in data/. Run the pipeline first.")
        audit_csv = candidates[-1]

    ts_tag = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = args.output or data_path(f"report_{ts_tag}.pdf")
    run_ts = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    print(f"[report_gen] Generating report from: {audit_csv}")
    out = generate_report(
        pipeline_result={},
        audit_csv=audit_csv,
        output_path=output_path,
        run_ts=run_ts,
    )
    print(f"[report_gen] Report saved to: {out}")
