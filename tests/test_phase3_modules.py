"""
tests/test_phase3_modules.py — ImmuneMethylTools Phase 3 Module Tests
======================================================================
Verifies the detection logic of all eight core/ detective modules
against data/mock_methylation.csv.

Run as pytest:
    pytest tests/test_phase3_modules.py -v

Run as standalone script (prints timestamped detections):
    python tests/test_phase3_modules.py
"""

import os
import sys
import tempfile
from datetime import datetime

import numpy as np
import pandas as pd
import pytest

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH  = os.path.join(REPO_ROOT, "data", "mock_methylation.csv")
sys.path.insert(0, REPO_ROOT)

from core.deconvolution import detect_lineage_shift, estimate_cell_fractions
from core.dmr_hunter import MIN_CPGS, find_dmrs
from core.ml_guard import N_TOP_CPGS, run_safe_model
from core.normalizer import check_confounding, robust_normalize
from core.qc_guard import audit_quality, detect_contamination
from core.repertoire_clonality import flag_clonal_artifacts, get_vdj_summary
from core.sample_audit import DUP_CORR_THRESH, detect_duplicates
from core.visuals import plot_beta_distribution, plot_pca, plot_qc_metrics


# ── Fixture ───────────────────────────────────────────────────────────────────

def load_data() -> pd.DataFrame:
    assert os.path.exists(CSV_PATH), (
        f"mock_methylation.csv not found at {CSV_PATH}. "
        "Run: python data/generate_mock_data.py"
    )
    return pd.read_csv(CSV_PATH)


# =============================================================================
# visuals.py
# =============================================================================

def test_plot_qc_metrics_creates_file():
    """plot_qc_metrics saves a non-empty PNG to the specified path."""
    df = load_data()
    with tempfile.TemporaryDirectory() as tmpdir:
        path   = os.path.join(tmpdir, "qc_metrics.png")
        result = plot_qc_metrics(df, save_path=path)
        assert result == path, "Return value should be the save path"
        assert os.path.exists(path), "PNG file was not created"
        assert os.path.getsize(path) > 0, "PNG file is empty"


def test_plot_beta_distribution_creates_file():
    """plot_beta_distribution saves a non-empty PNG to the specified path."""
    df = load_data()
    with tempfile.TemporaryDirectory() as tmpdir:
        path   = os.path.join(tmpdir, "beta_kde.png")
        result = plot_beta_distribution(df, save_path=path)
        assert result == path
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0


def test_plot_pca_creates_file():
    """plot_pca saves a non-empty PNG coloured by disease_label."""
    df = load_data()
    with tempfile.TemporaryDirectory() as tmpdir:
        path   = os.path.join(tmpdir, "pca.png")
        result = plot_pca(df, color_by="disease_label", save_path=path)
        assert result == path
        assert os.path.exists(path)
        assert os.path.getsize(path) > 0


# =============================================================================
# qc_guard.py
# =============================================================================

def test_audit_quality_excludes_bisulfite_failures():
    """
    S001 and S002 (non_cpg_meth_rate ≈ 5 %) must be excluded from the
    clean_samples_list; all other samples should pass coverage and conversion QC.
    """
    clean = audit_quality(load_data())
    assert "S001" not in clean, "S001 (bisulfite failure) should be excluded"
    assert "S002" not in clean, "S002 (bisulfite failure) should be excluded"


def test_audit_quality_clean_count():
    """
    41 samples total; S001 and S002 fail bisulfite QC → 39 clean samples.
    """
    clean = audit_quality(load_data())
    assert len(clean) == 39, (
        f"Expected 39 clean samples, got {len(clean)}"
    )


def test_detect_contamination_flags_s020():
    """
    S020 (50 % contamination mix) must be identified by the cohort-relative
    bimodality coefficient check AND its mean beta in the muddy range [0.40, 0.65].
    """
    flagged, _ = detect_contamination(load_data())
    assert "S020" in flagged, (
        f"S020 (contaminated) not detected. Flagged: {flagged}"
    )


def test_detect_contamination_report_structure():
    """Report DataFrame must contain all expected diagnostic columns."""
    _, report = detect_contamination(load_data())
    expected_cols = [
        "bimodality_coeff", "mean_beta",
        "low_bimodality", "muddy_mean", "contamination_flag",
    ]
    for col in expected_cols:
        assert col in report.columns, f"Missing column: {col}"


def test_detect_contamination_s020_has_lower_bc():
    """
    S020's bimodality coefficient must be below the cohort median, confirming
    the distribution collapse from contamination.
    """
    _, report = detect_contamination(load_data())
    s020_bc   = report.loc["S020", "bimodality_coeff"]
    cohort_bc = report["bimodality_coeff"].median()
    assert s020_bc < cohort_bc, (
        f"S020 BC ({s020_bc:.4f}) should be below cohort median ({cohort_bc:.4f})"
    )


# =============================================================================
# sample_audit.py
# =============================================================================

def test_detect_duplicates_flags_s010_sdup():
    """
    S010 and S_DUP (Pearson r ≈ 1.000) must appear as a flagged duplicate pair
    in the top-100 high-variance CpG correlation analysis.
    """
    result  = detect_duplicates(load_data())
    flagged = result[result["duplicate_flag"]]
    pairs   = set(
        frozenset({row.sample_a, row.sample_b}) for _, row in flagged.iterrows()
    )
    assert frozenset({"S010", "S_DUP"}) in pairs, (
        f"S010 ↔ S_DUP not flagged as duplicates. Flagged pairs: {pairs}"
    )


def test_detect_duplicates_only_one_pair_flagged():
    """
    Only the intentional duplicate (S010 ↔ S_DUP) should exceed r > 0.99;
    no spurious pairs should be flagged.
    """
    result = detect_duplicates(load_data())
    n_flagged = result["duplicate_flag"].sum()
    assert n_flagged == 1, (
        f"Expected exactly 1 flagged pair, got {n_flagged}. "
        f"Check DUP_CORR_THRESH ({DUP_CORR_THRESH}) or baseline noise level."
    )


def test_detect_duplicates_output_columns():
    """Result DataFrame must contain all expected columns."""
    result = detect_duplicates(load_data())
    for col in ["sample_a", "sample_b", "pearson_r", "duplicate_flag"]:
        assert col in result.columns, f"Missing column: {col}"


def test_detect_duplicates_sorted_descending():
    """Result must be sorted by pearson_r descending (S010/S_DUP pair first)."""
    result = detect_duplicates(load_data())
    assert result["pearson_r"].is_monotonic_decreasing, (
        "detect_duplicates result must be sorted by pearson_r descending"
    )


# =============================================================================
# normalizer.py
# =============================================================================

def test_check_confounding_detects_batch_disease():
    """
    Batch_01 is 81 % Case — strong association (Cramér's V > 0.5) must be
    flagged so the analyst knows naïve batch correction will absorb disease signal.
    """
    result = check_confounding(load_data(), "batch_id", "disease_label")
    assert result["confounded"] is True, (
        f"Batch × Disease confounding not detected. Cramér's V = {result['cramers_v']:.4f}"
    )
    assert result["cramers_v"] > 0.50, (
        f"Cramér's V = {result['cramers_v']:.4f} (expected > 0.50)"
    )


def test_check_confounding_result_keys():
    """Result dict must contain all expected keys."""
    result = check_confounding(load_data(), "batch_id", "disease_label")
    for key in ["chi2", "p_value", "cramers_v", "confounded", "contingency_table"]:
        assert key in result, f"Missing key: {key}"


def test_check_confounding_p_value_significant():
    """Chi-square p-value must be significant (< 0.05) for the known confound."""
    result = check_confounding(load_data(), "batch_id", "disease_label")
    assert result["p_value"] < 0.05, (
        f"Expected significant p-value, got p = {result['p_value']:.4f}"
    )


def test_robust_normalize_adds_beta_normalized_column():
    """robust_normalize must add a 'beta_normalized' column to the DataFrame."""
    df_norm = robust_normalize(load_data(), save_figure=False)
    assert "beta_normalized" in df_norm.columns, (
        "Column 'beta_normalized' not found after normalisation"
    )


def test_robust_normalize_centers_sample_medians():
    """
    After median-centring, each sample's median beta_normalized must be
    exactly zero (within floating-point tolerance).
    """
    df_norm       = robust_normalize(load_data(), save_figure=False)
    sample_medians = df_norm.groupby("sample_id")["beta_normalized"].median()
    worst          = sample_medians.abs().max()
    assert worst < 1e-8, (
        f"Median-centring failed: max |median| = {worst:.2e} (expected < 1e-8)"
    )


def test_robust_normalize_preserves_original_beta():
    """robust_normalize must not modify the original 'beta_value' column."""
    df      = load_data()
    df_norm = robust_normalize(df, save_figure=False)
    pd.testing.assert_series_equal(
        df["beta_value"].reset_index(drop=True),
        df_norm["beta_value"].reset_index(drop=True),
        check_names=False,
        obj="beta_value column after normalize",
    )


# =============================================================================
# repertoire_clonality.py
# =============================================================================

def test_flag_clonal_detects_p001_sample():
    """
    P001 (S001) is the injected clonal patient; their sample must appear
    in the flagged_samples list from flag_clonal_artifacts.
    """
    df = load_data()
    _, flagged = flag_clonal_artifacts(df)
    p001_samples = df[df["patient_id"] == "P001"]["sample_id"].unique().tolist()
    overlap = [s for s in p001_samples if s in flagged]
    assert overlap, (
        f"No P001 samples found in clonal flagged list. "
        f"P001 samples: {p001_samples}, flagged: {flagged}"
    )


def test_clonal_rows_meet_dual_criteria():
    """
    Every flagged row must simultaneously satisfy:
      is_vdj_region == True  AND  beta_value > 0.80  AND  fragment_length > 180 bp
    """
    clonal_rows, _ = flag_clonal_artifacts(load_data())
    assert len(clonal_rows) > 0, "No clonal rows detected"
    assert clonal_rows["is_vdj_region"].astype(bool).all(), "Non-VDJ row in clonal results"
    assert (clonal_rows["beta_value"] > 0.80).all(), "Clonal row with beta ≤ 0.80"
    assert (clonal_rows["fragment_length"] > 180).all(), "Clonal row with fragment ≤ 180 bp"


def test_vdj_summary_p001_highest_clonal_hits():
    """
    get_vdj_summary must rank P001 at the top (highest clonal_hits) with
    more than zero hits — confirming the clonal signal is captured.
    """
    summary = get_vdj_summary(load_data())
    assert summary.iloc[0]["patient_id"] == "P001", (
        f"Expected P001 at top of VDJ summary, got {summary.iloc[0]['patient_id']}"
    )
    assert summary.iloc[0]["clonal_hits"] > 0, "P001 clonal_hits should be > 0"


def test_vdj_summary_non_clonal_patients_zero_hits():
    """
    All patients other than P001 must have zero clonal hits — the artefact
    must not have leaked into the baseline.
    """
    summary = get_vdj_summary(load_data())
    others  = summary[summary["patient_id"] != "P001"]
    assert (others["clonal_hits"] == 0).all(), (
        f"Unexpected clonal hits in non-P001 patients:\n"
        f"{others[others['clonal_hits'] > 0][['patient_id','clonal_hits']]}"
    )


# =============================================================================
# deconvolution.py
# =============================================================================

def test_cell_fractions_sum_to_one():
    """
    All four fractions (B, T, Treg, other) must sum to 1.0 for every sample
    (within floating-point tolerance).
    """
    df    = load_data()
    fracs = estimate_cell_fractions(df)
    cols  = ["b_fraction", "t_fraction", "treg_fraction", "other_fraction"]
    row_sums = fracs[cols].sum(axis=1)
    worst    = (row_sums - 1.0).abs().max()
    # Fractions are rounded to 4 d.p.; max cumulative rounding error = 4 × 5e-5 = 2e-4
    assert worst < 1e-3, (
        f"Cell fractions do not sum to 1.0 (max deviation = {worst:.2e})"
    )


def test_cell_fractions_covers_all_samples():
    """estimate_cell_fractions must return one row per sample."""
    df    = load_data()
    fracs = estimate_cell_fractions(df)
    assert len(fracs) == df["sample_id"].nunique(), (
        f"Expected {df['sample_id'].nunique()} fraction rows, got {len(fracs)}"
    )


def test_cell_fractions_all_non_negative():
    """All fraction values must be in [0, 1]."""
    fracs = estimate_cell_fractions(load_data())
    cols  = ["b_fraction", "t_fraction", "treg_fraction", "other_fraction"]
    for col in cols:
        assert (fracs[col] >= 0).all(), f"Negative values in {col}"
        assert (fracs[col] <= 1).all(), f"Values > 1 in {col}"


def test_lineage_shift_output_structure():
    """detect_lineage_shift must return a DataFrame with all required columns."""
    shifts = detect_lineage_shift(load_data())
    for col in [
        "sample_id", "foxp3_mean_beta", "pax5_mean_beta",
        "treg_flag", "bcell_shift_flag", "any_lineage_flag",
    ]:
        assert col in shifts.columns, f"Missing column in lineage report: {col}"


def test_lineage_shift_covers_all_samples():
    """detect_lineage_shift must return exactly one row per sample."""
    df     = load_data()
    shifts = detect_lineage_shift(df)
    assert len(shifts) == df["sample_id"].nunique(), (
        f"Expected {df['sample_id'].nunique()} rows, got {len(shifts)}"
    )


# =============================================================================
# dmr_hunter.py
# =============================================================================

def test_dmr_hunter_safety_assertion_triggers():
    """
    Passing a DataFrame that contains samples NOT in clean_samples_list must
    raise AssertionError with 'SAFETY VIOLATION' in the message.
    """
    df            = load_data()
    clean_samples = audit_quality(df)  # excludes S001, S002
    # df still contains S001, S002 → should trigger the safety guard
    with pytest.raises(AssertionError, match="SAFETY VIOLATION"):
        find_dmrs(df, clean_samples)


def test_dmr_hunter_output_columns():
    """find_dmrs must return a DataFrame with all required columns."""
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = find_dmrs(df_clean, clean_samples)
    for col in [
        "window_id", "cpgs", "n_cpgs",
        "case_mean", "ctrl_mean", "delta_beta",
        "wilcoxon_stat", "p_value", "p_adj", "significant",
    ]:
        assert col in result.columns, f"Missing column: {col}"


def test_dmr_hunter_excludes_vdj_cpgs_from_windows():
    """
    VDJ-region CpGs must not appear in any window (clonal artefact exclusion).
    Clonal CpGs with high beta would otherwise inflate Case means and generate
    spurious DMR calls.
    """
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = find_dmrs(df_clean, clean_samples)
    vdj_cpgs      = set(df[df["is_vdj_region"].astype(bool)]["cpg_id"].unique())

    for window_str in result["cpgs"]:
        for cpg in window_str.split(","):
            assert cpg not in vdj_cpgs, (
                f"VDJ CpG {cpg} found in DMR window — clonal exclusion failed"
            )


def test_dmr_hunter_padj_bounded():
    """All BH-adjusted p-values must lie in [0, 1]."""
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = find_dmrs(df_clean, clean_samples)
    assert (result["p_adj"] >= 0.0).all(), "p_adj contains values < 0"
    assert (result["p_adj"] <= 1.0).all(), "p_adj contains values > 1"


def test_dmr_hunter_significant_windows_meet_all_criteria():
    """
    Every window marked 'significant=True' must pass all three DMR criteria:
    p_adj < 0.05, |ΔBeta| > 0.10, n_cpgs ≥ MIN_CPGS.
    """
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = find_dmrs(df_clean, clean_samples)
    sig           = result[result["significant"]]

    if len(sig):
        assert (sig["p_adj"] < 0.05).all(),              "Significant window with p_adj ≥ 0.05"
        assert (sig["delta_beta"].abs() > 0.10).all(),    "Significant window with |ΔBeta| ≤ 0.10"
        assert (sig["n_cpgs"] >= MIN_CPGS).all(),         "Significant window with n_cpgs < MIN_CPGS"


# =============================================================================
# ml_guard.py
# =============================================================================

def test_ml_guard_output_keys():
    """run_safe_model must return a dict with all required keys."""
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = run_safe_model(df_clean)
    for key in [
        "cv_results", "mean_auc", "std_auc",
        "mean_accuracy", "n_samples", "n_features",
    ]:
        assert key in result, f"Missing key in run_safe_model output: {key}"


def test_ml_guard_auc_in_valid_range():
    """mean_auc must be a valid probability: in [0, 1]."""
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = run_safe_model(df_clean)
    assert 0.0 <= result["mean_auc"] <= 1.0, (
        f"mean_auc = {result['mean_auc']:.4f} is outside [0, 1]"
    )


def test_ml_guard_n_samples_matches_clean_count():
    """n_samples must equal the number of clean samples passed in."""
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = run_safe_model(df_clean)
    assert result["n_samples"] == len(clean_samples), (
        f"n_samples={result['n_samples']} but len(clean_samples)={len(clean_samples)}"
    )


def test_ml_guard_n_features_bounded():
    """n_features must not exceed N_TOP_CPGS (the top-variance CpG cap)."""
    df            = load_data()
    clean_samples = audit_quality(df)
    df_clean      = df[df["sample_id"].isin(clean_samples)]
    result        = run_safe_model(df_clean)
    assert result["n_features"] <= N_TOP_CPGS, (
        f"n_features={result['n_features']} exceeds N_TOP_CPGS={N_TOP_CPGS}"
    )


# =============================================================================
# __main__ — timestamped standalone detection log
# =============================================================================

if __name__ == "__main__":
    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def log(status: str, module: str, test: str, detail: str):
        icon = "✓ PASS    " if status == "ok" else "✗ FAIL    "
        print(f"[{ts()}] [{module:15s}] {icon} | {test} | {detail}")

    def run(fn, module, test, detail_fn):
        try:
            fn()
            log("ok", module, test, detail_fn())
        except Exception as exc:
            log("fail", module, test, str(exc)[:120])

    print(f"[{ts()}] Loading {CSV_PATH}")
    df = load_data()
    print(f"[{ts()}] Shape: {df.shape}  |  Samples: {df['sample_id'].nunique()}")
    print()

    # ── visuals ───────────────────────────────────────────────────────────────
    with tempfile.TemporaryDirectory() as tmpdir:
        run(test_plot_qc_metrics_creates_file,       "VISUALS",  "plot_qc_metrics",        lambda: "PNG saved")
        run(test_plot_beta_distribution_creates_file,"VISUALS",  "plot_beta_distribution", lambda: "PNG saved")
        run(test_plot_pca_creates_file,              "VISUALS",  "plot_pca",               lambda: "PNG saved")

    # ── qc_guard ──────────────────────────────────────────────────────────────
    clean = audit_quality(df)
    run(test_audit_quality_excludes_bisulfite_failures, "QC_GUARD", "audit_quality/exclusions",
        lambda: f"S001, S002 excluded from n={len(clean)} clean samples")
    run(test_audit_quality_clean_count,                 "QC_GUARD", "audit_quality/count",
        lambda: f"clean_samples n={len(clean)}")
    flagged_c, report = detect_contamination(df)
    run(test_detect_contamination_flags_s020,           "QC_GUARD", "detect_contamination/S020",
        lambda: f"flagged={flagged_c}  S020 BC={report.loc['S020','bimodality_coeff']:.4f}")
    run(test_detect_contamination_report_structure,     "QC_GUARD", "detect_contamination/columns",
        lambda: "all columns present")
    run(test_detect_contamination_s020_has_lower_bc,    "QC_GUARD", "detect_contamination/BC_rank",
        lambda: f"S020 BC < cohort median ({report['bimodality_coeff'].median():.4f})")

    # ── sample_audit ──────────────────────────────────────────────────────────
    dup_result = detect_duplicates(df)
    run(test_detect_duplicates_flags_s010_sdup,   "SAMPLE_AUDIT", "detect_duplicates/pair",
        lambda: f"S010↔S_DUP r={dup_result.iloc[0].pearson_r:.6f}")
    run(test_detect_duplicates_only_one_pair_flagged, "SAMPLE_AUDIT", "detect_duplicates/count",
        lambda: f"n_flagged={dup_result['duplicate_flag'].sum()}")
    run(test_detect_duplicates_output_columns,    "SAMPLE_AUDIT", "detect_duplicates/columns",
        lambda: "all columns present")
    run(test_detect_duplicates_sorted_descending, "SAMPLE_AUDIT", "detect_duplicates/sort",
        lambda: "sorted by pearson_r descending")

    # ── normalizer ────────────────────────────────────────────────────────────
    conf = check_confounding(df, "batch_id", "disease_label")
    run(test_check_confounding_detects_batch_disease, "NORMALIZER", "check_confounding/detected",
        lambda: f"Cramér's V={conf['cramers_v']:.4f}  p={conf['p_value']:.4e}")
    run(test_check_confounding_result_keys,           "NORMALIZER", "check_confounding/keys",
        lambda: "all keys present")
    run(test_check_confounding_p_value_significant,   "NORMALIZER", "check_confounding/p_value",
        lambda: f"p={conf['p_value']:.4e} < 0.05")
    df_norm = robust_normalize(df, save_figure=False)
    run(test_robust_normalize_adds_beta_normalized_column, "NORMALIZER", "robust_normalize/column",
        lambda: "beta_normalized added")
    sample_meds = df_norm.groupby("sample_id")["beta_normalized"].median()
    run(test_robust_normalize_centers_sample_medians, "NORMALIZER", "robust_normalize/centers",
        lambda: f"max |median| = {sample_meds.abs().max():.2e}")
    run(test_robust_normalize_preserves_original_beta,"NORMALIZER", "robust_normalize/preserves",
        lambda: "beta_value unchanged")

    # ── repertoire_clonality ──────────────────────────────────────────────────
    clonal_rows, clonal_flagged = flag_clonal_artifacts(df)
    run(test_flag_clonal_detects_p001_sample,           "CLONALITY", "flag_clonal/P001",
        lambda: f"flagged_samples={clonal_flagged}  n_rows={len(clonal_rows)}")
    run(test_clonal_rows_meet_dual_criteria,            "CLONALITY", "flag_clonal/criteria",
        lambda: f"all rows: VDJ=True, β>0.80, frag>180bp")
    vdj_sum = get_vdj_summary(df)
    run(test_vdj_summary_p001_highest_clonal_hits,      "CLONALITY", "vdj_summary/P001_top",
        lambda: f"P001 clonal_hits={vdj_sum.iloc[0]['clonal_hits']}")
    run(test_vdj_summary_non_clonal_patients_zero_hits, "CLONALITY", "vdj_summary/others_zero",
        lambda: "all non-P001 patients: clonal_hits=0")

    # ── deconvolution ─────────────────────────────────────────────────────────
    fracs  = estimate_cell_fractions(df)
    run(test_cell_fractions_sum_to_one,        "DECONVOLVE", "estimate_fractions/sum",
        lambda: f"max |1-sum| = {(fracs[['b_fraction','t_fraction','treg_fraction','other_fraction']].sum(axis=1)-1.0).abs().max():.2e} (tol 1e-3)")
    run(test_cell_fractions_covers_all_samples,"DECONVOLVE", "estimate_fractions/count",
        lambda: f"n_rows={len(fracs)}")
    run(test_cell_fractions_all_non_negative,  "DECONVOLVE", "estimate_fractions/bounds",
        lambda: "all fractions in [0,1]")
    shifts = detect_lineage_shift(df)
    run(test_lineage_shift_output_structure,   "DECONVOLVE", "lineage_shift/columns",
        lambda: "all columns present")
    run(test_lineage_shift_covers_all_samples, "DECONVOLVE", "lineage_shift/count",
        lambda: f"n_rows={len(shifts)}")

    # ── dmr_hunter ────────────────────────────────────────────────────────────
    df_clean = df[df["sample_id"].isin(clean)]
    dmrs     = find_dmrs(df_clean, clean)
    sig_dmrs = dmrs[dmrs["significant"]]
    run(test_dmr_hunter_safety_assertion_triggers,       "DMR_HUNTER", "find_dmrs/safety",
        lambda: "AssertionError raised for dirty input")
    run(test_dmr_hunter_output_columns,                  "DMR_HUNTER", "find_dmrs/columns",
        lambda: "all columns present")
    run(test_dmr_hunter_excludes_vdj_cpgs_from_windows,  "DMR_HUNTER", "find_dmrs/vdj_excluded",
        lambda: f"VDJ CpGs excluded from all {len(dmrs)} windows")
    run(test_dmr_hunter_padj_bounded,                    "DMR_HUNTER", "find_dmrs/padj_range",
        lambda: f"p_adj in [0,1] across {len(dmrs)} windows")
    run(test_dmr_hunter_significant_windows_meet_all_criteria, "DMR_HUNTER", "find_dmrs/criteria",
        lambda: f"n_significant={len(sig_dmrs)}")

    # ── ml_guard ──────────────────────────────────────────────────────────────
    ml_result = run_safe_model(df_clean)
    run(test_ml_guard_output_keys,                "ML_GUARD", "run_safe_model/keys",
        lambda: "all keys present")
    run(test_ml_guard_auc_in_valid_range,         "ML_GUARD", "run_safe_model/auc_range",
        lambda: f"AUC={ml_result['mean_auc']:.4f} ± {ml_result['std_auc']:.4f}")
    run(test_ml_guard_n_samples_matches_clean_count, "ML_GUARD", "run_safe_model/n_samples",
        lambda: f"n_samples={ml_result['n_samples']} == len(clean)={len(clean)}")
    run(test_ml_guard_n_features_bounded,         "ML_GUARD", "run_safe_model/n_features",
        lambda: f"n_features={ml_result['n_features']} ≤ N_TOP_CPGS={N_TOP_CPGS}")

    print()
    print(f"[{ts()}] Done. Run full suite: pytest tests/test_phase3_modules.py -v")
