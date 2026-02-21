"""
tests/test_io_utils_and_pipeline.py — io_utils safe loader + pipeline tests
============================================================================
Covers the two modules added after Phase 3:

  core/io_utils.py  — project_root, data_path, load_methylation validation
  core/pipeline.py  — run_pipeline end-to-end results

The pipeline fixture is module-scoped so run_pipeline() executes exactly
once regardless of how many pipeline tests run.

Run as pytest:
    pytest tests/test_io_utils_and_pipeline.py -v

Run as standalone script (io_utils tests only — pipeline tests skipped for speed):
    python tests/test_io_utils_and_pipeline.py
"""

import os
import sys
import tempfile
from datetime import datetime

import pandas as pd
import pytest

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(REPO_ROOT, "data", "mock_methylation.csv")
sys.path.insert(0, REPO_ROOT)

from core.infrastructure.io_utils import REQUIRED_COLUMNS, data_path, load_methylation, project_root  # noqa: E402
from core.orchestration.pipeline import run_pipeline  # noqa: E402


# =============================================================================
# Helpers
# =============================================================================

def _minimal_valid_df() -> pd.DataFrame:
    """Two-row DataFrame that satisfies the full CLAUDE.md schema (16 columns)."""
    return pd.DataFrame({
        "sample_id": ["S001", "S001"],
        "patient_id": ["P001", "P001"],
        "batch_id": ["Batch_01", "Batch_01"],
        "age": [40, 40],
        "disease_label": ["Case", "Case"],
        "cpg_id": ["cg00000001", "cg00000002"],
        "beta_value": [0.2, 0.8],
        "depth": [30, 25],
        "fragment_length": [150, 160],
        "is_vdj_region": [False, False],
        "non_cpg_meth_rate": [0.004, 0.003],
        "sex": ["F", "F"],
        "is_x_chromosome": [False, True],
        "gc_content": [0.45, 0.55],
        "chrom": ["chr1", "chrX"],
        "pos": [1000000, 50000000],
    })


def _write_csv(df: pd.DataFrame, directory: str, name: str = "test.csv") -> str:
    path = os.path.join(directory, name)
    df.to_csv(path, index=False)
    return path


# =============================================================================
# io_utils — project_root / data_path
# =============================================================================

def test_project_root_is_parent_of_core():
    """project_root() must return the project root, not core/ itself."""
    root = project_root()
    assert os.path.isdir(root), \
        f"project_root() is not a directory: {root}"
    assert os.path.isdir(os.path.join(root, "core")), \
        f"project_root() has no core/ subdirectory: {root}"
    assert not root.rstrip(os.sep).endswith("core"), \
        "project_root() returned the core/ directory itself, not its parent"


def test_data_path_resolves_correctly():
    """data_path() must return an absolute path ending in data/<filename>."""
    result = data_path("foo.csv")
    assert os.path.isabs(result), \
        f"data_path() returned a relative path: {result}"
    assert result.endswith(os.path.join("data", "foo.csv")), \
        f"data_path('foo.csv') = {result} — expected to end with data/foo.csv"


# =============================================================================
# io_utils — load_methylation happy path
# =============================================================================

def test_load_methylation_valid_file():
    """Happy path: returns a DataFrame with the correct shape and all required columns."""
    assert os.path.exists(CSV_PATH), \
        "Generate mock data first: python data/generate_mock_data.py"
    df = load_methylation(CSV_PATH, verbose=False)
    assert isinstance(df, pd.DataFrame), "load_methylation must return a DataFrame"
    assert df.shape[0] > 0, "Loaded DataFrame is empty"
    for col in REQUIRED_COLUMNS:
        assert col in df.columns, f"Missing column after load: {col}"


# =============================================================================
# io_utils — load_methylation validation failures
# =============================================================================

def test_load_methylation_file_not_found():
    """Raises FileNotFoundError for a path that does not exist."""
    with pytest.raises(FileNotFoundError, match="file not found"):
        load_methylation("/nonexistent/path/data.csv", verbose=False)


def test_load_methylation_missing_column():
    """Raises ValueError when a required column is absent from the CSV."""
    with tempfile.TemporaryDirectory() as tmpdir:
        df = _minimal_valid_df().drop(columns=["beta_value"])
        path = _write_csv(df, tmpdir)
        with pytest.raises(ValueError, match="missing required column"):
            load_methylation(path, verbose=False)


def test_load_methylation_beta_out_of_range():
    """Raises ValueError when any beta_value exceeds 1."""
    with tempfile.TemporaryDirectory() as tmpdir:
        df = _minimal_valid_df()
        df.loc[0, "beta_value"] = 1.5
        path = _write_csv(df, tmpdir)
        with pytest.raises(ValueError, match="beta_value"):
            load_methylation(path, verbose=False)


def test_load_methylation_bad_disease_label():
    """Raises ValueError when disease_label contains an unrecognized value."""
    with tempfile.TemporaryDirectory() as tmpdir:
        df = _minimal_valid_df()
        df.loc[0, "disease_label"] = "Unknown"
        path = _write_csv(df, tmpdir)
        with pytest.raises(ValueError, match="disease_label"):
            load_methylation(path, verbose=False)


def test_load_methylation_negative_depth():
    """Raises ValueError when depth contains a negative value."""
    with tempfile.TemporaryDirectory() as tmpdir:
        df = _minimal_valid_df()
        df.loc[0, "depth"] = -1
        path = _write_csv(df, tmpdir)
        with pytest.raises(ValueError, match="depth"):
            load_methylation(path, verbose=False)


# =============================================================================
# pipeline — run_pipeline (module-scoped fixture: runs once for all tests below)
# =============================================================================

@pytest.fixture(scope="module")
def pipeline_result():
    """Execute run_pipeline once and share the result across all pipeline tests."""
    assert os.path.exists(CSV_PATH), \
        "Generate mock data first: python data/generate_mock_data.py"
    return run_pipeline(CSV_PATH, save_figures=False)


def test_pipeline_returns_expected_keys(pipeline_result):
    """run_pipeline must return a dict containing all documented output keys."""
    expected = [
        "clean_samples", "n_qc_failed", "n_contaminated", "n_sex_flagged", "n_deduped",
        "confounded", "cramers_v", "n_clonal_rows",
        "dmrs", "n_sig_dmrs", "mean_auc", "std_auc", "df_norm", "clean_csv",
    ]
    for key in expected:
        assert key in pipeline_result, f"Missing key in run_pipeline output: '{key}'"


def test_pipeline_clean_samples_count(pipeline_result):
    """
    101 total - 3 QC failures - 1 contaminated - 2 sex mixup - 1 duplicate = 94 clean samples.
    """
    n = len(pipeline_result["clean_samples"])
    assert n == 94, f"Expected 94 clean samples, got {n}"


def test_pipeline_n_sex_flagged(pipeline_result):
    """Exactly 2 samples must be flagged as sex-signal mismatches: S035 and S036."""
    assert pipeline_result["n_sex_flagged"] == 2, (
        f"Expected 2 sex mixups (S035, S036), got {pipeline_result['n_sex_flagged']}"
    )


def test_pipeline_n_qc_failed(pipeline_result):
    """Exactly 3 samples must fail bisulfite/depth QC: S001, S002, S030."""
    assert pipeline_result["n_qc_failed"] == 3, (
        f"Expected 3 QC failures (S001/S002 bisulfite, S030 depth), "
        f"got {pipeline_result['n_qc_failed']}"
    )


def test_pipeline_n_contaminated(pipeline_result):
    """Exactly 1 sample must be flagged as contaminated: S020."""
    assert pipeline_result["n_contaminated"] == 1, (
        f"Expected 1 contaminated sample (S020), "
        f"got {pipeline_result['n_contaminated']}"
    )


def test_pipeline_n_deduped(pipeline_result):
    """Exactly 1 sample must be dropped as a technical duplicate: S_DUP."""
    assert pipeline_result["n_deduped"] == 1, (
        f"Expected 1 duplicate dropped (S_DUP), "
        f"got {pipeline_result['n_deduped']}"
    )


def test_pipeline_confounded(pipeline_result):
    """Batch × disease confound must be detected with Cramér's V > 0.5."""
    assert pipeline_result["confounded"] is True, (
        f"Expected confounded=True; "
        f"Cramér's V = {pipeline_result['cramers_v']:.4f}"
    )
    assert pipeline_result["cramers_v"] > 0.5, (
        f"Cramér's V = {pipeline_result['cramers_v']:.4f}, expected > 0.5"
    )


def test_pipeline_exports_clean_csv(pipeline_result):
    """
    clean_methylation.csv must exist, contain 94 samples, and have all
    16 CLAUDE.md schema columns.
    """
    clean_csv = pipeline_result["clean_csv"]
    assert os.path.isfile(clean_csv), \
        f"clean_methylation.csv not found at {clean_csv}"
    df = pd.read_csv(clean_csv)
    assert df["sample_id"].nunique() == 94, (
        f"Expected 94 samples in clean_methylation.csv, "
        f"got {df['sample_id'].nunique()}"
    )
    for col in REQUIRED_COLUMNS:
        assert col in df.columns, \
            f"Column '{col}' missing from clean_methylation.csv"


# =============================================================================
# __main__ — standalone run (io_utils tests only; pipeline skipped for speed)
# =============================================================================

if __name__ == "__main__":
    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def run(fn, tag, detail_fn):
        try:
            fn()
            print(f"[{ts()}] [{tag:20s}] PASS | {detail_fn()}")
        except Exception as exc:
            print(f"[{ts()}] [{tag:20s}] FAIL | {str(exc)[:120]}")

    print(f"[{ts()}] io_utils tests (standalone — pipeline tests require pytest)")
    print()

    run(test_project_root_is_parent_of_core, "IO_UTILS/project_root",
        lambda: f"root={project_root()}")
    run(test_data_path_resolves_correctly, "IO_UTILS/data_path",
        lambda: f"data_path('x.csv')={data_path('x.csv')}")
    run(test_load_methylation_valid_file, "IO_UTILS/load_valid",
        lambda: f"loaded {CSV_PATH}")
    run(test_load_methylation_file_not_found, "IO_UTILS/load_missing",
        lambda: "FileNotFoundError raised correctly")
    run(test_load_methylation_missing_column, "IO_UTILS/load_no_col",
        lambda: "ValueError raised for missing column")
    run(test_load_methylation_beta_out_of_range, "IO_UTILS/load_bad_beta",
        lambda: "ValueError raised for beta > 1")
    run(test_load_methylation_bad_disease_label, "IO_UTILS/load_bad_label",
        lambda: "ValueError raised for unknown disease_label")
    run(test_load_methylation_negative_depth, "IO_UTILS/load_neg_depth",
        lambda: "ValueError raised for depth < 0")

    print()
    print(f"[{ts()}] Pipeline tests require pytest (module-scoped fixture):")
    print(f"[{ts()}]   pytest tests/test_io_utils_and_pipeline.py -v -k pipeline")
