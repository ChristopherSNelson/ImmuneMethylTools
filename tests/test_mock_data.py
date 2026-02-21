"""
tests/test_mock_data.py — ImmuneMethylTools Phase 2 Sanity Checks
==================================================================
Verifies that all five stumper artifacts are correctly embedded in
data/mock_methylation.csv.

Run as pytest:
    pytest tests/test_mock_data.py -v

Run as standalone script (prints timestamped detections):
    python tests/test_mock_data.py
"""

import os
from datetime import datetime

import pandas as pd
from scipy.stats import pearsonr

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(REPO_ROOT, "data", "mock_methylation.csv")

# ── Thresholds (mirror generate_mock_data.py constants) ───────────────────────
BATCH_SHIFT_MIN = 0.05   # minimum expected delta between Batch_01 Case vs Control
CLONAL_BETA_MIN = 0.80
CLONAL_FRAG_MIN = 195   # baseline clips at 220; clonal injected at 200–260 bp
BISULFITE_FAIL_THRESH = 0.02
BISULFITE_FAIL_COUNT = 2
DUP_CORR_MIN = 0.99
CONTAM_STD_MAX = 0.28   # contaminated sample should have lower std (collapsed bimodal)
CONTAM_SAMPLE = "S020"
REF_SAMPLE = "S005"


# ── Fixture ───────────────────────────────────────────────────────────────────

def load_data() -> pd.DataFrame:
    assert os.path.exists(CSV_PATH), (
        f"mock_methylation.csv not found at {CSV_PATH}. "
        "Run: python data/generate_mock_data.py"
    )
    return pd.read_csv(CSV_PATH)


# =============================================================================
# SCHEMA & COORDINATE TESTS
# =============================================================================

def test_chrom_and_pos_columns_exist():
    """chrom and pos columns must be present in the dataset."""
    df = load_data()
    assert "chrom" in df.columns, "Missing 'chrom' column"
    assert "pos" in df.columns, "Missing 'pos' column"


def test_total_cpg_count():
    """Dataset must contain 10,000 unique CpG sites."""
    df = load_data()
    n_cpgs = df["cpg_id"].nunique()
    assert n_cpgs == 10_000, f"Expected 10,000 CpGs, got {n_cpgs}"


def test_total_sample_count():
    """Dataset must contain 101 unique samples (100 patients + 1 duplicate)."""
    df = load_data()
    n_samples = df["sample_id"].nunique()
    assert n_samples == 101, f"Expected 101 samples, got {n_samples}"


def test_vdj_cpgs_in_real_coordinates():
    """Every VDJ CpG must fall within a GRCh38 VDJ locus interval."""
    from core.qc.repertoire_clonality import VDJ_LOCI_GRCH38

    df = load_data()
    vdj_cpgs = df[df["is_vdj_region"]][["cpg_id", "chrom", "pos"]].drop_duplicates("cpg_id")

    # Build interval lookup
    by_chrom = {}
    for _name, (chrom, start, end, _lineage) in VDJ_LOCI_GRCH38.items():
        by_chrom.setdefault(chrom, []).append((start, end))

    for _, row in vdj_cpgs.iterrows():
        intervals = by_chrom.get(row.chrom, [])
        in_any = any(s <= row.pos <= e for s, e in intervals)
        assert in_any, (
            f"VDJ CpG {row.cpg_id} at {row.chrom}:{row.pos} is not within any GRCh38 VDJ locus"
        )


def test_x_linked_cpgs_on_chrx():
    """All X-linked CpGs must have chrom == 'chrX'."""
    df = load_data()
    x_cpgs = df[df["is_x_chromosome"].astype(bool)]
    assert (x_cpgs["chrom"] == "chrX").all(), (
        "Some is_x_chromosome=True CpGs are not on chrX"
    )


# =============================================================================
# ARTIFACT 1 — Confounded Batch
# =============================================================================

def test_batch_case_confound():
    """
    Batch_01 Case samples must have a mean beta > Batch_01 Control mean beta
    by at least BATCH_SHIFT_MIN (0.05).
    """
    df = load_data()
    b1_case = df[(df["batch_id"] == "Batch_01") & (df["disease_label"] == "Case")]["beta_value"].mean()
    b1_ctrl = df[(df["batch_id"] == "Batch_01") & (df["disease_label"] == "Control")]["beta_value"].mean()
    delta = b1_case - b1_ctrl
    assert delta >= BATCH_SHIFT_MIN, (
        f"Artifact 1 FAIL: Batch_01 Case vs Control delta = {delta:.4f} "
        f"(expected >= {BATCH_SHIFT_MIN})"
    )


def test_batch1_case_enrichment():
    """
    At least 75% of Batch_01 samples should be Case (spec: 80%).
    """
    df = load_data()
    b1_samples = df[df["batch_id"] == "Batch_01"][["sample_id", "disease_label"]].drop_duplicates()
    case_frac = (b1_samples["disease_label"] == "Case").mean()
    assert case_frac >= 0.75, (
        f"Artifact 1 FAIL: Batch_01 Case fraction = {case_frac:.2f} (expected >= 0.75)"
    )


# =============================================================================
# ARTIFACT 2 — Clonal VDJ Artifact
# =============================================================================

def test_vdj_clonal_artifact_exists():
    """
    At least one patient must have VDJ-region CpGs with beta > 0.8 AND
    fragment_length > 195 bp simultaneously.  With the tightened baseline
    (Normal(150, 12), clipped at 220), background P(>195 bp) ≈ 0 %, so any
    hit here is the injected artifact, not noise.
    """
    df = load_data()
    vdj = df[df["is_vdj_region"]]
    clonal = vdj[(vdj["beta_value"] > CLONAL_BETA_MIN) & (vdj["fragment_length"] > CLONAL_FRAG_MIN)]
    assert len(clonal) > 0, (
        f"Artifact 2 FAIL: No VDJ rows with beta > {CLONAL_BETA_MIN} "
        f"AND fragment > {CLONAL_FRAG_MIN} bp"
    )


def test_vdj_clonal_is_single_case_patient():
    """
    The combined high-beta / long-fragment VDJ signal must come from exactly
    one patient (the injected clone) and that patient must be a Case.
    Multiple patients firing here means the baseline is leaking — not the artifact.
    """
    df = load_data()
    vdj = df[df["is_vdj_region"]]
    clonal = vdj[(vdj["beta_value"] > CLONAL_BETA_MIN) & (vdj["fragment_length"] > CLONAL_FRAG_MIN)]
    patients = clonal["patient_id"].unique()
    labels = clonal["disease_label"].unique()
    assert len(patients) == 1, (
        f"Artifact 2 FAIL: Expected exactly 1 clonal patient, "
        f"got {len(patients)}: {list(patients)}"
    )
    assert "Case" in labels, (
        f"Artifact 2 FAIL: Clonal patient is not a Case. Labels: {list(labels)}"
    )


# =============================================================================
# ARTIFACT 3 — Bisulfite Conversion Failure
# =============================================================================

def test_bisulfite_failure_count():
    """
    Exactly BISULFITE_FAIL_COUNT samples must have mean non_cpg_meth_rate > 2%.
    """
    df = load_data()
    sample_ncpg = df.groupby("sample_id")["non_cpg_meth_rate"].mean()
    failing = sample_ncpg[sample_ncpg > BISULFITE_FAIL_THRESH]
    assert len(failing) == BISULFITE_FAIL_COUNT, (
        f"Artifact 3 FAIL: Found {len(failing)} samples above 2% threshold, "
        f"expected {BISULFITE_FAIL_COUNT}. Samples: {list(failing.index)}"
    )


def test_bisulfite_failure_magnitude():
    """
    Failing samples should have non_cpg_meth_rate well above threshold (>3%),
    not just barely over — confirming genuine failure, not noise.
    """
    df = load_data()
    sample_ncpg = df.groupby("sample_id")["non_cpg_meth_rate"].mean()
    failing = sample_ncpg[sample_ncpg > BISULFITE_FAIL_THRESH]
    for sid, rate in failing.items():
        assert rate > 0.03, (
            f"Artifact 3 FAIL: Sample {sid} non_cpg_meth_rate = {rate:.4f} "
            "is barely above threshold — may be noise, not a true failure."
        )


# =============================================================================
# ARTIFACT 4 — Sample Duplication
# =============================================================================

def test_duplicate_sample_exists():
    """
    S_DUP must be present in the dataset.
    """
    df = load_data()
    assert "S_DUP" in df["sample_id"].values, (
        "Artifact 4 FAIL: Sample S_DUP not found in dataset."
    )


def test_duplicate_correlation():
    """
    S010 and S_DUP must have Pearson r > DUP_CORR_MIN (0.99) over all CpG beta values.
    """
    df = load_data()
    pivot = df.pivot_table(index="cpg_id", columns="sample_id", values="beta_value")
    assert "S010" in pivot.columns, "Artifact 4 FAIL: S010 not in pivot columns."
    assert "S_DUP" in pivot.columns, "Artifact 4 FAIL: S_DUP not in pivot columns."
    r, _ = pearsonr(pivot["S010"].fillna(0), pivot["S_DUP"].fillna(0))
    assert r >= DUP_CORR_MIN, (
        f"Artifact 4 FAIL: S010 vs S_DUP Pearson r = {r:.4f} "
        f"(expected >= {DUP_CORR_MIN})"
    )


# =============================================================================
# ARTIFACT 5 — Sample Contamination
# =============================================================================

def test_contamination_distribution_collapse():
    """
    The contaminated sample (S020) should have lower beta std than a clean
    reference (S005), indicating the bimodal distribution has been smeared.
    """
    df = load_data()
    std_contam = df[df["sample_id"] == CONTAM_SAMPLE]["beta_value"].std()
    std_ref = df[df["sample_id"] == REF_SAMPLE]["beta_value"].std()
    assert std_contam < std_ref, (
        f"Artifact 5 FAIL: Contaminated sample ({CONTAM_SAMPLE}) std = {std_contam:.4f} "
        f"is NOT less than reference ({REF_SAMPLE}) std = {std_ref:.4f}. "
        "Distribution collapse not detected."
    )


def test_contamination_mean_near_half():
    """
    The contaminated sample's mean beta should be closer to 0.5 than 0.4 or 0.6,
    reflecting the uniform contaminator pulling the distribution toward the center.
    """
    df = load_data()
    mean_contam = df[df["sample_id"] == CONTAM_SAMPLE]["beta_value"].mean()
    assert 0.40 <= mean_contam <= 0.65, (
        f"Artifact 5 FAIL: Contaminated sample mean beta = {mean_contam:.4f} "
        "not in expected muddy range [0.40, 0.65]."
    )


# =============================================================================
# ARTIFACT 6 — Low Coverage
# =============================================================================

def test_low_coverage_sample():
    """
    S030 must have mean depth close to 5x (Poisson λ=5), well below the
    10x quality threshold.
    """
    df = load_data()
    mean_depth = df[df["sample_id"] == "S030"]["depth"].mean()
    assert mean_depth < 10.0, (
        f"Artifact 6 FAIL: S030 mean depth = {mean_depth:.1f}x "
        "(expected < 10x for low-coverage failure)"
    )


# =============================================================================
# ARTIFACT 7 — Sex Metadata Mixup
# =============================================================================

def test_sex_mixup_artifact():
    """
    S035 must have female-like X-linked beta (~0.50, XCI signal) but sex='M';
    S036 must have male-like X-linked beta (~0.25) but sex='F'.
    This confirms both the injected XCI signal and the inverted metadata.
    """
    df = load_data()
    s35_x = df[(df["sample_id"] == "S035") & df["is_x_chromosome"].astype(bool)]["beta_value"].mean()
    s36_x = df[(df["sample_id"] == "S036") & df["is_x_chromosome"].astype(bool)]["beta_value"].mean()
    assert s35_x >= 0.40, (
        f"Artifact 7 FAIL: S035 X-beta = {s35_x:.3f} (expected female XCI ~0.50, >= 0.40)"
    )
    assert s36_x < 0.35, (
        f"Artifact 7 FAIL: S036 X-beta = {s36_x:.3f} (expected male ~0.25, < 0.35)"
    )
    assert df[df["sample_id"] == "S035"]["sex"].iloc[0] == "M", (
        "Artifact 7 FAIL: S035 metadata sex must be 'M' (true F, reported M)"
    )
    assert df[df["sample_id"] == "S036"]["sex"].iloc[0] == "F", (
        "Artifact 7 FAIL: S036 metadata sex must be 'F' (true M, reported F)"
    )


# =============================================================================
# STANDALONE __main__ — timestamped detection log
# =============================================================================

if __name__ == "__main__":
    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def log(status: str, artifact: str, detail: str):
        icon = "✓ DETECTED" if status == "ok" else "✗ MISSING "
        print(f"[{ts()}] [TEST] {icon} | {artifact} | {detail}")

    print(f"[{ts()}] [TEST] Loading {CSV_PATH}")
    df = load_data()
    print(f"[{ts()}] [TEST] Shape: {df.shape}  |  Samples: {df['sample_id'].nunique()}")
    print()

    # ── Artifact 1 ────────────────────────────────────────────────────────────
    b1_case = df[(df["batch_id"] == "Batch_01") & (df["disease_label"] == "Case")]["beta_value"].mean()
    b1_ctrl = df[(df["batch_id"] == "Batch_01") & (df["disease_label"] == "Control")]["beta_value"].mean()
    delta = b1_case - b1_ctrl
    b1_samples = df[df["batch_id"] == "Batch_01"][["sample_id", "disease_label"]].drop_duplicates()
    case_frac = (b1_samples["disease_label"] == "Case").mean()
    if delta >= BATCH_SHIFT_MIN and case_frac >= 0.75:
        log("ok", "Artifact 1 — Confounded Batch",
            f"Batch_01 Case mean β={b1_case:.3f}, Control mean β={b1_ctrl:.3f}, "
            f"Δ={delta:.3f}, Case fraction={case_frac:.0%}")
    else:
        log("fail", "Artifact 1 — Confounded Batch",
            f"Δ={delta:.3f} (need >={BATCH_SHIFT_MIN}), Case frac={case_frac:.0%} (need >=75%)")

    # ── Artifact 2 ────────────────────────────────────────────────────────────
    vdj = df[df["is_vdj_region"]]
    clonal = vdj[(vdj["beta_value"] > CLONAL_BETA_MIN) & (vdj["fragment_length"] > CLONAL_FRAG_MIN)]
    n_clone = len(clonal)
    patients = clonal["patient_id"].unique().tolist() if n_clone else []
    dis_labels = clonal["disease_label"].unique().tolist() if n_clone else []
    if n_clone > 0 and len(patients) == 1 and "Case" in dis_labels:
        log("ok", "Artifact 2 — Clonal VDJ",
            f"{n_clone} VDJ rows with β>{CLONAL_BETA_MIN} & frag>{CLONAL_FRAG_MIN}bp | "
            f"patient={patients[0]}, label={dis_labels[0]}")
    else:
        log("fail", "Artifact 2 — Clonal VDJ",
            f"rows={n_clone}, patients={patients}, labels={dis_labels} "
            f"(expected 1 Case patient)")

    # ── Artifact 3 ────────────────────────────────────────────────────────────
    sample_ncpg = df.groupby("sample_id")["non_cpg_meth_rate"].mean()
    failing = sample_ncpg[sample_ncpg > BISULFITE_FAIL_THRESH]
    if len(failing) == BISULFITE_FAIL_COUNT:
        log("ok", "Artifact 3 — Bisulfite Failure",
            f"{len(failing)} samples above 2% threshold: "
            + ", ".join(f"{s}={r:.4f}" for s, r in failing.items()))
    else:
        log("fail", "Artifact 3 — Bisulfite Failure",
            f"Found {len(failing)} failing samples (expected {BISULFITE_FAIL_COUNT}): "
            + str(list(failing.index)))

    # ── Artifact 4 ────────────────────────────────────────────────────────────
    if "S_DUP" in df["sample_id"].values:
        pivot = df.pivot_table(index="cpg_id", columns="sample_id", values="beta_value")
        r, _ = pearsonr(pivot["S010"].fillna(0), pivot["S_DUP"].fillna(0))
        if r >= DUP_CORR_MIN:
            log("ok", "Artifact 4 — Sample Duplication",
                f"S010 vs S_DUP Pearson r={r:.4f} (threshold >={DUP_CORR_MIN})")
        else:
            log("fail", "Artifact 4 — Sample Duplication",
                f"S010 vs S_DUP Pearson r={r:.4f} (need >={DUP_CORR_MIN})")
    else:
        log("fail", "Artifact 4 — Sample Duplication", "S_DUP not present in dataset")

    # ── Artifact 5 ────────────────────────────────────────────────────────────
    std_contam = df[df["sample_id"] == CONTAM_SAMPLE]["beta_value"].std()
    std_ref = df[df["sample_id"] == REF_SAMPLE]["beta_value"].std()
    mean_c = df[df["sample_id"] == CONTAM_SAMPLE]["beta_value"].mean()
    if std_contam < std_ref and 0.40 <= mean_c <= 0.65:
        log("ok", "Artifact 5 — Contamination",
            f"{CONTAM_SAMPLE} std={std_contam:.4f} < ref {REF_SAMPLE} std={std_ref:.4f}, "
            f"mean β={mean_c:.3f} (muddy range confirmed)")
    else:
        log("fail", "Artifact 5 — Contamination",
            f"{CONTAM_SAMPLE} std={std_contam:.4f}, ref std={std_ref:.4f}, mean β={mean_c:.3f}")

    print()
    print(f"[{ts()}] [TEST] Done. Run full suite: pytest tests/test_mock_data.py -v")
