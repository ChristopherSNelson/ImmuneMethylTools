"""
generate_mock_data.py — ImmuneMethylTools Lab Simulator
========================================================
Generates mock_methylation.csv with eight embedded "stumper" artifacts that
mimic real-world pitfalls in immune-cell WGBS / RRBS analysis:

  Artifact 1 — Confounded Batch:   Batch_01 enriched for Cases (+0.1 beta shift)
  Artifact 2 — Clonal Artifact:    VDJ locus, beta > 0.8, fragment > 180 bp
  Artifact 3 — Bisulfite Failure:  2 samples with non_cpg_meth_rate > 0.02
  Artifact 4 — Sample Duplication: 2 samples with Pearson r > 0.99
  Artifact 5 — Contamination:      1 sample with muddy beta (peak near 0.5)
  Artifact 6 — Low Coverage:       S030 depth forced to Poisson(lambda=5), mean ~5x
  Artifact 7 — Sex Metadata Mixup: S035 (true F) reported M; S036 (true M) reported F
  Artifact 8 — Lineage Anomaly:    S045/S046 FoxP3 proxy beta ~0.06 (Treg-enriched);
                                    S065/S066 PAX5 proxy beta ~0.65 (B-cell depleted)

In addition to the stumper artifacts, one genuine biological signal is injected:
  True DMR      — cg00003000-cg00003010: all Case samples +0.10 beta (autosomal; batch-independent)
  This gives dmr_hunter a real, non-artifactual DMR to detect after batch correction.

Outputs
-------
  data/mock_methylation.csv
  output/figures/qc_before_after.png   — Before/After visualization of each artifact
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless — no display needed
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.gridspec as gridspec  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy.stats import pearsonr  # noqa: E402

# ── Reproducibility ──────────────────────────────────────────────────────────
RNG = np.random.default_rng(seed=42)

# ── Study parameters ─────────────────────────────────────────────────────────
N_PATIENTS = 100         # 50 Case, 50 Control
N_CPGS = 10_000          # CpG sites per sample
N_X_CPGS = 600           # last N_X_CPGS of N_CPGS are X-linked
N_BATCHES = 3
CASE_LABEL = "Case"
CTRL_LABEL = "Control"

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "output", "figures")
OUT_CSV = os.path.join(os.path.dirname(__file__), "mock_methylation.csv")
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── GRCh38 VDJ locus coordinates (with 2 kb buffer) ─────────────────────────
VDJ_LOCI_GRCH38: dict[str, tuple[str, int, int, str]] = {
    "IGH":     ("chr14", 105_584_437, 106_881_844, "B-cell"),
    "IGK":     ("chr2",   88_855_361,  90_237_368, "B-cell"),
    "IGL":     ("chr22",  22_024_076,  22_924_913, "B-cell"),
    "TRA_TRD": ("chr14",  22_088_057,  23_023_075, "T-cell"),
    "TRB":     ("chr7",  142_297_011, 142_815_287, "T-cell"),
    "TRG":     ("chr7",   38_238_024,  38_370_055, "T-cell"),
}

# ── GRCh38 chromosome sizes (autosomes only) ────────────────────────────────
CHROM_SIZES_GRCH38: dict[str, int] = {
    "chr1":  248_956_422, "chr2":  242_193_529, "chr3":  198_295_559,
    "chr4":  190_214_555, "chr5":  181_538_259, "chr6":  170_805_979,
    "chr7":  159_345_973, "chr8":  145_138_636, "chr9":  138_394_717,
    "chr10": 133_797_422, "chr11": 135_086_622, "chr12": 133_275_309,
    "chr13": 114_364_328, "chr14": 107_043_718, "chr15": 101_991_189,
    "chr16":  90_338_345, "chr17":  83_257_441, "chr18":  80_373_285,
    "chr19":  58_617_616, "chr20":  64_444_167, "chr21":  46_709_983,
    "chr22":  50_818_468,
}
CHRX_SIZE = 156_040_895


# =============================================================================
# 0.  COORDINATE INFRASTRUCTURE
# =============================================================================

def _build_vdj_intervals() -> dict[str, list[tuple[int, int, str]]]:
    """Group VDJ loci by chromosome for interval checks."""
    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for name, (chrom, start, end, _lineage) in VDJ_LOCI_GRCH38.items():
        by_chrom.setdefault(chrom, []).append((start, end, name))
    return by_chrom


def _in_vdj(chrom: str, pos: int, vdj_by_chrom: dict) -> bool:
    """Check whether a (chrom, pos) falls within any VDJ locus."""
    for start, end, _name in vdj_by_chrom.get(chrom, []):
        if start <= pos <= end:
            return True
    return False


def assign_cpg_coordinates(n_cpgs: int, n_x_cpgs: int) -> dict[int, tuple[str, int]]:
    """
    Assign realistic GRCh38 (chrom, pos) to each CpG index 1..n_cpgs.

    Layout:
      - Last n_x_cpgs indices → chrX, evenly spaced
      - ~3% of autosomal indices → placed inside real VDJ loci (proportional to locus size)
      - Remaining autosomal → distributed across chr1-22 proportional to chromosome length
      - Reserved signal CpG ranges are placed on non-VDJ, non-X autosomes
    """
    coord_map: dict[int, tuple[str, int]] = {}
    vdj_by_chrom = _build_vdj_intervals()
    n_autosomal = n_cpgs - n_x_cpgs  # 9400

    # ── chrX CpGs (last n_x_cpgs indices) ──────────────────────────────────
    x_start_idx = n_cpgs - n_x_cpgs + 1  # 9401
    spacing = CHRX_SIZE // (n_x_cpgs + 1)
    for i, cpg_idx in enumerate(range(x_start_idx, n_cpgs + 1)):
        pos = (i + 1) * spacing
        coord_map[cpg_idx] = ("chrX", pos)

    # ── VDJ CpGs (~3% of autosomal) ───────────────────────────────────────
    n_vdj = int(0.03 * n_autosomal)  # ~282
    # Distribute proportionally by locus size
    loci = list(VDJ_LOCI_GRCH38.items())
    locus_sizes = np.array([end - start for _n, (_c, start, end, _l) in loci], dtype=float)
    locus_fracs = locus_sizes / locus_sizes.sum()
    locus_counts = np.round(locus_fracs * n_vdj).astype(int)
    # Fix rounding to hit exactly n_vdj
    diff = n_vdj - locus_counts.sum()
    locus_counts[np.argmax(locus_sizes)] += diff

    vdj_indices: list[int] = []
    vdj_coords: list[tuple[str, int]] = []
    # Collect available autosomal indices (1..n_autosomal), excluding reserved signal ranges
    reserved_ranges = set()
    for rng in [range(3000, 3011), range(1500, 1508), range(2000, 2006),
                range(1, 6), range(6, 11)]:
        reserved_ranges.update(rng)

    available_for_vdj = [i for i in range(1, n_autosomal + 1) if i not in reserved_ranges]
    chosen_vdj_indices = sorted(RNG.choice(available_for_vdj, size=n_vdj, replace=False))
    vdj_indices = list(chosen_vdj_indices)

    idx_cursor = 0
    for li, (name, (chrom, start, end, _lineage)) in enumerate(loci):
        count = int(locus_counts[li])
        for j in range(count):
            cpg_idx = vdj_indices[idx_cursor]
            # Place at evenly spaced positions within the locus
            pos = start + int((j + 1) * (end - start) / (count + 1))
            coord_map[cpg_idx] = (chrom, pos)
            idx_cursor += 1

    vdj_set = set(vdj_indices)

    # ── Place signal CpG clusters at fixed genomic positions ──────────────
    # Signal ranges must form tight clusters (< 1000 bp spacing) so that
    # the distance-based DMR caller can group them into testable regions.
    _SIGNAL_CLUSTERS: list[tuple[range, str, int, int]] = [
        # (cpg_index_range, chrom, start_pos, spacing_bp)
        (range(3000, 3011), "chr6",  30_000_000, 200),   # true bio DMR (11 CpGs)
        (range(1500, 1508), "chr11", 50_000_000, 400),   # borderline negative control
        (range(2000, 2006), "chr15", 40_000_000, 400),   # subtle negative control
    ]
    signal_placed = set()
    for idx_range, chrom, start_pos, spacing in _SIGNAL_CLUSTERS:
        for j, cpg_idx in enumerate(idx_range):
            coord_map[cpg_idx] = (chrom, start_pos + j * spacing)
            signal_placed.add(cpg_idx)

    # ── Remaining autosomal CpGs ──────────────────────────────────────────
    remaining_indices = [
        i for i in range(1, n_autosomal + 1)
        if i not in vdj_set and i not in signal_placed
    ]
    n_remaining = len(remaining_indices)

    # Distribute across chr1-22 proportional to chromosome length
    chroms = list(CHROM_SIZES_GRCH38.keys())
    sizes = np.array([CHROM_SIZES_GRCH38[c] for c in chroms], dtype=float)
    chrom_fracs = sizes / sizes.sum()
    chrom_counts = np.round(chrom_fracs * n_remaining).astype(int)
    # Fix rounding
    diff = n_remaining - chrom_counts.sum()
    chrom_counts[0] += diff

    idx_cursor = 0
    for ci, chrom in enumerate(chroms):
        count = int(chrom_counts[ci])
        chrom_size = CHROM_SIZES_GRCH38[chrom]
        intervals = vdj_by_chrom.get(chrom, [])

        for j in range(count):
            cpg_idx = remaining_indices[idx_cursor]
            # Generate a random position, avoiding VDJ intervals
            for _attempt in range(20):
                pos = int(RNG.integers(1_000_000, chrom_size - 1_000_000))
                if not any(s <= pos <= e for s, e, _n in intervals):
                    break
            coord_map[cpg_idx] = (chrom, pos)
            idx_cursor += 1

    return coord_map


# =============================================================================
# 1.  PATIENT / SAMPLE MANIFEST
# =============================================================================

def build_manifest() -> pd.DataFrame:
    """
    Create one row per (sample, CpG).  Assigns batch with confounded Case/Batch_01
    enrichment (Artifact 1 setup).
    """
    n_half = N_PATIENTS // 2
    patient_ids = [f"P{i:03d}" for i in range(1, N_PATIENTS + 1)]
    disease_lbls = ([CASE_LABEL] * n_half) + ([CTRL_LABEL] * n_half)
    ages = RNG.integers(25, 70, size=N_PATIENTS).tolist()

    # Artifact 1 — Batch_01 gets 80 % of Case patients
    case_patients = [p for p, d in zip(patient_ids, disease_lbls) if d == CASE_LABEL]
    ctrl_patients = [p for p, d in zip(patient_ids, disease_lbls) if d == CTRL_LABEL]

    n_case_b1 = int(0.80 * len(case_patients))      # 40 / 50
    n_ctrl_b1 = int(0.20 * len(ctrl_patients))       # 10 / 50

    # Remaining Case patients split evenly across Batch_02 and Batch_03
    n_case_remaining = len(case_patients) - n_case_b1  # 10
    n_case_b2 = n_case_remaining // 2                  # 5
    # Remaining Control patients split evenly across Batch_02 and Batch_03
    n_ctrl_remaining = len(ctrl_patients) - n_ctrl_b1  # 40
    ctrl_midpoint = n_ctrl_b1 + n_ctrl_remaining // 2  # 10 + 20 = 30

    batch_map = {}
    for p in case_patients[:n_case_b1]:
        batch_map[p] = "Batch_01"
    for p in case_patients[n_case_b1:n_case_b1 + n_case_b2]:
        batch_map[p] = "Batch_02"
    for p in case_patients[n_case_b1 + n_case_b2:]:
        batch_map[p] = "Batch_03"
    for p in ctrl_patients[:n_ctrl_b1]:
        batch_map[p] = "Batch_01"
    for p in ctrl_patients[n_ctrl_b1:ctrl_midpoint]:
        batch_map[p] = "Batch_02"
    for p in ctrl_patients[ctrl_midpoint:]:
        batch_map[p] = "Batch_03"

    age_map = dict(zip(patient_ids, ages))
    disease_map = dict(zip(patient_ids, disease_lbls))

    # Odd patient number → Female; even → Male
    sex_map = {
        f"P{i:03d}": ("F" if i % 2 == 1 else "M")
        for i in range(1, N_PATIENTS + 1)
    }

    # Pre-compute coordinate map once
    coord_map = assign_cpg_coordinates(N_CPGS, N_X_CPGS)

    rows = []
    sample_counter = 1
    for pid in patient_ids:
        sid = f"S{sample_counter:03d}"
        sample_counter += 1
        for cg_idx in range(1, N_CPGS + 1):
            chrom, pos = coord_map[cg_idx]
            rows.append({
                "sample_id": sid,
                "patient_id": pid,
                "batch_id": batch_map[pid],
                "age": age_map[pid],
                "disease_label": disease_map[pid],
                "cpg_id": f"cg{cg_idx:08d}",
                "sex": sex_map[pid],
                "is_x_chromosome": chrom == "chrX",
                "chrom": chrom,
                "pos": pos,
            })

    return pd.DataFrame(rows)


# =============================================================================
# 2.  BASELINE METHYLATION VALUES
# =============================================================================

def add_baseline_methylation(df: pd.DataFrame) -> pd.DataFrame:
    """
    Draw beta values from a bimodal distribution typical of CpG methylation:
      ~60 % fully methylated (beta ~ 0.85),  ~40 % fully unmethylated (beta ~ 0.10).
    Small Gaussian noise added per CpG site.

    Bisulfite intuition notes:
    - Fragment length baseline: Normal(150, 12).  P(>180 bp) ~ 0.6 %,  making
      the VDJ clonal artifact (>200 bp) a genuine outlier, not background noise.
    - VDJ-region baseline beta: set to low/unmethylated (~0.10) after the
      genome-wide draw.  Biologically, active VDJ loci in B/T cells must remain
      accessible (unmethylated) for recombination; hypermethylation there only
      arises AFTER a dominant clone silences the locus — i.e., the artifact.
    """
    n = len(df)
    # Site-level random component (same shape across samples for a given CpG)
    cpg_means = RNG.choice([0.10, 0.85], size=N_CPGS, p=[0.40, 0.60])
    cpg_mean_series = np.tile(cpg_means, N_PATIENTS)  # repeated for each sample

    noise = RNG.normal(0, 0.06, size=n)
    beta = np.clip(cpg_mean_series + noise, 0.0, 1.0)
    df["beta_value"] = beta

    # Sequencing depth: Negative Binomial ~ mean 30, overdispersion
    df["depth"] = RNG.negative_binomial(n=5, p=5 / (5 + 25), size=n).clip(1, None)

    # Fragment length: tight Normal(150, 12) — >180 bp is ~0.6 % of baseline.
    # This ensures the clonal VDJ signal (injected at 200-260 bp) is a clear outlier.
    df["fragment_length"] = RNG.normal(150, 12, size=n).astype(int).clip(80, 220)

    # Derive is_vdj_region from coordinate map (chrom/pos already in DataFrame)
    vdj_by_chrom = _build_vdj_intervals()
    # Compute VDJ status once per CpG (not per sample)
    cpg_coords = df[["cpg_id", "chrom", "pos"]].drop_duplicates("cpg_id")
    cpg_vdj = {
        row.cpg_id: _in_vdj(row.chrom, row.pos, vdj_by_chrom)
        for _, row in cpg_coords.iterrows()
    }
    df["is_vdj_region"] = df["cpg_id"].map(cpg_vdj)

    # Bisulfite intuition: VDJ baseline must be unmethylated.
    # Overwrite the bimodal draw for VDJ CpGs with a low-methylation distribution.
    vdj_mask = df["is_vdj_region"]
    df.loc[vdj_mask, "beta_value"] = RNG.normal(0.10, 0.05, size=vdj_mask.sum()).clip(0.0, 0.35)

    # Non-CpG methylation rate (bisulfite conversion proxy): ~ N(0.004, 0.001)
    df["non_cpg_meth_rate"] = RNG.normal(0.004, 0.001, size=n).clip(0, 1)

    # GC content: fraction of G+C dinucleotides around each CpG site.
    # Constant per CpG (identical across all samples for a given cpg_id).
    # Random Uniform(0.30, 0.70) in mock data — serves as a negative control
    # (no association with disease) and realistic structure for real-data runs.
    gc_map = {idx: round(RNG.uniform(0.30, 0.70), 2) for idx in range(1, N_CPGS + 1)}
    df["gc_content"] = df["cpg_id"].apply(
        lambda c: gc_map[int(c.lstrip("cg"))]
    )

    return df


# =============================================================================
# ARTIFACTS
# =============================================================================

def inject_artifact1_confounded_batch(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 1 — Confounded Batch
    Batch_01 Case samples receive a systematic +0.10 beta shift.
    This mimics a plate / reagent batch effect that is collinear with disease.
    """
    mask = (df["batch_id"] == "Batch_01") & (df["disease_label"] == CASE_LABEL)
    df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + 0.10).clip(0, 1)
    print(f"  [Artifact 1] +0.1 shift applied to {mask.sum()} rows "
          f"({df.loc[mask, 'sample_id'].nunique()} samples in Batch_01/Case)")
    return df


def inject_artifact2_clonal_vdj(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 2 — Clonal Expansion Artifact in VDJ Locus
    Pick one Case patient; for all their VDJ-region CpGs on chr14 (IGH locus)
    set beta > 0.8 and fragment_length > 180 bp.  This mimics a dominant clone
    in which VDJ loci are hypermethylated and the original fragment is long
    (compact chromatin).
    """
    case_patients = df.loc[df["disease_label"] == CASE_LABEL, "patient_id"].unique()
    # Skip index 0 (P001 -> S001) and index 1 (P002 -> S002): both carry the
    # bisulfite-failure artifact (inject_artifact3_bisulfite_failure targets
    # all_samples[:2]).  Selecting index 2 ensures the clonal patient passes
    # sample-level QC so the VDJ masking stage can be demonstrated.
    clonal_patient = case_patients[2]  # deterministic; avoids bisulfite-failure overlap

    mask = (
        (df["patient_id"] == clonal_patient)
        & (df["is_vdj_region"])
        & (df["chrom"] == "chr14")
    )
    n_affected = mask.sum()

    df.loc[mask, "beta_value"] = RNG.uniform(0.82, 0.97, size=n_affected)
    df.loc[mask, "fragment_length"] = RNG.integers(200, 260, size=n_affected)

    print(f"  [Artifact 2] Clonal VDJ artifact injected into patient {clonal_patient}: "
          f"{n_affected} CpG rows (beta > 0.8, fragment > 180 bp, chr14/IGH)")
    return df


def inject_artifact3_bisulfite_failure(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 3 — Incomplete Bisulfite Conversion
    Two samples receive non_cpg_meth_rate drawn from N(0.05, 0.01) — well above
    the 2 % failure threshold.  High non-CpG methylation = cytosines not converted
    = inflated beta values genome-wide.
    """
    all_samples = df["sample_id"].unique()
    # Choose 2 samples that are NOT already the clonal patient's sample
    bad_samples = all_samples[:2]  # deterministic; S001, S002 (both Case/Batch_01)

    for sid in bad_samples:
        mask = df["sample_id"] == sid
        df.loc[mask, "non_cpg_meth_rate"] = RNG.normal(0.05, 0.008, size=mask.sum()).clip(0.02, 1)
        # Bisulfite failure also artifactually inflates beta
        df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + RNG.normal(0.08, 0.02, size=mask.sum())).clip(0, 1)

    print(f"  [Artifact 3] Bisulfite failure injected into samples: {list(bad_samples)}")
    return df


def inject_artifact4_sample_duplication(df: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 4 — Technical Duplicate
    Clone sample S010's beta values into a new sample S_DUP with tiny noise so
    Pearson r > 0.99.  S_DUP gets a unique patient_id (P_DUP) so it does not
    inflate the CpG count attributed to S010's patient.
    """
    source_sid = "S010"
    dup_sid = "S_DUP"
    dup_pid = "P_DUP"

    source_rows = df[df["sample_id"] == source_sid].copy()
    source_rows["sample_id"] = dup_sid
    source_rows["patient_id"] = dup_pid
    # Add tiny noise (std = 0.002) to keep r > 0.99
    source_rows["beta_value"] = (
        source_rows["beta_value"] + RNG.normal(0, 0.002, size=len(source_rows))
    ).clip(0, 1)

    df = pd.concat([df, source_rows], ignore_index=True)

    # Verify
    orig = df.loc[df["sample_id"] == source_sid, "beta_value"].values
    clone = df.loc[df["sample_id"] == dup_sid, "beta_value"].values
    r, _ = pearsonr(orig, clone)
    print(f"  [Artifact 4] Duplicate {dup_sid} (patient {dup_pid}) cloned from {source_sid}: r = {r:.4f}")
    return df


def inject_artifact5_contamination(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 5 — Sample Contamination
    Pick one sample and collapse its bimodal beta distribution toward 0.5.
    Real-world: cross-contamination from another sample type smears the signal.
    We use a mixture: 50 % original + 50 % uniform(0.3, 0.7) noise.
    """
    contaminated_sid = "S020"
    mask = df["sample_id"] == contaminated_sid
    n = mask.sum()

    original_beta = df.loc[mask, "beta_value"].values
    contaminant = RNG.uniform(0.35, 0.65, size=n)      # peaks near 0.5
    mixed_beta = 0.50 * original_beta + 0.50 * contaminant
    df.loc[mask, "beta_value"] = mixed_beta.clip(0, 1)

    print(f"  [Artifact 5] Contamination injected into {contaminated_sid}: "
          f"mean beta shifted to {mixed_beta.mean():.3f}")
    return df


def inject_artifact6_low_depth(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 6 — Low Coverage (depth failure)
    Force sample S030's read depth to Poisson(lambda=5), yielding a mean of ~5x.
    Sites with mean depth < 10 reads have inflated binomial sampling variance;
    beta values from such sites are statistically unreliable and must be
    excluded by the depth filter in qc_guard.audit_quality().
    """
    mask = df["sample_id"] == "S030"
    df.loc[mask, "depth"] = RNG.poisson(lam=5, size=mask.sum())
    mean_depth = df.loc[mask, "depth"].mean()
    print(f"  [Artifact 6] Low-coverage failure injected into S030: "
          f"mean depth = {mean_depth:.1f}x (threshold: 10x)")
    return df


def inject_true_biological_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Inject a legitimate DMR that is NOT a batch effect or artifact.
    We pick a 11-CpG window and set a clean baseline (~0.35) then add
    +0.25 for Case samples only.  This avoids beta-ceiling clipping
    (which compresses the delta when baseline values are already high).

    Chosen region (cg00003000-cg00003010) is far from the proxy markers
    (FoxP3/PAX5/VDJ) and is autosomal (not X-linked), so it is unaffected
    by the XCI signal re-injection.

    Expected values after injection:
      Control: ~0.35 (low baseline, no shift)
      Case:    ~0.60 (baseline + 0.25)
      Raw delta: ~0.25

    After median-centering and per-sample median aggregation in the
    distance-based cluster DMR caller, the observed delta-beta is ~0.15,
    comfortably above DELTA_BETA_MIN = 0.10.

    Purpose: gives dmr_hunter at least one genuinely significant, non-artifactual
    DMR to find after batch correction, demonstrating the full detection pipeline.
    """
    target_cpgs = [f"cg{i:08d}" for i in range(3000, 3011)]
    all_signal_mask = df["cpg_id"].isin(target_cpgs)
    case_signal_mask = all_signal_mask & (df["disease_label"] == "Case")

    # Set a clean low baseline for ALL samples (Case + Control)
    n_all = int(all_signal_mask.sum())
    df.loc[all_signal_mask, "beta_value"] = RNG.normal(0.35, 0.04, size=n_all).clip(0.15, 0.55)

    # Add Case-only shift
    n_case = int(case_signal_mask.sum())
    df.loc[case_signal_mask, "beta_value"] = (
        df.loc[case_signal_mask, "beta_value"] + 0.25
    ).clip(0, 1)

    print(f"  [Signal Spike] Injected true biological signal into {len(target_cpgs)} CpGs "
          f"({n_case} Case rows; baseline=0.35 + Case shift=+0.25).")
    return df


def inject_borderline_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Inject a borderline sub-threshold biological signal at cg00001500-cg00001507
    (8 CpGs). Raw shift +0.09 to Case beta values; after median-centring the
    observed delta-beta should be ~0.08, just below the DELTA_BETA_MIN = 0.10 cutoff.

    Purpose: negative control confirming the pipeline does not over-call near
    the detection boundary.
    """
    target_cpgs = [f"cg{i:08d}" for i in range(1500, 1508)]
    mask = (df["disease_label"] == "Case") & (df["cpg_id"].isin(target_cpgs))
    df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + 0.09).clip(0, 1)
    print(f"  [Borderline Signal] Injected borderline signal into {len(target_cpgs)} CpGs "
          f"({mask.sum()} Case rows; +0.09 shift, expected delta-beta ~0.08).")
    return df


def inject_subtle_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Inject a clearly sub-threshold biological signal at cg00002000-cg00002005
    (6 CpGs). Raw shift +0.08 to Case beta values; after median-centring the
    observed delta-beta should be ~0.04, well below the DELTA_BETA_MIN = 0.10 cutoff.

    Purpose: negative control confirming the pipeline ignores weak signals.
    """
    target_cpgs = [f"cg{i:08d}" for i in range(2000, 2006)]
    mask = (df["disease_label"] == "Case") & (df["cpg_id"].isin(target_cpgs))
    df.loc[mask, "beta_value"] = (df.loc[mask, "beta_value"] + 0.08).clip(0, 1)
    print(f"  [Subtle Signal] Injected subtle signal into {len(target_cpgs)} CpGs "
          f"({mask.sum()} Case rows; +0.08 shift, expected delta-beta ~0.04).")
    return df


def inject_xci_signal(df: pd.DataFrame) -> pd.DataFrame:
    """
    Re-inject XCI-appropriate methylation signal for all X-linked CpGs.

    Called AFTER artifacts 1-6 so that batch shifts, bisulfite failures, and
    contamination do not corrupt the ground-truth X-linked beta used by the
    XCI guard.  Artifact 7 then swaps sex metadata for S035/S036 without
    touching beta values.

    Female (XX): XCI -> mean ~0.50 at X-linked sites (one X active, one silenced)
    Male   (XY): single active X -> unmethylated baseline, mean ~0.25
    """
    female_x_mask = df["is_x_chromosome"].astype(bool) & (df["sex"] == "F")
    df.loc[female_x_mask, "beta_value"] = (
        RNG.normal(0.50, 0.05, size=int(female_x_mask.sum())).clip(0.35, 0.65)
    )
    male_x_mask = df["is_x_chromosome"].astype(bool) & (df["sex"] == "M")
    df.loc[male_x_mask, "beta_value"] = (
        RNG.normal(0.25, 0.04, size=int(male_x_mask.sum())).clip(0.10, 0.33)
    )
    n_f = int(female_x_mask.sum())
    n_m = int(male_x_mask.sum())
    print(f"  [XCI Signal] X-linked betas re-injected: "
          f"{n_f} female rows (~0.50) | {n_m} male rows (~0.25)")
    return df


def inject_artifact8_lineage(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 8 — Lineage Composition Anomaly

    Injects two types of cell-fraction anomaly detectable by the deconvolution
    module (core/deconvolution.py):

    Treg-enriched samples — S045 (Case, P045) and S046 (Case, P046):
        FoxP3 proxy CpGs (cg00000001–cg00000005, chr1) are driven to low beta
        ~0.06, mimicking an active/open FoxP3 promoter — the epigenetic
        hallmark of Treg enrichment.  These betas fall clearly below both the
        male threshold (0.15) and the female threshold (0.30) used in
        detect_lineage_shift(), so both samples are flagged regardless of sex.

    B-cell-depleted samples — S065 (Control, P065) and S066 (Control, P066):
        PAX5 proxy CpGs (cg00000006–cg00000010, chr1) are driven to high beta
        ~0.65, mimicking PAX5 silencing — loss of B-cell epigenetic identity.
        These betas exceed the 0.50 PAX5_BCELL_SHIFT_THRESH and trigger a
        B-cell-shift flag in detect_lineage_shift().

    Both proxy panels are autosomal (chr1) in the mock CpG layout, so
    is_x_chromosome is False and no XCI correction is applied in
    estimate_cell_fractions().  The defensive foxp3_x / pax5_x guards in
    that function therefore evaluate to False — by design.

    Call order: must run AFTER inject_xci_signal() (XCI only touches
    is_x_chromosome rows and does not overwrite these chr1 proxies, so order
    is immaterial in practice, but placing it last keeps the intent explicit).
    """
    foxp3_cpgs = [f"cg{i:08d}" for i in range(1, 6)]
    pax5_cpgs  = [f"cg{i:08d}" for i in range(6, 11)]
    treg_samples  = ["S045", "S046"]
    bcell_samples = ["S065", "S066"]

    # ── Treg-enriched: FoxP3 hypomethylation ─────────────────────────────────
    treg_mask = df["sample_id"].isin(treg_samples) & df["cpg_id"].isin(foxp3_cpgs)
    n_treg = int(treg_mask.sum())
    df.loc[treg_mask, "beta_value"] = RNG.normal(0.06, 0.02, size=n_treg).clip(0.01, 0.12)

    # ── B-cell-depleted: PAX5 hypermethylation ────────────────────────────────
    bcell_mask = df["sample_id"].isin(bcell_samples) & df["cpg_id"].isin(pax5_cpgs)
    n_bcell = int(bcell_mask.sum())
    df.loc[bcell_mask, "beta_value"] = RNG.normal(0.65, 0.05, size=n_bcell).clip(0.52, 0.80)

    print(f"  [Artifact 8] Lineage anomalies injected:")
    print(f"    Treg-enriched   : {treg_samples} — FoxP3 proxy beta ~0.06 "
          f"({n_treg} rows; < M=0.15 / F=0.30 thresholds)")
    print(f"    B-cell depleted : {bcell_samples} — PAX5 proxy beta ~0.65 "
          f"({n_bcell} rows; > 0.50 threshold)")
    return df


def inject_artifact7_sex_mixup(df: pd.DataFrame) -> pd.DataFrame:
    """
    Artifact 7 — Sex Metadata Mixup

    S035 (truly female, P035): X-linked beta was injected as female (~0.50,
      XCI signal), but the 'sex' metadata column is overwritten to "M".
      Reported sex contradicts the observed X-linked methylation signal.

    S036 (truly male, P036): X-linked beta was injected as male (~0.25),
      but the 'sex' metadata column is overwritten to "F".
      Reported sex contradicts the observed X-linked methylation signal.

    detect_sex_mixups() (core/xci_guard.py) will flag both samples for exclusion.
    """
    df.loc[df["sample_id"] == "S035", "sex"] = "M"   # true F, reported M
    df.loc[df["sample_id"] == "S036", "sex"] = "F"   # true M, reported F
    n_s035 = (df["sample_id"] == "S035").sum()
    n_s036 = (df["sample_id"] == "S036").sum()
    print(f"  [Artifact 7] Sex metadata mixup injected: "
          f"S035 ({n_s035} rows, true F -> reported M), "
          f"S036 ({n_s036} rows, true M -> reported F)")
    return df


# =============================================================================
# 3.  BEFORE / AFTER VISUALIZATION
# =============================================================================

def plot_before_after(df_before: pd.DataFrame, df_after: pd.DataFrame) -> None:
    """
    Six-panel figure showing the fingerprint of each artifact before and after
    injection so developers can visually confirm the simulation worked.
    """
    fig = plt.figure(figsize=(20, 16))
    fig.suptitle(
        "ImmuneMethylTools — Mock Data: Before vs. After Artifact Injection",
        fontsize=14, fontweight="bold", y=0.98
    )
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.55, wspace=0.4)

    palette = {"Case": "#E74C3C", "Control": "#3498DB"}

    # ── Panel A: Batch x Disease (Artifact 1) ──────────────────────────────
    ax_a1 = fig.add_subplot(gs[0, 0])
    ax_a2 = fig.add_subplot(gs[0, 1])
    panels_a = [(ax_a1, df_before, "Before Artifact Injection"), (ax_a2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_a:
        sample_mean = (df_.groupby(["sample_id", "batch_id", "disease_label"])
                       ["beta_value"].mean().reset_index())
        sns.boxplot(data=sample_mean, x="batch_id", y="beta_value",
                    hue="disease_label", palette=palette, ax=ax,
                    linewidth=0.8, fliersize=2)
        # Annotate n per (batch, disease) group inside the top of each box column
        counts = sample_mean.groupby(["batch_id", "disease_label"]).size()
        batches = sorted(sample_mean["batch_id"].unique())
        hue_ord = sorted(sample_mean["disease_label"].unique())
        n_hue = len(hue_ord)
        bw = 0.8 / n_hue
        offsets = [(i - (n_hue - 1) / 2) * bw for i in range(n_hue)]
        ylo, yhi = ax.get_ylim()
        y_top = yhi - (yhi - ylo) * 0.02   # 2% below the top edge
        for bi, batch in enumerate(batches):
            for di, disease in enumerate(hue_ord):
                n = counts.get((batch, disease), 0)
                ax.text(bi + offsets[di], y_top, f"n={n}",
                        ha="center", va="top", fontsize=6)
        ax.set_ylim(ylo, yhi)   # lock ylim so text doesn't trigger autoscale
        ax.set_title(f"{title}: Batch x Disease\nmean beta", fontsize=9)
        ax.set_xlabel("Batch", fontsize=8)
        ax.set_ylabel("Mean beta", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, title_fontsize=6)

    # ── Panel B: VDJ Fragment Length (Artifact 2) ──────────────────────────
    ax_b1 = fig.add_subplot(gs[0, 2])
    ax_b2 = fig.add_subplot(gs[0, 3])
    panels_b = [(ax_b1, df_before, "Before Artifact Injection"), (ax_b2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_b:
        vdj_data = df_[df_["is_vdj_region"]]
        sns.scatterplot(data=vdj_data, x="fragment_length", y="beta_value",
                        hue="disease_label", palette=palette, ax=ax,
                        s=8, alpha=0.5, linewidth=0)
        ax.axvline(180, color="k", linestyle="--", linewidth=0.8, label="180 bp")
        ax.axhline(0.8, color="gray", linestyle="--", linewidth=0.8, label="beta=0.8")
        ax.set_title(f"{title}: VDJ Loci\nFragment vs Beta", fontsize=9)
        ax.set_xlabel("Fragment length (bp)", fontsize=8)
        ax.set_ylabel("Beta value", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=5)

    # ── Panel C: Non-CpG Meth Rate (Artifact 3) ────────────────────────────
    ax_c1 = fig.add_subplot(gs[1, 0])
    ax_c2 = fig.add_subplot(gs[1, 1])
    panels_c = [(ax_c1, df_before, "Before Artifact Injection"), (ax_c2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_c:
        sample_ncpg = df_.groupby("sample_id")["non_cpg_meth_rate"].mean().reset_index()
        ax.hist(sample_ncpg["non_cpg_meth_rate"], bins=30, range=(0, 0.07), color="#2ECC71", edgecolor="white")
        ax.axvline(0.02, color="red", linestyle="--", linewidth=1.2, label="2% threshold")
        ax.set_title(f"{title}: Non-CpG Meth Rate\n(bisulfite QC)", fontsize=9)
        ax.set_xlabel("Non-CpG meth rate", fontsize=8)
        ax.set_ylabel("# Samples", fontsize=8)
        ax.set_xlim(0, 0.07)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)

    # ── Panel D: Sample Correlation Heatmap (Artifact 4) ───────────────────
    ax_d1 = fig.add_subplot(gs[1, 2])
    ax_d2 = fig.add_subplot(gs[1, 3])
    panels_d = [(ax_d1, df_before, "Before Artifact Injection"), (ax_d2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_d:
        pivot = df_.pivot_table(index="cpg_id", columns="sample_id", values="beta_value")
        # Subset to a manageable number of samples for visibility
        cols = sorted(pivot.columns)[:12]
        if "S_DUP" in pivot.columns:
            cols = cols[:11] + ["S_DUP"]
        corr = pivot[cols].corr()
        sns.heatmap(corr, ax=ax, cmap="RdYlBu_r", vmin=0.7, vmax=1.0,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={"shrink": 0.7})
        subtitle = (
            "Pearson r, subset incl. duplicate pair" if "S_DUP" in pivot.columns
            else "Pearson r, 12-sample subset"
        )
        ax.set_title(f"{title}: Sample Correlation\n({subtitle})", fontsize=9)
        ax.tick_params(labelsize=5, rotation=45)

    # ── Panel E: Beta Distribution (Artifact 5) ────────────────────────────
    ax_e1 = fig.add_subplot(gs[2, 0:2])
    ax_e2 = fig.add_subplot(gs[2, 2:4])
    target_sid = "S020"
    ref_sid = "S005"
    panels_e = [(ax_e1, df_before, "Before Artifact Injection"), (ax_e2, df_after, "After Artifact Injection")]
    for ax, df_, title in panels_e:
        for sid, color, label in [
            (target_sid, "#E74C3C", f"{target_sid} (contaminated)"),
            (ref_sid, "#3498DB", f"{ref_sid} (reference)"),
        ]:
            sub = df_[df_["sample_id"] == sid]["beta_value"]
            if len(sub):
                ax.hist(sub, bins=50, alpha=0.55, color=color, label=label,
                        edgecolor="none", density=True)
        ax.set_title(f"{title}: Beta Distribution\n(contamination check)", fontsize=9)
        ax.set_xlabel("Beta value", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)

    out_path = os.path.join(FIGURES_DIR, "qc_before_after.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  [Visualization] Saved: {out_path}")


# =============================================================================
# 4.  MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("ImmuneMethylTools — Mock Data Generator")
    print("=" * 60)

    # ── Build manifest & baseline data ──────────────────────────────────────
    print("\n[1] Building patient manifest...")
    manifest = build_manifest()
    print(f"    Rows: {len(manifest):,}  |  Samples: {manifest['sample_id'].nunique()}"
          f"  |  CpGs per sample: {N_CPGS}")

    print("[2] Sampling baseline methylation...")
    df = add_baseline_methylation(manifest.copy())

    # ── Snapshot BEFORE artifacts ────────────────────────────────────────────
    df_before = df.copy()

    # ── Inject artifacts ─────────────────────────────────────────────────────
    print("\n[3] Injecting stumper artifacts...")
    df = inject_artifact1_confounded_batch(df)
    df = inject_artifact2_clonal_vdj(df)
    df = inject_artifact3_bisulfite_failure(df)
    df = inject_artifact5_contamination(df)
    df = inject_artifact6_low_depth(df)
    df = inject_true_biological_signal(df)
    df = inject_borderline_signal(df)
    df = inject_subtle_signal(df)
    df = inject_xci_signal(df)
    # Duplication AFTER XCI injection so S_DUP inherits the same X-linked
    # beta values as S010 — otherwise independent XCI re-injection would
    # destroy the high correlation on top-variance (sex-driven) CpGs.
    df = inject_artifact4_sample_duplication(df, manifest)
    df = inject_artifact7_sex_mixup(df)
    df = inject_artifact8_lineage(df)

    # ── Clip & round ─────────────────────────────────────────────────────────
    df["beta_value"] = df["beta_value"].clip(0.0, 1.0).round(4)
    df["non_cpg_meth_rate"] = df["non_cpg_meth_rate"].clip(0.0, 1.0).round(6)

    # ── Save CSV ─────────────────────────────────────────────────────────────
    print(f"\n[4] Saving CSV -> {OUT_CSV}")
    df.to_csv(OUT_CSV, index=False)
    print(f"    Final shape: {df.shape}  ({df['sample_id'].nunique()} samples)")

    # ── Summary statistics ────────────────────────────────────────────────────
    print("\n[5] Artifact Summary:")
    print(f"    Batch_01 Case mean beta   : "
          f"{df[(df.batch_id=='Batch_01') & (df.disease_label==CASE_LABEL)]['beta_value'].mean():.3f}")
    print(f"    Batch_01 Control mean beta: "
          f"{df[(df.batch_id=='Batch_01') & (df.disease_label==CTRL_LABEL)]['beta_value'].mean():.3f}")
    bis_fail = df.groupby("sample_id")["non_cpg_meth_rate"].mean()
    print(f"    Samples with non_cpg > 2% : "
          f"{(bis_fail > 0.02).sum()} -> {list(bis_fail[bis_fail > 0.02].index)}")
    if "S_DUP" in df["sample_id"].values:
        orig = df[df.sample_id == "S010"]["beta_value"].values
        clone = df[df.sample_id == "S_DUP"]["beta_value"].values
        r, _ = pearsonr(orig, clone)
        print(f"    S010 vs S_DUP Pearson r   : {r:.4f}")
    print(f"    S020 (contaminated) mean beta : "
          f"{df[df.sample_id == 'S020']['beta_value'].mean():.3f}")
    print(f"    S030 (low coverage) mean depth: "
          f"{df[df.sample_id == 'S030']['depth'].mean():.1f}x")
    s35_x_beta = df[(df.sample_id == "S035") & df["is_x_chromosome"].astype(bool)]["beta_value"].mean()
    s36_x_beta = df[(df.sample_id == "S036") & df["is_x_chromosome"].astype(bool)]["beta_value"].mean()
    print(f"    S035 reported sex='M', true X mean beta: {s35_x_beta:.3f} (expect ~0.50 — female XCI)")
    print(f"    S036 reported sex='F', true X mean beta: {s36_x_beta:.3f} (expect ~0.25 — male)")
    print(f"    is_x_chromosome: {df['is_x_chromosome'].sum()} X-linked rows  "
          f"({df['is_x_chromosome'].astype(bool).sum() // df['sample_id'].nunique()} per sample)")
    n_vdj = df["is_vdj_region"].sum()
    n_vdj_cpgs = df[df["is_vdj_region"]]["cpg_id"].nunique()
    print(f"    is_vdj_region: {n_vdj} VDJ rows  ({n_vdj_cpgs} unique CpGs)")
    print(f"    chrom values: {sorted(df['chrom'].unique())}")
    true_dmr_cpgs = [f"cg{i:08d}" for i in range(3000, 3011)]
    case_dmr_beta = df[(df["disease_label"] == "Case") & (df["cpg_id"].isin(true_dmr_cpgs))]["beta_value"].mean()
    ctrl_dmr_beta = df[(df["disease_label"] == "Control") & (df["cpg_id"].isin(true_dmr_cpgs))]["beta_value"].mean()
    print(f"    True DMR (cg3000-3010)  Case mean beta: {case_dmr_beta:.3f}  "
          f"Control mean beta: {ctrl_dmr_beta:.3f}  delta={case_dmr_beta - ctrl_dmr_beta:+.3f}")
    border_cpgs = [f"cg{i:08d}" for i in range(1500, 1508)]
    case_border = df[(df["disease_label"] == "Case") & (df["cpg_id"].isin(border_cpgs))]["beta_value"].mean()
    ctrl_border = df[(df["disease_label"] == "Control") & (df["cpg_id"].isin(border_cpgs))]["beta_value"].mean()
    print(f"    Borderline (cg1500-1507) Case mean beta: {case_border:.3f}  "
          f"Control mean beta: {ctrl_border:.3f}  delta={case_border - ctrl_border:+.3f}")
    subtle_cpgs = [f"cg{i:08d}" for i in range(2000, 2006)]
    case_subtle = df[(df["disease_label"] == "Case") & (df["cpg_id"].isin(subtle_cpgs))]["beta_value"].mean()
    ctrl_subtle = df[(df["disease_label"] == "Control") & (df["cpg_id"].isin(subtle_cpgs))]["beta_value"].mean()
    print(f"    Subtle (cg2000-2005)     Case mean beta: {case_subtle:.3f}  "
          f"Control mean beta: {ctrl_subtle:.3f}  delta={case_subtle - ctrl_subtle:+.3f}")
    foxp3_cpgs = [f"cg{i:08d}" for i in range(1, 6)]
    pax5_cpgs  = [f"cg{i:08d}" for i in range(6, 11)]
    for sid in ["S045", "S046"]:
        m = df[(df["sample_id"] == sid) & (df["cpg_id"].isin(foxp3_cpgs))]["beta_value"].mean()
        print(f"    {sid} (Treg-enriched)   FoxP3 mean beta: {m:.3f}  (expect ~0.06)")
    for sid in ["S065", "S066"]:
        m = df[(df["sample_id"] == sid) & (df["cpg_id"].isin(pax5_cpgs))]["beta_value"].mean()
        print(f"    {sid} (B-cell depleted) PAX5  mean beta: {m:.3f}  (expect ~0.65)")

    # ── Before/After visualization ────────────────────────────────────────────
    print("\n[6] Generating Before/After visualization...")
    # Align df_before to same columns for heatmap panel (no S_DUP row)
    plot_before_after(df_before, df)

    print("\n[Done] Mock data generation complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()
