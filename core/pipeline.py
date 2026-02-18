"""
core/pipeline.py — ImmuneMethylTools End-to-End Pipeline
=========================================================
Runs all artifact detector modules in the correct order, passing the
clean_samples list automatically from one stage to the next.

Stage order (CLAUDE.md architecture rule: QC before normalization before DMR):
  1a. QC Guard       — bisulfite / depth failures
  1b. QC Guard       — contamination detection
  2.  Sample Audit   — technical duplicate removal (retains sample_a, drops sample_b)
  2.5 QC Guard       — site-level depth filter (rows with depth < 5 removed from df_clean)
  3.  Repertoire     — clonal VDJ artifact flagging; Stage 3.5 masks VDJ-locus
                       beta_value to NaN in flagged samples before normalization

Output is written to output/ (logs, figures, audit CSVs, clean data, reports).
Input data (mock_methylation.csv) remains in data/.
  4.  Normalizer     — batch × disease confound check + median-centring
  5.  Deconvolution  — cell-fraction estimation + lineage shift detection
  6.  DMR Hunter     — sliding-window Wilcoxon DMR caller on normalized clean data
  7.  ML Guard       — ElasticNet + GroupKFold classification validation

The clean_samples list is updated after each filtering stage (1a, 1b, 2) and
propagated automatically to all downstream modules.
"""

import os
import sys
from datetime import datetime

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from deconvolution import detect_lineage_shift, estimate_cell_fractions  # noqa: E402
from dmr_hunter import find_dmrs  # noqa: E402
from io_utils import (  # noqa: E402
    Tee, append_flagged_samples, data_path, load_methylation, project_root, write_audit_log,
)
from ml_guard import run_safe_model  # noqa: E402
from normalizer import check_confounding, robust_normalize  # noqa: E402
from visuals import plot_exclusion_accounting, plot_volcano  # noqa: E402
from qc_guard import (  # noqa: E402
    BISULFITE_FAIL_THRESH,
    SITE_DEPTH_THRESH,
    SITE_LOW_DEPTH_SAMPLE_WARN,
    audit_quality,
    detect_contamination,
    filter_site_quality,
)
from repertoire_clonality import flag_clonal_artifacts, mask_clonal_vdj_sites  # noqa: E402
from sample_audit import detect_duplicates  # noqa: E402


# =============================================================================
# Public API
# =============================================================================


def run_pipeline(csv_path: str, save_figures: bool = True, save_report: bool = False) -> dict:
    """
    Execute the full ImmuneMethylTools artifact detection and analysis pipeline.

    Parameters
    ----------
    csv_path    : path to mock_methylation.csv (or real data matching the schema)
    save_figures: passed through to normalizer.robust_normalize
    save_report : if True, generate a PDF report at the end of the run

    Returns
    -------
    dict with keys:
        clean_samples   list[str]       final QC-passed, deduped sample IDs
        n_qc_failed     int             samples dropped by bisulfite/depth QC
        n_contaminated  int             samples dropped as contaminated
        n_deduped       int             samples dropped as technical duplicates
        confounded      bool            batch × disease Cramér's V > 0.5
        cramers_v       float
        n_clonal_rows   int             VDJ CpG rows bearing clonal signature
        dmrs            pd.DataFrame    full window table from DMR hunter
        n_sig_dmrs      int
        mean_auc        float
        std_auc         float
        df_norm         pd.DataFrame    clean, median-centred DataFrame
    """
    _base = project_root()
    _now = datetime.now()
    run_ts = _now.strftime("%Y-%m-%dT%H:%M:%S")
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _log = os.path.join(_base, "output", "logs", f"pipeline_{ts_tag}.log")
    _flag_csv = os.path.join(_base, "output", "flagged_samples.csv")
    _audit_csv = os.path.join(_base, "output", f"audit_log_pipeline_{ts_tag}.csv")
    os.makedirs(os.path.join(_base, "output", "logs"), exist_ok=True)

    audit_entries = []
    flagged_rows = []

    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def ae(module, sample_id, status, description, metric):
        return {
            "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            "module": module,
            "sample_id": sample_id,
            "status": status,
            "description": description,
            "metric": metric,
        }

    with Tee(_log):
        banner = "=" * 68
        print(f"\n{banner}")
        print(f"  ImmuneMethylTools Pipeline  —  {run_ts}")
        print(f"{banner}\n")

        # ── Load ───────────────────────────────────────────────────────────────
        print(f"[{ts()}] [PIPELINE] Loading {csv_path}")
        df = load_methylation(csv_path, verbose=False)
        all_samples = df["sample_id"].unique().tolist()
        n_total = len(all_samples)
        print(f"[{ts()}] [PIPELINE]           | Cohort size: n={n_total} samples")
        audit_entries.append(ae("PIPELINE", "cohort", "INFO", "Cohort loaded", f"n={n_total}"))

        # ── Stage 1a: bisulfite / depth QC ────────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 1a: QC Guard (bisulfite/depth) ──")
        clean_samples = audit_quality(df)
        qc_failed = [s for s in all_samples if s not in clean_samples]
        n_qc_failed = len(qc_failed)

        sample_stats = df.groupby("sample_id").agg(
            mean_ncpg=("non_cpg_meth_rate", "mean"),
            mean_depth=("depth", "mean"),
        )
        for sid in qc_failed:
            ncpg = float(sample_stats.loc[sid, "mean_ncpg"])
            depth = float(sample_stats.loc[sid, "mean_depth"])
            if ncpg > BISULFITE_FAIL_THRESH:
                flag_type = "bisulfite_failure"
                metric = f"non_cpg={ncpg:.3f}"
                reason = f"non_cpg_meth_rate={ncpg:.3f}"
            else:
                flag_type = "depth_failure"
                metric = f"mean_depth={depth:.1f}"
                reason = f"mean_depth={depth:.1f}"
            print(
                f"[{ts()}] [PIPELINE] DETECTED | {flag_type} | "
                f"sample={sid}  {reason}"
            )
            audit_entries.append(ae("QC_GUARD", sid, "DETECTED", flag_type, metric))
            flagged_rows.append({
                "run_timestamp": run_ts, "module": "qc_guard",
                "sample_id": sid, "flag_type": flag_type, "detail": reason,
            })

        print(
            f"[{ts()}] [PIPELINE]           | QC: {n_qc_failed} failed, "
            f"{len(clean_samples)} passed"
        )
        audit_entries.append(ae(
            "QC_GUARD", "cohort", "INFO",
            "Bisulfite/depth QC complete",
            f"passed={len(clean_samples)} failed={n_qc_failed}",
        ))

        # ── Stage 1b: contamination ────────────────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 1b: QC Guard (contamination) ──")
        contaminated, contam_report = detect_contamination(df)
        n_contaminated = 0
        for sid in contaminated:
            if sid in clean_samples:
                clean_samples.remove(sid)
                n_contaminated += 1
            row = contam_report.loc[sid]
            print(
                f"[{ts()}] [PIPELINE] DETECTED | contamination | "
                f"sample={sid}  BC={row.bimodality_coeff:.4f}  "
                f"mean_beta={row.mean_beta:.4f}"
            )
            audit_entries.append(ae(
                "QC_GUARD", sid, "DETECTED",
                "Contamination — low bimodality coefficient",
                f"BC={row.bimodality_coeff:.4f}",
            ))
            flagged_rows.append({
                "run_timestamp": run_ts, "module": "qc_guard",
                "sample_id": sid, "flag_type": "contamination",
                "detail": (
                    f"BC={row.bimodality_coeff:.4f} "
                    f"mean_beta={row.mean_beta:.4f}"
                ),
            })

        if n_contaminated == 0:
            print(f"[{ts()}] [PIPELINE]           | No contamination detected")
        print(
            f"[{ts()}] [PIPELINE]           | Contamination: {n_contaminated} flagged, "
            f"{len(clean_samples)} remaining"
        )
        audit_entries.append(ae(
            "QC_GUARD", "cohort", "INFO",
            "Contamination check complete",
            f"flagged={n_contaminated}",
        ))

        # ── Stage 2: duplicate detection ──────────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 2: Sample Audit (duplicates) ──")
        df_qc = df[df["sample_id"].isin(clean_samples)].copy()
        dup_result, ids_to_drop = detect_duplicates(df_qc)
        dup_pairs = dup_result[dup_result["duplicate_flag"]]
        deduped = set()

        for drop_sid in ids_to_drop:
            if drop_sid in clean_samples and drop_sid not in deduped:
                clean_samples.remove(drop_sid)
                deduped.add(drop_sid)

        for _, row in dup_pairs.iterrows():
            print(
                f"[{ts()}] [PIPELINE] DETECTED | duplicate pair | "
                f"{row.sample_a} ↔ {row.sample_b}  r={row.pearson_r:.6f}  "
                f"→ dropping {row.sample_b}"
            )
            audit_entries.append(ae(
                "SAMPLE_AUDIT", row.sample_b, "DETECTED",
                f"Technical duplicate of {row.sample_a} — dropped",
                f"r={row.pearson_r:.4f}",
            ))
            flagged_rows.append({
                "run_timestamp": run_ts, "module": "sample_audit",
                "sample_id": row.sample_b, "flag_type": "duplicate",
                "detail": f"r={row.pearson_r:.4f} paired with {row.sample_a} — dropped",
            })

        n_deduped = len(deduped)
        if n_deduped == 0:
            print(f"[{ts()}] [PIPELINE]           | No duplicates detected")
        print(
            f"[{ts()}] [PIPELINE]           | Dedup: {n_deduped} dropped, "
            f"{len(clean_samples)} remaining → final clean_samples list"
        )
        audit_entries.append(ae(
            "SAMPLE_AUDIT", "cohort", "INFO",
            "Duplicate check complete",
            f"dropped={n_deduped} remaining={len(clean_samples)}",
        ))

        # All downstream stages use df_clean filtered to the final clean_samples list
        df_clean = df[df["sample_id"].isin(clean_samples)].copy()

        # ── Exclusion accounting figure ────────────────────────────────────────
        if save_figures:
            excl_path = plot_exclusion_accounting(
                n_total,
                [
                    ("Bisulfite/Depth QC", n_qc_failed),
                    ("Contamination", n_contaminated),
                    ("Duplicate", n_deduped),
                ],
            )
            print(f"[{ts()}] [PIPELINE]           | Exclusion accounting figure → {excl_path}")
            audit_entries.append(ae(
                "PIPELINE", "cohort", "INFO",
                "Exclusion accounting figure saved", excl_path,
            ))

        # ── Stage 2.5: site-level depth QC ───────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 2.5: Site-level Depth QC ──")
        df_clean, site_stats = filter_site_quality(df_clean, min_depth=SITE_DEPTH_THRESH)
        pct_low = site_stats["pct_low"]
        print(
            f"[{ts()}] [PIPELINE]           | Site filter: "
            f"{site_stats['n_low']}/{site_stats['n_total']} rows below "
            f"{SITE_DEPTH_THRESH}x ({pct_low:.1f}%)"
        )
        audit_entries.append(ae(
            "QC_GUARD", "cohort", "INFO",
            "Site-level depth filtering complete",
            f"{pct_low:.1f}% sites < {SITE_DEPTH_THRESH}x dropped",
        ))
        for sid, pct_s in site_stats["per_sample_pct"].items():
            if pct_s > SITE_LOW_DEPTH_SAMPLE_WARN:
                print(
                    f"[{ts()}] [PIPELINE] DETECTED | high low-depth site rate | "
                    f"sample={sid}  {pct_s:.1f}% sites < {SITE_DEPTH_THRESH}x"
                )
                audit_entries.append(ae(
                    "QC_GUARD", sid, "DETECTED",
                    "High proportion of low-depth sites",
                    f"{pct_s:.1f}% sites < {SITE_DEPTH_THRESH}x",
                ))

        # ── Stage 3: clonality (informational) ────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 3: Repertoire Clonality (informational) ──")
        clonal_rows, clonal_samples = flag_clonal_artifacts(df_clean)
        n_clonal_rows = len(clonal_rows)

        if clonal_samples:
            for sid in clonal_samples:
                sub = clonal_rows[clonal_rows["sample_id"] == sid]
                print(
                    f"[{ts()}] [PIPELINE] DETECTED | clonal VDJ artifact | "
                    f"sample={sid}  n_rows={len(sub)}  "
                    f"mean_beta={sub['beta_value'].mean():.3f}  "
                    f"mean_frag={sub['fragment_length'].mean():.0f} bp"
                )
                audit_entries.append(ae(
                    "CLONALITY", sid, "DETECTED",
                    "Clonal VDJ artifact — VDJ CpGs masked (beta=NaN) before normalization and ML",
                    f"n_rows={len(sub)} mean_beta={sub['beta_value'].mean():.3f}",
                ))
        else:
            print(f"[{ts()}] [PIPELINE]           | No clonal VDJ artifacts detected")

        audit_entries.append(ae(
            "CLONALITY", "cohort", "INFO",
            "Clonality scan complete",
            f"n_clonal_rows={n_clonal_rows}",
        ))

        # ── Stage 3.5: mask clonal VDJ sites ──────────────────────────────────
        if clonal_samples:
            df_clean, n_masked = mask_clonal_vdj_sites(df_clean, clonal_samples)
            print(
                f"[{ts()}] [PIPELINE]           | VDJ sites masked "
                f"(beta->NaN) in {len(clonal_samples)} sample(s) | "
                f"n_sites={n_masked}"
            )
            audit_entries.append(ae(
                "CLONALITY", "cohort", "INFO",
                "VDJ-locus beta values masked to NaN in clonal sample(s)",
                f"n_masked_sites={n_masked}",
            ))

        # ── Stage 4: normalization ─────────────────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 4: Normalizer ──")
        confound = check_confounding(df_clean, "batch_id", "disease_label")
        v_status = "DETECTED" if confound["confounded"] else "INFO    "
        print(
            f"[{ts()}] [PIPELINE] {v_status} | batch × disease confound | "
            f"Cramér's V={confound['cramers_v']:.4f}  p={confound['p_value']:.4e}"
        )
        audit_entries.append(ae(
            "NORMALIZER", "cohort",
            "DETECTED" if confound["confounded"] else "INFO",
            "Batch × disease confound check"
            + (" — CONFOUNDED" if confound["confounded"] else " — OK"),
            f"V={confound['cramers_v']:.4f} p={confound['p_value']:.4e}",
        ))

        df_norm = robust_normalize(df_clean, save_figure=save_figures)
        post_std = df_norm.groupby("sample_id")["beta_normalized"].std()
        print(
            f"[{ts()}] [PIPELINE]           | Median-centring applied | "
            f"post-norm std [{post_std.min():.4f}, {post_std.max():.4f}]"
        )
        audit_entries.append(ae(
            "NORMALIZER", "cohort", "INFO",
            "Median-centring normalization applied",
            f"std_range=[{post_std.min():.4f},{post_std.max():.4f}]",
        ))

        # ── Stage 5: deconvolution (informational) ────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 5: Deconvolution (informational) ──")
        fracs = estimate_cell_fractions(df_clean)

        # Join disease_label for sorting and group summaries
        meta = df_clean[["sample_id", "disease_label"]].drop_duplicates("sample_id")
        fracs = fracs.merge(meta, on="sample_id")
        fracs["b_t_ratio"] = (fracs["b_fraction"] / fracs["t_fraction"].replace(0, float("nan"))).round(3)
        fracs_sorted = fracs.sort_values(["disease_label", "sample_id"]).reset_index(drop=True)

        mean_b = float(fracs["b_fraction"].mean())
        mean_t = float(fracs["t_fraction"].mean())
        mean_treg = float(fracs["treg_fraction"].mean())
        print(
            f"[{ts()}] [PIPELINE]           | Cohort cell fractions: "
            f"B={mean_b:.3f}  T={mean_t:.3f}  Treg={mean_treg:.3f}"
        )
        audit_entries.append(ae(
            "DECONVOLVE", "cohort", "INFO",
            "Cell fractions estimated",
            f"mean_B={mean_b:.3f} mean_T={mean_t:.3f} mean_Treg={mean_treg:.3f}",
        ))

        # Per-group means
        for label, grp in fracs_sorted.groupby("disease_label", sort=True):
            print(
                f"[{ts()}] [PIPELINE]           | {label:8s} mean: "
                f"B={grp['b_fraction'].mean():.3f}  "
                f"T={grp['t_fraction'].mean():.3f}  "
                f"Treg={grp['treg_fraction'].mean():.3f}  "
                f"B:T={grp['b_t_ratio'].mean():.2f}"
            )

        # Per-sample table sorted by Case / Control
        header = (
            f"\n{'sample_id':<10} {'label':<8} {'B':>6} {'T':>6} {'Treg':>6} {'Other':>6} {'B:T':>6}"
        )
        print(f"[{ts()}] [PIPELINE]           | Per-sample cell fractions:{header}")
        sep = "-" * 55
        print(f"[{ts()}] [PIPELINE]           | {sep}")
        for _, row in fracs_sorted.iterrows():
            print(
                f"[{ts()}] [PIPELINE]           | "
                f"{row.sample_id:<10} {row.disease_label:<8} "
                f"{row.b_fraction:>6.3f} {row.t_fraction:>6.3f} "
                f"{row.treg_fraction:>6.3f} {row.other_fraction:>6.3f} "
                f"{row.b_t_ratio:>6.2f}"
            )
        print(f"[{ts()}] [PIPELINE]           | {sep}")

        shifts = detect_lineage_shift(df_clean)
        ls_flagged = shifts[shifts["any_lineage_flag"]]
        if len(ls_flagged):
            for _, row in ls_flagged.iterrows():
                parts = []
                if row.treg_flag:
                    parts.append(f"FoxP3 β={row.foxp3_mean_beta:.3f}")
                if row.bcell_shift_flag:
                    parts.append(f"PAX5 β={row.pax5_mean_beta:.3f}")
                detail = " | ".join(parts)
                print(
                    f"[{ts()}] [PIPELINE] DETECTED | lineage shift | "
                    f"sample={row.sample_id}  {detail}"
                )
                audit_entries.append(ae(
                    "DECONVOLVE", row.sample_id, "DETECTED",
                    "Lineage shift detected", " ".join(parts),
                ))
        else:
            print(f"[{ts()}] [PIPELINE]           | No lineage shifts detected")

        # ── Stage 6: DMR Hunter ───────────────────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 6: DMR Hunter ──")
        dmrs = find_dmrs(df_norm, clean_samples, normalized_col="beta_normalized")
        sig_dmrs = dmrs[dmrs["significant"]]
        n_sig = len(sig_dmrs)
        print(
            f"[{ts()}] [PIPELINE]           | {n_sig} significant DMRs "
            f"of {len(dmrs)} windows tested"
        )
        audit_entries.append(ae(
            "DMR_HUNTER", "cohort",
            "DETECTED" if n_sig else "INFO",
            "Sliding-window DMR scan complete",
            f"n_sig={n_sig} of {len(dmrs)} windows",
        ))
        for _, row in sig_dmrs.head(5).iterrows():
            risk_tag = " ⚠ HIGH CLONALITY" if row.clonal_risk else ""
            print(
                f"[{ts()}] [PIPELINE] DETECTED | DMR{risk_tag} | "
                f"{row.window_id}  ΔBeta={row.delta_beta:+.4f}  "
                f"p_adj={row.p_adj:.3e}  n_cpgs={row.n_cpgs}  n_vdj_cpgs={row.n_vdj_cpgs}"
            )

        # Report clonal-risk significant DMRs separately for analyst review
        sig_clonal = sig_dmrs[sig_dmrs["clonal_risk"]] if "clonal_risk" in sig_dmrs.columns else pd.DataFrame()
        clonal_all = dmrs[dmrs["clonal_risk"]] if "clonal_risk" in dmrs.columns else pd.DataFrame()
        print(
            f"[{ts()}] [PIPELINE]           | VDJ clonal_risk windows: "
            f"{len(clonal_all)} total, {len(sig_clonal)} significant"
        )
        if save_figures:
            v1 = plot_volcano(dmrs, color_clonal_risk=False)
            v2 = plot_volcano(dmrs, color_clonal_risk=True)
            print(f"[{ts()}] [PIPELINE]           | Volcano (standard)      → {v1}")
            print(f"[{ts()}] [PIPELINE]           | Volcano (clonal risk)   → {v2}")
            audit_entries.append(ae("DMR_HUNTER", "cohort", "INFO", "Volcano plots saved",
                                    "volcano.png volcano_clonal_risk.png"))
        for _, row in sig_clonal.iterrows():
            print(
                f"[{ts()}] [PIPELINE] DETECTED | HIGH CLONALITY — sig DMR in VDJ locus | "
                f"{row.window_id}  ΔBeta={row.delta_beta:+.4f}  n_vdj_cpgs={row.n_vdj_cpgs}"
            )
            audit_entries.append(ae(
                "DMR_HUNTER", "cohort", "DETECTED",
                f"HIGH CLONALITY — significant DMR overlaps VDJ locus: {row.window_id}",
                f"delta_beta={row.delta_beta:+.4f} n_vdj_cpgs={row.n_vdj_cpgs}",
            ))

        # ── Stage 7: ML Guard ─────────────────────────────────────────────────
        print(f"\n[{ts()}] [PIPELINE] ── Stage 7: ML Guard ──")
        ml = run_safe_model(df_norm, feature_col="beta_normalized")
        print(
            f"[{ts()}] [PIPELINE]           | ElasticNet + GroupKFold | "
            f"AUC={ml['mean_auc']:.4f} ± {ml['std_auc']:.4f}  "
            f"Accuracy={ml['mean_accuracy']:.4f}  "
            f"n_samples={ml['n_samples']}  n_features={ml['n_features']}"
        )
        audit_entries.append(ae(
            "ML_GUARD", "cohort", "INFO",
            "ElasticNet GroupKFold CV complete",
            f"AUC={ml['mean_auc']:.4f}±{ml['std_auc']:.4f}",
        ))
        if ml["warning"]:
            print(f"[{ts()}] [PIPELINE] WARNING | ML | {ml['warning']}")
            audit_entries.append(ae(
                "ML_GUARD", "cohort", "DETECTED",
                "ML model warning", ml["warning"][:80],
            ))

        # ── Summary ───────────────────────────────────────────────────────────
        confound_label = (
            f"YES  (Cramér's V={confound['cramers_v']:.4f})"
            if confound["confounded"] else "no"
        )
        print(f"\n{banner}")
        print("  Pipeline Summary")
        print(f"{banner}")
        print(f"  Input samples              : {n_total}")
        print(f"  Bisulfite/depth failures   : {n_qc_failed}")
        print(f"  Contamination flagged      : {n_contaminated}")
        print(f"  Technical duplicates removed: {n_deduped}")
        print(f"  Final clean samples        : {len(clean_samples)}")
        print(f"  Batch confound             : {confound_label}")
        print(f"  Clonal VDJ rows (masked)   : {n_clonal_rows}")
        print(f"  Significant DMRs           : {n_sig}")
        print(f"  Classification AUC         : {ml['mean_auc']:.4f} ± {ml['std_auc']:.4f}")
        print("  Clean data export          : output/clean_methylation.csv")
        print(f"  Audit log                  : output/audit_log_pipeline_{ts_tag}.csv")
        print(f"  Run log                    : output/logs/pipeline_{ts_tag}.log")
        print(f"{banner}\n")

        # ── Export clean data ─────────────────────────────────────────────────
        clean_csv = os.path.join(_base, "output", "clean_methylation.csv")
        os.makedirs(os.path.dirname(clean_csv), exist_ok=True)
        df_clean.to_csv(clean_csv, index=False)
        print(f"[{ts()}] [PIPELINE] Clean data  → {clean_csv}")
        audit_entries.append(ae(
            "PIPELINE", "cohort", "INFO",
            "Clean methylation data exported",
            f"n_samples={len(clean_samples)} path=output/clean_methylation.csv",
        ))

        # ── Persist ───────────────────────────────────────────────────────────
        append_flagged_samples(flagged_rows, _flag_csv)
        write_audit_log(audit_entries, _audit_csv)
        print(f"[{ts()}] [PIPELINE] Audit log  → {_audit_csv}")
        print(f"[{ts()}] [PIPELINE] Flag log   → {_flag_csv}")
        print(f"[{ts()}] [PIPELINE] Run log    → {_log}")

        # ── Build result dict ─────────────────────────────────────────────────
        result = {
            "clean_samples": clean_samples,
            "n_total": n_total,
            "n_qc_failed": n_qc_failed,
            "n_contaminated": n_contaminated,
            "n_deduped": n_deduped,
            "confounded": confound["confounded"],
            "cramers_v": confound["cramers_v"],
            "n_clonal_rows": n_clonal_rows,
            "dmrs": dmrs,
            "n_sig_dmrs": n_sig,
            "mean_auc": ml["mean_auc"],
            "std_auc": ml["std_auc"],
            "df_norm": df_norm,
            "clean_csv": clean_csv,
            "audit_csv": _audit_csv,
            "run_ts": run_ts,
        }

        # ── Optional PDF report ───────────────────────────────────────────────
        if save_report:
            from report_gen import generate_report
            _report_path = os.path.join(_base, "output", f"report_{ts_tag}.pdf")
            generate_report(result, _audit_csv, _report_path, run_ts)
            result["report_path"] = _report_path
            print(f"[{ts()}] [PIPELINE] Report    → {_report_path}")

    return result


# =============================================================================
# __main__
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="ImmuneMethylTools -- full artifact detection pipeline."
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Generate a PDF report at the end of the run.",
    )
    parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Skip figure generation (faster for debugging).",
    )
    args = parser.parse_args()

    run_pipeline(
        data_path("mock_methylation.csv"),
        save_figures=not args.no_figures,
        save_report=args.report,
    )
