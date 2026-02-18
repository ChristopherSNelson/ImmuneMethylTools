"""
core/ml_guard.py — ImmuneMethylTools ML Validator
==================================================
Safe cross-validated classification to quantify separability of Case vs.
Control methylation profiles with explicit data-leakage prevention.

Biological intent
-----------------
After QC, normalization, and DMR identification, a key validation step is
asking: "Can a regularised model distinguish cases from controls in truly
held-out data, while respecting patient structure?"

Two critical safeguards are applied:

  1. ElasticNet regularisation — combines L1 (feature selection; handles
     correlated CpGs) and L2 (shrinkage; prevents weight explosion in high-dim
     methylation feature space).  l1_ratio = 0.5 balances sparsity and
     stability.

  2. GroupKFold cross-validation — ensures NO patient appears in both
     training and test folds.  Without this, if a patient has technical
     replicates or multiple samples, the model memorises their methylation
     profile and inflates AUC — a textbook data-leakage pattern in methylation
     studies.

Interpretation
--------------
  AUC ≈ 0.50  → no separability (model is guessing)
  AUC ≈ 0.65  → modest biological signal (investigate DMRs further)
  AUC > 0.80  → strong separation — FIRST confirm this is not a batch or
                QC artefact (re-run after removing confounded batches)
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# ── Parameters ────────────────────────────────────────────────────────────────
N_SPLITS = 5       # GroupKFold folds
L1_RATIO = 0.5     # ElasticNet mix: 0 = Ridge, 1 = Lasso
C_PARAM = 1.0     # inverse regularisation strength
MAX_ITER = 5000    # saga needs generous budget for high-dim methylation features
N_TOP_CPGS = 200     # restrict to top-variance CpGs (speed + stability)


# =============================================================================
# Public API
# =============================================================================


def run_safe_model(
    df: pd.DataFrame,
    feature_col: str = "beta_value",
    n_splits: int = N_SPLITS,
    min_site_depth: int = 5,
) -> dict:
    """
    ElasticNet LogisticRegression with GroupKFold cross-validation.

    Parameters
    ----------
    df             : long-format methylation DataFrame; should be QC-filtered
                     and normalized before calling.
    feature_col    : methylation column to use as features
                     ('beta_value' or 'beta_normalized')
    n_splits       : number of GroupKFold folds (default: 5)
    min_site_depth : per-row minimum read depth; rows below are excluded before
                     the pivot (default: 5 — matches SITE_DEPTH_THRESH in qc_guard)

    Returns
    -------
    dict with keys:
        cv_results    : raw sklearn cross_validate output dict
        mean_auc      : float — mean ROC-AUC across folds
        std_auc       : float — std of ROC-AUC across folds
        mean_accuracy : float — mean accuracy across folds
        n_samples     : int
        n_features    : int
        warning       : str or None (class imbalance / insufficient data)
    """
    # ── Exclude low-depth sites ────────────────────────────────────────────────
    if min_site_depth > 0:
        df = df[df["depth"] >= min_site_depth].copy()

    # ── Build Sample × CpG matrix ─────────────────────────────────────────────
    pivot = df.pivot_table(index="sample_id", columns="cpg_id", values=feature_col)
    pivot = pivot.fillna(pivot.mean())

    # Select top-variance CpGs
    cpg_var = pivot.var(axis=0).sort_values(ascending=False)
    top_cpgs = cpg_var.head(N_TOP_CPGS).index
    X = pivot[top_cpgs].values

    # ── Labels and groups ─────────────────────────────────────────────────────
    meta = (
        df[["sample_id", "disease_label", "patient_id"]]
        .drop_duplicates("sample_id")
        .set_index("sample_id")
    )
    meta = meta.loc[pivot.index]
    y = (meta["disease_label"] == "Case").astype(int).values
    groups = meta["patient_id"].values

    # ── Warnings ───────────────────────────────────────────────────────────────
    case_frac = float(y.mean())
    warning = None
    n_groups = len(np.unique(groups))
    if case_frac < 0.20 or case_frac > 0.80:
        warning = (
            f"Class imbalance: Case fraction={case_frac:.1%}.  "
            "Consider class_weight='balanced' or SMOTE oversampling."
        )
    if n_groups < n_splits:
        warning = (
            f"Insufficient patient groups ({n_groups}) for {n_splits}-fold "
            "GroupKFold.  Reduce n_splits."
        )
        n_splits = max(2, n_groups)

    # ── Pipeline ───────────────────────────────────────────────────────────────
    clf = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "model",
                # solver='saga' supports l1_ratio (ElasticNet) natively.
                # penalty= was deprecated in sklearn 1.8; l1_ratio alone is sufficient.
                # Keeping it here to avoid silent errors on older installs.
                LogisticRegression(
                    solver="saga",
                    penalty="elasticnet",  # here for old scikit defense
                    l1_ratio=L1_RATIO,
                    C=C_PARAM,
                    max_iter=MAX_ITER,
                    random_state=42,
                ),
            ),
        ]
    )

    cv = GroupKFold(n_splits=n_splits)
    scoring = {
        # Use string scorer "roc_auc" — handles predict_proba routing automatically
        # and is compatible across sklearn versions.
        "roc_auc": "roc_auc",
        "accuracy": "accuracy",
    }

    cv_results = cross_validate(
        clf, X, y, groups=groups, cv=cv,
        scoring=scoring, return_train_score=False,
        error_score=np.nan,
    )

    return {
        "cv_results": cv_results,
        "mean_auc": float(np.nanmean(cv_results["test_roc_auc"])),
        "std_auc": float(np.nanstd(cv_results["test_roc_auc"])),
        "mean_accuracy": float(np.nanmean(cv_results["test_accuracy"])),
        "n_samples": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "warning": warning,
    }


# =============================================================================
# __main__ — run safe model on clean, normalized mock data
# =============================================================================

if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))          # for io_utils
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # for core.*
    from core.normalizer import robust_normalize
    from core.qc_guard import audit_quality
    from io_utils import data_path, load_methylation, project_root, write_audit_log  # noqa: E402

    MODULE = "ML_GUARD"
    _now = datetime.now()
    ts_tag = _now.strftime("%Y%m%d_%H%M%S")
    _base = project_root()
    _audit_csv = os.path.join(_base, "output", f"audit_log_{MODULE}_{ts_tag}.csv")

    audit_entries = []

    def ts():
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def ae(sample_id, status, description, metric):
        return {
            "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
            "module": MODULE,
            "sample_id": sample_id,
            "status": status,
            "description": description,
            "metric": metric,
        }

    df = load_methylation(data_path("mock_methylation.csv"))

    # ── QC gate ────────────────────────────────────────────────────────────────
    clean_samples = audit_quality(df)
    df_clean = df[df["sample_id"].isin(clean_samples)].copy()
    print(f"[{ts()}] [ML_GUARD]           | Clean samples: n={len(clean_samples)}")
    audit_entries.append(ae(
        "cohort", "INFO", "Clean samples loaded for ML validation",
        f"n={len(clean_samples)}",
    ))

    # ── Normalize ──────────────────────────────────────────────────────────────
    df_norm = robust_normalize(df_clean, save_figure=False)
    audit_entries.append(ae(
        "cohort", "INFO", "Median-centring normalization applied",
        f"n_samples={len(clean_samples)}",
    ))

    # ── Run model ──────────────────────────────────────────────────────────────
    results = run_safe_model(df_norm, feature_col="beta_normalized")

    print(
        f"[{ts()}] [ML_GUARD] DETECTED | ElasticNet + GroupKFold CV | "
        f"AUC={results['mean_auc']:.4f} ± {results['std_auc']:.4f}  "
        f"Accuracy={results['mean_accuracy']:.4f}  "
        f"n_samples={results['n_samples']}  n_features={results['n_features']}"
    )
    audit_entries.append(ae(
        "cohort", "INFO", "ElasticNet + GroupKFold CV complete",
        f"AUC={results['mean_auc']:.4f}±{results['std_auc']:.4f}",
    ))

    if results["warning"]:
        print(f"[{ts()}] [ML_GUARD] WARNING | {results['warning']}")
        audit_entries.append(ae(
            "cohort", "DETECTED", "ML model warning",
            results["warning"][:80],
        ))

    per_fold = results["cv_results"]["test_roc_auc"]
    for i, auc in enumerate(per_fold, 1):
        auc_str = f"{auc:.4f}" if not np.isnan(auc) else "nan"
        print(f"[{ts()}] [ML_GUARD]           | Fold {i}: AUC={auc_str}")
        audit_entries.append(ae(
            "cohort", "INFO", f"GroupKFold fold {i} AUC", auc_str,
        ))

    write_audit_log(audit_entries, _audit_csv)
    print(f"[{ts()}] [ML_GUARD] Audit log → {_audit_csv}")
