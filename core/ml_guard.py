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

  1. Logit transform (beta → M-value) — raw beta values are bounded [0,1]
     with a bimodal distribution unsuitable for linear models.  log2(β/(1−β))
     maps them to ℝ with better statistical properties (Du et al. 2010).
     Applied inside the sklearn Pipeline so no leakage occurs across CV folds.

  2. ElasticNet regularisation — combines L1 (feature selection; handles
     correlated CpGs) and L2 (shrinkage; prevents weight explosion in high-dim
     methylation feature space).  l1_ratio = 0.5 balances sparsity and
     stability.

  3. GroupKFold cross-validation — ensures NO patient appears in both
     training and test folds.  Without this, if a patient has technical
     replicates or multiple samples, the model memorises their methylation
     profile and inflates AUC — a textbook data-leakage pattern in methylation
     studies.

Interpretation
--------------
  AUC ≈ 0.50  → no separability (model is guessing)
  AUC ≈ 0.65  → modest biological signal (investigate DMRs further)
  AUC > 0.80  → strong separation — FIRST confirm this is not a batch or
                QC artifact (re-run after removing confounded batches)
"""

import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler

# ── Parameters ────────────────────────────────────────────────────────────────
N_SPLITS = 5       # GroupKFold folds
L1_RATIO = 0.5     # ElasticNet mix: 0 = Ridge, 1 = Lasso
C_PARAM = 1.0     # inverse regularisation strength
MAX_ITER = 5000    # saga needs generous budget for high-dim methylation features
N_TOP_CPGS = 200   # restrict to top-variance CpGs (speed + stability)
LOGIT_EPSILON = 1e-4  # clip boundary before logit to avoid ±inf


# ── Logit transform (beta → M-value) ──────────────────────────────────────────

def _logit_transform(X: np.ndarray) -> np.ndarray:
    """Map beta values to M-values via logit: log2(β / (1 − β)).

    Beta values are clipped to [LOGIT_EPSILON, 1 − LOGIT_EPSILON] before
    transformation to prevent ±inf at the boundaries.  Du et al. (2010)
    showed M-values have better statistical properties for linear modelling
    compared to bounded beta values.
    """
    X_clipped = np.clip(X, LOGIT_EPSILON, 1 - LOGIT_EPSILON)
    return np.log2(X_clipped / (1 - X_clipped))


# =============================================================================
# Public API
# =============================================================================


def run_safe_model(
    df: pd.DataFrame,
    feature_col: str = "beta_value",
    n_splits: int = N_SPLITS,
    min_site_depth: int = 5,
    logit_transform: bool = True,
) -> dict:
    """
    ElasticNet LogisticRegression with GroupKFold cross-validation.

    Parameters
    ----------
    df              : long-format methylation DataFrame; should be QC-filtered
                      and normalized before calling.
    feature_col     : methylation column to use as features when
                      logit_transform=False ('beta_value' or 'beta_normalized').
                      Ignored when logit_transform=True — raw beta_value is
                      always used as the logit input because beta_normalized can
                      exceed [0,1] after median-centring, making logit undefined.
    n_splits        : number of GroupKFold folds (default: 5)
    min_site_depth  : per-row minimum read depth; rows below are excluded before
                      the pivot (default: 5 — matches SITE_DEPTH_THRESH in qc_guard)
    logit_transform : if True (default), apply log2(β/(1−β)) to convert beta
                      values to M-values before scaling.  M-values are unbounded
                      and have better statistical properties for linear models
                      (Du et al. 2010).

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
    # When logit_transform=True always pivot on raw beta_value — beta_normalized
    # can exceed [0,1] after median-centring, making logit undefined on it.
    pivot_col = "beta_value" if logit_transform else feature_col
    pivot = df.pivot_table(index="sample_id", columns="cpg_id", values=pivot_col)
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
    steps = []
    if logit_transform:
        # Convert beta → M-value before scaling.  Applied per CV fold inside
        # the Pipeline so no information leaks from test to train.
        steps.append(("logit", FunctionTransformer(_logit_transform)))
    steps.extend([
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
    ])
    clf = Pipeline(steps)

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

    # ── QC filter ──────────────────────────────────────────────────────────────
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
