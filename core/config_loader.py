"""
core/config_loader.py — ImmuneMethylTools Analysis Threshold Configuration
===========================================================================
Loads analysis thresholds from config.json at the project root.
All values fall back to hard-coded defaults if the file is absent or a key
is missing, so the pipeline runs identically with or without config.json.

Usage
-----
In pipeline.py::

    from config_loader import load_config
    cfg = load_config()                         # auto-discovers config.json
    cfg = load_config("path/to/custom.json")    # explicit path

In individual modules, functions accept threshold kwargs — pass values from
cfg to override without editing source files::

    audit_quality(df, bisulfite_thresh=cfg["qc"]["bisulfite_fail_thresh"])

Config file format (config.json)
---------------------------------
{
  "qc": {
    "bisulfite_fail_thresh":      0.01,
    "depth_fail_thresh":          10,
    "site_depth_thresh":          5,
    "site_low_depth_sample_warn": 20.0,
    "bc_sigma_thresh":            2.0,
    "contamination_mean_lo":      0.40,
    "contamination_mean_hi":      0.65
  },
  "duplicates": {
    "corr_thresh": 0.99
  },
  "clonality": {
    "beta_min":        0.80,
    "frag_sd_thresh":  3.0,
    "min_locus_hits":  3
  },
  "dmr": {
    "p_adj_thresh":   0.05,
    "delta_beta_min": 0.10,
    "chunk_size":     null
  },
  "ml": {
    "n_top_cpgs": 200,
    "l1_ratio":   0.5,
    "c_param":    1.0,
    "chunk_size": null
  }
}
"""

from __future__ import annotations

import json
import os

# ── Default thresholds (match module-level constants) ─────────────────────────
_DEFAULTS: dict[str, dict] = {
    "qc": {
        "bisulfite_fail_thresh":     0.01,   # non-CpG meth rate; > this → bisulfite failure
        "depth_fail_thresh":         10,     # mean sample depth; < this → low coverage
        "site_depth_thresh":         5,      # per-CpG row depth; < this → excluded
        "site_low_depth_sample_warn": 20.0,  # % low-depth rows per sample → log warning
        "bc_sigma_thresh":           2.0,    # BC z-score below cohort median → contamination
        "contamination_mean_lo":     0.40,   # muddy beta lower bound
        "contamination_mean_hi":     0.65,   # muddy beta upper bound
    },
    "duplicates": {
        "corr_thresh": 0.99,   # Pearson r ≥ this → technical duplicate
    },
    "clonality": {
        "beta_min": 0.80,        # VDJ beta above this → clonal hypermethylation
        "frag_sd_thresh": 3.0,   # fragment outlier: > sample_mean + this * sample_std
        "min_locus_hits": 3,     # require ≥ this many hits per locus to flag
    },
    "dmr": {
        "p_adj_thresh":   0.05,   # BH-corrected p-value cutoff
        "delta_beta_min": 0.10,   # minimum |ΔBeta| to call a DMR
        "chunk_size":     None,   # CpGs per pivot chunk; None = in-memory (set to e.g. 50_000 for EPIC)
        "covariate_cols": ["age", "sex"],  # OLS covariates; None or [] = Wilcoxon
    },
    "ml": {
        "n_top_cpgs": 200,    # restrict to top-variance CpGs before classification
        "l1_ratio":   0.5,    # ElasticNet mix (0 = Ridge, 1 = Lasso)
        "c_param":    1.0,    # inverse regularisation strength
        "chunk_size": None,   # CpGs per variance chunk; None = in-memory
    },
}


def load_config(path: str | None = None) -> dict[str, dict]:
    """
    Load analysis thresholds from a JSON config file.

    Parameters
    ----------
    path : str or None
        Path to a JSON config file.  If None, looks for ``config.json`` in the
        project root (two directory levels above this file).  If the file does
        not exist the function returns the built-in defaults silently.

    Returns
    -------
    dict
        Nested dict with sections ``"qc"``, ``"duplicates"``, ``"clonality"``,
        ``"dmr"``, ``"ml"``.  Each section contains threshold key–value pairs.
        Missing keys fall back to the defaults in ``_DEFAULTS``.
    """
    if path is None:
        path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "config.json",
        )

    # Start from a deep copy of defaults
    cfg: dict[str, dict] = {section: dict(values) for section, values in _DEFAULTS.items()}

    if os.path.isfile(path):
        with open(path) as f:
            user_cfg = json.load(f)
        for section, values in user_cfg.items():
            if section.startswith("_"):   # skip comment / metadata keys
                continue
            if section in cfg:
                cfg[section].update(values)
            else:
                cfg[section] = values

    return cfg
