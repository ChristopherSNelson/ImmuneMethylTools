# ImmuneMethylTools

**A DNA Methylation QC Pipeline for Autoimmune Disease Research.**

ImmuneMethylTools is a robust pipeline designed to identify and mitigate common wet-lab and bioinformatics artifacts in DNA methylation data, ensuring the integrity and reliability of downstream analyses. It provides comprehensive quality control, surgical data cleaning, and leverages advanced statistical and machine learning methods to isolate true biological signals.

## Key Features

*   **Artifact Detection & Mitigation:** Flags and handles issues like bisulfite failure, low coverage, sample contamination, duplicates, sex-metadata mixups, and clonal VDJ artifacts.
*   **Multi-Stage Workflow:** A structured pipeline (QC -> Normalization -> DMR -> ML) ensures data quality at every step.
*   **Configurable:** Easily adjust analysis thresholds via `config.json` without modifying code.
*   **Interactive Exploration:** Includes a Jupyter Notebook for an end-to-end interactive demonstration.

## Quick Start

1.  **Setup:** Run `bash setup.sh` and `python data/generate_mock_data.py`.
2.  **Run Pipeline:** Execute `python core/pipeline.py` to process data and generate outputs in the `output/` directory. Use `--report` for a full PDF report.
3.  **Explore:** Open the Jupyter Notebook: `jupyter notebook notebooks/ImmuneMethylTools_Validation.ipynb`.

---

**For a comprehensive guide, including detailed features, design choices, full setup, workflow, artifact map, and future work, please refer to [tldr_readme.md](tldr_readme.md).**
