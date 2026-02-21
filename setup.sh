#!/usr/bin/env bash
# setup.sh — Create virtual environment and install ImmuneMethylTools dependencies
set -euo pipefail

PYTHON=${PYTHON:-python3}
VENV_DIR="venv"

echo "==> ImmuneMethylTools Setup"
echo "    Python: $($PYTHON --version)"
echo "    Venv dir: $VENV_DIR"

# Create venv
if [ ! -d "$VENV_DIR" ]; then
    echo "==> Creating virtual environment..."
    $PYTHON -m venv "$VENV_DIR"
else
    echo "==> Virtual environment already exists, skipping creation."
fi

# Activate
source "$VENV_DIR/bin/activate"
echo "==> Activated venv: $(which python)"

# Upgrade pip
pip install --upgrade pip --quiet

# Install requirements
echo "==> Installing requirements..."
pip install -r requirements.txt

# Also install the package itself (picks up pyproject.toml dependencies)
echo "==> Installing package in editable mode..."
pip install -e . --quiet

echo ""
echo "Setup complete. Activate with:"
echo "  source $VENV_DIR/bin/activate"
echo ""
echo "Run mock data generation:"
echo "  python data/generate_mock_data.py"
echo ""
echo "Run tests:"
echo "  python -m pytest tests/ -v"
echo ""
echo "Run the full pipeline (with PDF report):"
echo "  python core/pipeline.py --report"
echo ""
echo "Open the demo notebook:"
echo "  jupyter notebook notebooks/ImmuneMethylTools_Validation.ipynb"
