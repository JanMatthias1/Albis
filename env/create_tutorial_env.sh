#!/usr/bin/env bash
#SBATCH --job-name=sim-app-env
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=sim-app-env-%j.log

set -euo pipefail

# --- Resolve paths -----------------------------------------------------------
# $SLURM_SUBMIT_DIR is where you ran sbatch from. If you submit from the repo
# root or env/ directory this resolves the project reliably. Override with
# SIM_APP_REPO_ROOT if your cluster copies the script to a spool directory.
SCRIPT_PATH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-${SCRIPT_PATH_DIR}}"

if [[ -f "${SUBMIT_DIR}/pyproject.toml" && -d "${SUBMIT_DIR}/sim_app" ]]; then
  REPO_ROOT="${SUBMIT_DIR}"
elif [[ -f "${SUBMIT_DIR}/../pyproject.toml" && -d "${SUBMIT_DIR}/../sim_app" ]]; then
  REPO_ROOT="$(cd "${SUBMIT_DIR}/.." && pwd)"
else
  REPO_ROOT="${SIM_APP_REPO_ROOT:-}"
fi

if [[ -z "${REPO_ROOT}" || ! -f "${REPO_ROOT}/pyproject.toml" ]]; then
  echo "Could not find sim_app repository root." >&2
  echo "Run from the repo root/env directory, or set SIM_APP_REPO_ROOT=/path/to/sim_app." >&2
  exit 1
fi

# --- Load conda --------------------------------------------------------------
# Compute nodes do not inherit your login-node module environment. Load conda
# here. Adjust the module name to whatever `module avail` shows on your cluster.
module load conda 2>/dev/null || true

if ! command -v conda >/dev/null 2>&1; then
  echo "conda is not on PATH." >&2
  echo "Load your cluster conda module first, or edit this script's module load line." >&2
  exit 1
fi

# --- Choose env name ---------------------------------------------------------
CONDA_ENV_NAME="${SIM_APP_CONDA_ENV:-sim-app-tutorial}"
PYTHON_VERSION="${SIM_APP_PYTHON_VERSION:-3.10}"
KERNEL_NAME="${SIM_APP_KERNEL_NAME:-sim-app-tutorial}"
KERNEL_DISPLAY_NAME="${SIM_APP_KERNEL_DISPLAY_NAME:-Python (sim-app tutorial)}"

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV_NAME}"; then
  echo "Conda environment already exists: ${CONDA_ENV_NAME}"
else
  echo "Creating conda environment: ${CONDA_ENV_NAME}"
  conda create -y -n "${CONDA_ENV_NAME}" "python=${PYTHON_VERSION}" pip
fi

conda activate "${CONDA_ENV_NAME}"

echo "Upgrading packaging tools"
python -m pip install --upgrade pip setuptools wheel

echo "Installing sim-app with notebook and plotting support"
python -m pip install -e "${REPO_ROOT}[tutorial,plot]"

echo "Verifying tutorial environment"
python - <<PY
import importlib.util
from pathlib import Path

required = ["anndata", "numpy", "scipy", "matplotlib", "notebook", "nbformat", "ipykernel"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    raise SystemExit("Missing Python package(s): " + ", ".join(missing))

import nbformat
import sim_app

notebook_path = Path("${REPO_ROOT}") / "tutorial" / "sim_app_tutorial.ipynb"
nbformat.read(notebook_path, as_version=4)

print(f"sim_app import OK: {sim_app.__version__}")
print(f"notebook file OK: {notebook_path}")
PY

if ! command -v jupyter-notebook >/dev/null 2>&1; then
  echo "jupyter-notebook is not available after installation." >&2
  exit 1
fi

echo "Jupyter command OK: $(command -v jupyter-notebook)"

echo "Registering Jupyter kernel: ${KERNEL_DISPLAY_NAME}"
python -m ipykernel install \
  --user \
  --name "${KERNEL_NAME}" \
  --display-name "${KERNEL_DISPLAY_NAME}"

HOSTNAME_VALUE="$(hostname -f 2>/dev/null || hostname)"

cat <<EOF

Environment ready.

Activate the tutorial environment with:
  conda activate "${CONDA_ENV_NAME}"

Start Jupyter Notebook with:
  jupyter-notebook --no-browser --ip=0.0.0.0 --port 8888

Then open:
  ${REPO_ROOT}/tutorial/sim_app_tutorial.ipynb

In Jupyter, select the kernel:
  ${KERNEL_DISPLAY_NAME}

If your local machine needs an SSH tunnel, run this from your local machine:
  ssh -N -L 8888:${HOSTNAME_VALUE}:8888 <user>@<cluster-login-host>

EOF
