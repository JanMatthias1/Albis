#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# create_sim_app_env.sh
#
# Creates a conda environment for the sim_app tutorial, installs sim_app
# with [tutorial,plot] extras, and registers a Jupyter kernel.
#
# Usage:
#   bash create_sim_app_env.sh
#
# Override defaults with env vars:
#   CONDA_ENV_NAME=my-env PYTHON_VERSION=3.11 bash create_sim_app_env.sh
# =============================================================================

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

CONDA_ENV_NAME="${CONDA_ENV_NAME:-sim-app-tutorial}"
PYTHON_VERSION="${PYTHON_VERSION:-3.10}"
KERNEL_NAME="${KERNEL_NAME:-sim-app-tutorial}"
KERNEL_DISPLAY_NAME="${KERNEL_DISPLAY_NAME:-Python (sim-app tutorial)}"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Run: module load conda (or your cluster's conda module)." >&2
    exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV_NAME}"; then
    echo "[setup] Environment already exists: ${CONDA_ENV_NAME} — skipping creation."
else
    echo "[setup] Creating conda environment: ${CONDA_ENV_NAME} (python=${PYTHON_VERSION})"
    conda create -y -n "${CONDA_ENV_NAME}" "python=${PYTHON_VERSION}" pip
fi

conda activate "${CONDA_ENV_NAME}"

echo "[setup] Upgrading pip/setuptools/wheel"
python -m pip install --upgrade pip setuptools wheel

echo "[setup] Installing sim_app with tutorial + plot extras"
python -m pip install -e "${REPO_ROOT}[tutorial,plot]"

echo "[setup] Ensuring jupyter notebook + ipykernel are present"
python -m pip install notebook ipykernel

echo "[setup] Registering Jupyter kernel: ${KERNEL_DISPLAY_NAME}"
python -m ipykernel install \
    --user \
    --name "${KERNEL_NAME}" \
    --display-name "${KERNEL_DISPLAY_NAME}"

echo "[setup] Verifying imports"
python -c "
import sim_app
import notebook
import ipykernel
print(f'sim_app {sim_app.__version__} OK')
print('notebook OK')
print('ipykernel OK')
"

echo ""
echo "[setup] Done. Environment: ${CONDA_ENV_NAME}"
echo ""
echo "Activate with:"
echo "  conda activate ${CONDA_ENV_NAME}"
echo ""
echo "Start notebook with:"
echo "  jupyter-notebook --no-browser --ip=0.0.0.0 --port=8888"