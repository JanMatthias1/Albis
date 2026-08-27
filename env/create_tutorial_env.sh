#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# create_tutorial_env.sh
#
# Creates a conda environment for the albis tutorial, installs albis
# with [tutorial,plot] extras, and registers a Jupyter kernel.
#
# Usage:
#   bash create_tutorial_env.sh
#
# Override defaults with env vars:
#   CONDA_ENV_NAME=my-env PYTHON_VERSION=3.11 bash create_tutorial_env.sh
# =============================================================================

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

CONDA_ENV_NAME="${CONDA_ENV_NAME:-sim-app-tutorial}"
PYTHON_VERSION="${PYTHON_VERSION:-3.10}"
KERNEL_NAME="${KERNEL_NAME:-sim-app-tutorial}"
KERNEL_DISPLAY_NAME="${KERNEL_DISPLAY_NAME:-Python (sim-app tutorial)}"
ENV_PREFIX="${REPO_ROOT}/env/${CONDA_ENV_NAME}"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found. Run: module load conda (or your cluster's conda module)." >&2
    exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ -d "${ENV_PREFIX}" ]]; then
    echo "[setup] Environment already exists at ${ENV_PREFIX} — skipping creation."
else
    echo "[setup] Creating conda environment at ${ENV_PREFIX} (python=${PYTHON_VERSION})"
    conda create -y --prefix "${ENV_PREFIX}" "python=${PYTHON_VERSION}" pip
fi

conda activate "${ENV_PREFIX}"

echo "[setup] Upgrading pip/setuptools/wheel"
python -m pip install --upgrade pip setuptools wheel

echo "[setup] Installing albis with tutorial + plot extras"
python -m pip install -e "${REPO_ROOT}[tutorial,plot]"

echo "[setup] Ensuring jupyter notebook + ipykernel are present"
python -m pip install notebook ipykernel

echo "[setup] Pinning jupyter_server to 2.18.2 (2.19+ omits the hostname URL needed for JHPCE portal access)"
python -m pip install "jupyter_server==2.18.2"

echo "[setup] Registering Jupyter kernel: ${KERNEL_DISPLAY_NAME}"
python -m ipykernel install \
    --user \
    --name "${KERNEL_NAME}" \
    --display-name "${KERNEL_DISPLAY_NAME}"

echo "[setup] Verifying imports"
python -c "
import albis
import notebook
import ipykernel
import jupyter_server
print(f'albis {albis.__version__} OK')
print('notebook OK')
print('ipykernel OK')
print(f'jupyter_server {jupyter_server.__version__} OK')
"

echo ""
echo "[setup] Done. Environment created at: ${ENV_PREFIX}"
echo ""
echo "Activate with:"
echo "  conda activate ${ENV_PREFIX}"
echo ""
echo "Start notebook with:"
echo "  jupyter-notebook --no-browser --ip=0.0.0.0 --port=8888"