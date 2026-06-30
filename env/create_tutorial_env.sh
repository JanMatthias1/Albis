#!/usr/bin/env bash
#SBATCH --job-name=sim-app-env
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=sim-app-env-%j.log

set -euo pipefail

# --- Resolve paths -----------------------------------------------------------
# $SLURM_SUBMIT_DIR is where you ran sbatch from. BASH_SOURCE is unreliable
# under SLURM because the script is copied to a spool dir before execution,
# so derive locations from the submit dir instead.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --- Load toolchain ----------------------------------------------------------
# Compute nodes do not inherit your login-node module environment. Load an
# appropriate Python here. Adjust the module name to whatever `module avail`
# shows on JHPCE.
module load conda 2>/dev/null || true
# module load python/3.11   # alternative if you want a bare python instead

PYTHON_BIN="${PYTHON:-python3}"

# --- Choose env location -----------------------------------------------------
# Default to scratch (home has quotas; scratch is the right place for a venv).
# Override by exporting SIM_APP_ENV_DIR before submitting.
ENV_DIR="${SIM_APP_ENV_DIR:-${SCRATCH:?SCRATCH is not set}/sim-app-tutorial-venv}"
KERNEL_NAME="${SIM_APP_KERNEL_NAME:-sim-app-tutorial}"
KERNEL_DISPLAY_NAME="${SIM_APP_KERNEL_DISPLAY_NAME:-Python (sim-app tutorial)}"

echo "Creating virtual environment at: ${ENV_DIR}"
"${PYTHON_BIN}" -m venv "${ENV_DIR}"

# shellcheck source=/dev/null
source "${ENV_DIR}/bin/activate"

echo "Upgrading packaging tools"
python -m pip install --upgrade pip setuptools wheel

echo "Installing sim-app with plotting support"
python -m pip install -e "${REPO_ROOT}/sim_app_package[plot]"

echo "Installing Jupyter tools"
python -m pip install jupyterlab notebook ipykernel

echo "Registering Jupyter kernel: ${KERNEL_DISPLAY_NAME}"
python -m ipykernel install \
  --user \
  --name "${KERNEL_NAME}" \
  --display-name "${KERNEL_DISPLAY_NAME}"

cat <<EOF

Environment ready.

Activate it with:
  source "${ENV_DIR}/bin/activate"

Open the tutorial with:
  jupyter lab "${REPO_ROOT}/tutorial/sim_app_tutorial.ipynb"

In Jupyter, select the kernel:
  ${KERNEL_DISPLAY_NAME}

EOF