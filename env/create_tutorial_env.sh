#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PYTHON_BIN="${PYTHON:-python3}"
ENV_DIR="${SIM_APP_ENV_DIR:-${SCRIPT_DIR}/.venv}"
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
