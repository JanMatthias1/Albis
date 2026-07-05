#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

PORT="${PORT:-8888}"
IP="${IP:-0.0.0.0}"
NOTEBOOK_DIR="${NOTEBOOK_DIR:-${REPO_ROOT}/tutorial}"
VENV_PATH="${VENV_PATH:-${REPO_ROOT}/.venv}"

cd "${REPO_ROOT}"

if [[ -n "${CONDA_ENV:-}" ]]; then
    if ! command -v conda >/dev/null 2>&1; then
        echo "CONDA_ENV is set to '${CONDA_ENV}', but conda is not on PATH." >&2
        exit 1
    fi
    CONDA_BASE="$(conda info --base)"
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV}"
elif [[ -d "${VENV_PATH}" ]]; then
    # shellcheck disable=SC1091
    source "${VENV_PATH}/bin/activate"
fi

python - <<'PY'
import importlib.util
import sys

required = ["anndata", "numpy", "scipy", "matplotlib", "notebook"]
missing = [name for name in required if importlib.util.find_spec(name) is None]

if missing:
    print("Missing Python package(s): " + ", ".join(missing), file=sys.stderr)
    print('Install them in this environment with: python -m pip install -e ".[tutorial,plot]"', file=sys.stderr)
    raise SystemExit(1)

import sim_app
print(f"Using sim_app {sim_app.__version__} from {sim_app.__file__}")
PY

if ! command -v jupyter-notebook >/dev/null 2>&1; then
    echo "jupyter-notebook is not on PATH after environment activation." >&2
    echo 'Install it with: python -m pip install -e ".[tutorial,plot]"' >&2
    exit 1
fi

HOSTNAME_VALUE="$(hostname -f 2>/dev/null || hostname)"

echo "Starting Jupyter Notebook"
echo "Host: ${HOSTNAME_VALUE}"
echo "Port: ${PORT}"
echo "Notebook directory: ${NOTEBOOK_DIR}"
echo
echo "If you need an SSH tunnel, run this from your local machine:"
echo "ssh -N -L ${PORT}:${HOSTNAME_VALUE}:${PORT} <user>@<cluster-login-host>"
echo

exec jupyter-notebook \
    --no-browser \
    --ip="${IP}" \
    --port="${PORT}" \
    --notebook-dir="${NOTEBOOK_DIR}"
