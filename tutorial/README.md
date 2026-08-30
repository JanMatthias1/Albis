# albis tutorial

This folder contains a package-focused tutorial for `albis`.

- `higher_level_api_tutorial.ipynb` shows the app-facing `generate_data()`
  workflow: generate one dataset, inspect the returned `AnnData`, plot
  coordinates, and save the result as `.h5ad`.
- `lower_level_api_tutorial.ipynb` shows the low-level simulator API
  (`simulate_3d_molecule_sphere_multires`/`_base`/`section_3d_molecule_sphere`),
  which exposes every simulation parameter as one plain dictionary you can
  pass with `**config` — including several expression-model parameters
  `generate_data()` doesn't expose at all. See
  [`../LOW_LEVEL_SIMULATOR.md`](../LOW_LEVEL_SIMULATOR.md) for the full
  written reference this notebook walks through.
- Both notebooks are intentionally small and use reduced simulation settings
  so they can run interactively during development.

Install `albis` with the notebook + plotting extras (from a clone of this
repo, or from PyPI):

```bash
pip install -e ".[tutorial,plot]"      # from a clone
# or: pip install "albis[tutorial,plot]"
```

Then launch Jupyter:

```bash
conda activate albis
jupyter-notebook --no-browser --ip=0.0.0.0 --port 8888
```

The default port is `8888`. Override it when needed:

```bash
jupyter-notebook --no-browser --ip=0.0.0.0 --port 8890
```

The paper, manuscript, and larger research examples live in the separate
`sim_paper` repository.
