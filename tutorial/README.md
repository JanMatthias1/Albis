# albis tutorial

This folder contains six tutorials for `albis`.

- `higher_level_api_tutorial.ipynb` shows the app-facing `generate_data()`
  workflow: generate one dataset, inspect the returned `AnnData`, plot
  coordinates, and save the result as `.h5ad`.
- `lower_level_api_tutorial.ipynb` shows the low-level simulator API
  (`simulate_3d_molecule_sphere_multires`/`_base`/`section_3d_molecule_sphere`),
  which exposes every simulation parameter as one plain dictionary you can
  pass with `**config` — including several expression-model parameters
  `generate_data()` doesn't expose at all. See
  [`../LOW_LEVEL_SIMULATOR.md`](../LOW_LEVEL_SIMULATOR.md) for the written
  reference to the one-shot `simulate_3d_molecule_sphere_multires` call this
  notebook walks through, and
  [`../LOW_LEVEL_SIMULATOR_STEPS.md`](../LOW_LEVEL_SIMULATOR_STEPS.md) for the
  `_base`/`section_3d_molecule_sphere` two-step pattern it also demonstrates.
- Both notebooks are intentionally small and use reduced simulation settings
  so they can run interactively during development.

- `xenium_cell_tutorial.ipynb` demonstrates cell-resolution simulation.
- `visiumHD_bin_tutorial.ipynb` demonstrates bin-resolution simulation.
- `visium_spot_tutorial.ipynb` demonstrates spot-resolution simulation.
- `strong_domain_mix_tutorial.ipynb` demonstrates stronger domain-specific
  cell-type mixtures across the three resolutions.

## Installation and launch

Install Albis, including plotting support:

```bash
pip install albis
```

Download the notebook you want to run. Open it in Jupyter using the same
Python environment; the notebooks import the installed package directly
and can run outside the repository.

If Jupyter is not already installed:

```bash
pip install notebook
jupyter notebook
```

Inside an existing notebook, use `%pip install albis` to install into the
active kernel, then restart the kernel if needed. Output files are saved
under `outputs/` relative to the notebook's working directory.

The paper, manuscript, and larger research examples live in the separate
`sim_paper` repository.
