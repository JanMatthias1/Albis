# sim-app tutorial

This folder contains a package-focused tutorial for `sim-app`.

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

Create the tutorial environment from the repository root:

```bash
bash env/create_tutorial_env.sh
```

On a SLURM cluster, submit the same setup script as a job:

```bash
sbatch env/create_tutorial_env.sh
```

The script creates a conda environment, installs `sim-app` with notebook and
plotting dependencies, verifies the tutorial notebook can be loaded, and
registers a Jupyter kernel named `sim-app-tutorial`.

The default conda environment name is `sim-app-tutorial`. Override it with
`SIM_APP_CONDA_ENV` if needed:

```bash
SIM_APP_CONDA_ENV=my-env bash env/create_tutorial_env.sh
```

After the environment is ready, launch Jupyter from the cluster like this:

```bash
conda activate sim-app-tutorial
jupyter-notebook --no-browser --ip=0.0.0.0 --port 8888
```

The default port is `8888`. Override it when needed:

```bash
jupyter-notebook --no-browser --ip=0.0.0.0 --port 8890
```

The paper, manuscript, and larger research examples live in the separate
`sim_paper` repository.
