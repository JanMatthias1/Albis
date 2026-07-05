# sim-app tutorial

This folder contains a package-focused tutorial for `sim-app`.

- `sim_app_tutorial.ipynb` shows the intended app workflow: generate one
  dataset, inspect the returned `AnnData`, plot coordinates, and save the
  result as `.h5ad`.
- The tutorial is intentionally small and uses reduced simulation settings so
  it can run interactively during development.

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
