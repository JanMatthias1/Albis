# sim-app tutorial

This folder contains a package-focused tutorial for `sim-app`.

- `sim_app_tutorial.ipynb` shows the intended app workflow: generate one
  dataset, inspect the returned `AnnData`, plot coordinates, and save the
  result as `.h5ad`.
- The tutorial is intentionally small and uses reduced simulation settings so
  it can run interactively during development.

From the repository root:

```bash
python -m pip install -e ".[plot]"
jupyter lab tutorial/sim_app_tutorial.ipynb
```

The paper, manuscript, and larger research examples live in the separate
`sim_paper` repository.
