# sim-app

Python package for the 3D transcript-level spatial transcriptomics simulator in this repository.

This packaging step keeps the original script in `../code/3d_simulation_sphere_new.py` unchanged and provides an importable package module:

```python
from sim_app import simulate_3d_molecule_sphere_multires
```

## Install for local development

From this directory:

```bash
python -m pip install -e .
```

## Minimal usage

```python
from sim_app import simulate_3d_molecule_sphere_multires

sim = simulate_3d_molecule_sphere_multires(
    n_cells=100,
    n_domains=4,
    n_cell_types=4,
    seed=2025,
)

adata_cell_obs = sim["adata_cell_obs"]
```

For helper functions, import from the module directly:

```python
from sim_app.simulation_sphere import sample_uniform_in_sphere
```

## Notes

- Package name for pip: `sim-app`
- Import package name in Python: `sim_app`
- Main module: `sim_app.simulation_sphere`
