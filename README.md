# sim-app

`sim-app` is a small Python package for generating synthetic 3D
spatial-transcriptomics data as app-ready
[`AnnData`](https://anndata.readthedocs.io/) objects.
It is structured as an installable package with a small tutorial workflow.

The simulator creates cells in a sphere, generates negative-binomial
gene-expression counts, expands those counts into individual transcript
locations, and aggregates transcripts to cell, bin, or spot observations.

```text
sim_app/
├── pyproject.toml
├── README.md
├── sim_app/
│   ├── api.py
│   ├── plotting.py
│   └── simulation_sphere.py
├── tests/
└── tutorial/
```

The public API is intentionally small:

```python
import sim_app

adata = sim_app.generate_data(...)
summary = sim_app.describe(adata)
figure = sim_app.plot(adata, ...)
path = sim_app.save(adata, "simulation.h5ad")
```

## Installation

For local development, install the package from the repository root:

```bash
python -m pip install -e .
```

Install optional static plotting support as well:

```bash
python -m pip install -e ".[plot]"
```

The PyPI package name is `sim-app`; the Python import name is `sim_app`.

## Development

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest
```

This repository now contains only the Python package. Manuscript, paper,
notebook, and tutorial material lives in the separate `sim_paper` repository.

## Tutorial

A package-focused app tutorial is available at
[`tutorial/sim_app_tutorial.ipynb`](tutorial/sim_app_tutorial.ipynb). It shows
how to generate a small dataset, inspect the returned `AnnData`, plot aligned
and unaligned coordinates, and save the result.

To run it from this environment:

```bash
bash env/create_tutorial_env.sh
conda activate sim-app-tutorial
jupyter-notebook --no-browser --ip=0.0.0.0 --port 8888
```

## Quick start

```python
import sim_app

adata = sim_app.generate_data(
    output="bin",
    slice_axis="Z",
    n_cells=1_000,
    n_slices=5,
    n_domains=4,
    seed=2025,
)

summary = sim_app.describe(adata)

figure = sim_app.plot(
    adata,
    view="2d",
    coordinates="aligned",
    color="domain_true",
)

sim_app.save(adata, "simulation_bins_z.h5ad")
```

## `generate_data`

```python
adata = sim_app.generate_data(output="bin", slice_axis="Z")
```

Parameters are usually passed as keyword arguments. They can also be omitted to
use `sim_app.DEFAULT_PARAMETERS`:

```python
adata = sim_app.generate_data(output="spot", slice_axis="X", n_cells=2_000)
```

If you build a configuration programmatically, passing a dictionary is still
supported:

```python
params = {"output": "spot", "slice_axis": "X", "n_cells": 2_000}
adata = sim_app.generate_data(params)
```

The function returns exactly **one** `AnnData` object. It does not generate
unrequested platforms or slice axes. For example, requesting Z-axis bins skips
cell-section outputs, spots, and X/Y bin outputs.

For a small tutorial or smoke-test dataset, use:

```python
adata = sim_app.example_data()
```

`example_data()` calls `generate_data()` with smaller defaults. Keyword
arguments can override those defaults:

```python
adata = sim_app.example_data(output="spot", slice_axis="Y", n_cells=500)
```

### Parameters

Parameters are grouped below in the order the simulator applies them: pick an
output, build the tissue sphere, carve it into domains, roughen the domain
boundaries, place cells and genes inside it, then cut, capture, and add batch
effects across slices.

#### Output selection

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `output` | `"bin"` | Observation type: `"cell"`, `"bin"`, or `"spot"`. |
| `slice_axis` | `"Z"` | Axis normal to the 2D slices: `"X"`, `"Y"`, or `"Z"`. |

#### 1. Tissue sphere

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `tissue_shape` | `"sphere"` | Tissue geometry. Only `"sphere"` is currently supported. |
| `sphere_radius_um` | `300.0` | Radius of the tissue sphere, in microns. |

#### 2. Spatial domains

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `n_domains` | `4` | Number of spatial domains: `n_domains - 1` angular wedges plus one central core. |
| `domain_layout` | `"core_wedges"` | Domain-generation strategy. Only `"core_wedges"` is currently supported. |
| `core_frac` | `0.35` | Core radius as a fraction of the sphere *radius* (not volume): a value of `0.55`, for example, gives a core spanning ~17% of the sphere's volume. |

#### 3. Domain boundary irregularity

By default the core and wedge boundaries are perfectly smooth and straight.
These parameters add two independent kinds of irregularity: a continuous
geometric warp of the boundary surfaces, and discrete relabeling ("fuzz") of
individual points near those boundaries.

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `core_bump_amp` | `0.12` | Fractional, direction-dependent bump applied to the core radius (`0.0` = perfect sphere). |
| `wedge_angle_amp_deg` | `12.0` | Angular wobble applied to wedge boundaries, in degrees (`0.0` = straight radial cuts). |
| `noise_terms` | `6` | Number of summed sinusoidal components in the smooth noise field driving `core_bump_amp` and `wedge_angle_amp_deg`. More terms produce smoother, less directional noise. |
| `noise_freq_range` | `(0.8, 2.2)` | Spatial-frequency range (cycles per micron) of that noise field. Higher frequencies produce finer-grained boundary texture. |
| `boundary_fuzz_width_deg` | `0.0` | Angular band, in degrees, around each wedge-wedge boundary in which points may be relabeled to the neighboring wedge. |
| `boundary_fuzz_flip_prob` | `0.0` | Probability that a point inside `boundary_fuzz_width_deg` of a boundary is relabeled to the neighboring wedge. |
| `core_fuzz_width_um` | `0.0` | Radial band, in microns, around the (already-warped) core boundary in which points may be relabeled across the core/wedge interface. |
| `core_fuzz_flip_prob` | `0.0` | Probability that a point inside `core_fuzz_width_um` of the core boundary is relabeled across the core/wedge interface. |

#### 4. Cells

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `n_cells` | `1000` | Number of simulated cells. |
| `allow_cell_overlap` | `False` | Permit overlapping cell spheres during placement. |
| `cell_radius_kwargs` | `None` | Optional low-level cell-radius distribution settings (`radius_dist`, `r_mean`, `r_sigma`, `r_min`, `r_max`). |
| `n_cell_types` | `4` | Number of cell types. |
| `domain_type_mix` | `None` | Optional `(n_domains, n_cell_types)` matrix; row `d` gives the cell-type probability distribution used to assign cells in domain `d`. |
| `marker_genes_per_type` | `80` | Number of marker genes assigned to each cell type. |
| `noise_gene_frac` | `0.10` | Fraction of the gene panel carrying no cell-type signal. |
| `shared_marker_frac` | `0.25` | Fraction of each cell type's markers that are shared with other cell types, rather than unique to it. |
| `inside_prob` | `0.95` | Target fraction of a cell's source molecules generated within its cell radius. |
| `assign_k` | `8` | Number of nearby cells considered when reassigning molecules to cells by containment for `output="cell"`. |

#### 5. Slicing, capture, and batch effects

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `n_slices` | `5` | Number of slices along the selected axis. |
| `capture_window_um` | `(500, 500)` | Bin/spot capture-window width and height, in microns. Crops each 2D slice before binning/spot-aggregation; applies only to `output="bin"` or `"spot"`. |
| `xenium_capture_window_um` | `(12000, 24000)` | Reserved for a future cell-level (Xenium-style) capture-window crop. Currently recorded in `adata.uns["captures"]` metadata only and **not yet applied** to filter or crop cells. |
| `bin_size_um` | `20.0` | Bin width for `output="bin"`. |
| `spot_spacing_um` | `100.0` | Spot-center spacing, in microns, for `output="spot"`. |
| `spot_radius_um` | `27.5` | Spot capture radius, in microns, for `output="spot"`. |
| `batch_sigma` | `0.15` | Standard deviation of the per-slice, per-gene log-fold-change applied as a simulated batch effect. |
| `unaligned_coordinates` | `True` | Also generate a randomly rotated/translated ("unaligned") copy of each slice's coordinates. |
| `max_deg` | `180.0` | Largest absolute per-slice in-plane rotation, in degrees, applied to unaligned coordinates. |
| `max_shift` | `200.0` | Largest absolute per-slice in-plane translation, in microns, applied to unaligned coordinates. |
| `include_truth` | `True` | Retain pre-batch-effect counts and ground-truth annotations in the returned object. |

#### Reproducibility

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `seed` | `2025` | Random seed for reproducibility. |
| `base_seed_unaligned` | `12345` | Base random seed for the per-slice unaligned-coordinate perturbation. |

Use `sim_app.DEFAULT_PARAMETERS` to inspect the complete supported parameter
set. Unsupported values fail early with a descriptive error.

## AnnData contract

All coordinates use microns. `adata.X` is a sparse count matrix containing the
final observed counts after slice-specific batch effects.

| Location | Contents |
| --- | --- |
| `adata.X` | Final observation-by-gene count matrix. |
| `adata.layers["counts_pre_batch"]` | Counts before batch effects; present when `include_truth=True`. |
| `adata.obs` | Slice IDs, capture/grid metadata, and available truth labels. |
| `adata.var` | Marker/noise-gene annotations. |
| `adata.obsm["spatial"]` | Canonical aligned 2D coordinates in the selected slice plane. |
| `adata.obsm["spatial_unaligned"]` | Per-slice rigidly transformed 2D coordinates. |
| `adata.obsm["spatial_3d"]` | Canonical 3D coordinates. |
| `adata.obsm["spatial_3d_unaligned"]` | 3D representation of the in-plane unaligned transforms. |
| `adata.obsm["cell_type_frac_true"]` | Source cell-type fractions for bin/spot observations. |
| `adata.obsm["domain_frac_true"]` | Source-domain fractions for bin/spot observations. |
| `adata.uns["batch_effect_factors"]` | Per-slice gene-wise multiplicative batch factors. |
| `adata.uns["sim_params"]` | Resolved simulation configuration and simulator metadata. |
| `adata.uns["output"]` | Selected platform and slice axis. |

`spatial_unaligned` is not shuffled data. Each slice receives one deterministic
random rigid transform: a rotation plus a translation. No scaling, shearing,
or molecule resimulation occurs. `spatial_3d_unaligned` applies that transform
to the corresponding plane while retaining the coordinate normal to the slice.

### Truth annotations

For `output="bin"` and `output="spot"`, composition fractions are computed
from each molecule's source cell type/domain. Molecules outside every cell can
still contribute to bin and spot counts. Cell output uses containment-based
assignment, so molecules can spill into another cell or be unassigned.

Set `include_truth=False` to remove `counts_pre_batch`, composition fractions,
and direct cell/domain truth labels from the returned object.

## Describe and Save

```python
summary = sim_app.describe(adata)
```

`describe()` returns a plain dictionary with the most important contents of the
object: observation and gene counts, total counts, output platform, slice axis,
available coordinate keys, layers, metadata keys, and truth annotations.

```python
path = sim_app.save(adata, "simulation_bins_z.h5ad")
```

`save()` writes the object with AnnData's `.h5ad` format and returns the path.
It is a convenience wrapper around `adata.write_h5ad(...)`.

## Plotting

```python
# Aligned 2D tissue view, colored by an observation annotation
sim_app.plot(adata, view="2d", coordinates="aligned", color="domain_true")

# Unaligned 3D view, colored by slice
sim_app.plot(adata, view="3d", coordinates="unaligned", color="slice_id")

# Expression of one gene in a single slice
sim_app.plot(adata, view="2d", color="G1", slice_id=0, point_size=8)
```

`sim_app.plot()` returns a Matplotlib `Figure`. It supports:

- `view="2d"` or `view="3d"`
- `coordinates="aligned"` or `coordinates="unaligned"`
- `color` set to an `obs` column or gene name
- `slice_id` to show one slice
- `point_size`, `alpha`, `max_points`, and `seed` rendering controls

Large datasets are deterministically downsampled to `max_points=50_000` by
default. Use `max_points=None` only when the output size is known to be safe.

## Performance guidance

Runtime and memory scale primarily with the number of generated molecules,
which in turn grows with `n_cells` and the simulated expression level.

- Start with `n_cells=1_000`, then increase gradually.
- Use a smaller `capture_window_um` or larger `bin_size_um` to reduce the
  number of bin observations.
- Request one output and one axis at a time through `generate_data()`.
- Use `max_points` in plotting rather than trying to render every bin/cell.

The default parameters are chosen for interactive experimentation, not for
reproducing the original large-scale standalone example.

## Low-level simulator

`sim_app.simulate_3d_molecule_sphere_multires(...)` remains available for
advanced workflows that need multiple modalities or axes in one call. It
returns a dictionary containing cell, bin, and spot `AnnData` objects.

The app-facing API should normally use `sim_app.generate_data(...)`, because
it builds only the selected modality and axis.

## Current scope

- Sphere tissue geometry only.
- Core-plus-wedge domain layout only.
- Static Matplotlib plotting only.
- File output is `.h5ad` via `sim_app.save(...)` or `adata.write_h5ad(...)`.

Future extensions can add ellipsoid/box geometries, additional domain layouts,
and interactive plotting without changing the high-level `generate_data()`
contract.
