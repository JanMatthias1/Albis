# sim-app

`sim-app` is a  Python package for generating synthetic 3D
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

Parameters are grouped below in the order the simulator applies them: select
an output, build the tissue sphere, place cells within it, assign each cell to
a spatial domain with irregular boundaries, assign cell types and generate
genes and molecules, then slice, capture, and apply batch effects.

#### Output selection

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `output` | `"bin"` | Observation type to generate: `"cell"`, `"bin"`, or `"spot"`. |
| `slice_axis` | `"Z"` | Axis normal to the 2D slice plane: `"X"`, `"Y"`, or `"Z"`. |

#### 1. Tissue sphere

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `tissue_shape` | `"sphere"` | Overall tissue geometry; only `"sphere"` is currently supported. |
| `sphere_radius_um` | `300.0` | Radius of the tissue sphere, in microns. |

#### 2. Cell placement

Cell centroids and radii are sampled within the tissue sphere before domains
are assigned; domain membership is determined afterward, from each cell's
position.

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `n_cells` | `1000` | Number of simulated cells. |
| `allow_cell_overlap` | `False` | Whether overlapping cell spheres are permitted during placement. |
| `cell_radius_kwargs` | `None` | Optional low-level cell-radius distribution settings (`radius_dist`, `r_mean`, `r_sigma`, `r_min`, `r_max`). |

#### 3. Spatial domains

Domain membership is evaluated per cell, using each cell's position within
the sphere.

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `n_domains` | `4` | Number of spatial domains: `n_domains - 1` angular wedges surrounding one central core. |
| `domain_layout` | `"core_wedges"` | Domain-generation strategy; only `"core_wedges"` is currently supported. |
| `core_frac` | `0.35` | Core radius as a fraction of the sphere radius, not its volume. For example, `0.55` yields a core occupying roughly 17% of the sphere's volume. |

#### 4. Domain boundary irregularity

By default, core and wedge boundaries are smooth and straight. These
parameters introduce two independent forms of irregularity: a continuous
geometric warp of the boundary surfaces, and discrete relabeling ("fuzz") of
individual points near those boundaries.

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `core_bump_amp` | `0.12` | Fractional, direction-dependent perturbation of the core radius; `0.0` yields a perfect sphere. |
| `wedge_angle_amp_deg` | `12.0` | Angular perturbation applied to wedge boundaries, in degrees; `0.0` yields straight radial cuts. |
| `noise_terms` | `6` | Number of sinusoidal components summed to construct the smooth noise field underlying `core_bump_amp` and `wedge_angle_amp_deg`. Higher values yield smoother, less directional noise. |
| `noise_freq_range` | `(0.8, 2.2)` | Spatial-frequency range (cycles per micron) of that noise field. Higher frequencies yield finer-grained boundary texture. |
| `boundary_fuzz_width_deg` | `0.0` | Angular band, in degrees, around each wedge-wedge boundary within which points may be relabeled to the neighboring wedge. |
| `boundary_fuzz_flip_prob` | `0.0` | Probability that a point within `boundary_fuzz_width_deg` of a boundary is relabeled to the neighboring wedge. |
| `core_fuzz_width_um` | `0.0` | Radial band, in microns, around the core boundary within which points may be relabeled across the core/wedge interface. |
| `core_fuzz_flip_prob` | `0.0` | Probability that a point within `core_fuzz_width_um` of the core boundary is relabeled across the core/wedge interface. |

#### 5. Cell types, genes, and molecules

Once each cell has a domain, it is assigned a cell type, a gene panel is
built, and gene-expression counts are expanded into individual transcript
locations.

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `n_cell_types` | `4` | Number of cell types. |
| `domain_type_mix` | `None` | Optional `(n_domains, n_cell_types)` composition matrix specifying which cell types occur in each domain, and in what proportions. Each row is renormalized to a probability distribution and used to draw the cell type of every cell assigned to that domain. Defaults to a small built-in 4x4 example composition, or a uniform mix when `n_domains`/`n_cell_types` differ from 4/4. |
| `marker_genes_per_type` | `80` | Number of marker genes assigned to each cell type. |
| `noise_gene_frac` | `0.10` | Fraction of the gene panel carrying no cell-type signal. |
| `shared_marker_frac` | `0.25` | Fraction of each cell type's markers shared with other cell types, rather than unique to it. |
| `inside_prob` | `0.95` | Target fraction of a cell's transcripts generated within its cell radius. |
| `assign_k` | `8` | Number of nearest cells considered when reassigning transcripts to cells by containment, for `output="cell"`. |

#### 6. Slicing, capture, and batch effects

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `n_slices` | `5` | Number of slices generated along the selected axis. |
| `capture_window_um` | `(500, 500)` | Width and height, in microns, of the rectangular capture window applied to each 2D slice before binning or spot aggregation. Applies only to `output="bin"` or `"spot"`. |
| `bin_size_um` | `20.0` | Bin width, in microns, for `output="bin"`. |
| `spot_spacing_um` | `100.0` | Center-to-center spacing between spots, in microns, for `output="spot"`. |
| `spot_radius_um` | `27.5` | Capture radius of each spot, in microns, for `output="spot"`. |
| `batch_sigma` | `0.15` | Standard deviation of the per-slice, per-gene log-fold-change applied to simulate batch effects across slices. |
| `unaligned_coordinates` | `True` | Whether to additionally generate a randomly rotated and translated ("unaligned") copy of each slice's coordinates. |
| `max_deg` | `180.0` | Maximum absolute per-slice in-plane rotation, in degrees, applied when generating unaligned coordinates. |
| `max_shift` | `200.0` | Maximum absolute per-slice in-plane translation, in microns, applied when generating unaligned coordinates. |
| `include_truth` | `True` | Whether to retain pre-batch-effect counts and ground-truth annotations in the returned object. |

#### Reproducibility

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `seed` | `2025` | Random seed governing cell placement, gene expression, and batch effects. |
| `base_seed_unaligned` | `12345` | Random seed governing the per-slice rotation and translation used to generate unaligned coordinates. |

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

`generate_data()` always builds exactly one modality/axis pair, so generating
several combinations that way re-simulates the tissue sphere, cells, domains,
and gene panel from scratch for each call. If you need multiple
resolutions — several modalities, several slice axes, or both — from the same
underlying tissue, call the lower-level `sim_app.simulate_3d_molecule_sphere_multires(...)`
directly instead. It builds the shared simulation once and returns every
requested modality/axis combination from it.

```python
sim = sim_app.simulate_3d_molecule_sphere_multires(
    sphere_R_um=300.0,
    n_cells=1_000,
    n_domains=4,
    n_slices=5,
    capture_window_um=(300.0, 300.0),
    bin_size_um=30.0,
    output_modalities=("bin", "spot"),
    slice_axes=("X", "Z"),
    seed=2025,
)

bin_adata_x = sim["bin_adatas"]["X"]
spot_adata_z = sim["spot_adatas"]["Z"]

sim_app.save(bin_adata_x, "bin_x.h5ad")
sim_app.save(spot_adata_z, "spot_z.h5ad")
```

`output_modalities` and `slice_axes` each default to every supported value
(`{"cell", "bin", "spot"}` and `("X", "Y", "Z")`) when omitted, so a bare call
with no filters generates everything at once.

`adata_cell_sectioned`, `bin_adatas`, and `spot_adatas` are each `{axis:
AnnData}` dicts, since a separate object is built per requested slicing axis —
that's what `sim["bin_adatas"]["X"]` above is indexing into. The dictionary
also includes `meta` (simulation parameters and summary statistics) and two
single, pre-sectioning cell-level objects: `adata_cell_true` (ideal
gene counts from the expression model) and `adata_cell_obs` (the same cells'
counts after molecules are sampled in 3D and reassigned to nearby cells,
which lets some counts spill into neighboring cells or go unassigned).
Comparing the two quantifies that spillover noise:

```python
import numpy as np

true_counts = np.asarray(sim["adata_cell_true"].X.sum(axis=1)).ravel()
obs_counts = np.asarray(sim["adata_cell_obs"].X.sum(axis=1)).ravel()

correlation = np.corrcoef(true_counts, obs_counts)[0, 1]
print(f"true vs. observed per-cell total-count correlation: {correlation:.3f}")
# true vs. observed per-cell total-count correlation: 0.999
```

`simulate_3d_molecule_sphere_multires` accepts the same capture-window
parameter names as `generate_data()` (`capture_window_um`,
`capture_window_center_um`, `xenium_capture_window_um`), along with several
lower-level, expression-model parameters not exposed through `generate_data()`
— see the function definition in `sim_app/simulation_sphere.py` for the full
signature.

For everyday use, prefer `sim_app.generate_data(...)`: it wraps this function
and returns exactly one `AnnData` object for the requested modality and axis,
without building the others.
