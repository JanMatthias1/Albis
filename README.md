# Albis: simulation of multi-resolution and multi-dimensional spatial transcriptomics data

[![PyPI](https://img.shields.io/pypi/v/albis.svg)](https://pypi.org/project/albis/)
[![Python](https://img.shields.io/pypi/pyversions/albis.svg)](https://pypi.org/project/albis/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**Albis** is a Python package for simulating 3D spatial-transcriptomics data
with known ground truth. It builds a spherical tissue with irregular spatial
domains, places cells, draws negative-binomial gene-expression counts,
expands them into individual transcript locations, then slices the tissue and
aggregates transcripts into **cells**, **bins**, or **spots**. Every output is
an analysis-ready [`AnnData`](https://anndata.readthedocs.io/) object carrying
2D and 3D coordinates, domain and cell-type labels, and per-slice batch
effects, so it can be used to benchmark clustering, deconvolution, and slice
alignment methods.

## Installation

Albis is on PyPI: **<https://pypi.org/project/albis/>**

```bash
pip install albis
```

```python
import albis as ab
```

Requires Python 3.10 or newer. Plotting support is included.

## Two levels of the simulator

Albis has two entry points. Both run the same model; they differ in how much
of it you control.

| | High level | Low level |
| --- | --- | --- |
| Function | `ab.generate_data(...)` | `ab.simulate_3d_molecule_sphere_multires(...)` |
| Use it for | Quick, small datasets; learning the package | Custom or production-scale simulations |
| Parameters | A curated subset with tutorial-scale defaults | Every model parameter |
| Output per call | One `AnnData` (one modality, one slice axis) | Several `AnnData` objects from the same tissue (any mix of cell/bin/spot and X/Y/Z) |
| Documentation | [HIGH_LEVEL_SIMULATOR.md](HIGH_LEVEL_SIMULATOR.md) | [LOW_LEVEL_SIMULATOR.md](LOW_LEVEL_SIMULATOR.md) |
| Tutorial | [`higher_level_api_tutorial.ipynb`](tutorial/higher_level_api_tutorial.ipynb) | [`lower_level_api_tutorial.ipynb`](tutorial/lower_level_api_tutorial.ipynb) |

### 1. High level: `generate_data()`

Start here. Choose an output type and slice axis, and optionally change a few
parameters:

```python
import albis as ab

adata = ab.generate_data(
    output="bin",        # "cell", "bin", or "spot"
    slice_axis="Z",      # "X", "Y", or "Z"
    n_cells=1_000,
    n_slices=5,
    n_domains=4,
    seed=2025,
)

print(ab.describe(adata))
ab.plot(adata, view="2d", coordinates="aligned", color="domain_true")
ab.save(adata, "simulation_bins_z.h5ad")
```

`ab.example_data()` gives a smaller dataset for a first try.

**→ Full guide:** [HIGH_LEVEL_SIMULATOR.md](HIGH_LEVEL_SIMULATOR.md) covers
all parameters, the `AnnData` output layout, empty bins/spots, plotting, and
performance.

### 2. Low level: `simulate_3d_molecule_sphere_multires()`

Use this when you need parameters `generate_data()` does not expose, such as
the NB dispersion (`theta`, `theta_jitter`), baseline expression
(`base_gene_lognormal`), marker fold-changes, noise-gene scale, or per-domain
size factors. It also builds the 3D tissue once and returns several
modalities and slice axes from that same tissue:

```python
sim = ab.simulate_3d_molecule_sphere_multires(
    sphere_R_um=300.0,
    n_cells=1_000,
    n_slices=5,
    capture_window_um=(300.0, 300.0),
    output_modalities=("cell", "bin", "spot"),
    slice_axes=("Z",),
    theta=2.0,
    seed=2025,
)

cell_z = sim["adata_cell_sectioned"]["Z"]
bin_z  = sim["bin_adatas"]["Z"]
spot_z = sim["spot_adatas"]["Z"]
```

**→ Full guide:** [LOW_LEVEL_SIMULATOR.md](LOW_LEVEL_SIMULATOR.md) has the
parameter reference, grouped to match the model description.

#### Platform-like configurations: Visium HD, Visium, and Xenium

We tuned low-level configurations so that simulated counts match real
platform data. These are the configurations used in the manuscript, and each
has a tutorial with the full parameter set:

| Output | Approximates | Tutorial |
| --- | --- | --- |
| **Bins** (8 µm and 16 µm) | 10x Visium HD | [`visiumHD_bin_tutorial.ipynb`](tutorial/visiumHD_bin_tutorial.ipynb) |
| **Spots** | 10x Visium | [`visium_spot_tutorial.ipynb`](tutorial/visium_spot_tutorial.ipynb) |
| **Cells** | 10x Xenium | [`xenium_cell_tutorial.ipynb`](tutorial/xenium_cell_tutorial.ipynb) |

[`strong_domain_mix_tutorial.ipynb`](tutorial/strong_domain_mix_tutorial.ipynb)
starts from these configurations and swaps in a cell-type composition with
sharper differences between domains. We use it for the spatial-domain
recovery benchmark.

## Tutorials

All notebooks are in [`tutorial/`](tutorial/); see the
[tutorial guide](tutorial/README.md). They import the installed package, so
you can download one and run it without cloning the repository.

## License

MIT. See [LICENSE](LICENSE).
