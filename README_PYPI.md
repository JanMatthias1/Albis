# Albis: simulation of multi-resolution and multi-dimensional spatial transcriptomics data

**Albis** is a Python package for simulating spatial-transcriptomics data
at cell, bin, and spot resolutions. It builds a three-dimensional tissue
model, generates gene-expression counts and transcript locations, and
aggregates them into observations across tissue sections.

The resulting datasets include spatial coordinates and known simulation
ground truth, making them useful for developing and evaluating methods
such as spatial clustering, cell-type composition estimation,
and slice alignment. Outputs are returned as
[AnnData](https://anndata.readthedocs.io/) objects for use in Python analysis
workflows and can be saved as `.h5ad` files.

## Features

- **Multiple resolutions:** generate cell-level observations, regular bins,
  or spatial spots with configurable sizes and spacing.
- **Multiple sectioning directions:** slice a 3D tissue model along the X,
  Y, or Z axis, retaining both 2D and 3D spatial coordinates.
- **Configurable tissue composition:** control spatial domains, cell types,
  and cell-type mixtures within domains.
- **Technical variation:** introduce slice-specific batch effects and
  coordinate rotations and translations for alignment experiments.
- **Ground-truth annotations:** retain domain and cell-type labels and
  composition information for evaluating downstream analyses.
- **Reproducible workflows:** set random seeds, inspect simulation summaries,
  and visualize results with built-in plotting functions.

The current simulator uses spherical tissue geometry with a central domain
and surrounding wedge-shaped domains. See the GitHub documentation for
the supported parameters and their interpretation.

## Installation

```bash
pip install albis
```

Requires Python 3.10 or newer. Plotting dependencies are included in the
standard installation.

## Quick start

Generate a small example, inspect it, plot its spatial domains, and save it:

```python
import albis as ab

adata = ab.example_data(seed=2025)
print(ab.describe(adata))

figure = ab.plot(adata, view="2d", color="domain_true", slice_id=1)
figure.savefig("albis_example.png", dpi=150, bbox_inches="tight")

ab.save(adata, "albis_example.h5ad")
```

Use `ab.generate_data()` to configure a simulation, including its output
resolution (`"cell"`, `"bin"`, or `"spot"`), sectioning axis, number of
slices, and tissue composition. Lower-level functions provide access to
the tissue-generation and sectioning steps.

## Documentation and tutorials

Full documentation and worked examples are maintained on GitHub:

- [Overview and API guide](https://github.com/JanMatthias1/Albis#readme)
- [Tutorial notebooks](https://github.com/JanMatthias1/Albis/tree/HEAD/tutorial)
- [Low-level simulator documentation](https://github.com/JanMatthias1/Albis/blob/HEAD/LOW_LEVEL_SIMULATOR.md)
- [Source code and issue tracker](https://github.com/JanMatthias1/Albis)

## License

Albis is distributed under the
[MIT license](https://github.com/JanMatthias1/Albis/blob/HEAD/LICENSE).
