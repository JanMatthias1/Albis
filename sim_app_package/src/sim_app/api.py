"""Public API for generating one app-ready AnnData object."""

from collections.abc import Mapping
from copy import deepcopy

from .simulation_sphere import _to_serializable, simulate_3d_molecule_sphere_multires


DEFAULT_PARAMETERS = {
    "output": "bin",
    "slice_axis": "Z",
    "n_slices": 5,
    "include_truth": True,
    "tissue_shape": "sphere",
    "sphere_radius_um": 300.0,
    "n_domains": 4,
    "domain_layout": "core_wedges",
    "core_frac": 0.35,
    "n_cell_types": 4,
    "n_cells": 1_000,
    "allow_cell_overlap": False,
    "cell_radius_kwargs": None,
    "domain_type_mix": None,
    "capture_window_um": (500.0, 500.0),
    "capture_window_center_um": (0.0, 0.0),
    "bin_size_um": 20.0,
    "spot_spacing_um": 100.0,
    "spot_radius_um": 27.5,
    "marker_genes_per_type": 80,
    "noise_gene_frac": 0.10,
    "shared_marker_frac": 0.25,
    "inside_prob": 0.95,
    "assign_k": 8,
    "batch_sigma": 0.15,
    "unaligned_coordinates": True,
    "max_deg": 180.0,
    "max_shift": 200.0,
    "seed": 2025,
}


def _resolve_parameters(parameters, overrides):
    if parameters is None:
        supplied = {}
    elif isinstance(parameters, Mapping):
        supplied = dict(parameters)
    else:
        raise TypeError("parameters must be a mapping or None.")

    supplied.update(overrides)
    unknown = set(supplied) - set(DEFAULT_PARAMETERS)
    if unknown:
        names = ", ".join(sorted(unknown))
        raise ValueError(f"Unsupported generate_data parameter(s): {names}.")

    config = deepcopy(DEFAULT_PARAMETERS)
    config.update(supplied)
    return config


def _validate_parameters(config):
    output = str(config["output"]).lower()
    if output not in {"cell", "bin", "spot"}:
        raise ValueError("output must be one of: 'cell', 'bin', 'spot'.")
    config["output"] = output

    axis = str(config["slice_axis"]).upper()
    if axis not in {"X", "Y", "Z"}:
        raise ValueError("slice_axis must be one of: 'X', 'Y', 'Z'.")
    config["slice_axis"] = axis

    if config["tissue_shape"] != "sphere":
        raise ValueError("Only tissue_shape='sphere' is supported in this release.")
    if config["domain_layout"] != "core_wedges":
        raise ValueError("Only domain_layout='core_wedges' is supported in this release.")

    for key in ("n_slices", "n_domains", "n_cell_types", "n_cells", "marker_genes_per_type", "assign_k"):
        if int(config[key]) <= 0:
            raise ValueError(f"{key} must be positive.")

    for key in ("sphere_radius_um", "bin_size_um", "spot_spacing_um", "spot_radius_um", "inside_prob"):
        if float(config[key]) <= 0:
            raise ValueError(f"{key} must be positive.")

    if not 0 < float(config["inside_prob"]) < 1:
        raise ValueError("inside_prob must be between 0 and 1.")

    for key in ("capture_window_um", "capture_window_center_um"):
        value = config[key]
        if len(value) != 2:
            raise ValueError(f"{key} must contain exactly two values.")
    if any(float(value) <= 0 for value in config["capture_window_um"]):
        raise ValueError("capture_window_um values must be positive.")


def _remove_truth_annotations(adata):
    adata.layers.pop("counts_pre_batch", None)
    for key in ("cell_type_frac_true", "domain_frac_true"):
        if key in adata.obsm:
            del adata.obsm[key]
    for key in ("cell_type_true", "domain_true", "cell_radius", "cell_sigma_for_95pct"):
        if key in adata.obs:
            del adata.obs[key]


def generate_data(parameters=None, /, **overrides):
    """Generate one selected, app-ready :class:`anndata.AnnData` object.

    Parameters are usually supplied as keyword arguments, for example
    ``generate_data(output="bin", slice_axis="Z")``. A mapping may also be
    passed when callers need to build the configuration programmatically. The
    initial public API supports a spherical tissue with the existing
    ``core_wedges`` domain layout.

    The returned object always contains ``obsm["spatial_3d"]`` and aligned
    2D ``obsm["spatial"]`` coordinates. It includes
    ``obsm["spatial_unaligned"]`` unless ``unaligned_coordinates`` is false.
    """
    config = _resolve_parameters(parameters, overrides)
    _validate_parameters(config)

    simulation = simulate_3d_molecule_sphere_multires(
        sphere_R_um=float(config["sphere_radius_um"]),
        n_domains=int(config["n_domains"]),
        core_frac=float(config["core_frac"]),
        n_cells=int(config["n_cells"]),
        allow_cell_overlap=bool(config["allow_cell_overlap"]),
        cell_radius_kwargs=config["cell_radius_kwargs"],
        n_cell_types=int(config["n_cell_types"]),
        domain_type_mix=config["domain_type_mix"],
        marker_genes_per_type=int(config["marker_genes_per_type"]),
        noise_gene_frac=float(config["noise_gene_frac"]),
        shared_marker_frac=float(config["shared_marker_frac"]),
        inside_prob=float(config["inside_prob"]),
        assign_k=int(config["assign_k"]),
        n_slices=int(config["n_slices"]),
        batch_sigma=float(config["batch_sigma"]),
        visium_capture_size_um=tuple(config["capture_window_um"]),
        visium_capture_center_um=tuple(config["capture_window_center_um"]),
        bin_size_um=float(config["bin_size_um"]),
        spot_spacing_um=float(config["spot_spacing_um"]),
        spot_radius_um=float(config["spot_radius_um"]),
        max_deg=float(config["max_deg"]),
        max_shift=float(config["max_shift"]),
        seed=int(config["seed"]),
        output_modalities=(config["output"],),
        slice_axes=(config["slice_axis"],),
    )

    output_key = {
        "cell": "adata_cell_sectioned",
        "bin": "bin_adatas",
        "spot": "spot_adatas",
    }[config["output"]]
    adata = simulation[output_key][config["slice_axis"]]

    if not config["unaligned_coordinates"]:
        adata.obsm.pop("spatial_unaligned", None)
        adata.obsm.pop("spatial_3d_unaligned", None)
        adata.uns.pop("rigid_perturb_inplane", None)
    if not config["include_truth"]:
        _remove_truth_annotations(adata)

    sim_params = _to_serializable(simulation["meta"])
    sim_params["generate_data"] = _to_serializable(config)
    adata.uns["sim_params"] = sim_params
    adata.uns["output"] = {
        "platform": config["output"],
        "slice_axis": config["slice_axis"],
    }
    return adata
