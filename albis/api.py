"""Public API for generating one app-ready AnnData object."""

from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path

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
    "core_bump_amp": 0.12,
    "wedge_angle_amp_deg": 12.0,
    "noise_terms": 6,
    "noise_freq_range": (0.8, 2.2),
    "boundary_fuzz_width_deg": 0.0,
    "boundary_fuzz_flip_prob": 0.0,
    "core_fuzz_width_um": 0.0,
    "core_fuzz_flip_prob": 0.0,
    "n_cell_types": 4,
    "n_cells": 1_000,
    "allow_cell_overlap": False,
    "cell_radius_kwargs": None,
    "domain_type_mix": None,
    "xenium_capture_window_um": (12_000.0, 24_000.0),
    "capture_window_um": "platform",
    "capture_window_center_um": (0.0, 0.0),
    "bin_size_um": 16.0,
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
    "base_seed_unaligned": 12345,
    "sync_unaligned_seed": False,
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

    for key in ("n_slices", "n_domains", "n_cell_types", "n_cells", "marker_genes_per_type", "assign_k", "noise_terms"):
        if int(config[key]) <= 0:
            raise ValueError(f"{key} must be positive.")

    for key in ("sphere_radius_um", "bin_size_um", "spot_spacing_um", "spot_radius_um", "inside_prob"):
        if float(config[key]) <= 0:
            raise ValueError(f"{key} must be positive.")
    for key in ("core_bump_amp", "wedge_angle_amp_deg", "boundary_fuzz_width_deg", "core_fuzz_width_um"):
        if float(config[key]) < 0:
            raise ValueError(f"{key} must be nonnegative.")

    if not 0 < float(config["inside_prob"]) < 1:
        raise ValueError("inside_prob must be between 0 and 1.")
    for key in ("boundary_fuzz_flip_prob", "core_fuzz_flip_prob"):
        if not 0 <= float(config[key]) <= 1:
            raise ValueError(f"{key} must be between 0 and 1.")

    for key in ("xenium_capture_window_um", "capture_window_center_um", "noise_freq_range"):
        value = config[key]
        if len(value) != 2:
            raise ValueError(f"{key} must contain exactly two values.")
    for key in ("xenium_capture_window_um", "noise_freq_range"):
        if any(float(value) <= 0 for value in config[key]):
            raise ValueError(f"{key} values must be positive.")

    # capture_window_um: "platform" (per-modality real window), False (no crop),
    # or an explicit (width, height) pair applied to every modality.
    cw = config["capture_window_um"]
    if cw is None or (isinstance(cw, str) and cw.lower() == "platform"):
        config["capture_window_um"] = "platform"
    elif cw is False:
        config["capture_window_um"] = False
    elif isinstance(cw, str) or not hasattr(cw, "__len__") or len(cw) != 2 or any(float(v) <= 0 for v in cw):
        raise ValueError(
            "capture_window_um must be 'platform', False, or a (width, height) pair of positive numbers."
        )
    else:
        config["capture_window_um"] = (float(cw[0]), float(cw[1]))


def _remove_truth_annotations(adata):
    adata.layers.pop("counts_pre_batch", None)
    for key in ("cell_type_frac_true", "domain_frac_true"):
        if key in adata.obsm:
            del adata.obsm[key]
    for key in ("cell_type_true", "domain_true", "cell_radius", "cell_sigma_for_95pct"):
        if key in adata.obs:
            del adata.obs[key]


def _matrix_total(X):
    total = X.sum()
    if hasattr(total, "item"):
        total = total.item()
    return int(total)


def _matrix_nnz(X):
    if hasattr(X, "nnz"):
        return int(X.nnz)
    return int((X != 0).sum())


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

    ``capture_window_um`` controls the capture area cropped from each slice:

    * ``"platform"`` (default) -- the real slide area for the chosen ``output``:
      6.5 x 6.5 mm for ``"bin"``/``"spot"`` (Visium / Visium HD),
      ``xenium_capture_window_um`` (12 x 24 mm) for ``"cell"`` (Xenium).
    * ``False`` -- no crop; bin/spot grids span the molecule bounding box, with
      no empty border.
    * a ``(width, height)`` pair -- that window, applied to any modality.

    When the tissue is smaller than the capture window, ``"bin"``/``"spot"``
    outputs include all-zero observations around the tissue. That is expected.
    They are flagged with ``obs["is_empty"]`` and labelled
    ``domain_true``/``cell_type_true`` = ``"unassigned"``; drop them with
    ``adata[~adata.obs["is_empty"]]`` or any per-observation minimum-count QC
    filter.
    """
    config = _resolve_parameters(parameters, overrides)
    _validate_parameters(config)

    # Resolve the capture window per output modality. "platform" -> the real
    # slide area (Visium/VisiumHD 6.5 mm for bin/spot, Xenium 12x24 mm for
    # cell); False -> no crop; an explicit pair -> that window for any modality.
    cw = config["capture_window_um"]
    if cw is False:
        bin_spot_window = False
        cell_window = False
    elif cw == "platform":
        bin_spot_window = (6500.0, 6500.0)
        cell_window = tuple(config["xenium_capture_window_um"])
    else:
        bin_spot_window = tuple(cw)
        cell_window = tuple(cw)

    simulation = simulate_3d_molecule_sphere_multires(
        sphere_R_um=float(config["sphere_radius_um"]),
        xenium_capture_window_um=cell_window,
        n_domains=int(config["n_domains"]),
        core_frac=float(config["core_frac"]),
        core_bump_amp=float(config["core_bump_amp"]),
        wedge_angle_amp_deg=float(config["wedge_angle_amp_deg"]),
        noise_terms=int(config["noise_terms"]),
        noise_freq_range=tuple(config["noise_freq_range"]),
        boundary_fuzz_width_deg=float(config["boundary_fuzz_width_deg"]),
        boundary_fuzz_flip_prob=float(config["boundary_fuzz_flip_prob"]),
        core_fuzz_width_um=float(config["core_fuzz_width_um"]),
        core_fuzz_flip_prob=float(config["core_fuzz_flip_prob"]),
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
        capture_window_um=bin_spot_window,
        capture_window_center_um=tuple(config["capture_window_center_um"]),
        bin_size_um=float(config["bin_size_um"]),
        spot_spacing_um=float(config["spot_spacing_um"]),
        spot_radius_um=float(config["spot_radius_um"]),
        max_deg=float(config["max_deg"]),
        max_shift=float(config["max_shift"]),
        base_seed_unaligned=int(config["base_seed_unaligned"]),
        sync_unaligned_seed=bool(config["sync_unaligned_seed"]),
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


def describe(adata):
    """Return a compact summary of an app-ready :class:`anndata.AnnData` object."""
    output = dict(adata.uns.get("output", {}))
    slice_ids = []
    if "slice_id" in adata.obs:
        slice_ids = sorted(int(value) for value in adata.obs["slice_id"].unique())

    truth_keys = []
    for key in ("counts_pre_batch",):
        if key in adata.layers:
            truth_keys.append(f"layers['{key}']")
    for key in ("cell_type_true", "domain_true", "cell_radius", "cell_sigma_for_95pct"):
        if key in adata.obs:
            truth_keys.append(f"obs['{key}']")
    for key in ("cell_type_frac_true", "domain_frac_true"):
        if key in adata.obsm:
            truth_keys.append(f"obsm['{key}']")

    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "total_counts": _matrix_total(adata.X),
        "nonzero_counts": _matrix_nnz(adata.X),
        "platform": output.get("platform"),
        "slice_axis": output.get("slice_axis"),
        "n_slices": len(slice_ids),
        "slice_ids": slice_ids,
        "obs_columns": list(map(str, adata.obs.columns)),
        "var_columns": list(map(str, adata.var.columns)),
        "layers": list(map(str, adata.layers.keys())),
        "obsm_keys": list(map(str, adata.obsm.keys())),
        "uns_keys": list(map(str, adata.uns.keys())),
        "truth_annotations": truth_keys,
    }


def save(adata, path):
    """Write an AnnData object to ``path`` and return the file path."""
    output_path = Path(path)
    adata.write_h5ad(output_path)
    return output_path


def example_data(**overrides):
    """Generate a small example dataset for tutorials, tests, and quick checks."""
    config = {
        "output": "bin",
        "slice_axis": "Z",
        "n_cells": 250,
        "n_slices": 3,
        "sphere_radius_um": 200.0,
        "capture_window_um": (300.0, 300.0),
        "bin_size_um": 30.0,
        "seed": 2025,
    }
    config.update(overrides)
    return generate_data(**config)
