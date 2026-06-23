"""Static visualizations for app-ready simulation outputs."""

import numpy as np


def _load_matplotlib():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise ImportError(
            "Plotting requires the optional dependency. Install it with "
            "`python -m pip install 'sim-app[plot]'`."
        ) from exc
    return plt, Line2D


def _select_indices(adata, slice_id, max_points, seed):
    indices = np.arange(adata.n_obs)
    if slice_id is not None:
        if "slice_id" not in adata.obs:
            raise ValueError("slice_id filtering requires adata.obs['slice_id'].")
        indices = indices[np.asarray(adata.obs["slice_id"], dtype=int) == int(slice_id)]

    if max_points is not None:
        max_points = int(max_points)
        if max_points <= 0:
            raise ValueError("max_points must be positive or None.")
        if indices.size > max_points:
            rng = np.random.default_rng(seed)
            indices = np.sort(rng.choice(indices, size=max_points, replace=False))
    return indices


def _color_values(adata, color, indices):
    if color is None:
        return None, None
    if color in adata.obs:
        return np.asarray(adata.obs[color])[indices], "obs"
    if color in adata.var_names:
        values = adata[indices, color].X
        if hasattr(values, "toarray"):
            values = values.toarray()
        return np.asarray(values).reshape(-1), "gene"
    raise ValueError(f"color must be an obs column or gene name. Received: {color!r}.")


def plot(
    adata,
    *,
    view="2d",
    coordinates="aligned",
    color=None,
    slice_id=None,
    point_size=4.0,
    max_points=50_000,
    seed=0,
    alpha=0.8,
):
    """Plot an aligned or unaligned simulation view.

    Parameters
    ----------
    adata
        An AnnData object returned by :func:`sim_app.generate_data`.
    view
        ``"2d"`` uses slice-plane coordinates; ``"3d"`` uses tissue
        coordinates.
    coordinates
        ``"aligned"`` or ``"unaligned"``. Unaligned 3D coordinates preserve
        the slicing coordinate and apply each slice's rigid in-plane transform.
    color
        An ``adata.obs`` column or a gene name. Omit for a neutral scatter.
    slice_id
        Optionally restrict the plot to one slice.
    point_size, max_points, seed, alpha
        Scatter rendering controls. Downsampling is deterministic for a seed.
    """
    view = str(view).lower()
    coordinates = str(coordinates).lower()
    if view not in {"2d", "3d"}:
        raise ValueError("view must be '2d' or '3d'.")
    if coordinates not in {"aligned", "unaligned"}:
        raise ValueError("coordinates must be 'aligned' or 'unaligned'.")

    coordinate_key = {
        ("2d", "aligned"): "spatial",
        ("2d", "unaligned"): "spatial_unaligned",
        ("3d", "aligned"): "spatial_3d",
        ("3d", "unaligned"): "spatial_3d_unaligned",
    }[(view, coordinates)]
    if coordinate_key not in adata.obsm:
        raise ValueError(f"adata.obsm['{coordinate_key}'] is not available.")

    indices = _select_indices(adata, slice_id, max_points, seed)
    if indices.size == 0:
        raise ValueError("No observations remain after applying the plot filters.")
    coords = np.asarray(adata.obsm[coordinate_key])[indices]
    expected_dims = 2 if view == "2d" else 3
    if coords.ndim != 2 or coords.shape[1] != expected_dims:
        raise ValueError(f"adata.obsm['{coordinate_key}'] must have {expected_dims} columns.")

    plt, Line2D = _load_matplotlib()
    if view == "3d":
        figure = plt.figure()
        axis = figure.add_subplot(projection="3d")
    else:
        figure, axis = plt.subplots()

    values, color_kind = _color_values(adata, color, indices)
    scatter_kwargs = {"s": float(point_size), "alpha": float(alpha), "linewidths": 0}
    if values is None:
        if view == "3d":
            axis.scatter(coords[:, 0], coords[:, 1], coords[:, 2], color="tab:blue", **scatter_kwargs)
        else:
            axis.scatter(coords[:, 0], coords[:, 1], color="tab:blue", **scatter_kwargs)
    elif np.issubdtype(np.asarray(values).dtype, np.number):
        if view == "3d":
            artist = axis.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=values, cmap="viridis", **scatter_kwargs)
        else:
            artist = axis.scatter(coords[:, 0], coords[:, 1], c=values, cmap="viridis", **scatter_kwargs)
        figure.colorbar(artist, ax=axis, label=color)
    else:
        labels, inverse = np.unique(np.asarray(values, dtype=str), return_inverse=True)
        cmap = plt.get_cmap("tab20", max(len(labels), 1))
        colors = cmap(inverse)
        if view == "3d":
            axis.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, **scatter_kwargs)
        else:
            axis.scatter(coords[:, 0], coords[:, 1], c=colors, **scatter_kwargs)
        if len(labels) <= 20:
            handles = [
                Line2D([], [], marker="o", linestyle="", color=cmap(index), label=label)
                for index, label in enumerate(labels)
            ]
            axis.legend(handles=handles, title=color, bbox_to_anchor=(1.02, 1), loc="upper left")

    if view == "3d":
        axis.set_xlabel("X (µm)")
        axis.set_ylabel("Y (µm)")
        axis.set_zlabel("Z (µm)")
    else:
        plane_dims = {"X": ("Y", "Z"), "Y": ("X", "Z"), "Z": ("X", "Y")}
        slice_axis = str(adata.obs["slice_axis"].iloc[0]).upper() if "slice_axis" in adata.obs else "Z"
        dim0, dim1 = plane_dims.get(slice_axis, ("X", "Y"))
        axis.set_xlabel(f"{dim0} (µm)")
        axis.set_ylabel(f"{dim1} (µm)")
        axis.set_aspect("equal", adjustable="box")

    title = f"{view.upper()} {coordinates} coordinates"
    if color_kind is not None:
        title += f" — {color}"
    if slice_id is not None:
        title += f" — slice {slice_id}"
    axis.set_title(title)
    figure.tight_layout()
    return figure
