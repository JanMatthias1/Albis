#!/usr/bin/env python3
"""
3D transcript-level Spatial Transcriptomics Simulator on a SPHERE (AnnData).

Core idea:
  1) Simulate CELLS in a 3D sphere: centroid, radius, domain, cell type
  2) Simulate per-cell per-gene counts with NB, with marker genes + noise genes
  3) Expand counts into RNA MOLECULES (transcripts): for each molecule sample a 3D coordinate
     from a Gaussian centered at the cell centroid, with sigma chosen so 95% molecules lie
     inside the cell boundary (sphere of radius R_cell)
  4) Build multiple resolutions from the SAME molecule truth:
     - cell-level (Xenium-like): assign each molecule to whichever cell boundary contains it
     - bin-level (VisiumHD-like): voxel/bin counts (and 2D slices as needed)
     - spot-level (Visium-like): larger capture regions
  5) Cut sphere into 10 slices along x, y, z (separately) and add slice-specific batch effects
     (SpatialMNN-like, implemented as gene-specific multiplicative effects per slice).

Notes:
- We allow ~5% molecules to fall outside the source cell boundary by NOT doing rejection sampling.
- Spillover happens automatically because we assign molecules by "observed location":
    the cell that contains a molecule gets the count, regardless of source cell.
- We drop molecules that do not fall inside ANY cell boundary (no "unassigned" label).

-------------------------------
(1) Domain difficulty ↑:
    - Core domain boundary becomes "bumpy" via smooth noise in 3D.
    - Wedge domain boundaries become wavy/irregular via smooth angular noise.
    - (Optional) extra "boundary fuzz" flips a small fraction of labels near boundaries.

(2) Different platform capture windows (same tissue sphere):
    - Tissue sphere diameter defaults to 12 mm => sphere_R_um = 6000 µm.
    - Xenium-like capture window default: 12×24 mm (crops cell-level output)
    - Visium/VisiumHD capture window default: 6.5×6.5 mm (crops bins/spots)
    - Spots: radius 27.5 µm (55 µm diam) + spacing 100 µm (unchanged).

(3) Coordinate storage:
    - obsm["spatial_3d"]: aligned 3D coords
    - obsm["spatial"]: aligned 2D coords (plane for each axis)
    - obsm["spatial_unaligned"]: per-slice rigidly perturbed 2D coords
"""

import os
import warnings
import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from scipy.stats import chi2
import anndata as ad


# -------------------------
#   Math helpers
# -------------------------

def mean_disp_to_nbinom_params(mu, theta):
    """
    NB2 parameterization:
      Var = mu + mu^2 / theta

    Sampling form (NumPy):
      r = theta
      p = theta / (theta + mu)
      X ~ NegBin(r, p)
    """
    r = theta
    p = theta / (theta + mu)
    return r, p


def sigma_from_radius_inside_prob(radius, inside_prob=0.95, dim=3):
    """
    Choose isotropic Gaussian sigma so that:
      P(||X - mu|| <= radius) = inside_prob
    when X ~ N(mu, sigma^2 I_dim).
    """
    c = chi2.ppf(inside_prob, df=dim)
    return radius / np.sqrt(c + 1e-12)


def _to_serializable(x):
    if isinstance(x, dict):
        return {k: _to_serializable(v) for k, v in x.items()}
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (list, tuple)):
        return [_to_serializable(v) for v in x]
    return x



#   Smooth noise utilities (for irregular domains)

def smooth_noise_3d(points_xyz, n_terms=6, freq_range=(0.6, 2.0), seed=0):
    """
    Smooth-ish random field by summing sinusoids:
      f(p) = sum_k a_k * sin(w_k · p + b_k)
    Returned noise is normalized to approx [-1, 1].

    points_xyz: (N,3) in same units (µm). To control smoothness, tune freq_range.
    """
    rng = np.random.default_rng(seed)
    P = np.asarray(points_xyz, dtype=np.float64)
    N = P.shape[0]
    out = np.zeros(N, dtype=np.float64)

    for _ in range(n_terms):
        # random direction for wavevector
        v = rng.normal(size=3)
        v /= (np.linalg.norm(v) + 1e-12)

        # frequency in "1 / (units)"
        f = rng.uniform(freq_range[0], freq_range[1])
        w = v * f

        phase = rng.uniform(0, 2 * np.pi)
        amp = rng.uniform(0.4, 1.0)

        out += amp * np.sin(P @ w + phase)

    # normalize to roughly [-1,1]
    m = np.max(np.abs(out)) + 1e-12
    return (out / m).astype(np.float32)


def smooth_noise_on_unit_sphere(unit_dirs, n_terms=6, freq_range=(1.0, 3.0), seed=0):
    """
    Same as smooth_noise_3d but used on direction vectors (unit sphere).
    This makes core boundary bumps depend on direction, not distance.
    """
    return smooth_noise_3d(unit_dirs, n_terms=n_terms, freq_range=freq_range, seed=seed)



#   Geometry / sampling

def sample_uniform_in_sphere(n, R=6000.0, center=(0.0, 0.0, 0.0), rng=None):
    """
    Uniformly sample points inside a 3D ball of radius R.
    Units: µm (by default).
    """
    if rng is None:
        rng = np.random.default_rng(0)

    v = rng.normal(size=(n, 3))
    v /= (np.linalg.norm(v, axis=1, keepdims=True) + 1e-12)

    u = rng.uniform(0.0, 1.0, size=n)
    r = R * np.cbrt(u)

    pts = v * r[:, None]
    pts[:, 0] += center[0]
    pts[:, 1] += center[1]
    pts[:, 2] += center[2]
    return pts


#   Noisy domains

def sphere_domain_labels_noisy(
    coords,
    center=(0.0, 0.0, 0.0),
    n_domains=4,
    core_frac=0.35,
    # noise controls
    core_bump_amp=0.12,          # fractional bump on core radius (0.0 -> perfect sphere)
    wedge_angle_amp_deg=12.0,    # angular wobble of wedge boundaries (degrees)
    noise_terms=6,
    noise_freq_range=(0.8, 2.2),
    seed=0,
    # optional wedge boundary fuzz (OFF by default)
    boundary_fuzz_width_deg=0.0,  # e.g. 6.0
    boundary_fuzz_flip_prob=0.0,  # e.g. 0.10
    # optional core↔wedge boundary fuzz (OFF by default)
    core_fuzz_width_um=0.0,       # e.g. 120.0 (µm band thickness around r_core)
    core_fuzz_flip_prob=0.0,      # e.g. 0.10
):
    """
    "3 wedges + core", but boundaries are irregular.

    - Core boundary: r_core(direction) = r_core0 * (1 + core_bump_amp * noise(direction))
    - Wedge boundaries: phi' = phi + angle_amp * noise(point)

    boundary_fuzz_* flips some points near wedge↔wedge boundaries (outside core).
    core_fuzz_* flips some points in a radial band around the core boundary (core↔wedge interface).
    """
    assert n_domains >= 2
    rng = np.random.default_rng(seed)

    xyz = np.asarray(coords, dtype=np.float32)
    x = xyz[:, 0] - center[0]
    y = xyz[:, 1] - center[1]
    z = xyz[:, 2] - center[2]
    r = np.sqrt(x * x + y * y + z * z)

    R_est = np.max(r) + 1e-12
    r_core0 = core_frac * R_est

    # phi in [0, 2pi)
    phi = np.arctan2(y, x)
    phi = np.where(phi < 0, phi + 2.0 * np.pi, phi)

    # noise for wedge angle perturbation (depends on position)
    noise_pos = smooth_noise_3d(
        points_xyz=np.stack([x, y, z], axis=1),
        n_terms=noise_terms,
        freq_range=noise_freq_range,
        seed=seed + 11,
    )

    angle_amp = np.deg2rad(wedge_angle_amp_deg)
    phi_warp = np.mod(phi + angle_amp * noise_pos, 2.0 * np.pi)

    n_wedges = n_domains - 1
    wedge_width = 2.0 * np.pi / n_wedges

    domain = np.floor(phi_warp / wedge_width).astype(int)
    domain = np.clip(domain, 0, n_wedges - 1)

    # noisy core boundary depends on direction
    unit = np.stack([x, y, z], axis=1)
    unit /= (np.linalg.norm(unit, axis=1, keepdims=True) + 1e-12)

    noise_dir = smooth_noise_on_unit_sphere(
        unit_dirs=unit,
        n_terms=noise_terms,
        freq_range=(1.0, 3.0),
        seed=seed + 23,
    )
    r_core = r_core0 * (1.0 + core_bump_amp * noise_dir)
    r_core = np.clip(r_core, 0.05 * R_est, 0.95 * R_est)

    core_label = n_domains - 1
    domain[r <= r_core] = core_label

    # core↔wedge boundary fuzz

    if core_fuzz_width_um > 0 and core_fuzz_flip_prob > 0:
        band = np.abs(r - r_core) <= float(core_fuzz_width_um)
        if np.any(band):
            band_idx = np.where(band)[0]
            flip = rng.uniform(size=band_idx.size) < core_fuzz_flip_prob
            flip_idx = band_idx[flip]

            if flip_idx.size:
                # native wedge based on warped angle
                native_wedge = np.floor(phi_warp[flip_idx] / wedge_width).astype(int)
                native_wedge = np.clip(native_wedge, 0, n_wedges - 1)

                is_core = domain[flip_idx] == core_label

                # core -> wedge
                if np.any(is_core):
                    fi = flip_idx[is_core]
                    domain[fi] = native_wedge[is_core]

                # wedge -> core
                if np.any(~is_core):
                    fi = flip_idx[~is_core]
                    domain[fi] = core_label

    # -----------------------------------------
    # Optional: wedge↔wedge boundary fuzz (ONLY outside core)
    # -----------------------------------------
    if boundary_fuzz_width_deg > 0 and boundary_fuzz_flip_prob > 0:
        width = np.deg2rad(boundary_fuzz_width_deg)
        frac = np.mod(phi_warp, wedge_width)
        dist_to_edge = np.minimum(frac, wedge_width - frac)

        outside_core = r > r_core
        near = outside_core & (dist_to_edge <= width)
        if np.any(near):
            near_idx = np.where(near)[0]
            flip = rng.uniform(size=near_idx.size) < boundary_fuzz_flip_prob
            flip_idx = near_idx[flip]
            if flip_idx.size:
                go_left = frac[flip_idx] < (wedge_width / 2)
                dom_new = domain[flip_idx].copy()
                dom_new[go_left] = (dom_new[go_left] - 1) % n_wedges
                dom_new[~go_left] = (dom_new[~go_left] + 1) % n_wedges
                domain[flip_idx] = dom_new

    return domain


# -------------------------
#   Cell placement (centroids + radii)
# -------------------------

def sample_cell_radii(
    n_cells,
    radius_dist="lognormal",
    r_mean=6.0,
    r_sigma=0.35,
    r_min=3.0,
    r_max=12.0,
    rng=None,
):
    if rng is None:
        rng = np.random.default_rng(0)

    if radius_dist == "lognormal":
        radii = np.exp(rng.normal(loc=np.log(max(r_mean, 1e-6)), scale=r_sigma, size=n_cells))
    elif radius_dist == "normal":
        radii = rng.normal(loc=r_mean, scale=r_sigma, size=n_cells)
    else:
        raise ValueError("radius_dist must be 'lognormal' or 'normal'.")

    radii = np.clip(radii, r_min, r_max)
    return radii


def sample_nonoverlapping_cells_in_sphere(
    n_cells,
    sphere_R=6000.0,
    center=(0.0, 0.0, 0.0),
    radius_kwargs=None,
    allow_overlap=False,
    max_attempts=50_000_000,
    rebuild_every=2000,
    seed=0,
):
    rng = np.random.default_rng(seed)
    if radius_kwargs is None:
        radius_kwargs = dict(radius_dist="lognormal", r_mean=6.0, r_sigma=0.35, r_min=3.0, r_max=12.0)

    centers = []
    radii = []
    tree = None
    r_max_current = 0.0
    radii_pool = sample_cell_radii(n_cells, rng=rng, **radius_kwargs)

    accepted, attempts = 0, 0
    while accepted < n_cells and attempts < max_attempts:
        attempts += 1
        c = sample_uniform_in_sphere(1, R=sphere_R, center=center, rng=rng)[0]
        r = float(radii_pool[accepted])

        if np.linalg.norm(c - np.array(center)) + r > sphere_R:
            continue

        if (not allow_overlap) and (tree is not None) and (accepted > 0):
            idxs = tree.query_ball_point(c, r + r_max_current)
            if idxs:
                neigh_centers = np.asarray([centers[i] for i in idxs])
                neigh_radii = np.asarray([radii[i] for i in idxs])
                d = np.linalg.norm(neigh_centers - c[None, :], axis=1)
                if np.any(d < (neigh_radii + r)):
                    continue

        centers.append(c)
        radii.append(r)
        accepted += 1
        r_max_current = max(r_max_current, r)

        if (accepted % rebuild_every) == 0:
            tree = cKDTree(np.asarray(centers))

    if accepted < n_cells:
        raise RuntimeError(
            f"Could only place {accepted}/{n_cells} cells after {attempts} attempts. "
            f"Try fewer cells, allow_overlap=True, smaller radii, or larger sphere_R."
        )

    centers = np.asarray(centers)
    radii = np.asarray(radii)
    tree = cKDTree(centers)
    return centers, radii, tree


# -------------------------
#   Genes: marker sets + noise genes
# -------------------------

def make_marker_and_noise_genes(
    n_cell_types=8,
    marker_genes_per_type=80,
    noise_gene_frac=0.10,
    shared_marker_frac=0.25,
    seed=0,
):
    rng = np.random.default_rng(seed)

    n_shared = int(round(shared_marker_frac * marker_genes_per_type))
    n_unique = marker_genes_per_type - n_shared

    shared_genes = list(range(n_shared))

    marker_sets = []
    next_idx = n_shared
    for _t in range(n_cell_types):
        uniq = list(range(next_idx, next_idx + n_unique))
        next_idx += n_unique
        marker_sets.append(set(shared_genes + uniq))

    n_marker_total = next_idx

    n_noise = int(np.ceil((noise_gene_frac / max(1e-9, (1.0 - noise_gene_frac))) * n_marker_total))
    noise_start = n_marker_total
    noise_end = n_marker_total + n_noise
    noise_idx = set(range(noise_start, noise_end))

    G = noise_end
    gene_names = [f"G{g+1}" for g in range(G)]

    gene_class = np.array(["noise"] * G, dtype=object)
    for t, ms in enumerate(marker_sets):
        for g in ms:
            gene_class[g] = f"marker_type{t+1}"
    for g in shared_genes:
        gene_class[g] = "marker_shared"

    gene_info = dict(
        n_genes=G,
        n_marker_total=n_marker_total,
        n_noise=n_noise,
        shared_genes=shared_genes,
        marker_genes_per_type=marker_genes_per_type,
        noise_gene_frac=noise_gene_frac,
        shared_marker_frac=shared_marker_frac,
    )
    return gene_names, marker_sets, noise_idx, gene_class, gene_info


# -------------------------
#   Cell type sampling by domain mixture
# -------------------------

def sample_cell_types_by_domain(domain_ids, domain_type_mix, rng):
    domain_ids = np.asarray(domain_ids)
    n_domains, n_cell_types = domain_type_mix.shape
    cell_type = np.empty(domain_ids.shape[0], dtype=int)
    for d in range(n_domains):
        idx = np.where(domain_ids == d)[0]
        if idx.size:
            cell_type[idx] = rng.choice(n_cell_types, size=idx.size, p=domain_type_mix[d])
    return cell_type


# -------------------------
#   Simulate per-cell gene counts (NB)
# -------------------------

def simulate_cell_gene_counts_nb(
    cell_type_ids,
    domain_ids,
    marker_sets,
    noise_idx,
    base_gene_lognormal=(1.0, 0.7),
    marker_foldchange=3.5,
    shared_marker_foldchange=2.5,
    noise_scale=0.9,
    cell_size_lognormal=(0.0, 0.35),
    domain_size_factors=None,
    theta=25.0,
    theta_jitter=2.0,
    seed=0,
    sparse_X=True,
):
    rng = np.random.default_rng(seed)
    cell_type_ids = np.asarray(cell_type_ids)
    domain_ids = np.asarray(domain_ids)
    n_cells = cell_type_ids.shape[0]
    n_cell_types = len(marker_sets)

    all_marker = set().union(*marker_sets)
    n_genes = max(max(all_marker, default=-1), max(noise_idx, default=-1)) + 1

    log_mu, log_sigma = base_gene_lognormal
    base_gene = np.exp(rng.normal(loc=log_mu, scale=log_sigma, size=n_genes))

    size_mu, size_sigma = cell_size_lognormal
    size_factor = np.exp(rng.normal(loc=size_mu, scale=size_sigma, size=n_cells))

    n_domains = int(np.max(domain_ids)) + 1
    if domain_size_factors is None:
        domain_size_factors = np.ones(n_domains, dtype=float)
    else:
        domain_size_factors = np.asarray(domain_size_factors, dtype=float)
        assert domain_size_factors.shape[0] == n_domains

    if theta_jitter > 0:
        theta_g = np.maximum(1e-3, rng.normal(loc=theta, scale=theta_jitter, size=n_genes))
    else:
        theta_g = np.full(n_genes, theta, dtype=float)

    shared_markers = set.intersection(*marker_sets) if n_cell_types > 1 else set(marker_sets[0])

    multipliers = np.ones((n_cell_types, n_genes), dtype=float)
    for g in noise_idx:
        multipliers[:, g] *= noise_scale
    for g in shared_markers:
        multipliers[:, g] *= shared_marker_foldchange
    for t in range(n_cell_types):
        for g in marker_sets[t]:
            if g in shared_markers:
                continue
            multipliers[t, g] *= marker_foldchange

    rows, cols, data = [], [], []
    for c in range(n_cells):
        t = cell_type_ids[c]
        d = domain_ids[c]
        mu = base_gene * size_factor[c] * domain_size_factors[d] * multipliers[t, :]

        r_nb, p_nb = mean_disp_to_nbinom_params(mu, theta_g)
        r_nb = np.clip(r_nb, 1e-6, 1e6)
        p_nb = np.clip(p_nb, 1e-9, 1 - 1e-9)
        x = rng.negative_binomial(r_nb, p_nb)

        nz = np.where(x > 0)[0]
        if nz.size:
            rows.append(np.full(nz.size, c, dtype=np.int32))
            cols.append(nz.astype(np.int32))
            data.append(x[nz].astype(np.int32))

    if rows:
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        data = np.concatenate(data)
        X = csr_matrix((data, (rows, cols)), shape=(n_cells, n_genes))
    else:
        X = csr_matrix((n_cells, n_genes), dtype=np.int32)

    return X if sparse_X else X.toarray(), base_gene, size_factor, theta_g, multipliers, shared_markers


# -------------------------
#   Molecules: sample coordinates + observed assignment (spillover)
# -------------------------

def sample_molecule_coords_for_cell(n_mols, centroid, radius, inside_prob=0.95, dim=3, rng=None):
    if rng is None:
        rng = np.random.default_rng(0)
    sigma = sigma_from_radius_inside_prob(radius, inside_prob=inside_prob, dim=dim)
    coords = rng.normal(loc=0.0, scale=sigma, size=(n_mols, dim)) + centroid[None, :]
    return coords, sigma


def assign_points_to_cells_by_containment(points_xyz, cell_centers, cell_radii, tree=None, k=8):
    points_xyz = np.asarray(points_xyz)
    if tree is None:
        tree = cKDTree(cell_centers)

    dists, idxs = tree.query(points_xyz, k=min(k, cell_centers.shape[0]))
    if idxs.ndim == 1:
        idxs = idxs[:, None]
        dists = dists[:, None]

    assigned = np.full(points_xyz.shape[0], -1, dtype=np.int32)
    for j in range(idxs.shape[1]):
        cand = idxs[:, j]
        ok = (assigned < 0) & (dists[:, j] <= cell_radii[cand])
        assigned[ok] = cand[ok]
    return assigned


def aggregate_molecules_to_cell_level(gene_ids, assigned_cell_ids, n_cells, n_genes):
    gene_ids = np.asarray(gene_ids, dtype=np.int32)
    assigned_cell_ids = np.asarray(assigned_cell_ids, dtype=np.int32)

    mask = assigned_cell_ids >= 0
    if not np.any(mask):
        return csr_matrix((n_cells, n_genes), dtype=np.int32)

    rows = assigned_cell_ids[mask]
    cols = gene_ids[mask]
    key = rows.astype(np.int64) * np.int64(n_genes) + cols.astype(np.int64)
    uniq, cnt = np.unique(key, return_counts=True)
    r = (uniq // n_genes).astype(np.int32)
    c = (uniq % n_genes).astype(np.int32)
    return csr_matrix((cnt.astype(np.int32), (r, c)), shape=(n_cells, n_genes))


# -------------------------
#   Slices + batch effect
# -------------------------

def slice_id_along_axis(x, R, n_slices=10):
    edges = np.linspace(-R, R, n_slices + 1)
    sid = np.searchsorted(edges, x, side="right") - 1
    sid = np.where((sid >= 0) & (sid < n_slices), sid, -1)
    return sid.astype(np.int32), edges


def apply_batch_effect_per_slice_sparse(X_csr, slice_ids, n_slices, batch_sigma=0.15, seed=0):
    """
    factor_{s,g} = exp(N(0, batch_sigma^2))
    then x -> Poisson(x * factor_{s,g})
    """
    rng = np.random.default_rng(seed)
    X = X_csr.tocsr(copy=True)
    n_obs, n_genes = X.shape

    slice_ids = np.asarray(slice_ids).astype(int)
    factors = {}

    for s in range(n_slices):
        rows = np.where(slice_ids == s)[0]
        if rows.size == 0:
            continue

        logfc = rng.normal(loc=0.0, scale=batch_sigma, size=n_genes)
        fac = np.exp(logfc).astype(np.float32)
        factors[int(s)] = fac

        for r in rows:
            start, end = X.indptr[r], X.indptr[r + 1]
            if start == end:
                continue
            cols = X.indices[start:end]
            dat = X.data[start:end].astype(np.float32)
            lam = dat * fac[cols]
            X.data[start:end] = rng.poisson(lam).astype(np.int32)

    X.eliminate_zeros()
    return X, factors



#   Rigid perturbation

def random_rigid_inplane_2d(
    coords2d: np.ndarray,
    max_deg: float = 10.0,
    max_shift: float = 50.0,
    rng: np.random.Generator | None = None,
):
    """
    Random rigid transform in 2D: rotation + translation.

    coords2d: (n,2)
    Returns: (coords_new, R2x2, t2,)
    """
    if rng is None:
        rng = np.random.default_rng()

    P = np.asarray(coords2d, dtype=float)
    if P.ndim != 2 or P.shape[1] != 2:
        raise ValueError(f"coords2d must be (n,2). Got {P.shape}.")

    theta = np.deg2rad(rng.uniform(-max_deg, max_deg))
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s],
                  [s,  c]], dtype=float)

    t = rng.uniform(-max_shift, max_shift, size=(2,)).astype(float)

    out = P @ R.T + t[None, :]
    return out.astype(np.float32), R, t



def axis_plane_dims(axis_letter):
    """
    For a given slicing axis, define the 2D plane dims:
      axis X => plane YZ => dims (1,2)
      axis Y => plane XZ => dims (0,2)
      axis Z => plane XY => dims (0,1)
    """
    ax = axis_letter.upper()
    if ax == "X":
        return (1, 2)
    if ax == "Y":
        return (0, 2)
    if ax == "Z":
        return (0, 1)
    raise ValueError("axis_letter must be X/Y/Z")


def add_spatial_keys_for_axis(
    adata,
    axis_letter,
    slice_key="slice_id",
    spatial3d_key_candidates=("spatial_3d", "spatial"),
    spatial2d_key="spatial",
    unaligned_key="spatial_unaligned",
    base_seed=12345,
    max_deg=180.0,
    max_shift=200.0,
):
    """
    Adds/overwrites:
      - obsm['spatial_3d']: aligned 3D coords
      - obsm['spatial']   : aligned 2D coords for this axis plane
      - obsm['spatial_unaligned']: per-slice rigid(2D) coords (rotation+translation)

    Also stores per-slice R/t and centroid deltas in:
      uns['rigid_perturb_inplane'][unaligned_key][slice_id]
    """
    # --- get aligned 3D coords ---
    coords3d = None
    for k in spatial3d_key_candidates:
        if k in adata.obsm:
            coords3d = np.asarray(adata.obsm[k])
            break
    if coords3d is None:
        raise KeyError(f"Need one of {spatial3d_key_candidates} in adata.obsm.")

    if coords3d.ndim != 2 or coords3d.shape[1] < 3:
        raise ValueError(f"Expected (n,>=3) coords for 3D. Got {coords3d.shape}.")

    coords3d = coords3d[:, :3].astype(np.float32)
    adata.obsm["spatial_3d"] = coords3d

    # --- choose correct plane for this cutting axis ---
    dims = axis_plane_dims(axis_letter)  # X->(1,2), Y->(0,2), Z->(0,1)
    coords2d = coords3d[:, dims].astype(np.float32)
    adata.obsm[spatial2d_key] = coords2d

    # --- per-slice rigid perturbation in THAT plane ---
    if slice_key not in adata.obs.columns:
        raise KeyError(f"Expected adata.obs['{slice_key}'].")

    sids = np.asarray(adata.obs[slice_key]).astype(int)
    out2d = coords2d.copy()

    adata.uns.setdefault("rigid_perturb_inplane", {})
    adata.uns["rigid_perturb_inplane"].setdefault(unaligned_key, {})

    for sid in np.unique(sids):
        if sid < 0:
            continue  # don't perturb invalid slices
        idx = np.where(sids == sid)[0]
        if idx.size == 0:
            continue

        rng = np.random.default_rng(base_seed + int(sid))

        before_mean = out2d[idx].mean(axis=0)
        coords_new, R, t = random_rigid_inplane_2d(
            out2d[idx],
            max_deg=max_deg,
            max_shift=max_shift,
            rng=rng,
        )
        after_mean = coords_new.mean(axis=0)

        out2d[idx] = coords_new

        adata.uns["rigid_perturb_inplane"][unaligned_key][str(int(sid))] = {
            "max_deg": float(max_deg),
            "max_shift": float(max_shift),
            "seed": int(base_seed + int(sid)),
            "R": R.tolist(),
            "t": t.tolist(),  # shape (2,)
            "mean_before": before_mean.tolist(),
            "mean_after": after_mean.tolist(),
            "mean_delta": (after_mean - before_mean).tolist(),
        }

    adata.obsm[unaligned_key] = out2d
    spatial3d_unaligned = coords3d.copy()
    spatial3d_unaligned[:, dims] = out2d
    adata.obsm["spatial_3d_unaligned"] = spatial3d_unaligned


# -------------------------
#   Bin / Spot aggregation with capture windows
# -------------------------

def _mix_fraction(obs_id, labels, n_classes, n_obs):
    key = obs_id.astype(np.int64) * np.int64(n_classes) + labels.astype(np.int64)
    uniq, cnt = np.unique(key, return_counts=True)
    r = (uniq // n_classes).astype(np.int32)
    c = (uniq % n_classes).astype(np.int32)
    mat = csr_matrix((cnt.astype(np.int32), (r, c)), shape=(n_obs, n_classes))
    s = np.asarray(mat.sum(axis=1)).ravel()
    frac = np.zeros((n_obs, n_classes), dtype=np.float32)
    nz = s > 0
    if np.any(nz):
        frac[nz] = (mat[nz].toarray() / s[nz, None]).astype(np.float32)
    return frac


def _empty_obs_mask(frac):
    """True for observations with no source molecules (all-zero fraction row).

    For bin/spot output these are grid cells that fall inside the capture
    window but outside the tissue. They are still returned -- expected, and
    dropped downstream with a per-observation min-count QC filter -- just
    flagged and labelled honestly rather than silently mislabelled.
    """
    return np.asarray(frac).sum(axis=1) == 0


def _labels_from_fracs(frac, names, empty_label="unassigned"):
    """Per-row argmax label, but rows that sum to 0 get ``empty_label``.

    ``np.argmax`` of an all-zero row returns 0, which would silently label
    every empty bin/spot as the first domain / first cell type. Those
    observations have no ground-truth class, so they are marked explicitly.
    Nothing is dropped here -- callers filter empty observations themselves.
    """
    frac = np.asarray(frac)
    labels = np.asarray(names, dtype=object)[np.argmax(frac, axis=1)]
    labels[_empty_obs_mask(frac)] = empty_label
    return labels.astype(str)


def _capture_window_from_center(center_xy, size_xy):
    """
    Return (min0, max0, min1, max1) for a 2D window.
    """
    cx, cy = center_xy
    sx, sy = size_xy
    min0 = cx - sx / 2.0
    max0 = cx + sx / 2.0
    min1 = cy - sy / 2.0
    max1 = cy + sy / 2.0
    return float(min0), float(max0), float(min1), float(max1)


def aggregate_molecules_to_grid_bins_2d_slices_window(
    mol_xyz,
    mol_gene,
    mol_src_celltype,
    mol_src_domain,
    sphere_R_um,
    axis="z",
    n_slices=10,
    bin_size_um=8.0,
    # capture window on the 2D plane
    window_center=(0.0, 0.0),
    window_size=(6500.0, 6500.0),
    n_cell_types=4,
    n_domains=4,
    center=(0.0, 0.0, 0.0),
):
    """
    2D bins per slice, within a rectangular capture window (VisiumHD-like).

    ``window_size=None`` disables the capture crop: the bin grid then spans
    exactly the molecule bounding box in the slice plane, so no molecule is
    dropped and there is no empty border around the tissue. (Empty bins
    *inside* the tissue footprint, from gaps between sparsely packed cells,
    are unaffected -- that is not a cropping effect.)
    """
    mol_xyz = np.asarray(mol_xyz)
    mol_gene = np.asarray(mol_gene, dtype=np.int32)
    mol_src_celltype = np.asarray(mol_src_celltype, dtype=np.int32)
    mol_src_domain = np.asarray(mol_src_domain, dtype=np.int32)

    if axis.lower() == "z":
        sd = 2; pd = (0, 1); plane_labels = ("X", "Y")
    elif axis.lower() == "x":
        sd = 0; pd = (1, 2); plane_labels = ("Y", "Z")
    elif axis.lower() == "y":
        sd = 1; pd = (0, 2); plane_labels = ("X", "Z")
    else:
        raise ValueError("axis must be one of {'x','y','z'}")

    # slice assignment by axis coordinate
    sid, edges = slice_id_along_axis(mol_xyz[:, sd] - center[sd], R=sphere_R_um, n_slices=n_slices)
    keep = sid >= 0
    mol_xyz = mol_xyz[keep]
    mol_gene = mol_gene[keep]
    mol_src_celltype = mol_src_celltype[keep]
    mol_src_domain = mol_src_domain[keep]
    sid = sid[keep]

    # apply plane capture window crop; window_size=None -> span the molecule
    # bounding box instead (no crop, no empty border)
    p0 = mol_xyz[:, pd[0]]
    p1 = mol_xyz[:, pd[1]]
    if window_size is None:
        if p0.size:
            min0, max0 = float(p0.min()), float(p0.max()) + 1e-6
            min1, max1 = float(p1.min()), float(p1.max()) + 1e-6
        else:
            min0 = max0 = min1 = max1 = 0.0
        eff_center = ((min0 + max0) / 2.0, (min1 + max1) / 2.0)
        eff_size = (max0 - min0, max1 - min1)
    else:
        min0, max0, min1, max1 = _capture_window_from_center(window_center, window_size)
        inwin = (p0 >= min0) & (p0 <= max0) & (p1 >= min1) & (p1 <= max1)
        mol_xyz = mol_xyz[inwin]
        mol_gene = mol_gene[inwin]
        mol_src_celltype = mol_src_celltype[inwin]
        mol_src_domain = mol_src_domain[inwin]
        sid = sid[inwin]
        p0 = p0[inwin]
        p1 = p1[inwin]
        eff_center = tuple(window_center)
        eff_size = tuple(window_size)

    # define bin grid over window
    nx = int(np.ceil((max0 - min0) / bin_size_um))
    ny = int(np.ceil((max1 - min1) / bin_size_um))

    ix = np.floor((p0 - min0) / bin_size_um).astype(np.int32)
    iy = np.floor((p1 - min1) / bin_size_um).astype(np.int32)

    inside = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    mol_xyz = mol_xyz[inside]
    mol_gene = mol_gene[inside]
    mol_src_celltype = mol_src_celltype[inside]
    mol_src_domain = mol_src_domain[inside]
    sid = sid[inside]
    ix = ix[inside]
    iy = iy[inside]

    bins_per_slice = nx * ny
    obs_id = sid * bins_per_slice + ix * ny + iy
    n_obs = n_slices * bins_per_slice

    n_genes = int(np.max(mol_gene)) + 1 if mol_gene.size else 0

    key = obs_id.astype(np.int64) * np.int64(n_genes) + mol_gene.astype(np.int64)
    uniq, cnt = np.unique(key, return_counts=True)
    rows = (uniq // n_genes).astype(np.int32)
    cols = (uniq % n_genes).astype(np.int32)
    X = csr_matrix((cnt.astype(np.int32), (rows, cols)), shape=(n_obs, n_genes))

    ct_frac = _mix_fraction(obs_id, mol_src_celltype, n_cell_types, n_obs)
    dom_frac = _mix_fraction(obs_id, mol_src_domain, n_domains, n_obs)

    # spatial centers
    c0 = min0 + (np.arange(nx) + 0.5) * bin_size_um
    c1 = min1 + (np.arange(ny) + 0.5) * bin_size_um
    slice_centers = 0.5 * (edges[:-1] + edges[1:]) + center[sd]

    spatial3d = np.zeros((n_obs, 3), dtype=np.float32)
    for s in range(n_slices):
        base = s * bins_per_slice
        gx = np.repeat(np.arange(nx), ny)
        gy = np.tile(np.arange(ny), nx)
        if axis.lower() == "z":
            spatial3d[base:base + bins_per_slice, 0] = c0[gx]
            spatial3d[base:base + bins_per_slice, 1] = c1[gy]
            spatial3d[base:base + bins_per_slice, 2] = slice_centers[s]
        elif axis.lower() == "x":
            spatial3d[base:base + bins_per_slice, 0] = slice_centers[s]
            spatial3d[base:base + bins_per_slice, 1] = c0[gx]
            spatial3d[base:base + bins_per_slice, 2] = c1[gy]
        else:
            spatial3d[base:base + bins_per_slice, 0] = c0[gx]
            spatial3d[base:base + bins_per_slice, 1] = slice_centers[s]
            spatial3d[base:base + bins_per_slice, 2] = c1[gy]

    obs = {
        "slice_axis": np.array([axis.upper()] * n_obs),
        "slice_id": np.repeat(np.arange(n_slices), bins_per_slice).astype(int),
        "bin_size_um": np.array([bin_size_um] * n_obs, dtype=float),
        "window_center0": np.array([eff_center[0]] * n_obs, dtype=float),
        "window_center1": np.array([eff_center[1]] * n_obs, dtype=float),
        "window_size0": np.array([eff_size[0]] * n_obs, dtype=float),
        "window_size1": np.array([eff_size[1]] * n_obs, dtype=float),
        "plane_dim0": np.array([plane_labels[0]] * n_obs),
        "plane_dim1": np.array([plane_labels[1]] * n_obs),
    }
    return X, obs, spatial3d, ct_frac, dom_frac


def aggregate_molecules_to_spots_2d_slices_window(
    mol_xyz,
    mol_gene,
    mol_src_celltype,
    mol_src_domain,
    sphere_R_um,
    axis="z",
    n_slices=10,
    spot_spacing_um=100.0,
    spot_radius_um=27.5,
    # capture window on the 2D plane
    window_center=(0.0, 0.0),
    window_size=(6500.0, 6500.0),
    n_cell_types=4,
    n_domains=4,
    center=(0.0, 0.0, 0.0),
):
    """
    Spot-level per slice within capture window (Visium-like).
    Assign each molecule to nearest spot center, keep if within spot_radius_um.

    ``window_size=None`` disables the capture crop: the spot lattice then spans
    exactly the molecule bounding box in the slice plane, so no molecule is
    dropped for being outside a slide and there is no empty border.
    """
    mol_xyz = np.asarray(mol_xyz)
    mol_gene = np.asarray(mol_gene, dtype=np.int32)
    mol_src_celltype = np.asarray(mol_src_celltype, dtype=np.int32)
    mol_src_domain = np.asarray(mol_src_domain, dtype=np.int32)

    if axis.lower() == "z":
        sd = 2; pd = (0, 1); plane_labels = ("X", "Y")
    elif axis.lower() == "x":
        sd = 0; pd = (1, 2); plane_labels = ("Y", "Z")
    elif axis.lower() == "y":
        sd = 1; pd = (0, 2); plane_labels = ("X", "Z")
    else:
        raise ValueError("axis must be one of {'x','y','z'}")

    sid, edges = slice_id_along_axis(mol_xyz[:, sd] - center[sd], R=sphere_R_um, n_slices=n_slices)
    keep = sid >= 0
    mol_xyz = mol_xyz[keep]
    mol_gene = mol_gene[keep]
    mol_src_celltype = mol_src_celltype[keep]
    mol_src_domain = mol_src_domain[keep]
    sid = sid[keep]

    # capture window; window_size=None -> span the molecule bounding box (no crop)
    p0 = mol_xyz[:, pd[0]]
    p1 = mol_xyz[:, pd[1]]
    if window_size is None:
        if p0.size:
            min0, max0 = float(p0.min()), float(p0.max())
            min1, max1 = float(p1.min()), float(p1.max())
        else:
            min0 = max0 = min1 = max1 = 0.0
        eff_center = ((min0 + max0) / 2.0, (min1 + max1) / 2.0)
        eff_size = (max0 - min0, max1 - min1)
    else:
        min0, max0, min1, max1 = _capture_window_from_center(window_center, window_size)
        inwin = (p0 >= min0) & (p0 <= max0) & (p1 >= min1) & (p1 <= max1)
        mol_xyz = mol_xyz[inwin]
        mol_gene = mol_gene[inwin]
        mol_src_celltype = mol_src_celltype[inwin]
        mol_src_domain = mol_src_domain[inwin]
        sid = sid[inwin]
        p0 = p0[inwin]
        p1 = p1[inwin]
        eff_center = tuple(window_center)
        eff_size = tuple(window_size)

    # spot center grid over window
    grid0 = np.arange(min0, max0 + 1e-6, spot_spacing_um, dtype=np.float32)
    grid1 = np.arange(min1, max1 + 1e-6, spot_spacing_um, dtype=np.float32)
    nx, ny = grid0.size, grid1.size
    n_spots_plane = nx * ny

    ix = np.rint((p0 - grid0[0]) / spot_spacing_um).astype(np.int32)
    iy = np.rint((p1 - grid1[0]) / spot_spacing_um).astype(np.int32)
    inside = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)

    mol_xyz = mol_xyz[inside]
    mol_gene = mol_gene[inside]
    mol_src_celltype = mol_src_celltype[inside]
    mol_src_domain = mol_src_domain[inside]
    sid = sid[inside]
    p0 = p0[inside]
    p1 = p1[inside]
    ix = ix[inside]
    iy = iy[inside]

    cx = grid0[ix]
    cy = grid1[iy]
    dist2 = (p0 - cx) ** 2 + (p1 - cy) ** 2
    ok = dist2 <= (spot_radius_um ** 2)

    mol_xyz = mol_xyz[ok]
    mol_gene = mol_gene[ok]
    mol_src_celltype = mol_src_celltype[ok]
    mol_src_domain = mol_src_domain[ok]
    sid = sid[ok]
    ix = ix[ok]
    iy = iy[ok]

    spot_id_plane = ix * ny + iy
    obs_id = sid * n_spots_plane + spot_id_plane
    n_obs = n_slices * n_spots_plane

    n_genes = int(np.max(mol_gene)) + 1 if mol_gene.size else 0

    key = obs_id.astype(np.int64) * np.int64(n_genes) + mol_gene.astype(np.int64)
    uniq, cnt = np.unique(key, return_counts=True)
    rows = (uniq // n_genes).astype(np.int32)
    cols = (uniq % n_genes).astype(np.int32)
    X = csr_matrix((cnt.astype(np.int32), (rows, cols)), shape=(n_obs, n_genes))

    ct_frac = _mix_fraction(obs_id, mol_src_celltype, n_cell_types, n_obs)
    dom_frac = _mix_fraction(obs_id, mol_src_domain, n_domains, n_obs)

    slice_centers = 0.5 * (edges[:-1] + edges[1:]) + center[sd]

    spatial3d = np.zeros((n_obs, 3), dtype=np.float32)
    for s in range(n_slices):
        base = s * n_spots_plane
        gx = np.repeat(np.arange(nx), ny)
        gy = np.tile(np.arange(ny), nx)
        if axis.lower() == "z":
            spatial3d[base:base + n_spots_plane, 0] = grid0[gx]
            spatial3d[base:base + n_spots_plane, 1] = grid1[gy]
            spatial3d[base:base + n_spots_plane, 2] = slice_centers[s]
        elif axis.lower() == "x":
            spatial3d[base:base + n_spots_plane, 0] = slice_centers[s]
            spatial3d[base:base + n_spots_plane, 1] = grid0[gx]
            spatial3d[base:base + n_spots_plane, 2] = grid1[gy]
        else:
            spatial3d[base:base + n_spots_plane, 0] = grid0[gx]
            spatial3d[base:base + n_spots_plane, 1] = slice_centers[s]
            spatial3d[base:base + n_spots_plane, 2] = grid1[gy]

    obs = {
        "slice_axis": np.array([axis.upper()] * n_obs),
        "slice_id": np.repeat(np.arange(n_slices), n_spots_plane).astype(int),
        "spot_spacing_um": np.array([spot_spacing_um] * n_obs, dtype=float),
        "spot_radius_um": np.array([spot_radius_um] * n_obs, dtype=float),
        "window_center0": np.array([eff_center[0]] * n_obs, dtype=float),
        "window_center1": np.array([eff_center[1]] * n_obs, dtype=float),
        "window_size0": np.array([eff_size[0]] * n_obs, dtype=float),
        "window_size1": np.array([eff_size[1]] * n_obs, dtype=float),
        "plane_dim0": np.array([plane_labels[0]] * n_obs),
        "plane_dim1": np.array([plane_labels[1]] * n_obs),
    }
    return X, obs, spatial3d, ct_frac, dom_frac


# -------------------------
#   Slice ids for cell centroids
# -------------------------

def add_slice_ids_to_cells(adata_cell, sphere_R_um, n_slices=10, center=(0.0, 0.0, 0.0)):
    xyz = np.asarray(adata_cell.obsm["spatial"])
    sx, _ = slice_id_along_axis(xyz[:, 0] - center[0], R=sphere_R_um, n_slices=n_slices)
    sy, _ = slice_id_along_axis(xyz[:, 1] - center[1], R=sphere_R_um, n_slices=n_slices)
    sz, _ = slice_id_along_axis(xyz[:, 2] - center[2], R=sphere_R_um, n_slices=n_slices)
    adata_cell.obs["slice_id_X"] = sx
    adata_cell.obs["slice_id_Y"] = sy
    adata_cell.obs["slice_id_Z"] = sz


def make_cell_sectioned_with_batch(adata_cell_obs, axis="Z", n_slices=10, batch_sigma=0.15, seed=0):
    axis = axis.upper()
    key = f"slice_id_{axis}"
    assert key in adata_cell_obs.obs.columns, f"Missing {key} in cell obs."

    sec = adata_cell_obs.copy()
    sec.obs["slice_axis"] = axis
    sec.obs["slice_id"] = sec.obs[key].astype(int)
    sec.layers["counts_pre_batch"] = sec.X.copy()

    X_be, factors = apply_batch_effect_per_slice_sparse(
        sec.X.tocsr(), slice_ids=sec.obs["slice_id"].values, n_slices=n_slices,
        batch_sigma=batch_sigma, seed=seed
    )
    sec.X = X_be
    sec.uns["batch_effect_factors"] = _to_serializable({str(int(k)): v.tolist() for k, v in factors.items()})
    return sec


#   Main simulator

def simulate_3d_molecule_sphere_multires(
    # tissue
    sphere_R_um=6000.0,               # 12 mm diameter default
    center=(0.0, 0.0, 0.0),

    # capture windows -- crop cell (xenium) and bin/spot (visium) output to the
    # slide area; pass False for a window to disable that crop
    xenium_capture_window_um=(12000.0, 24000.0),   # 12×24 mm, crops cells
    capture_window_um=(6500.0, 6500.0),            # 6.5×6.5 mm, crops bins/spots
    capture_window_center_um=(0.0, 0.0),           # can shift to include domains

    # domains
    n_domains=4,
    core_frac=0.35,
    # domain irregularity knobs
    core_bump_amp=0.12,
    wedge_angle_amp_deg=12.0,
    noise_terms=6,
    noise_freq_range=(0.8, 2.2),
    boundary_fuzz_width_deg=0.0,
    boundary_fuzz_flip_prob=0.0,
    core_fuzz_width_um=0.0,
    core_fuzz_flip_prob=0.0,

    # cells
    n_cells=20000,
    allow_cell_overlap=False,
    cell_radius_kwargs=None,

    # cell types
    n_cell_types=4,
    domain_type_mix=None,  # (n_domains, n_cell_types)

    # genes
    marker_genes_per_type=80,
    noise_gene_frac=0.10,
    shared_marker_frac=0.25,

    # expression model (NB)
    base_gene_lognormal=(0.7, 0.7),
    marker_foldchange=3.5,
    shared_marker_foldchange=2.5,
    noise_scale=0.9,
    cell_size_lognormal=(0.0, 0.35),
    domain_size_factors=None,
    theta=25.0,
    theta_jitter=2.0,

    # molecule placement
    inside_prob=0.95,
    assign_k=8,

    # slicing + batch effect
    n_slices=10,
    batch_sigma=0.15,

    # bin-level (VisiumHD-like)
    bin_size_um=8.0,

    # spot-level (Visium-like)
    spot_spacing_um=100.0,
    spot_radius_um=27.5,

    # coordinate perturbation for unaligned
    max_deg=180.0,
    max_shift=200.0,
    base_seed_unaligned=12345,
    sync_unaligned_seed=False,

    # general
    seed=2025,
    sparse_X=True,

    # optional output selection
    output_modalities=None,
    slice_axes=None,
    _section=True,
):
    if _section:
        base = simulate_3d_molecule_sphere_base(
            sphere_R_um=sphere_R_um,
            center=center,
            xenium_capture_window_um=xenium_capture_window_um,
            capture_window_um=capture_window_um,
            capture_window_center_um=capture_window_center_um,
            n_domains=n_domains,
            core_frac=core_frac,
            core_bump_amp=core_bump_amp,
            wedge_angle_amp_deg=wedge_angle_amp_deg,
            noise_terms=noise_terms,
            noise_freq_range=noise_freq_range,
            boundary_fuzz_width_deg=boundary_fuzz_width_deg,
            boundary_fuzz_flip_prob=boundary_fuzz_flip_prob,
            core_fuzz_width_um=core_fuzz_width_um,
            core_fuzz_flip_prob=core_fuzz_flip_prob,
            n_cells=n_cells,
            allow_cell_overlap=allow_cell_overlap,
            cell_radius_kwargs=cell_radius_kwargs,
            n_cell_types=n_cell_types,
            domain_type_mix=domain_type_mix,
            marker_genes_per_type=marker_genes_per_type,
            noise_gene_frac=noise_gene_frac,
            shared_marker_frac=shared_marker_frac,
            base_gene_lognormal=base_gene_lognormal,
            marker_foldchange=marker_foldchange,
            shared_marker_foldchange=shared_marker_foldchange,
            noise_scale=noise_scale,
            cell_size_lognormal=cell_size_lognormal,
            domain_size_factors=domain_size_factors,
            theta=theta,
            theta_jitter=theta_jitter,
            inside_prob=inside_prob,
            assign_k=assign_k,
            seed=seed,
            sparse_X=sparse_X,
            output_modalities=output_modalities,
        )
        return section_3d_molecule_sphere(
            base,
            n_slices=n_slices,
            batch_sigma=batch_sigma,
            capture_window_um=capture_window_um,
            capture_window_center_um=capture_window_center_um,
            xenium_capture_window_um=xenium_capture_window_um,
            bin_size_um=bin_size_um,
            spot_spacing_um=spot_spacing_um,
            spot_radius_um=spot_radius_um,
            max_deg=max_deg,
            max_shift=max_shift,
            base_seed_unaligned=base_seed_unaligned,
            sync_unaligned_seed=sync_unaligned_seed,
            output_modalities=output_modalities,
            slice_axes=slice_axes,
        )

    rng = np.random.default_rng(seed)

    valid_modalities = {"cell", "bin", "spot"}
    if output_modalities is None:
        output_modalities = valid_modalities
    elif isinstance(output_modalities, str):
        output_modalities = {output_modalities.lower()}
    else:
        output_modalities = {str(modality).lower() for modality in output_modalities}
    invalid_modalities = output_modalities - valid_modalities
    if invalid_modalities:
        names = ", ".join(sorted(invalid_modalities))
        raise ValueError(f"Unsupported output modalities: {names}.")

    if slice_axes is None:
        slice_axes = ("X", "Y", "Z")
    elif isinstance(slice_axes, str):
        slice_axes = (slice_axes,)
    output_axes = tuple(dict.fromkeys(str(axis).upper() for axis in slice_axes))
    if not output_axes or any(axis not in {"X", "Y", "Z"} for axis in output_axes):
        raise ValueError("slice_axes must contain one or more values from 'X', 'Y', 'Z'.")

    generate_cells = "cell" in output_modalities
    generate_bins = "bin" in output_modalities
    generate_spots = "spot" in output_modalities

    # ---- domain->type mix ----
    if domain_type_mix is None:
        dm = np.array([
            [0.70, 0.20, 0.07, 0.03],
            [0.10, 0.70, 0.15, 0.05],
            [0.08, 0.12, 0.70, 0.10],
            [0.25, 0.25, 0.25, 0.25],
        ], dtype=float)
        if (n_domains, n_cell_types) != dm.shape:
            warnings.warn(
                f"domain_type_mix was not provided and (n_domains={n_domains}, "
                f"n_cell_types={n_cell_types}) does not match the built-in "
                f"{dm.shape} example composition, so every domain will be assigned "
                "the exact same uniform cell-type distribution (no domain-specific "
                "enrichment at all). Pass an explicit domain_type_mix of shape "
                f"({n_domains}, {n_cell_types}) if domains should differ in "
                "cell-type composition.",
                stacklevel=2,
            )
            dm = np.tile(1.0 / n_cell_types, (n_domains, n_cell_types))
        domain_type_mix = dm
    else:
        domain_type_mix = np.asarray(domain_type_mix, dtype=float)
        assert domain_type_mix.shape == (n_domains, n_cell_types)
        domain_type_mix = domain_type_mix / domain_type_mix.sum(axis=1, keepdims=True)

    # ---- cells: centroids + radii ----
    centers, radii, tree = sample_nonoverlapping_cells_in_sphere(
        n_cells=n_cells,
        sphere_R=sphere_R_um,
        center=center,
        radius_kwargs=cell_radius_kwargs,
        allow_overlap=allow_cell_overlap,
        seed=seed,
    )

    # ---- domains/types per cell (NOISY domains) ----
    domain_ids = sphere_domain_labels_noisy(
        centers,
        center=center,
        n_domains=n_domains,
        core_frac=core_frac,
        core_bump_amp=core_bump_amp,
        wedge_angle_amp_deg=wedge_angle_amp_deg,
        noise_terms=noise_terms,
        noise_freq_range=noise_freq_range,
        seed=seed + 100,
        boundary_fuzz_width_deg=boundary_fuzz_width_deg,
        boundary_fuzz_flip_prob=boundary_fuzz_flip_prob,
        core_fuzz_width_um=core_fuzz_width_um,
        core_fuzz_flip_prob=core_fuzz_flip_prob
    )
    cell_type_ids = sample_cell_types_by_domain(domain_ids, domain_type_mix, rng=rng)

    # ---- genes ----
    gene_names, marker_sets, noise_idx, gene_class, gene_info = make_marker_and_noise_genes(
        n_cell_types=n_cell_types,
        marker_genes_per_type=marker_genes_per_type,
        noise_gene_frac=noise_gene_frac,
        shared_marker_frac=shared_marker_frac,
        seed=seed,
    )
    n_genes = len(gene_names)

    # ---- true counts ----
    X_cell_true, base_gene, size_factor, theta_g, multipliers, shared_markers = simulate_cell_gene_counts_nb(
        cell_type_ids=cell_type_ids,
        domain_ids=domain_ids,
        marker_sets=marker_sets,
        noise_idx=noise_idx,
        base_gene_lognormal=base_gene_lognormal,
        marker_foldchange=marker_foldchange,
        shared_marker_foldchange=shared_marker_foldchange,
        noise_scale=noise_scale,
        cell_size_lognormal=cell_size_lognormal,
        domain_size_factors=domain_size_factors,
        theta=theta,
        theta_jitter=theta_jitter,
        seed=seed + 1,
        sparse_X=sparse_X,
    )
    X_true_csr = X_cell_true.tocsr()

    # ---- molecules ----
    #
    # We keep TWO streams:
    #   (A) cell-level stream: only molecules assigned to some cell (assigned >= 0)
    #       -> used for X_cell_obs (spillover into OTHER cells is allowed by containment)
    #   (B) full stream: ALL molecules (assigned may be -1)
    #       -> used for bins/spots (so ambient/outside-all-cells molecules still contribute)
    #
    mol_xyz_all = []
    mol_gene_all = []
    mol_assigned_cell_all = []

    mol_xyz_all_full = []
    mol_gene_all_full = []
    mol_src_celltype_all_full = []
    mol_src_domain_all_full = []

    sigmas = np.zeros(n_cells, dtype=np.float32)

    n_mols_generated_total = 0
    n_mols_assigned_total = 0
    n_mols_unassigned_total = 0

    for c in range(n_cells):
        start, end = X_true_csr.indptr[c], X_true_csr.indptr[c + 1]
        if start == end:
            continue

        cols = X_true_csr.indices[start:end]
        cnts = X_true_csr.data[start:end].astype(np.int32)
        n_mols = int(cnts.sum())
        if n_mols <= 0:
            continue

        gene_ids = np.repeat(cols, cnts)
        xyz, sigma = sample_molecule_coords_for_cell(
            n_mols=n_mols,
            centroid=centers[c],
            radius=radii[c],
            inside_prob=inside_prob,
            dim=3,
            rng=rng,
        )
        sigmas[c] = float(sigma)

        n_mols_generated_total += xyz.shape[0]

        # ---- (B) FULL stream: keep everything for bins/spots ----
        if generate_bins or generate_spots:
            mol_xyz_all_full.append(xyz.astype(np.float32))
            mol_gene_all_full.append(gene_ids.astype(np.int32))
            mol_src_celltype_all_full.append(np.full(xyz.shape[0], cell_type_ids[c], dtype=np.int32))
            mol_src_domain_all_full.append(np.full(xyz.shape[0], domain_ids[c], dtype=np.int32))

        # ---- (A) CELL stream: only those that land in some cell ----
        if not generate_cells:
            continue

        assigned = assign_points_to_cells_by_containment(
            xyz, cell_centers=centers, cell_radii=radii, tree=tree, k=assign_k
        )
        n_mols_unassigned_total += int(np.sum(assigned < 0))
        keep = assigned >= 0
        if not np.any(keep):
            continue

        xyz_keep = xyz[keep]
        gene_keep = gene_ids[keep]
        assigned_keep = assigned[keep]

        n_mols_assigned_total += xyz_keep.shape[0]

        mol_xyz_all.append(xyz_keep.astype(np.float32))
        mol_gene_all.append(gene_keep.astype(np.int32))
        mol_assigned_cell_all.append(assigned_keep.astype(np.int32))

    # Build global arrays
    # (A) assigned-to-a-cell molecules (for cell-level observed)
    if mol_xyz_all:
        mol_xyz = np.vstack(mol_xyz_all)
        mol_gene = np.concatenate(mol_gene_all)
        mol_assigned_cell = np.concatenate(mol_assigned_cell_all)
    else:
        mol_xyz = np.zeros((0, 3), dtype=np.float32)
        mol_gene = np.zeros((0,), dtype=np.int32)
        mol_assigned_cell = np.zeros((0,), dtype=np.int32)

    # (B) full molecules (for bins/spots)
    if mol_xyz_all_full:
        mol_xyz_full = np.vstack(mol_xyz_all_full)
        mol_gene_full = np.concatenate(mol_gene_all_full)
        mol_src_celltype_full = np.concatenate(mol_src_celltype_all_full)
        mol_src_domain_full = np.concatenate(mol_src_domain_all_full)
    else:
        mol_xyz_full = np.zeros((0, 3), dtype=np.float32)
        mol_gene_full = np.zeros((0,), dtype=np.int32)
        mol_src_celltype_full = np.zeros((0,), dtype=np.int32)
        mol_src_domain_full = np.zeros((0,), dtype=np.int32)

    # ---- observed cell-level counts (spillover allowed by observed containment) ----
    if generate_cells:
        X_cell_obs = aggregate_molecules_to_cell_level(
            gene_ids=mol_gene, assigned_cell_ids=mol_assigned_cell, n_cells=n_cells, n_genes=n_genes
        )
    else:
        X_cell_obs = X_cell_true.copy()

    # ---- AnnData: cell true + cell observed ----
    dom_names = [f"D{d}" for d in range(n_domains)]
    ct_names = [f"type{t+1}" for t in range(n_cell_types)]

    cell_obs = {
        "domain_true": np.array([dom_names[d] for d in domain_ids]),
        "cell_type_true": np.array([ct_names[t] for t in cell_type_ids]),
        "cell_radius": radii.astype(np.float32),
        "cell_sigma_for_95pct": sigmas.astype(np.float32),
    }

    var = {
        "gene_class": gene_class,
        "is_noise": np.array([g in noise_idx for g in range(n_genes)], dtype=bool),
        "is_shared_marker": np.array([g in set(shared_markers) for g in range(n_genes)], dtype=bool),
    }

    adata_cell_true = ad.AnnData(X=X_cell_true, obs=cell_obs, var=var)
    adata_cell_true.var_names = gene_names
    adata_cell_true.obsm["spatial"] = centers.astype(np.float32)  # aligned 3D for base cell object

    adata_cell_obs = ad.AnnData(X=X_cell_obs, obs=cell_obs, var=var)
    adata_cell_obs.var_names = gene_names
    adata_cell_obs.obsm["spatial"] = centers.astype(np.float32)  # aligned 3D for base cell object

    if not _section:
        meta = dict(
            units="microns",
            sphere_R_um=float(sphere_R_um),
            center=center,
            n_domains=n_domains,
            core_frac=float(core_frac),
            domain_irregularity=dict(
                core_bump_amp=float(core_bump_amp),
                wedge_angle_amp_deg=float(wedge_angle_amp_deg),
                noise_terms=int(noise_terms),
                noise_freq_range=noise_freq_range,
                boundary_fuzz_width_deg=float(boundary_fuzz_width_deg),
                boundary_fuzz_flip_prob=float(boundary_fuzz_flip_prob),
                core_fuzz_width_um=float(core_fuzz_width_um),
                core_fuzz_flip_prob=float(core_fuzz_flip_prob),
            ),
            captures=dict(
                xenium_capture_window_um=xenium_capture_window_um,
                capture_window_um=capture_window_um,
                capture_window_center_um=capture_window_center_um,
            ),
            n_cells=int(n_cells),
            n_cell_types=int(n_cell_types),
            domain_type_mix=domain_type_mix,
            gene_info=gene_info,
            inside_prob=float(inside_prob),
            assign_k=int(assign_k),
            seed=int(seed),
            output_modalities=sorted(output_modalities),
            n_molecules_generated_total=int(n_mols_generated_total),
            n_molecules_assigned_total=(int(n_mols_assigned_total) if generate_cells else None),
            n_molecules_unassigned_total=(int(n_mols_unassigned_total) if generate_cells else None),
            n_molecules_total=(int(n_mols_assigned_total) if generate_cells else None),
        )
        adata_cell_true.uns["sim_params"] = _to_serializable(meta)
        adata_cell_obs.uns["sim_params"] = _to_serializable(meta)
        return dict(
            adata_cell_true=adata_cell_true,
            adata_cell_obs=adata_cell_obs,
            molecules=dict(
                assigned_xyz=mol_xyz,
                assigned_gene=mol_gene,
                assigned_cell=mol_assigned_cell,
                full_xyz=mol_xyz_full,
                full_gene=mol_gene_full,
                full_src_celltype=mol_src_celltype_full,
                full_src_domain=mol_src_domain_full,
            ),
            meta=meta,
        )

    # slice ids for cells along X/Y/Z
    add_slice_ids_to_cells(adata_cell_obs, sphere_R_um=sphere_R_um, n_slices=n_slices, center=center)
    add_slice_ids_to_cells(adata_cell_true, sphere_R_um=sphere_R_um, n_slices=n_slices, center=center)

    # per-modality seed offset for the unaligned-coordinate perturbation -- see
    # sync_unaligned_seed docstring note in section_3d_molecule_sphere.
    modality_offset = {"cell": 0, "bin": 0, "spot": 0} if sync_unaligned_seed else \
        {"cell": 10_000, "bin": 20_000, "spot": 30_000}

    # ---- cell-level: sectioned versions with batch effects ----
    cell_sectioned = {}
    for AX in (output_axes if generate_cells else ()):
        sec = make_cell_sectioned_with_batch(
            adata_cell_obs, axis=AX, n_slices=n_slices, batch_sigma=batch_sigma, seed=seed + 1000 + ord(AX)
        )
        add_spatial_keys_for_axis(
            sec,
            axis_letter=AX,
            slice_key="slice_id",
            unaligned_key="spatial_unaligned",
            base_seed=base_seed_unaligned + modality_offset["cell"] + ord(AX),
            max_deg=max_deg,
            max_shift=max_shift,
        )
        cell_sectioned[AX] = sec

    # ---- bin/spot AnnData per axis (cropped to Visium window) ----
    # IMPORTANT CHANGE:
    #   bins/spots use the FULL molecule stream (mol_xyz_full, mol_gene_full, ...)
    bin_adatas = {}
    spot_adatas = {}

    for ax in (tuple(axis.lower() for axis in output_axes) if generate_bins else ()):
        # bins (VisiumHD) within 6.5×6.5mm window
        Xb, obs_b, spatial3d_b, ct_frac_b, dom_frac_b = aggregate_molecules_to_grid_bins_2d_slices_window(
            mol_xyz=mol_xyz_full,
            mol_gene=mol_gene_full,
            mol_src_celltype=mol_src_celltype_full,
            mol_src_domain=mol_src_domain_full,
            sphere_R_um=sphere_R_um,
            axis=ax,
            n_slices=n_slices,
            bin_size_um=bin_size_um,
            window_center=capture_window_center_um,
            window_size=capture_window_um,
            n_cell_types=n_cell_types,
            n_domains=n_domains,
            center=center,
        )
        adb = ad.AnnData(X=Xb, obs=obs_b, var=var)
        adb.var_names = gene_names
        adb.layers["counts_pre_batch"] = adb.X.copy()
        adb.obsm["spatial"] = spatial3d_b  # temporarily 3D, will be re-mapped below
        adb.obsm["cell_type_frac_true"] = ct_frac_b
        adb.obsm["domain_frac_true"] = dom_frac_b
        adb.obs["cell_type_true"] = _labels_from_fracs(ct_frac_b, ct_names)
        adb.obs["domain_true"] = _labels_from_fracs(dom_frac_b, dom_names)
        adb.obs["is_empty"] = _empty_obs_mask(dom_frac_b)

        # batch effects for bins
        Xb_be, factors = apply_batch_effect_per_slice_sparse(
            adb.X.tocsr(),
            slice_ids=adb.obs["slice_id"].values,
            n_slices=n_slices,
            batch_sigma=batch_sigma,
            seed=seed + 2000 + ord(ax),
        )
        adb.X = Xb_be
        adb.uns["batch_effect_factors"] = _to_serializable({str(int(k)): v.tolist() for k, v in factors.items()})

        add_spatial_keys_for_axis(
            adb,
            axis_letter=ax.upper(),
            slice_key="slice_id",
            unaligned_key="spatial_unaligned",
            base_seed=base_seed_unaligned + modality_offset["bin"] + ord(ax.upper() if sync_unaligned_seed else ax),
            max_deg=max_deg,
            max_shift=max_shift,
        )
        bin_adatas[ax.upper()] = adb

    for ax in (tuple(axis.lower() for axis in output_axes) if generate_spots else ()):
        # spots (Visium) within 6.5×6.5mm window
        Xs, obs_s, spatial3d_s, ct_frac_s, dom_frac_s = aggregate_molecules_to_spots_2d_slices_window(
            mol_xyz=mol_xyz_full,
            mol_gene=mol_gene_full,
            mol_src_celltype=mol_src_celltype_full,
            mol_src_domain=mol_src_domain_full,
            sphere_R_um=sphere_R_um,
            axis=ax,
            n_slices=n_slices,
            spot_spacing_um=spot_spacing_um,
            spot_radius_um=spot_radius_um,
            window_center=capture_window_center_um,
            window_size=capture_window_um,
            n_cell_types=n_cell_types,
            n_domains=n_domains,
            center=center,
        )
        ads = ad.AnnData(X=Xs, obs=obs_s, var=var)
        ads.var_names = gene_names
        ads.layers["counts_pre_batch"] = ads.X.copy()
        ads.obsm["spatial"] = spatial3d_s  # temporarily 3D, will be re-mapped below
        ads.obsm["cell_type_frac_true"] = ct_frac_s
        ads.obsm["domain_frac_true"] = dom_frac_s
        ads.obs["cell_type_true"] = _labels_from_fracs(ct_frac_s, ct_names)
        ads.obs["domain_true"] = _labels_from_fracs(dom_frac_s, dom_names)
        ads.obs["is_empty"] = _empty_obs_mask(dom_frac_s)

        # batch effects for spots
        Xs_be, factors_s = apply_batch_effect_per_slice_sparse(
            ads.X.tocsr(),
            slice_ids=ads.obs["slice_id"].values,
            n_slices=n_slices,
            batch_sigma=batch_sigma,
            seed=seed + 3000 + ord(ax),
        )
        ads.X = Xs_be
        ads.uns["batch_effect_factors"] = _to_serializable({str(int(k)): v.tolist() for k, v in factors_s.items()})

        add_spatial_keys_for_axis(
            ads,
            axis_letter=ax.upper(),
            slice_key="slice_id",
            unaligned_key="spatial_unaligned",
            base_seed=base_seed_unaligned + modality_offset["spot"] + ord(ax.upper() if sync_unaligned_seed else ax),
            max_deg=max_deg,
            max_shift=max_shift,
        )
        spot_adatas[ax.upper()] = ads

    # ---- meta ----
    meta = dict(
        units="microns",
        sphere_R_um=float(sphere_R_um),
        center=center,
        n_domains=n_domains,
        core_frac=float(core_frac),
        domain_irregularity=dict(
            core_bump_amp=float(core_bump_amp),
            wedge_angle_amp_deg=float(wedge_angle_amp_deg),
            noise_terms=int(noise_terms),
            noise_freq_range=noise_freq_range,
            boundary_fuzz_width_deg=float(boundary_fuzz_width_deg),
            boundary_fuzz_flip_prob=float(boundary_fuzz_flip_prob),
            core_fuzz_width_um=float(core_fuzz_width_um),
            core_fuzz_flip_prob=float(core_fuzz_flip_prob),
        ),
        captures=dict(
            xenium_capture_window_um=xenium_capture_window_um,
            capture_window_um=capture_window_um,
            capture_window_center_um=capture_window_center_um,
        ),
        n_cells=int(n_cells),
        n_cell_types=int(n_cell_types),
        domain_type_mix=domain_type_mix,
        gene_info=gene_info,
        inside_prob=float(inside_prob),
        assign_k=int(assign_k),
        n_slices=int(n_slices),
        batch_sigma=float(batch_sigma),
        bin_size_um=float(bin_size_um),
        spot_spacing_um=float(spot_spacing_um),
        spot_radius_um=float(spot_radius_um),
        unaligned_xy=dict(
            max_deg=float(max_deg),
            max_shift=float(max_shift),
            base_seed_unaligned=int(base_seed_unaligned),
        ),
        seed=int(seed),

        # molecule accounting
        n_molecules_generated_total=int(n_mols_generated_total),
        n_molecules_assigned_total=(int(n_mols_assigned_total) if generate_cells else None),
        n_molecules_unassigned_total=(int(n_mols_unassigned_total) if generate_cells else None),

        # Keep backward-compatible key name
        n_molecules_total=(int(n_mols_assigned_total) if generate_cells else None),
    )

    adata_cell_true.uns["sim_params"] = _to_serializable(meta)
    adata_cell_obs.uns["sim_params"] = _to_serializable(meta)

    out = dict(
        adata_cell_true=adata_cell_true,
        adata_cell_obs=adata_cell_obs,
        adata_cell_sectioned=cell_sectioned,  # X/Y/Z sectioned (with batch) + aligned/unaligned 2D
        bin_adatas=bin_adatas,                # X/Y/Z + capture crop + aligned/unaligned 2D
        spot_adatas=spot_adatas,              # X/Y/Z + capture crop + aligned/unaligned 2D
        meta=meta,
    )
    return out



def simulate_3d_molecule_sphere_base(**kwargs):
    """Generate the intact 3D sphere before slicing, capture, or batch effects.

    The returned dictionary contains `adata_cell_true`, `adata_cell_obs`, `meta`,
    and a `molecules` dictionary with the transcript coordinates needed for
    later bin/spot sectioning.
    """
    kwargs["_section"] = False
    return simulate_3d_molecule_sphere_multires(**kwargs)


def section_3d_molecule_sphere(
    base,
    n_slices=10,
    batch_sigma=0.15,
    capture_window_um=None,
    capture_window_center_um=None,
    xenium_capture_window_um=None,
    bin_size_um=8.0,
    spot_spacing_um=100.0,
    spot_radius_um=27.5,
    max_deg=180.0,
    max_shift=200.0,
    base_seed_unaligned=12345,
    sync_unaligned_seed=False,
    output_modalities=None,
    slice_axes=None,
):
    """Slice the intact sphere and apply per-slice capture + batch effects.

    Capture windows (all in the slice plane, centred on
    ``capture_window_center_um``):

    * ``capture_window_um`` -- bin/spot crop. ``None`` inherits the value
      recorded by step 1; ``False`` disables the crop (grid spans the molecule
      bounding box, no empty border); a ``(w, h)`` pair is used as-is.
    * ``xenium_capture_window_um`` -- cell crop. Same ``None`` / ``False`` /
      ``(w, h)`` semantics. Cells whose in-plane centroid falls outside the
      window are dropped, matching a real Xenium capture area.

    NOTE: when the tissue is smaller than the capture window, bin/spot outputs
    contain all-zero observations around the tissue. They are flagged with
    ``obs["is_empty"]`` and labelled ``domain_true``/``cell_type_true`` =
    ``"unassigned"``; filter them downstream with ``adata[~adata.obs["is_empty"]]``
    or a per-observation min-count QC.

    ``sync_unaligned_seed``: the per-slice unaligned-coordinate perturbation
    (rotation + translation) is seeded as ``base_seed_unaligned + modality_offset
    + ord(axis) + slice_id``, where ``modality_offset`` is normally a different
    fixed constant per modality (10_000 cell / 20_000 bin / 30_000 spot) so that
    a single call requesting multiple modalities at once doesn't give them all
    the identical perturbation. Pass ``sync_unaligned_seed=True`` to use
    ``modality_offset=0`` for all three instead, so cell/bin/spot draw the SAME
    per-slice perturbation when generated with matching ``base_seed_unaligned``/
    ``max_deg``/``max_shift`` (even across separate single-modality calls, e.g.
    one ``generate_simulation_noisy.py --modality X`` run per modality) --
    useful for isolating how each modality's own alignment method resolves an
    identical starting misalignment. Default False preserves prior behavior.
    """
    modality_offset = {"cell": 0, "bin": 0, "spot": 0} if sync_unaligned_seed else \
        {"cell": 10_000, "bin": 20_000, "spot": 30_000}
    valid_modalities = {"cell", "bin", "spot"}
    if output_modalities is None:
        output_modalities = valid_modalities
    elif isinstance(output_modalities, str):
        output_modalities = {output_modalities.lower()}
    else:
        output_modalities = {str(modality).lower() for modality in output_modalities}
    invalid_modalities = output_modalities - valid_modalities
    if invalid_modalities:
        names = ", ".join(sorted(invalid_modalities))
        raise ValueError(f"Unsupported output modalities: {names}.")

    if slice_axes is None:
        slice_axes = ("X", "Y", "Z")
    elif isinstance(slice_axes, str):
        slice_axes = (slice_axes,)
    output_axes = tuple(dict.fromkeys(str(axis).upper() for axis in slice_axes))
    if not output_axes or any(axis not in {"X", "Y", "Z"} for axis in output_axes):
        raise ValueError("slice_axes must contain one or more values from 'X', 'Y', 'Z'.")

    generate_cells = "cell" in output_modalities
    generate_bins = "bin" in output_modalities
    generate_spots = "spot" in output_modalities

    meta = dict(base["meta"])
    captures = dict(meta.get("captures", {}))
    sphere_R_um = float(meta["sphere_R_um"])
    center = tuple(meta["center"])
    n_domains = int(meta["n_domains"])
    n_cell_types = int(meta["n_cell_types"])
    seed = int(meta["seed"])

    if capture_window_um is None:
        inherited = captures.get("capture_window_um", (6500.0, 6500.0))
        capture_window_um = None if inherited in (None, False) else tuple(inherited)
    elif capture_window_um is False:
        capture_window_um = None  # aggregators: None -> bounding-box grid, no crop
    if capture_window_center_um is None:
        capture_window_center_um = tuple(captures.get("capture_window_center_um", (0.0, 0.0)))

    if xenium_capture_window_um is None:
        inherited = captures.get("xenium_capture_window_um", (12000.0, 24000.0))
        if inherited is False:
            xenium_capture_window_um = False
        elif inherited is None:
            xenium_capture_window_um = (12000.0, 24000.0)
        else:
            xenium_capture_window_um = tuple(inherited)
    elif xenium_capture_window_um is not False:
        xenium_capture_window_um = tuple(xenium_capture_window_um)

    adata_cell_true = base["adata_cell_true"].copy()
    adata_cell_obs = base["adata_cell_obs"].copy()
    add_slice_ids_to_cells(adata_cell_obs, sphere_R_um=sphere_R_um, n_slices=n_slices, center=center)
    add_slice_ids_to_cells(adata_cell_true, sphere_R_um=sphere_R_um, n_slices=n_slices, center=center)

    gene_names = list(adata_cell_true.var_names)
    var = adata_cell_true.var.copy()
    dom_names = [f"D{d}" for d in range(n_domains)]
    ct_names = [f"type{t+1}" for t in range(n_cell_types)]
    molecules = base["molecules"]

    cell_sectioned = {}
    for AX in (output_axes if generate_cells else ()):
        sec = make_cell_sectioned_with_batch(
            adata_cell_obs,
            axis=AX,
            n_slices=n_slices,
            batch_sigma=batch_sigma,
            seed=seed + 1000 + ord(AX),
        )
        add_spatial_keys_for_axis(
            sec,
            axis_letter=AX,
            slice_key="slice_id",
            unaligned_key="spatial_unaligned",
            base_seed=base_seed_unaligned + modality_offset["cell"] + ord(AX),
            max_deg=max_deg,
            max_shift=max_shift,
        )
        if xenium_capture_window_um is not False:
            # crop cells to the (Xenium) capture area, on the aligned in-plane
            # coords -- same plane-dim order as obsm['spatial']
            xmin0, xmax0, xmin1, xmax1 = _capture_window_from_center(
                capture_window_center_um, xenium_capture_window_um
            )
            xy = np.asarray(sec.obsm["spatial"])
            in_win = (
                (xy[:, 0] >= xmin0) & (xy[:, 0] <= xmax0)
                & (xy[:, 1] >= xmin1) & (xy[:, 1] <= xmax1)
            )
            if not in_win.all():
                sec = sec[in_win].copy()
        cell_sectioned[AX] = sec

    if (generate_bins or generate_spots) and molecules["full_gene"].size == 0:
        raise ValueError(
            "This base sphere does not include the full molecule stream needed "
            "for bin/spot outputs. Generate it with output_modalities including "
            "'bin' or 'spot'."
        )

    bin_adatas = {}
    for ax in (tuple(axis.lower() for axis in output_axes) if generate_bins else ()):
        Xb, obs_b, spatial3d_b, ct_frac_b, dom_frac_b = aggregate_molecules_to_grid_bins_2d_slices_window(
            mol_xyz=molecules["full_xyz"],
            mol_gene=molecules["full_gene"],
            mol_src_celltype=molecules["full_src_celltype"],
            mol_src_domain=molecules["full_src_domain"],
            sphere_R_um=sphere_R_um,
            axis=ax,
            n_slices=n_slices,
            bin_size_um=bin_size_um,
            window_center=capture_window_center_um,
            window_size=capture_window_um,
            n_cell_types=n_cell_types,
            n_domains=n_domains,
            center=center,
        )
        adb = ad.AnnData(X=Xb, obs=obs_b, var=var)
        adb.var_names = gene_names
        adb.layers["counts_pre_batch"] = adb.X.copy()
        adb.obsm["spatial"] = spatial3d_b
        adb.obsm["cell_type_frac_true"] = ct_frac_b
        adb.obsm["domain_frac_true"] = dom_frac_b
        adb.obs["cell_type_true"] = _labels_from_fracs(ct_frac_b, ct_names)
        adb.obs["domain_true"] = _labels_from_fracs(dom_frac_b, dom_names)
        adb.obs["is_empty"] = _empty_obs_mask(dom_frac_b)

        Xb_be, factors = apply_batch_effect_per_slice_sparse(
            adb.X.tocsr(),
            slice_ids=adb.obs["slice_id"].values,
            n_slices=n_slices,
            batch_sigma=batch_sigma,
            seed=seed + 2000 + ord(ax),
        )
        adb.X = Xb_be
        adb.uns["batch_effect_factors"] = _to_serializable({str(int(k)): v.tolist() for k, v in factors.items()})
        add_spatial_keys_for_axis(
            adb,
            axis_letter=ax.upper(),
            slice_key="slice_id",
            unaligned_key="spatial_unaligned",
            base_seed=base_seed_unaligned + modality_offset["bin"] + ord(ax.upper() if sync_unaligned_seed else ax),
            max_deg=max_deg,
            max_shift=max_shift,
        )
        bin_adatas[ax.upper()] = adb

    spot_adatas = {}
    for ax in (tuple(axis.lower() for axis in output_axes) if generate_spots else ()):
        Xs, obs_s, spatial3d_s, ct_frac_s, dom_frac_s = aggregate_molecules_to_spots_2d_slices_window(
            mol_xyz=molecules["full_xyz"],
            mol_gene=molecules["full_gene"],
            mol_src_celltype=molecules["full_src_celltype"],
            mol_src_domain=molecules["full_src_domain"],
            sphere_R_um=sphere_R_um,
            axis=ax,
            n_slices=n_slices,
            spot_spacing_um=spot_spacing_um,
            spot_radius_um=spot_radius_um,
            window_center=capture_window_center_um,
            window_size=capture_window_um,
            n_cell_types=n_cell_types,
            n_domains=n_domains,
            center=center,
        )
        ads = ad.AnnData(X=Xs, obs=obs_s, var=var)
        ads.var_names = gene_names
        ads.layers["counts_pre_batch"] = ads.X.copy()
        ads.obsm["spatial"] = spatial3d_s
        ads.obsm["cell_type_frac_true"] = ct_frac_s
        ads.obsm["domain_frac_true"] = dom_frac_s
        ads.obs["cell_type_true"] = _labels_from_fracs(ct_frac_s, ct_names)
        ads.obs["domain_true"] = _labels_from_fracs(dom_frac_s, dom_names)
        ads.obs["is_empty"] = _empty_obs_mask(dom_frac_s)

        Xs_be, factors_s = apply_batch_effect_per_slice_sparse(
            ads.X.tocsr(),
            slice_ids=ads.obs["slice_id"].values,
            n_slices=n_slices,
            batch_sigma=batch_sigma,
            seed=seed + 3000 + ord(ax),
        )
        ads.X = Xs_be
        ads.uns["batch_effect_factors"] = _to_serializable({str(int(k)): v.tolist() for k, v in factors_s.items()})
        add_spatial_keys_for_axis(
            ads,
            axis_letter=ax.upper(),
            slice_key="slice_id",
            unaligned_key="spatial_unaligned",
            base_seed=base_seed_unaligned + modality_offset["spot"] + ord(ax.upper() if sync_unaligned_seed else ax),
            max_deg=max_deg,
            max_shift=max_shift,
        )
        spot_adatas[ax.upper()] = ads

    meta.update(
        captures=dict(
            captures,
            capture_window_um=capture_window_um,
            capture_window_center_um=capture_window_center_um,
            xenium_capture_window_um=xenium_capture_window_um,
        ),
        n_slices=int(n_slices),
        batch_sigma=float(batch_sigma),
        bin_size_um=float(bin_size_um),
        spot_spacing_um=float(spot_spacing_um),
        spot_radius_um=float(spot_radius_um),
        unaligned_xy=dict(
            max_deg=float(max_deg),
            max_shift=float(max_shift),
            base_seed_unaligned=int(base_seed_unaligned),
        ),
    )
    adata_cell_true.uns["sim_params"] = _to_serializable(meta)
    adata_cell_obs.uns["sim_params"] = _to_serializable(meta)

    return dict(
        adata_cell_true=adata_cell_true,
        adata_cell_obs=adata_cell_obs,
        adata_cell_sectioned=cell_sectioned,
        bin_adatas=bin_adatas,
        spot_adatas=spot_adatas,
        meta=meta,
    )



#   Example main

if __name__ == "__main__":
    outdir = "/dcs04/hicks/data/multi-sample-alignment-benchmark/code/simulation/sim_out_molecule3d_fixed"
    os.makedirs(outdir, exist_ok=True)

    # cell radii still in microns
    cell_radius_kwargs = dict(
        radius_dist="lognormal",
        r_mean=7.5,
        r_sigma=0.28,
        r_min=4.0,
        r_max=14.0,
    )

    domain_type_mix = np.array([
        # D0: slightly enriched for type1/type2
        [0.18, 0.18, 0.13, 0.12, 0.11, 0.10, 0.09, 0.09],

        # D1: slightly enriched for type3/type4
        [0.11, 0.12, 0.18, 0.18, 0.13, 0.10, 0.09, 0.09],

        # D2: slightly enriched for type5/type6
        [0.10, 0.11, 0.12, 0.13, 0.18, 0.18, 0.09, 0.09],

        # D3: slightly enriched for type7/type8
        [0.10, 0.10, 0.11, 0.12, 0.13, 0.13, 0.16, 0.15],

        # D4: mild gradient-ish mix (still balanced)
        [0.14, 0.13, 0.12, 0.11, 0.12, 0.13, 0.13, 0.12],

        # D5: near-uniform “mixed” domain
        [0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125],
    ], dtype=float)

    sim = simulate_3d_molecule_sphere_multires(
        # tissue: 12 mm diameter
        sphere_R_um=6000.0,

        # capture windows
        xenium_capture_window_um=(12000.0, 24000.0),  # 12×24 mm
        capture_window_um=(6500.0, 6500.0),           # 6.5×6.5 mm
        capture_window_center_um=(0.0, 0.0),

        # domains: noisier boundaries
        n_domains=6,
        core_frac=0.55,
        core_bump_amp=0.25,
        wedge_angle_amp_deg=25.0,
        noise_terms=16,
        noise_freq_range=(3.0,6.0),
        boundary_fuzz_width_deg=6.0,
        boundary_fuzz_flip_prob=0.15,
        core_fuzz_width_um=300.0,
        core_fuzz_flip_prob=0.25,

        # cells
        n_cells=600_000,
        cell_radius_kwargs=cell_radius_kwargs,
        allow_cell_overlap=False,

        # cell types
        n_cell_types=8,
        domain_type_mix=domain_type_mix,

        # slices + batch effect
        n_slices=10,
        batch_sigma=0.22,

        # VisiumHD bins / Visium spots
        bin_size_um=8.0,
        spot_spacing_um=100.0,
        spot_radius_um=27.5,

        # unaligned perturbation
        max_deg=270.0,
        max_shift=3000.0,
        base_seed_unaligned=12345,

        seed=2025,
    )

    ad_true = sim["adata_cell_true"]
    ad_obs = sim["adata_cell_obs"]

    # ---- Save base cell-level ----
    ad_true.write_h5ad(os.path.join(outdir, "cell_true.h5ad"), compression="gzip")
    ad_obs.write_h5ad(os.path.join(outdir, "cell_obs.h5ad"), compression="gzip")

    # ---- Save sectioned cell-level (X/Y/Z) with batch effects + aligned/unaligned coords ----
    for AX, sec in sim["adata_cell_sectioned"].items():
        sec.write_h5ad(os.path.join(outdir, f"cell_obs_sectioned_{AX}.h5ad"), compression="gzip")

    # ---- Save bin/spot all-slices (cropped to Visium window) ----
    for ax, adb in sim["bin_adatas"].items():
        adb.write_h5ad(os.path.join(outdir, f"bins_{ax}.h5ad"), compression="gzip")

    for ax, ads in sim["spot_adatas"].items():
        ads.write_h5ad(os.path.join(outdir, f"spots_{ax}.h5ad"), compression="gzip")

    print("Done. Outputs saved to:", outdir)
    print("Total molecules kept (assigned to some cell):", sim["meta"]["n_molecules_total"])
    print("Capture settings:", sim["meta"]["captures"])
    print("Domain irregularity:", sim["meta"]["domain_irregularity"])
