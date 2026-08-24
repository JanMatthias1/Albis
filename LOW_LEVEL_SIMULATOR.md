# Low-level simulator reference

`sim_app.generate_data(...)` (documented in the main [README](README.md)) is a
curated, tutorial-scale entry point: it always builds exactly one
modality/axis pair, and only exposes the parameters listed in its
`DEFAULT_PARAMETERS` dict (`sim_app/api.py`) — anything not in that dict
raises `ValueError: Unsupported generate_data parameter(s)` if you try to
pass it.

The functions documented here — `simulate_3d_molecule_sphere_multires`,
`simulate_3d_molecule_sphere_base`, and `section_3d_molecule_sphere` — are
what `generate_data()` itself calls internally, with every parameter exposed
directly. This is the API used to generate the larger, production-scale
datasets (e.g. the manuscript benchmark data), and it's the only way to reach
several parameters `generate_data()` doesn't expose at all (see
[Parameters exclusive to the low-level API](#parameters-exclusive-to-the-low-level-api)
below).

## One-shot vs. two-step

`generate_data()` always re-simulates the tissue sphere, cells, domains, and
gene panel from scratch for every call. If you need multiple resolutions —
several modalities, several slice axes, or both — from the same underlying
tissue, use one of the two calling patterns below instead.

**One-shot** — build the shared 3D sphere once, then slice/capture/apply
batch effects for every requested modality/axis combination, all in one call:

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
```

**Two-step** — split tissue/cell/gene/molecule generation from
slicing/aggregation, so you can reuse the same simulated tissue across
multiple sectioning configurations without regenerating cells and genes each
time:

```python
base = sim_app.simulate_3d_molecule_sphere_base(
    sphere_R_um=300.0,
    n_cells=1_000,
    n_domains=4,
    output_modalities=("bin", "spot"),
    seed=2025,
)

sim = sim_app.section_3d_molecule_sphere(
    base,
    n_slices=5,
    capture_window_um=(300.0, 300.0),
    bin_size_um=30.0,
    output_modalities=("bin", "spot"),
    slice_axes=("X", "Z"),
)

bin_adata_x = sim["bin_adatas"]["X"]
spot_adata_z = sim["spot_adatas"]["Z"]
```

`simulate_3d_molecule_sphere_base(**kwargs)` is literally
`simulate_3d_molecule_sphere_multires(**kwargs, _section=False)` — it stops
right after tissue/cell/domain/gene/molecule generation and returns that
checkpoint (`adata_cell_true`, `adata_cell_obs`, `meta`, and a `molecules`
dict with transcript coordinates and source labels) instead of proceeding to
slicing. `section_3d_molecule_sphere(base, ...)` picks up from there. To make
bin or spot outputs later, generate `base` with `output_modalities` including
`"bin"`/`"spot"` so the full molecule stream is retained.

`output_modalities` and `slice_axes` each default to every supported value
(`{"cell", "bin", "spot"}` and `("X", "Y", "Z")`) when omitted, so a bare
call with no filters generates everything at once. `adata_cell_sectioned`,
`bin_adatas`, and `spot_adatas` are each `{axis: AnnData}` dicts, one object
per requested slicing axis.

## Full parameter reference

Parameters are grouped and ordered to match the manuscript's Methods section,
paragraph by paragraph. The **Step** column says which call each parameter
belongs to: **1** = `simulate_3d_molecule_sphere_base(...)` only, **2** =
`section_3d_molecule_sphere(...)` only, **both** = accepted by both (needed
by both if you use the two-step pattern — see the `capture_window_um` note
below for why this matters). Defaults are each function's own defaults —
several differ from `generate_data()`'s smaller tutorial-scale defaults
(noted where relevant).

### 3D tissue geometry and cell placement

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `sphere_R_um` | 1 | `6000.0` | Radius of the tissue sphere, in microns (`R`). (`generate_data()`'s `sphere_radius_um` defaults to `300.0`.) |
| `center` | 1 | `(0.0, 0.0, 0.0)` | Tissue center, in microns (`c_0`). Not exposed via `generate_data()`. |
| `n_cells` | 1 | `20000` | Number of simulated cells. (`generate_data()` defaults to `1000`.) |
| `allow_cell_overlap` | 1 | `False` | Whether overlapping cell spheres are permitted during placement. |
| `cell_radius_kwargs` | 1 | `None` | Cell radius `r_i` ~ truncated log-normal; settings are `radius_dist`, `r_mean`, `r_sigma`, `r_min`, `r_max`. `None` uses `radius_dist="lognormal", r_mean=6.0, r_sigma=0.35, r_min=3.0, r_max=12.0`. |

### Spatial domains with irregular boundaries

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `n_domains` | 1 | `4` | Number of spatial domains: `n_domains - 1` angular wedges surrounding one central core. |
| `core_frac` | 1 | `0.35` | Core radius as a fraction of the sphere radius. |
| `core_bump_amp` | 1 | `0.12` | Fractional, direction-dependent perturbation of the core radius; `0.0` yields a perfect sphere. |
| `wedge_angle_amp_deg` | 1 | `12.0` | Angular perturbation applied to wedge boundaries, in degrees; `0.0` yields straight radial cuts. |
| `noise_terms` | 1 | `6` | Number of sinusoidal components summed to build the smooth noise field underlying `core_bump_amp`/`wedge_angle_amp_deg`. |
| `noise_freq_range` | 1 | `(0.8, 2.2)` | Spatial-frequency range (cycles/micron) of that noise field. |
| `boundary_fuzz_width_deg` | 1 | `0.0` | Angular band, in degrees, around each wedge-wedge boundary within which points may be relabeled to the neighboring wedge. Off by default. |
| `boundary_fuzz_flip_prob` | 1 | `0.0` | Probability a point within `boundary_fuzz_width_deg` gets relabeled. |
| `core_fuzz_width_um` | 1 | `0.0` | Radial band, in microns, around the core boundary within which points may be relabeled across the core/wedge interface. Off by default. |
| `core_fuzz_flip_prob` | 1 | `0.0` | Probability a point within `core_fuzz_width_um` gets relabeled. |

### Cell type composition by domain

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `n_cell_types` | 1 | `4` | Number of cell types (`K`). |
| `domain_type_mix` | 1 | `None` | `(n_domains, n_cell_types)` composition matrix `pi_d`. `None` uses a built-in example composition **only** when `n_domains=4, n_cell_types=4` exactly; for any other shape it silently falls back to a **fully uniform** mix (a `UserWarning` is now raised in this case — see below). |

> **`domain_type_mix` fallback now warns.** If you change `n_domains`/
> `n_cell_types` away from `4`/`4` and don't supply `domain_type_mix`
> yourself, `simulate_3d_molecule_sphere_multires` raises a `UserWarning`
> explaining that every domain will get an identical, fully uniform
> cell-type distribution — no domain-specific enrichment at all. See the
> main README's [section 5 note](README.md#5-cell-types-genes-and-molecules)
> for the full explanation.

### Gene programs and per-cell gene expression

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `marker_genes_per_type` | 1 | `80` | Number of marker genes assigned to each cell type. |
| `noise_gene_frac` | 1 | `0.10` | Fraction of the gene panel carrying no cell-type signal. |
| `shared_marker_frac` | 1 | `0.25` | Fraction of each cell type's markers shared with other cell types, rather than unique to it. |
| `base_gene_lognormal` | 1 | `(0.7, 0.7)` | `(mu, sigma)` of the log-normal distribution each gene's baseline expression level `lambda_g` is drawn from. Not exposed via `generate_data()`. |
| `marker_foldchange` | 1 | `3.5` | Expression multiplier `m_{t,g}` applied to a cell type's own (unique) marker genes. Not exposed via `generate_data()`. |
| `shared_marker_foldchange` | 1 | `2.5` | Expression multiplier `m_{t,g}` applied to markers shared across cell types (between `marker_foldchange` and baseline). Not exposed via `generate_data()`. |
| `noise_scale` | 1 | `0.9` | Multiplier `m_{t,g}` applied to non-marker "noise gene" expression, a separate background-noise lever. Not exposed via `generate_data()`. |
| `cell_size_lognormal` | 1 | `(0.0, 0.35)` | `(mu, sigma)` of the log-normal distribution each cell's size factor `s_i` is drawn from. Not exposed via `generate_data()`. |
| `domain_size_factors` | 1 | `None` | Optional length-`n_domains` array of per-domain expression scale factors (`delta_{d_i}` in the count model). `None` means all domains get a factor of `1.0` — i.e. **no domain-level expression shift is applied**; every domain-driven signal comes from `domain_type_mix` instead. Not exposed via `generate_data()`. Passing e.g. `[1.0, 1.0, 1.0, 1.0, 1.0, 1.3]` would make domain 5 run 30% higher overall expression regardless of cell type. |
| `theta` | 1 | `25.0` | Center of the per-gene negative-binomial dispersion `theta_g` (`Var = mean + mean^2/theta`); lower values give noisier, more overdispersed counts, higher values approach Poisson. Not exposed via `generate_data()`. |
| `theta_jitter` | 1 | `2.0` | Spreads `theta` per gene via `Normal(theta, theta_jitter)`, so dispersion varies gene-to-gene instead of being fixed. Not exposed via `generate_data()`. |

### Transcript-level molecule simulation

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `inside_prob` | 1 | `0.95` | Target fraction of a cell's transcripts generated within its own cell radius `r_i`. |

### Cell-level observation and spillover

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `assign_k` | 1 | `8` | Number of nearest cells considered when reassigning transcripts to cells by containment, for `output_modalities` including `"cell"`. |

### Platform-specific capture windows

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `xenium_capture_window_um` | 1 | `(12000.0, 24000.0)` | Conceptual Xenium-like capture region (12x24 mm); recorded in output metadata but not used to crop/window cell-level data. Only accepted by step 1 — `section_3d_molecule_sphere` has no equivalent parameter. |
| `capture_window_um` | both | `(6500.0, 6500.0)` | Width/height, in microns, of the capture window used for bin/spot aggregation. (`generate_data()` defaults to `(500.0, 500.0)`.) Accepted by **both** calls: step 1 records it into `base["meta"]`, and step 2 actually uses it for aggregation — if omitted in step 2 (`None`), it inherits whatever was recorded in step 1. **If you split into two steps, pass this in step 1 too** (or explicitly in both), or step 2 silently inherits step 1's own default of `(6500.0, 6500.0)` instead of whatever you intended. |
| `capture_window_center_um` | both | `(0.0, 0.0)` | Center of the capture window, in microns. Same both-steps/inheritance behavior as `capture_window_um`. |
| `bin_size_um` | 2 | `8.0` | Bin width, in microns, for bin-level output. (`generate_data()` defaults to `20.0`.) Does **not** inherit from step 1 if omitted in step 2 — always falls back to this function's own default. |
| `spot_spacing_um` | 2 | `100.0` | Center-to-center spacing between spots, in microns. |
| `spot_radius_um` | 2 | `27.5` | Capture radius of each spot, in microns. |

### Slicing, batch effects, and unaligned coordinates

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `n_slices` | 2 | `10` | Number of slices generated along each requested axis. (`generate_data()` defaults to `5`.) |
| `slice_axes` | 2 | `None` (-> `("X","Y","Z")`) | Which axes to slice along. `generate_data()` only builds one (`slice_axis`). |
| `batch_sigma` | 2 | `0.15` | Standard deviation of the per-slice, per-gene log-fold-change (`epsilon_{s,g}`) applied to simulate batch effects across slices. |
| `max_deg` | 2 | `180.0` | Maximum absolute per-slice in-plane rotation, in degrees, applied when generating unaligned coordinates. Sampled uniformly on `[-max_deg, +max_deg]`, independently per slice. |
| `max_shift` | 2 | `200.0` | Maximum absolute per-slice in-plane translation, in microns, applied when generating unaligned coordinates. Sampled uniformly on `[-max_shift, +max_shift]` independently in x and y, so the resultant 2D displacement can exceed `max_shift` (up to `max_shift * sqrt(2)`). |
| `base_seed_unaligned` | 2 | `12345` | Base seed for the per-slice unaligned-coordinate RNG (actual seed per slice is `base_seed_unaligned + slice_id`). |

### Simulation outputs

| Parameter | Step | Default | Meaning |
| --- | :-: | ---: | --- |
| `output_modalities` | both | `None` (-> all) | Which of `{"cell", "bin", "spot"}` to generate. `generate_data()` only builds one (`output`). Accepted by **both** calls: step 1 needs it to know which molecule streams to retain for later sectioning, step 2 needs it to know what to actually build. |
| `sparse_X` | 1 | `True` | Whether count matrices are stored as sparse (`scipy.sparse`) or dense arrays. Not exposed via `generate_data()`. |
| `seed` | 1 | `2025` | Random seed. |

## Parameters exclusive to the low-level API

These have no equivalent in `generate_data()`'s `DEFAULT_PARAMETERS` at all —
reachable only through `simulate_3d_molecule_sphere_multires`/`_base`/
`section_3d_molecule_sphere`:

`center`, `base_gene_lognormal`, `marker_foldchange`,
`shared_marker_foldchange`, `noise_scale`, `cell_size_lognormal`,
`domain_size_factors`, `theta`, `theta_jitter`, `sparse_X`,
`output_modalities` (multi-modality in one call), `slice_axes` (multi-axis in
one call).
