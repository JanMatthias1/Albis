import unittest

import numpy as np

import albis as ab


_SMALL = dict(slice_axis="Z", n_cells=300, n_slices=1, sphere_radius_um=250.0, seed=7)


class PublicApiTest(unittest.TestCase):
    def test_public_api_exports_expected_symbols(self):
        self.assertEqual(ab.__version__, "0.1.0")
        self.assertTrue(callable(ab.generate_data))
        self.assertTrue(callable(ab.example_data))
        self.assertTrue(callable(ab.describe))
        self.assertTrue(callable(ab.save))
        self.assertTrue(callable(ab.plot))
        self.assertTrue(callable(ab.simulate_3d_molecule_sphere_base))
        self.assertTrue(callable(ab.section_3d_molecule_sphere))

    def test_generate_data_rejects_unknown_parameters(self):
        with self.assertRaisesRegex(ValueError, "Unsupported generate_data parameter"):
            ab.generate_data(unknown_parameter=True)

    def test_capture_window_rejects_bad_values(self):
        with self.assertRaisesRegex(ValueError, "capture_window_um"):
            ab.generate_data(output="bin", capture_window_um=(0.0, 100.0))
        with self.assertRaisesRegex(ValueError, "capture_window_um"):
            ab.generate_data(output="bin", capture_window_um="visium")

    def test_platform_default_bins_cover_full_visium_window(self):
        adata = ab.generate_data(output="bin", **_SMALL)
        self.assertEqual(adata.obs["window_size0"].iloc[0], 6500.0)
        self.assertEqual(adata.obs["window_size1"].iloc[0], 6500.0)

    def test_capture_window_false_shrinks_grid_to_tissue(self):
        full = ab.generate_data(output="bin", **_SMALL)
        nocrop = ab.generate_data(output="bin", capture_window_um=False, **_SMALL)
        self.assertLess(nocrop.obs["window_size0"].iloc[0], 6500.0)
        self.assertLess(nocrop.n_obs, full.n_obs)

    def test_xenium_window_crops_cells(self):
        keep_all = ab.generate_data(output="cell", capture_window_um=False, **_SMALL)
        cropped = ab.generate_data(output="cell", capture_window_um=(150.0, 150.0), **_SMALL)
        self.assertLess(cropped.n_obs, keep_all.n_obs)
        xy = np.asarray(cropped.obsm["spatial"])
        self.assertTrue(np.all(np.abs(xy) <= 75.0 + 1e-6))

    def test_empty_bins_are_flagged_and_not_mislabelled(self):
        # 250 um tissue inside the 6.5 mm platform window -> most bins are empty.
        adata = ab.generate_data(output="bin", bin_size_um=60.0, **_SMALL)
        empty = adata.obs["is_empty"].to_numpy(dtype=bool)
        x_empty = np.asarray(adata.layers["counts_pre_batch"].sum(axis=1)).ravel() == 0
        self.assertTrue(empty.any() and not empty.all())
        np.testing.assert_array_equal(empty, x_empty)
        # empty observations carry no ground-truth class ...
        self.assertEqual(set(adata.obs.loc[empty, "domain_true"]), {"unassigned"})
        self.assertEqual(set(adata.obs.loc[empty, "cell_type_true"]), {"unassigned"})
        # ... and no non-empty observation is labelled "unassigned"
        self.assertNotIn("unassigned", set(adata.obs.loc[~empty, "domain_true"]))
        self.assertNotIn("unassigned", set(adata.obs.loc[~empty, "cell_type_true"]))

    def test_cell_output_has_no_is_empty_flag(self):
        adata = ab.generate_data(output="cell", **_SMALL)
        self.assertNotIn("is_empty", adata.obs)
        self.assertNotIn("unassigned", set(adata.obs["domain_true"]))


if __name__ == "__main__":
    unittest.main()
