import unittest

import albis as ab


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


if __name__ == "__main__":
    unittest.main()
