import unittest

import sim_app


class PublicApiTest(unittest.TestCase):
    def test_public_api_exports_expected_symbols(self):
        self.assertEqual(sim_app.__version__, "0.1.0")
        self.assertTrue(callable(sim_app.generate_data))
        self.assertTrue(callable(sim_app.example_data))
        self.assertTrue(callable(sim_app.describe))
        self.assertTrue(callable(sim_app.save))
        self.assertTrue(callable(sim_app.plot))

    def test_generate_data_rejects_unknown_parameters(self):
        with self.assertRaisesRegex(ValueError, "Unsupported generate_data parameter"):
            sim_app.generate_data(unknown_parameter=True)


if __name__ == "__main__":
    unittest.main()
