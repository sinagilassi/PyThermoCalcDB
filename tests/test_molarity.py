import unittest

import numpy as np

from pythermocalcdb.compositions.molarity import _calc_molarities


class TestMolarity(unittest.TestCase):
    def test_scalar_volume_returns_component_moles_shape(self):
        result = _calc_molarities(
            component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
            solution_volume=2.0,
        )

        np.testing.assert_allclose(result, np.array([[0.5, 1.0], [1.5, 2.0]]))
        self.assertEqual(result.shape, (2, 2))

    def test_matching_volume_shape_is_allowed(self):
        result = _calc_molarities(
            component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
            solution_volume=np.array([[1.0, 2.0], [3.0, 4.0]]),
        )

        np.testing.assert_allclose(result, np.ones((2, 2)))

    def test_state_volume_shape_is_allowed(self):
        result = _calc_molarities(
            component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
            solution_volume=np.array([1.0, 2.0]),
        )

        np.testing.assert_allclose(result, np.array([[1.0, 2.0], [1.5, 2.0]]))
        self.assertEqual(result.shape, (2, 2))

    def test_non_broadcastable_volume_shape_raises_validation_error(self):
        with self.assertRaisesRegex(ValueError, "one entry per state"):
            _calc_molarities(
                component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
                solution_volume=np.array([1.0, 2.0, 3.0]),
            )

    def test_column_state_volume_shape_is_allowed(self):
        result = _calc_molarities(
            component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
            solution_volume=np.array([[1.0], [2.0]]),
        )

        np.testing.assert_allclose(result, np.array([[1.0, 2.0], [1.5, 2.0]]))
        self.assertEqual(result.shape, (2, 2))


if __name__ == "__main__":
    unittest.main()
