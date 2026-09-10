import math
import unittest

import numpy as np

from pythermocalcdb.conversions.extensive_intensive import (
    _calc_molar_property_to_total,
    _calc_specific_property_to_total,
    _calc_total_to_molar_property,
    _calc_total_to_specific_property,
    calc_molar_property_to_total,
    calc_specific_property_to_total,
    calc_total_to_molar_property,
    calc_total_to_specific_property,
    molar_property_to_total,
    specific_property_to_total,
    total_to_molar_property,
    total_to_specific_property,
)


class TestExtensiveIntensiveConversions(unittest.TestCase):
    def test_public_calc_wrappers_delegate_to_core_equations(self):
        self.assertTrue(math.isclose(calc_molar_property_to_total(2.0, 10.0), 20.0))
        self.assertTrue(math.isclose(calc_total_to_molar_property(20.0, 2.0), 10.0))
        self.assertTrue(math.isclose(calc_specific_property_to_total(3.0, 4.0), 12.0))
        self.assertTrue(math.isclose(calc_total_to_specific_property(12.0, 3.0), 4.0))

    def test_legacy_public_names_alias_calc_wrappers(self):
        self.assertIs(molar_property_to_total, calc_molar_property_to_total)
        self.assertIs(total_to_molar_property, calc_total_to_molar_property)
        self.assertIs(specific_property_to_total, calc_specific_property_to_total)
        self.assertIs(total_to_specific_property, calc_total_to_specific_property)

    def test_core_scalar_inputs_return_float(self):
        result = _calc_molar_property_to_total(2, 10)
        self.assertIsInstance(result, float)
        self.assertTrue(math.isclose(result, 20.0))

    def test_core_sequence_inputs_return_float_array(self):
        result = _calc_total_to_molar_property([10, 20, 30], [1, 2, 3])
        np.testing.assert_allclose(result, np.array([10.0, 10.0, 10.0]))
        self.assertEqual(result.dtype, np.float64)

    def test_core_2d_inputs_support_broadcasting(self):
        total = np.array([[10.0, 20.0], [30.0, 40.0]])
        moles = np.array([[1.0], [10.0]])
        np.testing.assert_allclose(
            _calc_total_to_molar_property(total, moles),
            np.array([[10.0, 20.0], [3.0, 4.0]]),
        )

    def test_core_mass_specific_round_trip_with_negative_property(self):
        total = _calc_specific_property_to_total([1.0, 2.0], [-5.0, 10.0])
        np.testing.assert_allclose(total, np.array([-5.0, 20.0]))
        np.testing.assert_allclose(
            _calc_total_to_specific_property(total, [1.0, 2.0]),
            np.array([-5.0, 10.0]),
        )

    def test_invalid_amount_mass_and_shapes_raise(self):
        with self.assertRaises(ValueError):
            _calc_molar_property_to_total(0.0, 10.0)
        with self.assertRaises(ValueError):
            _calc_specific_property_to_total(-1.0, 10.0)
        with self.assertRaises(ValueError):
            _calc_total_to_specific_property([1.0, 2.0], [1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            _calc_total_to_molar_property(float("nan"), 1.0)


if __name__ == "__main__":
    unittest.main()
