import math
import unittest

import numpy as np

from pythermocalcdb.conversions.property_basis import (
    _calc_molar_to_mass_specific,
    _calc_mass_specific_to_molar,
    calc_molar_to_mass_specific,
    calc_mass_specific_to_molar,
    calc_molar_cp_to_mass_cp,
    calc_mass_cp_to_molar_cp,
    molar_to_mass_specific,
    mass_specific_to_molar,
    molar_cp_to_mass_cp,
    mass_cp_to_molar_cp,
)


class TestPropertyBasisConversions(unittest.TestCase):
    def test_public_calc_wrappers_delegate_to_core_equations(self):
        self.assertTrue(math.isclose(calc_molar_to_mass_specific(10.0, 2.0), 5.0))
        self.assertTrue(math.isclose(calc_mass_specific_to_molar(5.0, 2.0), 10.0))

    def test_legacy_public_names_alias_calc_wrappers(self):
        self.assertIs(molar_to_mass_specific, calc_molar_to_mass_specific)
        self.assertIs(mass_specific_to_molar, calc_mass_specific_to_molar)
        self.assertIs(molar_cp_to_mass_cp, calc_molar_cp_to_mass_cp)
        self.assertIs(mass_cp_to_molar_cp, calc_mass_cp_to_molar_cp)

    def test_core_scalar_inputs_return_float(self):
        result = _calc_molar_to_mass_specific(10, 2)
        self.assertIsInstance(result, float)
        self.assertTrue(math.isclose(result, 5.0))

    def test_core_sequence_inputs_return_float_array(self):
        result = _calc_mass_specific_to_molar([5, 10, 15], [2, 2, 2])
        np.testing.assert_allclose(result, np.array([10.0, 20.0, 30.0]))
        self.assertEqual(result.dtype, np.float64)

    def test_core_2d_inputs_support_broadcasting(self):
        molar_property = np.array([[10.0, 20.0], [30.0, 40.0]])
        molecular_weight = np.array([[2.0], [10.0]])
        np.testing.assert_allclose(
            _calc_molar_to_mass_specific(molar_property, molecular_weight),
            np.array([[5.0, 10.0], [3.0, 4.0]]),
        )

    def test_core_round_trip_with_negative_property(self):
        mass_specific = _calc_molar_to_mass_specific([-10.0, 20.0], [2.0, 4.0])
        np.testing.assert_allclose(mass_specific, np.array([-5.0, 5.0]))
        np.testing.assert_allclose(
            _calc_mass_specific_to_molar(mass_specific, [2.0, 4.0]),
            np.array([-10.0, 20.0]),
        )

    def test_invalid_molecular_weight_and_shapes_raise(self):
        with self.assertRaises(ValueError):
            _calc_molar_to_mass_specific(10.0, 0.0)
        with self.assertRaises(ValueError):
            _calc_mass_specific_to_molar(10.0, -1.0)
        with self.assertRaises(ValueError):
            _calc_molar_to_mass_specific([1.0, 2.0], [1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            _calc_mass_specific_to_molar(float("nan"), 1.0)


if __name__ == "__main__":
    unittest.main()
