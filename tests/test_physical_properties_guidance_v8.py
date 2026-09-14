import math
import unittest

import numpy as np
from pythermodb_settings.models import AnnotatedValue

from pythermocalcdb.critical import (
    calc_acentric_factor_from_reduced_vapor_pressure,
    calc_critical_compressibility_factor,
)
from pythermocalcdb.critical.core import (
    _calc_acentric_factor_from_vapor_pressure,
    _calc_critical_compressibility_factor,
    _calc_critical_volume_from_zc,
)
from pythermocalcdb.heat_capacity import (
    calc_enthalpy_change_from_cp_polynomial,
    calc_entropy_change_from_cp_polynomial,
    calc_heat_capacity_polynomial,
)
from pythermocalcdb.heat_capacity.core import _calc_heat_capacity_polynomial
from pythermocalcdb.phase_change import calc_heat_of_vaporization_watson
from pythermocalcdb.phase_change.core import _calc_heat_of_vaporization_watson
from pythermocalcdb.solubility import calc_solubility_parameter_from_hvap_volume
from pythermocalcdb.solubility.core import (
    _calc_internal_energy_of_vaporization,
    _calc_solubility_parameter_from_hvap_density,
)
from pythermocalcdb.states import calc_reduced_temperature
from pythermocalcdb.states.core import (
    _calc_reduced_pressure,
    _calc_reduced_temperature,
    _calc_reduced_volume,
)


class TestPhysicalPropertiesGuidanceV8(unittest.TestCase):
    def test_reduced_properties_scalar_array_and_invalid(self):
        self.assertEqual(_calc_reduced_temperature(300.0, 600.0), 0.5)
        self.assertTrue(np.allclose(
            _calc_reduced_pressure([1.0e5, 2.0e5], 1.0e6),
            [0.1, 0.2],
        ))
        volumes = _calc_reduced_volume([[1.0, 2.0], [3.0, 4.0]], [2.0, 4.0])
        self.assertTrue(np.allclose(volumes, [[0.5, 0.5], [1.5, 1.0]]))
        result = calc_reduced_temperature(300.0, 600.0)
        self.assertIsInstance(result, AnnotatedValue)
        self.assertEqual(result.value, 0.5)
        self.assertIsNone(result.unit)
        with self.assertRaises(ValueError):
            _calc_reduced_temperature(300.0, 0.0)

    def test_critical_compressibility_and_acentric_factor(self):
        zc = _calc_critical_compressibility_factor(5.0e6, 2.0e-4, 400.0)
        expected = 5.0e6 * 2.0e-4 / (8.31446261815324 * 400.0)
        self.assertTrue(math.isclose(zc, expected, rel_tol=1e-12))
        self.assertTrue(math.isclose(
            _calc_critical_volume_from_zc(zc, 400.0, 5.0e6),
            2.0e-4,
            rel_tol=1e-12,
        ))
        zc_public = calc_critical_compressibility_factor(5.0e6, 2.0e-4, 400.0)
        self.assertIsInstance(zc_public, AnnotatedValue)
        self.assertTrue(math.isclose(zc_public.value, expected))
        self.assertIsNone(zc_public.unit)
        omega = calc_acentric_factor_from_reduced_vapor_pressure(0.1)
        self.assertIsInstance(omega, AnnotatedValue)
        self.assertEqual(omega.value, 0.0)
        self.assertTrue(np.allclose(
            _calc_acentric_factor_from_vapor_pressure([1.0e5, 1.0e4], 1.0e6),
            [0.0, 1.0],
        ))
        with self.assertRaises(ValueError):
            _calc_acentric_factor_from_vapor_pressure(0.0, 1.0e6)

    def test_heat_capacity_polynomial_and_integrals(self):
        cp = calc_heat_capacity_polynomial(300.0, 10.0, 0.1)
        self.assertIsInstance(cp, AnnotatedValue)
        self.assertEqual(cp.value, 40.0)
        self.assertEqual(cp.unit, "J/(mol.K)")
        self.assertTrue(np.allclose(
            _calc_heat_capacity_polynomial([300.0, 400.0], 10.0, 0.1),
            [40.0, 50.0],
        ))
        dh = calc_enthalpy_change_from_cp_polynomial(300.0, 400.0, 10.0, 0.1)
        self.assertIsInstance(dh, AnnotatedValue)
        self.assertEqual(dh.value, 10.0 * 100.0 + 0.1 / 2.0 * (400.0**2 - 300.0**2))
        self.assertEqual(dh.unit, "J/mol")
        ds = calc_entropy_change_from_cp_polynomial(300.0, 400.0, 10.0, 0.1)
        self.assertIsInstance(ds, AnnotatedValue)
        self.assertTrue(math.isclose(
            ds.value,
            10.0 * math.log(400.0 / 300.0) + 0.1 * 100.0,
            rel_tol=1e-12,
        ))
        with self.assertRaises(ValueError):
            _calc_heat_capacity_polynomial(0.0, 10.0)

    def test_phase_change_watson_alias(self):
        expected = 40000.0 * ((1.0 - 350.0 / 500.0) / (1.0 - 300.0 / 500.0)) ** 0.38
        self.assertTrue(math.isclose(
            _calc_heat_of_vaporization_watson(40000.0, 350.0, 300.0, 500.0),
            expected,
            rel_tol=1e-12,
        ))
        result = calc_heat_of_vaporization_watson(40000.0, 350.0, 300.0, 500.0)
        self.assertIsInstance(result, AnnotatedValue)
        self.assertTrue(math.isclose(result.value, expected, rel_tol=1e-12))
        self.assertEqual(result.unit, "J/mol")

    def test_hildebrand_solubility_helpers(self):
        du = _calc_internal_energy_of_vaporization(40000.0, 300.0)
        expected_du = 40000.0 - 8.31446261815324 * 300.0
        self.assertTrue(math.isclose(du, expected_du, rel_tol=1e-12))
        result = calc_solubility_parameter_from_hvap_volume(40000.0, 300.0, 1.0e-4)
        self.assertIsInstance(result, AnnotatedValue)
        self.assertTrue(math.isclose(result.value, math.sqrt(expected_du / 1.0e-4), rel_tol=1e-12))
        self.assertEqual(result.unit, "Pa^0.5")
        self.assertTrue(math.isclose(
            _calc_solubility_parameter_from_hvap_density(40000.0, 300.0, 1.0e4),
            math.sqrt(expected_du * 1.0e4),
            rel_tol=1e-12,
        ))
        with self.assertRaises(ValueError):
            _calc_internal_energy_of_vaporization(1000.0, 300.0)


if __name__ == "__main__":
    unittest.main()
