import math
import unittest

import numpy as np
from pythermodb_settings.models import AnnotatedValue

from pythermocalcdb.density import (
    calc_density_from_molar_volume,
    calc_gas_density_from_compressibility,
    calc_multiphase_density_from_volume_fractions,
    calc_rackett_constant_from_acentric_factor,
    calc_saturated_liquid_molar_volume_rackett,
)
from pythermocalcdb.density.core import (
    _calc_density_from_molar_volume,
    _calc_gas_density_from_compressibility,
    _calc_multiphase_density_from_volume_fractions,
    _calc_rackett_constant_from_acentric_factor,
    _calc_saturated_liquid_molar_volume_rackett,
)
from pythermocalcdb.reactions import (
    calc_component_moles_from_reaction_extent,
    calc_component_moles_from_reaction_extents,
)
from pythermocalcdb.reactions.core import (
    _calc_component_moles_from_reaction_extent,
    _calc_component_moles_from_reaction_extents,
)
from pythermocalcdb.thermo import (
    calc_activity_coefficient_from_fugacity,
    calc_cp_departure_from_eos_derivatives,
    calc_dimensionless_enthalpy_departure,
    calc_dimensionless_entropy_departure,
    calc_enthalpy_departure,
    calc_enthalpy_from_ideal_and_departure,
    calc_entropy_departure,
    calc_entropy_from_ideal_and_departure,
    calc_fugacity_coefficient,
    calc_fugacity_from_coefficient,
    calc_phase_equilibrium_fugacity_residual,
)
from pythermocalcdb.thermo.core import (
    _calc_cp_departure_from_eos_derivatives,
    _calc_enthalpy_departure,
    _calc_fugacity_coefficient,
)
from pythermocalcdb.transport import (
    calc_liquid_mixture_viscosity_log_rule,
    calc_viscosity_exponential_correlation,
)
from pythermocalcdb.transport.core import (
    _calc_liquid_mixture_viscosity_log_rule,
    _calc_viscosity_exponential_correlation,
)


class TestDwsimGuidanceV9(unittest.TestCase):
    def test_density_definitions_and_rackett(self):
        self.assertTrue(math.isclose(_calc_density_from_molar_volume(0.018, 18.0e-6), 1000.0, rel_tol=1e-12))
        result = calc_density_from_molar_volume(0.018, 18.0e-6)
        self.assertIsInstance(result, AnnotatedValue)
        self.assertTrue(math.isclose(result.value, 1000.0, rel_tol=1e-12))
        self.assertEqual(result.unit, "kg/m3")

        rho_z = _calc_gas_density_from_compressibility(0.02897, 101325.0, 300.0, 1.0)
        expected_rho_z = 0.02897 * 101325.0 / (8.314462618 * 300.0)
        self.assertTrue(math.isclose(rho_z, expected_rho_z, rel_tol=1e-12))
        expected_rho_z_public = 0.02897 * 101325.0 / (8.31446261815324 * 300.0)
        self.assertTrue(math.isclose(
            calc_gas_density_from_compressibility(0.02897, 101325.0, 300.0, 1.0).value,
            expected_rho_z_public,
            rel_tol=1e-12,
        ))

        self.assertEqual(_calc_multiphase_density_from_volume_fractions([800.0, 2.0], [0.75, 0.25]), 600.5)
        self.assertEqual(calc_multiphase_density_from_volume_fractions([800.0, 2.0], [0.75, 0.25]).value, 600.5)

        zra = _calc_rackett_constant_from_acentric_factor(0.2)
        self.assertTrue(math.isclose(zra, 0.2956 - 0.08775 * 0.2, rel_tol=1e-12))
        self.assertTrue(math.isclose(calc_rackett_constant_from_acentric_factor(0.2).value, zra, rel_tol=1e-12))
        volume = _calc_saturated_liquid_molar_volume_rackett(300.0, 500.0, 5.0e6, zra)
        self.assertGreater(volume, 0.0)
        tr = 300.0 / 500.0
        expected_volume_public = (8.31446261815324 * 500.0 / 5.0e6) * zra ** (1.0 + (1.0 - tr) ** (2.0 / 7.0))
        self.assertTrue(math.isclose(
            calc_saturated_liquid_molar_volume_rackett(300.0, 500.0, 5.0e6, zra).value,
            expected_volume_public,
            rel_tol=1e-12,
        ))
        with self.assertRaises(ValueError):
            _calc_saturated_liquid_molar_volume_rackett(500.0, 500.0, 5.0e6, zra)

    def test_reaction_extent_stoichiometry(self):
        n = _calc_component_moles_from_reaction_extent([2.0, 1.0, 0.0], [-1.0, -0.5, 1.0], 1.0)
        self.assertTrue(np.allclose(n, [1.0, 0.5, 1.0]))
        result = calc_component_moles_from_reaction_extent([2.0, 1.0, 0.0], [-1.0, -0.5, 1.0], 1.0)
        self.assertIsInstance(result, AnnotatedValue)
        self.assertEqual(result.value, [1.0, 0.5, 1.0])
        n_multi = _calc_component_moles_from_reaction_extents([2.0, 1.0, 0.0], [[-1.0, -0.5, 1.0]], [1.0])
        self.assertTrue(np.allclose(n_multi, [1.0, 0.5, 1.0]))
        self.assertEqual(calc_component_moles_from_reaction_extents([2.0, 1.0, 0.0], [[-1.0, -0.5, 1.0]], [1.0]).value, [1.0, 0.5, 1.0])
        with self.assertRaises(ValueError):
            _calc_component_moles_from_reaction_extent([0.1], [-1.0], 1.0)

    def test_departure_fugacity_and_activity_helpers(self):
        self.assertEqual(_calc_enthalpy_departure(1200.0, 1000.0), 200.0)
        self.assertEqual(calc_enthalpy_departure(1200.0, 1000.0).value, 200.0)
        self.assertEqual(calc_entropy_departure(12.0, 10.0).value, 2.0)
        self.assertTrue(math.isclose(calc_dimensionless_enthalpy_departure(2494.338785445972, 300.0).value, 1.0, rel_tol=1e-12))
        self.assertTrue(math.isclose(calc_dimensionless_entropy_departure(8.31446261815324).value, 1.0, rel_tol=1e-12))
        self.assertEqual(calc_enthalpy_from_ideal_and_departure(1000.0, 200.0).value, 1200.0)
        self.assertEqual(calc_entropy_from_ideal_and_departure(10.0, 2.0).value, 12.0)
        cp_dep = _calc_cp_departure_from_eos_derivatives(300.0, 0.2, 10.0, -1000.0)
        expected_cp_dep = 300.0 * 0.2 - 300.0 * 10.0**2 / -1000.0 - 8.314462618
        self.assertTrue(math.isclose(cp_dep, expected_cp_dep, rel_tol=1e-12))
        expected_cp_dep_public = 300.0 * 0.2 - 300.0 * 10.0**2 / -1000.0 - 8.31446261815324
        self.assertTrue(math.isclose(calc_cp_departure_from_eos_derivatives(300.0, 0.2, 10.0, -1000.0).value, expected_cp_dep_public, rel_tol=1e-12))

        self.assertEqual(_calc_fugacity_coefficient(5.0e4, 0.5, 1.0e5), 1.0)
        self.assertEqual(calc_fugacity_coefficient(5.0e4, 0.5, 1.0e5), 1.0)
        self.assertEqual(calc_fugacity_from_coefficient(1.0, 0.5, 1.0e5), 5.0e4)
        self.assertEqual(calc_phase_equilibrium_fugacity_residual(5.0e4, 4.0e4), 1.0e4)
        self.assertEqual(calc_activity_coefficient_from_fugacity(5.0e4, 0.5, 1.0e5), 1.0)

    def test_viscosity_rules(self):
        eta = _calc_liquid_mixture_viscosity_log_rule([0.25, 0.75], [1.0e-3, 4.0e-3])
        expected = math.exp(0.25 * math.log(1.0e-3) + 0.75 * math.log(4.0e-3))
        self.assertTrue(math.isclose(eta, expected, rel_tol=1e-12))
        self.assertTrue(math.isclose(calc_liquid_mixture_viscosity_log_rule([0.25, 0.75], [1.0e-3, 4.0e-3]), expected, rel_tol=1e-12))
        corr = _calc_viscosity_exponential_correlation(300.0, 1.0, 2.0, 0.0, 0.0, 1.0)
        self.assertTrue(math.isclose(corr, math.exp(1.0 + 2.0 / 300.0), rel_tol=1e-12))
        self.assertTrue(math.isclose(calc_viscosity_exponential_correlation(300.0, 1.0, 2.0, 0.0, 0.0, 1.0), corr, rel_tol=1e-12))
        with self.assertRaises(ValueError):
            _calc_liquid_mixture_viscosity_log_rule([0.2, 0.2], [1.0, 2.0])


if __name__ == "__main__":
    unittest.main()





