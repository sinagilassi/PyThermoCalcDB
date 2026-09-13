import math
import unittest

from pythermocalcdb.compositions import calc_ideal_osmotic_pressure
from pythermocalcdb.compositions.core.osmotic import _calc_ideal_osmotic_pressure
from pythermocalcdb.mixtures import (
    calc_binary_partial_molar_properties,
    calc_excess_gibbs_energy_from_activity_coefficients,
    calc_excess_property,
    calc_gibbs_duhem_residual,
    calc_molar_property_from_partial_molar_properties,
    calc_total_property_from_partial_molar_properties,
    check_gibbs_duhem_consistency,
)
from pythermocalcdb.mixtures.core.excess import (
    _calc_excess_entropy_from_gibbs_enthalpy,
)
from pythermocalcdb.thermo import (
    calc_cp_from_cv_general,
    calc_cp_minus_cv_general,
    calc_cv_from_cp_general,
    calc_enthalpy_vaporization_watson,
    calc_isentropic_compressibility,
    calc_joule_thomson_coefficient_from_alpha,
    calc_speed_of_sound,
    calc_speed_of_sound_from_isentropic_compressibility,
    calc_transition_enthalpy_from_constant_delta_cp,
    calc_transition_enthalpy_from_cp_integral,
)
from pythermocalcdb.thermo.core.derivatives import (
    _calc_isothermal_compressibility_from_density,
)


class TestMissingLowLevelGuidanceV6(unittest.TestCase):
    def test_derivative_primitives(self):
        self.assertTrue(math.isclose(
            _calc_isothermal_compressibility_from_density(1000.0, 1.0e-6),
            1.0e-9,
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(
            calc_isentropic_compressibility(1.0e-9, 20.0, 30.0),
            2.0e-9 / 3.0,
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(
            calc_joule_thomson_coefficient_from_alpha(300.0, 0.024, 30.0, 1.0e-3),
            0.024 * (0.3 - 1.0) / 30.0,
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(
            calc_speed_of_sound_from_isentropic_compressibility(1000.0, 1.0e-9),
            math.sqrt(1.0 / (1000.0 * 1.0e-9)),
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(
            calc_speed_of_sound(1000.0, 1.0e-9, 30.0, 20.0),
            math.sqrt(30.0 / (1000.0 * 1.0e-9 * 20.0)),
            rel_tol=1e-12,
        ))

    def test_general_heat_capacity_relation(self):
        delta = calc_cp_minus_cv_general(300.0, 1.0e-4, 1.0e-3, 1.0e-9)
        self.assertTrue(math.isclose(delta, 30.0, rel_tol=1e-12))
        self.assertTrue(math.isclose(
            calc_cv_from_cp_general(80.0, 300.0, 1.0e-4, 1.0e-3, 1.0e-9),
            50.0,
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(
            calc_cp_from_cv_general(50.0, 300.0, 1.0e-4, 1.0e-3, 1.0e-9),
            80.0,
            rel_tol=1e-12,
        ))

    def test_phase_change_corrections(self):
        expected = 40000.0 * ((1.0 - 350.0 / 500.0) / (1.0 - 300.0 / 500.0)) ** 0.38
        self.assertTrue(math.isclose(
            calc_enthalpy_vaporization_watson(40000.0, 300.0, 350.0, 500.0),
            expected,
            rel_tol=1e-12,
        ))
        self.assertEqual(
            calc_transition_enthalpy_from_constant_delta_cp(10000.0, 300.0, 350.0, 10.0),
            10500.0,
        )
        self.assertEqual(
            calc_transition_enthalpy_from_cp_integral(10000.0, 500.0),
            10500.0,
        )

    def test_partial_molar_and_excess_primitives(self):
        self.assertEqual(
            calc_total_property_from_partial_molar_properties([1.0, 2.0], [10.0, 20.0]),
            50.0,
        )
        self.assertEqual(
            calc_molar_property_from_partial_molar_properties([0.25, 0.75], [10.0, 20.0]),
            17.5,
        )
        self.assertEqual(calc_binary_partial_molar_properties(15.0, 0.25, 4.0), (18.0, 14.0))
        self.assertEqual(calc_excess_property(12.0, 10.0), 2.0)
        self.assertTrue(math.isclose(
            calc_excess_gibbs_energy_from_activity_coefficients([0.5, 0.5], [2.0, 1.0], 300.0),
            8.314462618 * 300.0 * 0.5 * math.log(2.0),
            rel_tol=1e-12,
        ))
        self.assertEqual(_calc_excess_entropy_from_gibbs_enthalpy(100.0, 250.0, 300.0), 0.5)
        self.assertEqual(calc_gibbs_duhem_residual([0.25, 0.75], [3.0, -1.0]), 0.0)
        self.assertTrue(check_gibbs_duhem_consistency([0.25, 0.75], [3.0, -1.0]))

    def test_ideal_osmotic_pressure(self):
        self.assertTrue(math.isclose(
            _calc_ideal_osmotic_pressure(1000.0, 300.0),
            1000.0 * 8.314462618 * 300.0,
            rel_tol=1e-12,
        ))
        result = calc_ideal_osmotic_pressure(1000.0, 300.0, vant_hoff_factor=2.0)
        self.assertTrue(math.isclose(result.value, 2.0 * 1000.0 * 8.314462618 * 300.0, rel_tol=1e-12))
        self.assertEqual(result.unit, "Pa")


if __name__ == "__main__":
    unittest.main()
