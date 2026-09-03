import math
import unittest

from pythermodb_settings.models import CustomProp, Pressure, Temperature

from pythermocalcdb.compositions import (
    calc_charge_balance,
    calc_ionic_strength_molality,
    calc_ionic_strength_molarity,
    calc_normality,
    check_electroneutrality,
)
from pythermocalcdb.mixtures import (
    calc_ideal_gibbs_energy_of_mixing,
    calc_ideal_molar_entropy_of_mixing,
    calc_ideal_molar_gibbs_energy_of_mixing,
    calc_mixture_molecular_weight_from_mass_fractions,
    calc_mixture_molecular_weight_from_mole_fractions,
)
from pythermocalcdb.reactions import (
    calc_dlnK_dT,
    calc_equilibrium_constant,
    calc_equilibrium_constant_at_temperature,
    calc_log_equilibrium_constant,
    calc_log_reaction_quotient,
    calc_reaction_entropy_std,
    calc_reaction_entropy_std_from_enthalpy_gibbs,
    calc_reaction_gibbs_energy,
    calc_reaction_quotient,
)
from pythermocalcdb.thermo import (
    calc_chemical_potential_from_activity,
    calc_gas_density_from_z,
    calc_gas_molar_volume_from_z,
    calc_heat_capacity_ratio,
    calc_ideal_gas_cp_from_cv,
    calc_ideal_gas_cv_from_cp,
    calc_ideal_gas_density,
    calc_ideal_gas_molar_volume,
)


class TestLowLevelThermoExtensions(unittest.TestCase):
    def test_mixture_molecular_weight(self):
        self.assertEqual(
            calc_mixture_molecular_weight_from_mole_fractions([1.0], [18.015]),
            18.015,
        )
        result = calc_mixture_molecular_weight_from_mass_fractions(
            [0.25, 0.75], [18.015, 46.07]
        )
        self.assertTrue(math.isclose(result, 33.15990810567847, rel_tol=1e-12))
        keyed = calc_mixture_molecular_weight_from_mole_fractions(
            {"water": 0.25, "ethanol": 0.75},
            {"water": 18.015, "ethanol": 46.07},
        )
        self.assertTrue(math.isclose(keyed, 39.05625, rel_tol=1e-12))

    def test_reaction_entropy(self):
        self.assertTrue(math.isclose(
            calc_reaction_entropy_std([-1.0, -0.5, 1.0], [130.68, 205.15, 213.74]),
            213.74 - 130.68 - 0.5 * 205.15,
            rel_tol=1e-12,
        ))
        result = calc_reaction_entropy_std_from_enthalpy_gibbs(
            -10000.0,
            -20000.0,
            Temperature(value=298.15, unit="K"),
        )
        self.assertTrue(math.isclose(result, 10000.0 / 298.15, rel_tol=1e-12))

    def test_equilibrium_primitives(self):
        temperature = Temperature(value=298.15, unit="K")
        self.assertEqual(calc_log_equilibrium_constant(0.0, temperature), 0.0)
        self.assertEqual(calc_equilibrium_constant(0.0, temperature), 1.0)
        self.assertEqual(calc_log_reaction_quotient([-1.0, 1.0], [1.0, 1.0]), 0.0)
        self.assertEqual(calc_reaction_quotient([-1.0, 1.0], [1.0, 1.0]), 1.0)

        delta_g_std = -5000.0
        ln_k = calc_log_equilibrium_constant(delta_g_std, temperature)
        delta_g = calc_reaction_gibbs_energy(
            delta_g_std,
            temperature,
            log_reaction_quotient=ln_k,
        )
        self.assertTrue(math.isclose(delta_g, 0.0, abs_tol=1e-12))

    def test_vant_hoff_primitives(self):
        t1 = Temperature(value=298.15, unit="K")
        t2 = Temperature(value=350.0, unit="K")
        delta_h = 40000.0

        dlnk_dt = calc_dlnK_dT(delta_h, t1)
        self.assertTrue(math.isclose(
            dlnk_dt,
            delta_h / (8.314462618 * 298.15 ** 2),
            rel_tol=1e-12,
        ))

        k2 = calc_equilibrium_constant_at_temperature(2.0, delta_h, t1, t2)
        expected = 2.0 * math.exp(-(delta_h / 8.314462618) * (1.0 / 350.0 - 1.0 / 298.15))
        self.assertTrue(math.isclose(k2, expected, rel_tol=1e-12))

    def test_chemical_potential(self):
        result = calc_chemical_potential_from_activity(
            -1000.0,
            1.0,
            Temperature(value=298.15, unit="K"),
        )
        self.assertEqual(result, -1000.0)

    def test_ideal_mixing(self):
        temperature = Temperature(value=300.0, unit="K")
        self.assertEqual(calc_ideal_molar_entropy_of_mixing([1.0]), 0.0)
        self.assertEqual(calc_ideal_molar_gibbs_energy_of_mixing([1.0], temperature), 0.0)
        total_g = calc_ideal_gibbs_energy_of_mixing(2.0, [0.5, 0.5], temperature)
        self.assertTrue(total_g < 0.0)

    def test_ionic_strength_and_charge_balance(self):
        self.assertEqual(calc_ionic_strength_molality([0.1, 0.1], [1.0, -1.0]), 0.1)
        self.assertEqual(calc_ionic_strength_molarity({"Na+": 0.1, "Cl-": 0.1}, {"Na+": 1.0, "Cl-": -1.0}), 0.1)
        self.assertEqual(calc_charge_balance([0.1, 0.1], [1.0, -1.0]), 0.0)
        self.assertTrue(check_electroneutrality([0.1, 0.1], [1.0, -1.0]))

    def test_z_based_gas_properties_reduce_to_ideal_gas_at_z_one(self):
        temperature = Temperature(value=300.0, unit="K")
        pressure = Pressure(value=101325.0, unit="Pa")
        molecular_weight = CustomProp(value=0.02897, unit="kg/mol")

        ideal_density = calc_ideal_gas_density(pressure, molecular_weight, temperature)
        z_density = calc_gas_density_from_z(pressure, molecular_weight, temperature, 1.0)
        ideal_volume = calc_ideal_gas_molar_volume(temperature, pressure)
        z_volume = calc_gas_molar_volume_from_z(temperature, pressure, 1.0)

        self.assertIsNotNone(ideal_density)
        self.assertIsNotNone(ideal_volume)
        self.assertTrue(math.isclose(z_density.value, ideal_density.value, rel_tol=1e-12))
        self.assertTrue(math.isclose(z_volume.value, ideal_volume.value, rel_tol=1e-12))

    def test_heat_capacity_helpers(self):
        cp = 29.1
        cv = calc_ideal_gas_cv_from_cp(cp)
        self.assertTrue(math.isclose(cp - cv, 8.31446261815324, rel_tol=1e-12))
        self.assertTrue(math.isclose(calc_ideal_gas_cp_from_cv(cv), cp, rel_tol=1e-12))
        self.assertTrue(math.isclose(calc_heat_capacity_ratio(cp, cv), cp / cv, rel_tol=1e-12))

    def test_normality(self):
        self.assertEqual(calc_normality(0.25, 2.0), 0.5)
        with self.assertRaises(ValueError):
            calc_normality(0.25, 0.0)


if __name__ == "__main__":
    unittest.main()
