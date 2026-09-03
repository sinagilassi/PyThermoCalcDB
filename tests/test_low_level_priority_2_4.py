import math
import unittest

from pythermodb_settings.models import CustomProp, Pressure, Temperature

from pythermocalcdb.basics.mixtures import (
    calc_mixture_molecular_weight_from_mass_fractions_1,
)
from pythermocalcdb.conversions import (
    mass_specific_to_molar,
    molar_cp_to_mass_cp,
    molar_property_to_total,
    molar_to_mass_specific,
    specific_property_to_total,
    total_to_molar_property,
    total_to_specific_property,
)
from pythermocalcdb.mixtures import (
    calc_ideal_mixture_density,
    calc_ideal_mixture_heat_capacity,
    calc_volume_fractions,
    mass_fraction_to_volume_fraction,
)
from pythermocalcdb.thermo import (
    calc_gibbs_energy,
    calc_gibbs_energy_change,
    calc_helmholtz_energy,
    calc_ideal_gas_density,
    calc_ideal_gas_internal_energy,
    calc_ideal_gas_molar_volume,
    calc_internal_energy,
    density_to_specific_volume,
    specific_volume_to_density,
)


class TestLowLevelPriorityTwoToFour(unittest.TestCase):
    def test_property_basis_round_trip(self):
        cp_mass = molar_to_mass_specific(75.3, 0.01801528)
        self.assertTrue(math.isclose(cp_mass, 4179.785160152936, rel_tol=1e-12))
        self.assertTrue(math.isclose(mass_specific_to_molar(cp_mass, 0.01801528), 75.3, rel_tol=1e-12))
        self.assertTrue(math.isclose(molar_cp_to_mass_cp(75.3, 0.01801528), cp_mass, rel_tol=1e-12))

    def test_extensive_intensive_round_trip(self):
        total = molar_property_to_total(2.0, 100.0)
        self.assertEqual(total, 200.0)
        self.assertEqual(total_to_molar_property(total, 2.0), 100.0)
        total_mass = specific_property_to_total(4.0, 10.0)
        self.assertEqual(total_mass, 40.0)
        self.assertEqual(total_to_specific_property(total_mass, 4.0), 10.0)

    def test_specific_volume_round_trip(self):
        specific_volume = density_to_specific_volume(1000.0)
        self.assertEqual(specific_volume, 0.001)
        self.assertEqual(specific_volume_to_density(specific_volume), 1000.0)

    def test_thermo_identities(self):
        temperature = Temperature(value=300.0, unit="K")
        self.assertEqual(calc_gibbs_energy(10000.0, temperature, 20.0), 4000.0)
        self.assertEqual(calc_gibbs_energy_change(5000.0, 10.0, temperature), 2000.0)
        self.assertEqual(calc_internal_energy(1000.0, 100000.0, 0.002), 800.0)
        self.assertTrue(math.isclose(calc_ideal_gas_internal_energy(10000.0, temperature), 7505.661214554028))
        self.assertEqual(calc_helmholtz_energy(8000.0, temperature, 10.0), 5000.0)

        # NOTE: With output_temperature_unit=None, temperature.value is used as-is.
        c_temperature = Temperature(value=26.85, unit="C")
        self.assertEqual(calc_gibbs_energy(10000.0, c_temperature, 20.0), 9463.0)
        self.assertEqual(calc_gibbs_energy_change(5000.0, 10.0, c_temperature), 4731.5)
        self.assertTrue(math.isclose(
            calc_ideal_gas_internal_energy(10000.0, c_temperature),
            10000.0 - 8.31446261815324 * 26.85,
        ))
        self.assertEqual(calc_helmholtz_energy(8000.0, c_temperature, 10.0), 7731.5)

        # NOTE: Explicit output_temperature_unit converts temperature before T*S or R*T.
        self.assertTrue(math.isclose(
            calc_gibbs_energy(
                10000.0,
                c_temperature,
                20.0,
                output_temperature_unit="K",
            ),
            4000.0,
        ))
        self.assertTrue(math.isclose(
            calc_gibbs_energy_change(
                5000.0,
                10.0,
                c_temperature,
                output_temperature_unit="K",
            ),
            2000.0,
        ))
        self.assertTrue(math.isclose(
            calc_ideal_gas_internal_energy(
                10000.0,
                c_temperature,
                output_temperature_unit="K",
            ),
            7505.661214554028,
        ))
        self.assertTrue(math.isclose(
            calc_helmholtz_energy(
                8000.0,
                c_temperature,
                10.0,
                output_temperature_unit="K",
            ),
            5000.0,
        ))
        with self.assertRaises(ValueError):
            calc_gibbs_energy(8000.0, Temperature(value=-300.0, unit="K"), 10.0)
        with self.assertRaises(ValueError):
            calc_ideal_gas_internal_energy(8000.0, Temperature(value=-300.0, unit="K"))
        with self.assertRaises(ValueError):
            calc_helmholtz_energy(8000.0, Temperature(value=-300.0, unit="K"), 10.0)

    def test_ideal_gas_properties(self):
        temperature = Temperature(value=300.0, unit="K")
        pressure = Pressure(value=101325.0, unit="Pa")
        molecular_weight = CustomProp(value=0.02897, unit="kg/mol")
        density = calc_ideal_gas_density(pressure, molecular_weight, temperature)
        molar_volume = calc_ideal_gas_molar_volume(temperature, pressure)
        self.assertIsNotNone(density)
        self.assertIsNotNone(molar_volume)
        self.assertTrue(math.isclose(density.value, 1.1768189899172703, rel_tol=1e-12))
        self.assertTrue(math.isclose(molar_volume.value, 0.024617209824287906, rel_tol=1e-12))

    def test_ideal_mixture_rules(self):
        self.assertTrue(math.isclose(calc_ideal_mixture_heat_capacity([0.25, 0.75], [20.0, 40.0]), 35.0))
        self.assertTrue(math.isclose(calc_ideal_mixture_density([0.25, 0.75], [800.0, 1000.0]), 941.1764705882354))
        self.assertEqual(calc_volume_fractions({"a": 2.0, "b": 3.0}), {"a": 0.4, "b": 0.6})
        result = mass_fraction_to_volume_fraction([0.25, 0.75], [800.0, 1000.0])
        self.assertTrue(math.isclose(sum(result), 1.0, rel_tol=1e-12))
        self.assertTrue(math.isclose(result[0], 0.29411764705882354, rel_tol=1e-12))

    def test_mixture_molecular_weight_from_mass_fractions(self):
        result = calc_mixture_molecular_weight_from_mass_fractions_1([0.25, 0.75], [18.015, 46.07])
        self.assertTrue(math.isclose(result, 33.15990810567847, rel_tol=1e-12))


if __name__ == "__main__":
    unittest.main()






