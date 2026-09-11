import math
import unittest

import numpy as np

from pythermocalcdb.mixtures import (
    calc_ideal_mixture_density,
    calc_ideal_mixture_heat_capacity,
    calc_mixture_molecular_weight_from_mass_fractions,
    calc_mixture_molecular_weight_from_mole_fractions,
    mass_fraction_to_volume_fraction,
)
from pythermocalcdb.mixtures.core import (
    _calc_ideal_mixture_density,
    _calc_ideal_mixture_density_from_props,
    _calc_ideal_entropy_of_mixing_from_props,
    _calc_ideal_gibbs_energy_of_mixing_from_props,
    _calc_ideal_mixture_heat_capacity,
    _calc_ideal_mixture_heat_capacity_from_props,
    _calc_ideal_molar_gibbs_energy_of_mixing_from_props,
    _calc_ideal_molar_entropy_of_mixing,
    _calc_ideal_molar_entropy_of_mixing_from_props,
    _calc_mass_fraction_to_volume_fraction,
    _calc_mass_fraction_to_volume_fraction_from_props,
    _calc_mixture_molecular_weight_from_mass_fractions_from_props,
    _calc_mixture_molecular_weight_from_mass_fractions,
    _calc_mixture_molecular_weight_from_mole_fractions_from_props,
    _calc_mixture_molecular_weight_from_mole_fractions,
    _calc_volume_fractions,
    _calc_volume_fractions_from_props,
)
from pythermodb_settings.models import CustomProp, Temperature


class TestMixtureCoreApi(unittest.TestCase):
    def test_core_mixture_rules_support_2d_state_arrays(self):
        fractions = np.array([[0.25, 0.75], [0.5, 0.5]])
        properties = np.array([[20.0, 40.0], [10.0, 30.0]])

        np.testing.assert_allclose(
            _calc_ideal_mixture_heat_capacity(fractions, properties),
            np.array([35.0, 20.0]),
        )
        np.testing.assert_allclose(
            _calc_ideal_mixture_density(fractions, np.array([[800.0, 1000.0], [10.0, 20.0]])),
            np.array([941.1764705882354, 13.333333333333334]),
        )
        np.testing.assert_allclose(
            _calc_mixture_molecular_weight_from_mole_fractions(fractions, properties),
            np.array([35.0, 20.0]),
        )
        np.testing.assert_allclose(
            _calc_mixture_molecular_weight_from_mass_fractions(fractions, properties),
            np.array([32.0, 15.0]),
        )

    def test_core_volume_fraction_rules_support_2d_state_arrays(self):
        np.testing.assert_allclose(
            _calc_volume_fractions(np.array([[1.0, 3.0], [2.0, 2.0]])),
            np.array([[0.25, 0.75], [0.5, 0.5]]),
        )
        np.testing.assert_allclose(
            _calc_mass_fraction_to_volume_fraction(
                np.array([[0.25, 0.75], [0.5, 0.5]]),
                np.array([[800.0, 1000.0], [10.0, 20.0]]),
            ),
            np.array([[0.29411764705882354, 0.7058823529411765], [2 / 3, 1 / 3]]),
        )

    def test_core_entropy_uses_zero_fraction_limit(self):
        self.assertTrue(math.isclose(
            _calc_ideal_molar_entropy_of_mixing([0.0, 1.0]),
            0.0,
            abs_tol=1e-12,
        ))

    def test_core_props_adapters_delegate_unit_and_mapping_work(self):
        mole_fractions = {
            "a": CustomProp(value=0.5, unit=""),
            "b": CustomProp(value=0.5, unit=""),
        }
        mass_fractions = {
            "a": CustomProp(value=0.25, unit=""),
            "b": CustomProp(value=0.75, unit=""),
        }
        densities = {
            "a": CustomProp(value=800.0, unit="kg/m^3"),
            "b": CustomProp(value=1000.0, unit="kg/m^3"),
        }
        properties = {
            "a": CustomProp(value=20.0, unit="J/mol/K"),
            "b": CustomProp(value=40.0, unit="J/mol/K"),
        }
        total_moles = CustomProp(value=2.0, unit="mol")

        self.assertTrue(math.isclose(
            _calc_ideal_mixture_density_from_props(mass_fractions, densities),
            941.1764705882354,
        ))
        self.assertTrue(math.isclose(
            _calc_ideal_mixture_heat_capacity_from_props(mass_fractions, properties),
            35.0,
        ))
        self.assertTrue(math.isclose(
            _calc_mixture_molecular_weight_from_mole_fractions_from_props(mole_fractions, properties),
            30.0,
        ))
        self.assertTrue(math.isclose(
            _calc_mixture_molecular_weight_from_mass_fractions_from_props(mass_fractions, properties),
            32.0,
        ))
        self.assertTrue(math.isclose(
            _calc_ideal_molar_entropy_of_mixing_from_props(mole_fractions),
            5.763146321643829,
        ))
        self.assertTrue(math.isclose(
            _calc_ideal_entropy_of_mixing_from_props(total_moles, mole_fractions),
            11.526292643287658,
        ))
        self.assertTrue(math.isclose(
            _calc_ideal_molar_gibbs_energy_of_mixing_from_props(
                mole_fractions,
                Temperature(value=300.0, unit="K"),
            ),
            -1728.9438964931485,
        ))
        self.assertTrue(math.isclose(
            _calc_ideal_gibbs_energy_of_mixing_from_props(
                total_moles,
                mole_fractions,
                Temperature(value=300.0, unit="K"),
            ),
            -3457.887792986297,
        ))
        self.assertEqual(
            _calc_volume_fractions_from_props({
                "a": CustomProp(value=1.0, unit="L"),
                "b": CustomProp(value=3.0, unit="L"),
            }),
            {"a": 0.25, "b": 0.75},
        )
        volume_fractions = _calc_mass_fraction_to_volume_fraction_from_props(
            mass_fractions,
            densities,
        )
        self.assertTrue(math.isclose(volume_fractions["a"], 0.29411764705882354))
        self.assertTrue(math.isclose(volume_fractions["b"], 0.7058823529411765))

    def test_core_props_adapters_reject_numeric_mappings(self):
        with self.assertRaises(TypeError):
            _calc_ideal_mixture_density_from_props({"a": 0.5}, {"a": 1000.0})
        with self.assertRaises(TypeError):
            _calc_ideal_mixture_heat_capacity_from_props({"a": 1.0}, {"a": 20.0})
        with self.assertRaises(TypeError):
            _calc_mixture_molecular_weight_from_mole_fractions_from_props({"a": 1.0}, {"a": 20.0})
        with self.assertRaises(TypeError):
            _calc_volume_fractions_from_props({"a": 1.0})
        with self.assertRaises(TypeError):
            _calc_ideal_entropy_of_mixing_from_props(1.0, {"a": CustomProp(value=1.0, unit="")})

    def test_public_mapping_inputs_validate_matching_keys(self):
        with self.assertRaises(ValueError):
            calc_ideal_mixture_heat_capacity({"a": 0.5, "b": 0.5}, {"a": 10.0, "c": 20.0})
        with self.assertRaises(ValueError):
            calc_ideal_mixture_density({"a": 0.5, "b": 0.5}, {"a": 10.0, "c": 20.0})
        with self.assertRaises(ValueError):
            calc_mixture_molecular_weight_from_mole_fractions(
                {"a": 0.5, "b": 0.5},
                {"a": 10.0, "c": 20.0},
            )
        with self.assertRaises(ValueError):
            calc_mixture_molecular_weight_from_mass_fractions(
                {"a": 0.5, "b": 0.5},
                {"a": 10.0, "c": 20.0},
            )
        with self.assertRaises(ValueError):
            mass_fraction_to_volume_fraction({"a": 0.5, "b": 0.5}, {"a": 10.0, "c": 20.0})


if __name__ == "__main__":
    unittest.main()
