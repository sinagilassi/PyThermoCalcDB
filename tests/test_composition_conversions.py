import math
import unittest

from pythermodb_settings.models import Component, CustomProp

from pythermocalcdb.compositions.conversions import (
    mass_concentration_to_molarity,
    mass_fraction_to_molality,
    mass_fraction_to_mole_fraction,
    mass_fraction_to_ppm,
    mass_fraction_to_weight_percent,
    molality_to_molarity,
    molality_to_mass_fraction,
    molality_to_mole_fraction,
    molarity_to_mass_concentration,
    molarity_to_mass_fraction,
    molarity_to_molality,
    molarities_to_molalities,
    mole_fraction_to_mass_fraction,
    mole_fraction_to_molality,
    mole_fraction_to_mole_percent,
    mole_fraction_to_ppb,
    ppm_mass_to_mass_fraction,
    ppb_mole_to_mole_fraction,
    weight_percent_to_mass_fraction,
)


class TestCompositionConversions(unittest.TestCase):
    def assert_close_sequence(self, values, expected):
        self.assertEqual(len(values), len(expected))
        for value, expected_value in zip(values, expected):
            self.assertTrue(math.isclose(value, expected_value, rel_tol=1e-12))

    def test_mole_fraction_to_mass_fraction_sequence(self):
        result = mole_fraction_to_mass_fraction([0.5, 0.5], [0.018015, 0.04607])
        self.assert_close_sequence(result, [0.28111102442069125, 0.7188889755793086])

    def test_mass_fraction_to_mole_fraction_mapping(self):
        result = mass_fraction_to_mole_fraction(
            {"water": 0.25, "ethanol": 0.75},
            {"water": 0.018015, "ethanol": 0.04607},
        )
        self.assertTrue(math.isclose(sum(result.values()), 1.0, rel_tol=1e-12))
        self.assertTrue(math.isclose(result["water"], 0.46017080357588774, rel_tol=1e-12))

    def test_molarity_molality_round_trip(self):
        molarity = 1.0
        molecular_weight = 0.05844
        density = 1.02
        molality = molarity_to_molality(molarity, molecular_weight, density)
        self.assertTrue(math.isclose(
            molality_to_molarity(molality, molecular_weight, density),
            molarity,
            rel_tol=1e-12,
        ))

    def test_multisolute_molarity_to_molality_mapping(self):
        result = molarities_to_molalities(
            {"a": 1.0, "b": 0.5},
            {"a": 0.04, "b": 0.06},
            1.05,
        )
        self.assertTrue(math.isclose(result["a"], 1.0204081632653061, rel_tol=1e-12))
        self.assertTrue(math.isclose(result["b"], 0.5102040816326531, rel_tol=1e-12))

    def test_molality_to_mole_fraction_includes_solvent(self):
        result = molality_to_mole_fraction({"NaCl": 1.0}, 0.01801528)
        self.assertTrue(math.isclose(sum(result.values()), 1.0, rel_tol=1e-12))
        self.assertTrue(math.isclose(result["NaCl"], 0.01769647308240796, rel_tol=1e-12))
        self.assertIn("solvent", result)

    def test_scalar_conversion_round_trips(self):
        self.assertTrue(math.isclose(
            mass_fraction_to_molality(molality_to_mass_fraction(2.0, 0.05844), 0.05844),
            2.0,
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(
            mass_concentration_to_molarity(molarity_to_mass_concentration(3.0, 0.02), 0.02),
            3.0,
            rel_tol=1e-12,
        ))
        self.assertTrue(math.isclose(weight_percent_to_mass_fraction(mass_fraction_to_weight_percent(0.12)), 0.12))
        self.assertTrue(math.isclose(ppm_mass_to_mass_fraction(mass_fraction_to_ppm(1e-4)), 1e-4))
        self.assertTrue(math.isclose(ppb_mole_to_mole_fraction(mole_fraction_to_ppb(1e-6)), 1e-6))
        self.assertTrue(math.isclose(mole_fraction_to_molality(0.02, 0.98, 0.01801528), 1.1328252053426935))
        self.assertTrue(math.isclose(mole_fraction_to_mole_percent(0.25), 25.0))

    def test_component_amounts_with_custom_prop_molecular_weights(self):
        result = mole_fraction_to_mass_fraction(
            {"water": 0.5, "ethanol": 0.5},
            {
                "water": CustomProp(value=18.015, unit="g/mol"),
                "ethanol": CustomProp(value=46.07, unit="g/mol"),
            },
            output_molecular_weight_unit="kg/mol",
        )
        self.assertTrue(math.isclose(sum(result.values()), 1.0, rel_tol=1e-12))
        self.assertTrue(math.isclose(result["water"], 0.28111102442069125, rel_tol=1e-12))

    def test_scalar_custom_prop_units(self):
        result = molarity_to_mass_fraction(
            molarity=CustomProp(value=1.0, unit="mol/L"),
            molecular_weight=CustomProp(value=58.44, unit="g/mol"),
            solution_density=CustomProp(value=1020.0, unit="kg/m^3"),
            output_molarity_unit="mol/L",
            output_molecular_weight_unit="kg/mol",
            output_solution_density_unit="kg/L",
        )
        self.assertTrue(math.isclose(result, 0.05729411764705882, rel_tol=1e-12))

    def test_component_mapping_orders_conversion_results(self):
        components = [
            Component(name="ethanol", formula="C2H6O", state="l"),
            Component(name="water", formula="H2O", state="l"),
        ]
        result = mole_fraction_to_mass_fraction(
            {"water": 0.5, "ethanol": 0.5},
            {"water": 18.015, "ethanol": 46.07},
            components=components,
            component_key="Formula",
        )
        self.assertEqual(list(result), ["C2H6O", "H2O"])


    def test_custom_prop_values_are_used_as_is_without_output_units(self):
        result = molarity_to_mass_fraction(
            molarity=CustomProp(value=1.0, unit="mol/L"),
            molecular_weight=CustomProp(value=0.05844, unit="kg/mol"),
            solution_density=CustomProp(value=1.02, unit="kg/L"),
        )
        self.assertTrue(math.isclose(result, 0.05729411764705882, rel_tol=1e-12))

    def test_invalid_inputs_raise(self):
        with self.assertRaises(ValueError):
            mole_fraction_to_mass_fraction([0.5, -0.5], [0.018, 0.046])
        with self.assertRaises(ValueError):
            molarity_to_molality(30.0, 0.05844, 1.0)


if __name__ == "__main__":
    unittest.main()


