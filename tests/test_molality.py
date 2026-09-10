import unittest

import numpy as np

from pythermocalcdb.compositions.molality import (
    _calc_molalities,
    calc_component_molalities_from_props,
    calc_molalities,
    calc_molalities_from_sequence,
    calc_molalities_from_mapping,
)
from pythermodb_settings.models import Component, CustomProp


class TestMolality(unittest.TestCase):
    def test_scalar_mass_returns_component_moles_shape(self):
        result = _calc_molalities(
            component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
            solvent_mass=2.0,
        )

        np.testing.assert_allclose(result, np.array([[0.5, 1.0], [1.5, 2.0]]))
        self.assertEqual(result.shape, (2, 2))

    def test_sequence_can_return_list(self):
        result = _calc_molalities(
            component_moles=[1.0, 2.0],
            solvent_mass=2.0,
            as_list=True,
        )

        self.assertEqual(result, [0.5, 1.0])
        self.assertIsInstance(result, list)

    def test_annotated_sequence_can_return_list(self):
        result = calc_molalities(
            component_moles=[1.0, 2.0],
            solvent_mass=2.0,
            as_list=True,
        )

        self.assertEqual(result.value, [0.5, 1.0])
        self.assertIsInstance(result.value, list)
        self.assertEqual(result.name, "molality")

    def test_sequence_api_returns_list(self):
        result = calc_molalities_from_sequence(
            component_moles=[1.0, 2.0],
            solvent_mass=2.0,
        )

        self.assertEqual(result.value, [0.5, 1.0])
        self.assertIsInstance(result.value, list)
        self.assertEqual(result.name, "molality")

    def test_matching_mass_shape_is_allowed(self):
        result = _calc_molalities(
            component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
            solvent_mass=np.array([[1.0, 2.0], [3.0, 4.0]]),
        )

        np.testing.assert_allclose(result, np.ones((2, 2)))

    def test_state_mass_shape_is_allowed(self):
        result = _calc_molalities(
            component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
            solvent_mass=np.array([1.0, 2.0]),
        )

        np.testing.assert_allclose(result, np.array([[1.0, 2.0], [1.5, 2.0]]))
        self.assertEqual(result.shape, (2, 2))

    def test_non_broadcastable_mass_shape_raises_validation_error(self):
        with self.assertRaisesRegex(ValueError, "one entry per state"):
            _calc_molalities(
                component_moles=np.array([[1.0, 2.0], [3.0, 4.0]]),
                solvent_mass=np.array([1.0, 2.0, 3.0]),
            )

    def test_non_positive_mass_raises_validation_error(self):
        with self.assertRaisesRegex(ValueError, "greater than zero"):
            _calc_molalities(
                component_moles=np.array([1.0, 2.0]),
                solvent_mass=0.0,
            )

    def test_annotated_mapping_api_matches_molarity_pattern(self):
        result = calc_molalities_from_mapping(
            component_moles={"A": 2.0, "B": 3.0},
            solvent_mass=10.0,
        )

        self.assertEqual(result.value, {"A": 0.2, "B": 0.3})
        self.assertEqual(result.name, "molality")

    def test_component_props_api_matches_molarity_pattern(self):
        components = [
            Component(name="carbon dioxide", formula="CO2", state="g"),
            Component(name="methane", formula="CH4", state="g"),
            Component(name="oxygen", formula="O2", state="g"),
        ]

        result = calc_component_molalities_from_props(
            component_moles={
                "A": CustomProp(value=2.0, unit="mol"),
                "B": CustomProp(value=3.0, unit="mol"),
            },
            solvent_mass=CustomProp(value=10.0, unit="kg"),
        )

        self.assertEqual(result.value, {"A": 0.2, "B": 0.3})
        self.assertEqual(result.unit, "mol/kg")

        mapped = calc_component_molalities_from_props(
            component_moles={
                "oxygen": CustomProp(value=2.0, unit="mol"),
                "carbon dioxide": CustomProp(value=1.0, unit="mol"),
                "methane": CustomProp(value=1.0, unit="mol"),
            },
            solvent_mass=CustomProp(value=2.0, unit="kg"),
            components=components,
            component_key="Formula-State",
            case_sensitive=False,
            sort_by_components_order=True,
        )

        self.assertEqual(mapped.value, {"CO2-g": 0.5, "CH4-g": 0.5, "O2-g": 1.0})


if __name__ == "__main__":
    unittest.main()
