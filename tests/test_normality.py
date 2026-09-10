import unittest

import numpy as np

from pythermocalcdb.compositions.normality import (
    _calc_normality,
    _calc_normality_from_props,
    calc_normality,
    calc_normality_from_sequence,
    calc_normality_from_props,
)
from pythermodb_settings.models import CustomProp


class TestNormality(unittest.TestCase):
    def test_calc_normality_returns_raw_float(self):
        self.assertEqual(_calc_normality(0.25, 2.0), 0.5)

    def test_calc_normality_returns_array_for_array_like_inputs(self):
        result = _calc_normality([0.25, 0.5], 2.0)

        np.testing.assert_allclose(result, np.array([0.5, 1.0]))
        self.assertIsInstance(result, np.ndarray)

    def test_calc_normality_can_return_list(self):
        result = _calc_normality([0.25, 0.5], 2.0, as_list=True)

        self.assertEqual(result, [0.5, 1.0])
        self.assertIsInstance(result, list)

    def test_calc_normality_rejects_non_positive_inputs(self):
        with self.assertRaises(ValueError):
            _calc_normality(0.0, 2.0)

        with self.assertRaises(ValueError):
            _calc_normality(0.25, 0.0)

    def test_annotated_normality_matches_composition_api_pattern(self):
        result = calc_normality(0.25, 2.0)

        self.assertEqual(result.value, 0.5)
        self.assertEqual(result.name, "normality")
        self.assertIsNone(result.unit)

    def test_annotated_normality_returns_array_for_array_like_inputs(self):
        result = calc_normality([0.25, 0.5], 2.0)

        np.testing.assert_allclose(result.value, np.array([0.5, 1.0]))
        self.assertIsInstance(result.value, np.ndarray)

    def test_annotated_normality_from_sequence_returns_list(self):
        result = calc_normality_from_sequence([0.25, 0.5], 2.0)

        self.assertEqual(result.value, [0.5, 1.0])
        self.assertIsInstance(result.value, list)
        self.assertEqual(result.name, "normality")

    def test_annotated_normality_allows_explicit_numeric_unit(self):
        result = calc_normality(0.25, 2.0, unit="eq/L")

        self.assertEqual(result.value, 0.5)
        self.assertEqual(result.unit, "eq/L")

    def test_calc_normality_from_props_returns_raw_float(self):
        result = _calc_normality_from_props(
            molarity=CustomProp(value=0.25, unit="mol/L"),
            equivalence_factor=2.0,
            output_unit="eq/L",
        )

        self.assertEqual(result, 0.5)

    def test_annotated_normality_from_props_sets_output_unit(self):
        result = calc_normality_from_props(
            molarity=CustomProp(value=0.25, unit="mol/L"),
            equivalence_factor=2.0,
            output_unit="eq/L",
        )

        self.assertEqual(result.value, 0.5)
        self.assertEqual(result.unit, "eq/L")

    def test_annotated_normality_from_props_rejects_mismatched_unit(self):
        with self.assertRaisesRegex(ValueError, "Mismatch"):
            calc_normality_from_props(
                molarity=CustomProp(value=0.25, unit="mol/L"),
                equivalence_factor=2.0,
                output_unit="eq/L",
                unit="eq/m3",
            )

    def test_custom_prop_molarity_uses_output_unit_denominator(self):
        result = calc_normality_from_props(
            molarity=CustomProp(value=0.25, unit="mol/L"),
            equivalence_factor=2.0,
            output_unit="eq/m3",
        )

        self.assertEqual(result.value, 500.0)
        self.assertEqual(result.unit, "eq/m3")


if __name__ == "__main__":
    unittest.main()
