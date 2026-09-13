import math
import unittest

import numpy as np
from pythermodb_settings.models import CustomProp, Temperature

from pythermocalcdb.compositions import (
    calc_osmolality,
    calc_osmolality_from_mapping,
    calc_osmolality_from_props,
    calc_osmolarity,
    calc_osmolarity_from_mapping,
    calc_osmolarity_from_props,
)
from pythermocalcdb.compositions.core.osmotic import _calc_osmolarity
from pythermocalcdb.mixtures import calc_ideal_enthalpy_of_mixing
from pythermocalcdb.reactions import (
    classify_reaction_direction,
    classify_reaction_direction_from_logs,
)
from pythermocalcdb.thermo import (
    calc_clapeyron_slope,
    calc_enthalpy_of_sublimation,
    calc_ideal_gas_isentropic_temperature,
    calc_isothermal_compressibility,
    calc_joule_thomson_coefficient,
    calc_phase_transition_entropy,
    calc_thermal_expansion_coefficient,
)
from pythermocalcdb.thermo.core.derivatives import (
    _calc_joule_thomson_coefficient,
)
from pythermocalcdb.thermo.core.phase_change import _calc_clapeyron_slope
from pythermocalcdb.thermo.vapor_pressure import (
    calc_enthalpy_vaporization_clausius_clapeyron,
    calc_log_vapor_pressure_ratio_clausius_clapeyron,
    calc_vapor_pressure_clausius_clapeyron,
)


class TestRemainingCalculationsV3(unittest.TestCase):
    def test_phase_change_reference_values(self):
        self.assertEqual(calc_phase_transition_entropy(6000.0, 300.0), 20.0)
        self.assertEqual(calc_enthalpy_of_sublimation(10000.0, 40000.0), 50000.0)
        self.assertTrue(math.isclose(
            calc_clapeyron_slope(40000.0, Temperature(value=373.15, unit="K"), 0.030),
            40000.0 / (373.15 * 0.030),
            rel_tol=1e-12,
        ))

    def test_phase_change_core_vectorized_and_invalid(self):
        result = _calc_clapeyron_slope([40000.0, 20000.0], 300.0, [0.02, 0.01])
        self.assertTrue(np.allclose(result, [6666.666666666667, 6666.666666666667]))
        with self.assertRaises(ValueError):
            _calc_clapeyron_slope(40000.0, 300.0, 0.0)

    def test_derivative_property_reference_values(self):
        self.assertEqual(calc_thermal_expansion_coefficient(1.0, 0.003), 0.003)
        self.assertEqual(calc_isothermal_compressibility(1.0, -1e-9), 1e-9)
        self.assertTrue(math.isclose(
            calc_joule_thomson_coefficient(300.0, 0.024, 30.0, 1.0e-4),
            (300.0 * 1.0e-4 - 0.024) / 30.0,
            rel_tol=1e-12,
        ))

    def test_derivative_core_vectorized_and_invalid(self):
        result = _calc_joule_thomson_coefficient(
            [300.0, 400.0],
            [0.024, 0.030],
            30.0,
            [1.0e-4, 1.0e-4],
        )
        self.assertTrue(np.allclose(result, [0.0002, 0.0003333333333333334]))
        with self.assertRaises(ValueError):
            _calc_joule_thomson_coefficient(300.0, 0.024, 0.0, 1.0e-4)

    def test_clausius_clapeyron_helpers(self):
        t1 = Temperature(value=350.0, unit="K")
        t2 = Temperature(value=360.0, unit="K")
        dh = 40000.0
        ln_ratio = calc_log_vapor_pressure_ratio_clausius_clapeyron(dh, t1, t2)
        expected = -(dh / 8.314462618) * (1.0 / 360.0 - 1.0 / 350.0)
        self.assertTrue(math.isclose(ln_ratio, expected, rel_tol=1e-12))
        p2 = calc_vapor_pressure_clausius_clapeyron(100000.0, dh, t1, t2)
        self.assertTrue(math.isclose(p2, 100000.0 * math.exp(expected), rel_tol=1e-12))
        inferred = calc_enthalpy_vaporization_clausius_clapeyron(100000.0, p2, t1, t2)
        self.assertTrue(math.isclose(inferred, dh, rel_tol=1e-12))

    def test_reaction_direction_classifier(self):
        self.assertEqual(classify_reaction_direction_from_logs(0.0, 1.0), "forward")
        self.assertEqual(classify_reaction_direction_from_logs(1.0, 1.0), "equilibrium")
        self.assertEqual(classify_reaction_direction_from_logs(2.0, 1.0), "reverse")
        self.assertEqual(classify_reaction_direction(0.5, 1.0), "forward")
        self.assertEqual(classify_reaction_direction(2.0, 1.0), "reverse")

    def test_osmotic_primitives(self):
        self.assertEqual(calc_osmolarity([0.1, 0.1], unit="mol/L").value, 0.2)
        self.assertEqual(calc_osmolarity_from_mapping({"Na+": 0.1, "Cl-": 0.1}).value, 0.2)
        self.assertEqual(calc_osmolality([0.2, 0.3], unit="mol/kg").value, 0.5)
        self.assertEqual(calc_osmolality_from_mapping({"a": 0.2, "b": 0.3}).value, 0.5)
        self.assertTrue(np.allclose(_calc_osmolarity([[0.1, 0.1], [0.2, 0.3]]), [0.2, 0.5]))

    def test_osmotic_props(self):
        def convert(value, from_unit, to_unit):
            if from_unit == "mmol/L" and to_unit == "mol/L":
                return value / 1000.0
            if from_unit == "mmol/kg" and to_unit == "mol/kg":
                return value / 1000.0
            return value

        osmolarity = calc_osmolarity_from_props(
            {
                "Na+": CustomProp(value=100.0, unit="mmol/L"),
                "Cl-": CustomProp(value=100.0, unit="mmol/L"),
            },
            output_unit="mol/L",
            unit_conversion_fn=convert,
        )
        self.assertEqual(osmolarity.value, 0.2)
        osmolality = calc_osmolality_from_props(
            {
                "a": CustomProp(value=200.0, unit="mmol/kg"),
                "b": CustomProp(value=300.0, unit="mmol/kg"),
            },
            output_unit="mol/kg",
            unit_conversion_fn=convert,
        )
        self.assertEqual(osmolality.value, 0.5)
        converted = calc_osmolarity_from_props(
            {"a": CustomProp(value=100.0, unit="mmol/L")},
            output_unit="mol/L",
            unit_conversion_fn=convert,
        )
        self.assertEqual(converted.value, convert(100.0, "mmol/L", "mol/L"))

    def test_low_priority_symmetry_helpers(self):
        self.assertEqual(calc_ideal_enthalpy_of_mixing([0.25, 0.75]), 0.0)
        self.assertTrue(math.isclose(
            calc_ideal_gas_isentropic_temperature(300.0, 100000.0, 200000.0, 1.4),
            300.0 * (2.0 ** ((1.4 - 1.0) / 1.4)),
            rel_tol=1e-12,
        ))


if __name__ == "__main__":
    unittest.main()
