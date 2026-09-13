import unittest

from pythermodb_settings.models import CustomProp, Temperature

from pythermocalcdb.thermo import (
    calc_heat_capacity_ratio,
    calc_helmholtz_energy,
    calc_ideal_gas_cp_from_cv,
    calc_ideal_gas_cv_from_cp,
)
from pythermocalcdb.thermo.core.heat_capacity import (
    _calc_heat_capacity_ratio,
    _calc_heat_capacity_ratio_from_props,
    _calc_ideal_gas_cp_from_cv,
    _calc_ideal_gas_cv_from_cp,
    _calc_ideal_gas_cv_from_cp_from_props,
)
from pythermocalcdb.thermo.core.helmholtz import (
    _calc_helmholtz_energy,
    _calc_helmholtz_energy_from_props,
)


class TestThermoCoreApi(unittest.TestCase):
    def test_public_heat_capacity_wrappers_delegate_to_core_adapters(self):
        cp = 29.0
        cv = calc_ideal_gas_cv_from_cp(cp)

        self.assertAlmostEqual(cv, 20.68553738184676)
        self.assertAlmostEqual(calc_ideal_gas_cp_from_cv(cv), cp)
        self.assertAlmostEqual(calc_heat_capacity_ratio(cp, cv), cp / cv)

    def test_heat_capacity_core_numeric_contract(self):
        cv = _calc_ideal_gas_cv_from_cp([29.0, 31.0], 8.0)
        self.assertEqual(cv.tolist(), [21.0, 23.0])

        cp = _calc_ideal_gas_cp_from_cv([[20.0, 21.0], [22.0, 23.0]], [[8.0], [9.0]])
        self.assertEqual(cp.tolist(), [[28.0, 29.0], [31.0, 32.0]])

        ratio = _calc_heat_capacity_ratio([28.0, 30.0], [20.0, 25.0])
        self.assertEqual(ratio.tolist(), [1.4, 1.2])

        with self.assertRaises(ValueError):
            _calc_ideal_gas_cv_from_cp(8.0, 8.0)
        with self.assertRaises(ValueError):
            _calc_heat_capacity_ratio([28.0, 30.0], [20.0, 25.0, 26.0])

    def test_heat_capacity_props_adapter_contract(self):
        cp = CustomProp(value=0.029, unit="kJ/mol.K")
        r = CustomProp(value=8.31446261815324, unit="J/mol.K")

        self.assertAlmostEqual(
            _calc_ideal_gas_cv_from_cp_from_props(
                cp,
                r,
                output_heat_capacity_unit="J/mol.K",
            ),
            20.68553738184676,
        )

        with self.assertRaises(TypeError):
            _calc_heat_capacity_ratio_from_props(29.0, CustomProp(value=20.0, unit="J/mol.K"))

    def test_public_helmholtz_wrapper_delegates_to_core_adapter(self):
        temperature = Temperature(value=300.0, unit="K")

        self.assertEqual(calc_helmholtz_energy(8000.0, temperature, 10.0), 5000.0)

    def test_helmholtz_core_numeric_contract(self):
        self.assertEqual(_calc_helmholtz_energy(8000.0, 300.0, 10.0), 5000.0)

        vector = _calc_helmholtz_energy([8000.0, 9000.0], 300.0, [10.0, 12.0])
        self.assertEqual(vector.tolist(), [5000.0, 5400.0])

        matrix = _calc_helmholtz_energy(
            [[8000.0, 9000.0], [7000.0, 8500.0]],
            [[300.0], [250.0]],
            [[10.0, 12.0], [8.0, 10.0]],
        )
        self.assertEqual(matrix.tolist(), [[5000.0, 5400.0], [5000.0, 6000.0]])

        with self.assertRaises(ValueError):
            _calc_helmholtz_energy([1.0, 2.0], [1.0, 2.0, 3.0], [1.0, 2.0])
        with self.assertRaises(ValueError):
            _calc_helmholtz_energy(float("inf"), 300.0, 10.0)

    def test_helmholtz_props_adapter_contract(self):
        temperature = Temperature(value=300.0, unit="K")
        internal_energy = CustomProp(value=8.0, unit="kJ/mol")
        entropy = CustomProp(value=10.0, unit="J/mol.K")

        self.assertEqual(
            _calc_helmholtz_energy_from_props(
                internal_energy,
                temperature,
                entropy,
                output_internal_energy_unit="J/mol",
            ),
            5000.0,
        )

        with self.assertRaises(TypeError):
            _calc_helmholtz_energy_from_props(8.0, temperature, entropy)


if __name__ == "__main__":
    unittest.main()
