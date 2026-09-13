import unittest

from pythermodb_settings.models import CustomProp, Temperature

from pythermocalcdb.thermo import calc_gibbs_energy, calc_gibbs_energy_change
from pythermocalcdb.thermo.core.gibbs import (
    _calc_gibbs_energy,
    _calc_gibbs_energy_change,
    _calc_gibbs_energy_from_props,
)


class TestGibbsCoreApi(unittest.TestCase):
    def test_public_gibbs_wrappers_delegate_to_core_adapters(self):
        temperature = Temperature(value=300.0, unit="K")

        self.assertEqual(calc_gibbs_energy(10000.0, temperature, 20.0), 4000.0)
        self.assertEqual(
            calc_gibbs_energy_change(5000.0, 10.0, temperature),
            2000.0,
        )

    def test_core_numeric_contract(self):
        self.assertEqual(_calc_gibbs_energy(10000.0, 300.0, 20.0), 4000.0)

        vector = _calc_gibbs_energy([10000.0, 12000.0], 300.0, [20.0, 25.0])
        self.assertEqual(vector.tolist(), [4000.0, 4500.0])

        matrix = _calc_gibbs_energy(
            [[10000.0, 12000.0], [8000.0, 9000.0]],
            [[300.0], [250.0]],
            [[20.0, 25.0], [10.0, 12.0]],
        )
        self.assertEqual(
            matrix.tolist(),
            [[4000.0, 4500.0], [5500.0, 6000.0]],
        )

        change = _calc_gibbs_energy_change(
            [5000.0, 6000.0],
            [10.0, 12.0],
            300.0,
        )
        self.assertEqual(change.tolist(), [2000.0, 2400.0])

    def test_core_numeric_validation(self):
        with self.assertRaises(ValueError):
            _calc_gibbs_energy([1.0, 2.0], [1.0, 2.0, 3.0], [1.0, 2.0])
        with self.assertRaises(ValueError):
            _calc_gibbs_energy(float("nan"), 300.0, 20.0)

    def test_props_adapter_contract(self):
        temperature = Temperature(value=300.0, unit="K")
        enthalpy = CustomProp(value=10.0, unit="kJ/mol")
        entropy = CustomProp(value=20.0, unit="J/mol.K")

        self.assertEqual(
            _calc_gibbs_energy_from_props(
                enthalpy,
                temperature,
                entropy,
                output_enthalpy_unit="J/mol",
            ),
            4000.0,
        )

        with self.assertRaises(TypeError):
            _calc_gibbs_energy_from_props(10.0, temperature, entropy)


if __name__ == "__main__":
    unittest.main()
