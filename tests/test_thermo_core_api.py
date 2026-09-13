import unittest

from pythermodb_settings.models import CustomProp, Temperature

from pythermocalcdb.thermo import (
    calc_heat_capacity_ratio,
    calc_helmholtz_energy,
    calc_ideal_gas_internal_energy,
    calc_ideal_gas_cp_from_cv,
    calc_ideal_gas_cv_from_cp,
    calc_internal_energy,
    density_to_specific_volume,
    specific_volume_to_density,
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
from pythermocalcdb.thermo.core.internal_energy import (
    _calc_ideal_gas_internal_energy,
    _calc_ideal_gas_internal_energy_from_props,
    _calc_internal_energy,
    _calc_internal_energy_from_props,
)
from pythermocalcdb.thermo.core.specific_volume import (
    _calc_density_to_specific_volume,
    _calc_density_to_specific_volume_from_props,
    _calc_specific_volume_to_density,
    _calc_specific_volume_to_density_from_props,
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

    def test_public_internal_energy_wrappers_delegate_to_core_adapters(self):
        temperature = Temperature(value=300.0, unit="K")

        self.assertEqual(calc_internal_energy(1000.0, 100000.0, 0.002), 800.0)
        self.assertAlmostEqual(
            calc_ideal_gas_internal_energy(10000.0, temperature),
            7505.661214554028,
        )

    def test_internal_energy_core_numeric_contract(self):
        self.assertEqual(_calc_internal_energy(1000.0, 100000.0, 0.002), 800.0)

        vector = _calc_internal_energy(
            [1000.0, 1200.0],
            [100000.0, 200000.0],
            [0.002, 0.001],
        )
        self.assertEqual(vector.tolist(), [800.0, 1000.0])

        matrix = _calc_internal_energy(
            [[1000.0, 1200.0], [900.0, 1100.0]],
            [[100000.0], [200000.0]],
            [[0.002, 0.001], [0.001, 0.002]],
        )
        self.assertEqual(matrix.tolist(), [[800.0, 1100.0], [700.0, 700.0]])

        ideal = _calc_ideal_gas_internal_energy([10000.0, 12000.0], 300.0, 8.0)
        self.assertEqual(ideal.tolist(), [7600.0, 9600.0])

        with self.assertRaises(ValueError):
            _calc_internal_energy(1000.0, -100000.0, 0.002)
        with self.assertRaises(ValueError):
            _calc_ideal_gas_internal_energy(
                [1.0, 2.0],
                [300.0, 301.0, 302.0],
                8.0,
            )

    def test_internal_energy_props_adapter_contract(self):
        temperature = Temperature(value=300.0, unit="K")
        enthalpy = CustomProp(value=1000.0, unit="J")
        pressure = CustomProp(value=100000.0, unit="Pa")
        volume = CustomProp(value=0.002, unit="m^3")
        molar_enthalpy = CustomProp(value=10.0, unit="kJ/mol")

        self.assertEqual(
            _calc_internal_energy_from_props(enthalpy, pressure, volume),
            800.0,
        )
        self.assertAlmostEqual(
            _calc_ideal_gas_internal_energy_from_props(
                molar_enthalpy,
                temperature,
                output_molar_enthalpy_unit="J/mol",
            ),
            7505.661214554028,
        )

        with self.assertRaises(TypeError):
            _calc_internal_energy_from_props(1000.0, pressure, volume)

    def test_public_specific_volume_wrappers_delegate_to_core_adapters(self):
        self.assertEqual(density_to_specific_volume(1000.0), 0.001)
        self.assertEqual(specific_volume_to_density(0.001), 1000.0)

    def test_specific_volume_core_numeric_contract(self):
        specific_volume = _calc_density_to_specific_volume([1000.0, 800.0])
        self.assertEqual(specific_volume.tolist(), [0.001, 0.00125])

        density = _calc_specific_volume_to_density(
            [[0.001, 0.002], [0.004, 0.005]]
        )
        self.assertEqual(density.tolist(), [[1000.0, 500.0], [250.0, 200.0]])

        with self.assertRaises(ValueError):
            _calc_density_to_specific_volume(0.0)
        with self.assertRaises(ValueError):
            _calc_specific_volume_to_density(float("nan"))

    def test_specific_volume_props_adapter_contract(self):
        density = CustomProp(value=1000.0, unit="kg/m^3")
        specific_volume = CustomProp(value=0.001, unit="m^3/kg")

        self.assertEqual(
            _calc_density_to_specific_volume_from_props(density),
            0.001,
        )
        self.assertEqual(
            _calc_specific_volume_to_density_from_props(specific_volume),
            1000.0,
        )

        with self.assertRaises(TypeError):
            _calc_density_to_specific_volume_from_props(1000.0)


if __name__ == "__main__":
    unittest.main()
