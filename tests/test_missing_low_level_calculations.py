import math
import unittest

import numpy as np
from pythermodb_settings.models import CustomProp, Pressure, Temperature

from pythermocalcdb.flows import (
    calc_additive_liquid_volumetric_flow_rate,
    calc_concentration_from_molar_flow_rate,
    calc_enthalpy_flow_rate,
    calc_flowing_heat_capacity,
    calc_gas_molar_flow_rate_from_z,
    calc_gas_volumetric_flow_rate_from_z,
    calc_ideal_gas_molar_flow_rate,
    calc_ideal_gas_volumetric_flow_rate,
    calc_molar_flow_rate_from_concentration,
)
from pythermocalcdb.flows.core.energy import _calc_enthalpy_flow_rate
from pythermocalcdb.mixtures import (
    calc_additive_liquid_volume,
    calc_total_heat_capacity,
)
from pythermocalcdb.mixtures.core.heat_capacity import _calc_total_heat_capacity
from pythermocalcdb.reactions import (
    calc_reaction_enthalpy_from_constant_delta_cp,
    calc_reaction_heat_capacity_change,
    calc_reaction_heat_rate,
    calc_reaction_volumetric_heat_source,
)
from pythermocalcdb.reactions.core.energetics import (
    _calc_reaction_heat_capacity_change,
    _calc_reaction_volumetric_heat_source,
)
from pythermocalcdb.thermo import (
    calc_gas_pressure_from_z,
    calc_gas_volume_from_z,
    calc_ideal_gas_pressure,
    calc_ideal_gas_volume,
)
from pythermocalcdb.thermo.core.density import _calc_ideal_gas_volume


class MissingLowLevelCalculationTests(unittest.TestCase):
    def test_total_and_flowing_heat_capacity_examples(self):
        self.assertTrue(math.isclose(calc_total_heat_capacity([1.0, 2.0], [10.0, 20.0]), 50.0))
        self.assertTrue(math.isclose(calc_flowing_heat_capacity([1.0, 2.0], [10.0, 20.0]), 50.0))

        batched = _calc_total_heat_capacity([[1.0, 2.0], [2.0, 3.0]], [[10.0, 20.0], [4.0, 5.0]])
        np.testing.assert_allclose(batched, np.array([50.0, 23.0]))

    def test_reaction_heat_capacity_and_kirchhoff_examples(self):
        self.assertTrue(math.isclose(calc_reaction_heat_capacity_change([-1.0, -2.0, 1.0], [10.0, 20.0, 70.0]), 20.0))

        enthalpy = calc_reaction_enthalpy_from_constant_delta_cp(
            -100000.0,
            20.0,
            Temperature(value=398.15, unit="K"),
            Temperature(value=298.15, unit="K"),
        )
        self.assertTrue(math.isclose(enthalpy, -98000.0))

        batched = _calc_reaction_heat_capacity_change([[-1.0, 1.0], [-2.0, 1.0]], [[10.0, 30.0], [20.0, 70.0]])
        np.testing.assert_allclose(batched, np.array([20.0, 30.0]))

    def test_ideal_and_z_gas_state_round_trips(self):
        moles = CustomProp(value=1.0, unit="mol")
        temperature = Temperature(value=298.15, unit="K")
        pressure = Pressure(value=101325.0, unit="Pa")

        volume = calc_ideal_gas_volume(moles, temperature, pressure)
        self.assertTrue(math.isclose(volume.value, 0.024465403697038125, rel_tol=1e-12))

        pressure_roundtrip = calc_ideal_gas_pressure(moles, temperature, volume)
        self.assertTrue(math.isclose(pressure_roundtrip.value, pressure.value, rel_tol=1e-12))

        z_volume = calc_gas_volume_from_z(moles, temperature, pressure, 0.9)
        self.assertTrue(math.isclose(z_volume.value, 0.9 * volume.value, rel_tol=1e-12))

        z_pressure = calc_gas_pressure_from_z(moles, temperature, z_volume, 0.9)
        self.assertTrue(math.isclose(z_pressure.value, pressure.value, rel_tol=1e-12))

        np.testing.assert_allclose(
            _calc_ideal_gas_volume([1.0, 2.0], 298.15, 101325.0),
            np.array([volume.value, 2.0 * volume.value]),
        )

    def test_gas_flow_round_trips(self):
        temperature = Temperature(value=298.15, unit="K")
        pressure = Pressure(value=101325.0, unit="Pa")

        q = calc_ideal_gas_volumetric_flow_rate(1.0, temperature, pressure)
        f = calc_ideal_gas_molar_flow_rate(q, temperature, pressure)
        self.assertTrue(math.isclose(f.value, 1.0, rel_tol=1e-12))

        q_z = calc_gas_volumetric_flow_rate_from_z(1.0, temperature, pressure, 0.9)
        f_z = calc_gas_molar_flow_rate_from_z(q_z, temperature, pressure, 0.9)
        self.assertTrue(math.isclose(f_z.value, 1.0, rel_tol=1e-12))

    def test_concentration_enthalpy_and_liquid_helpers(self):
        flow = calc_molar_flow_rate_from_concentration(2.0, 3.0)
        self.assertTrue(math.isclose(flow.value, 6.0))

        concentration = calc_concentration_from_molar_flow_rate(flow, 3.0)
        self.assertTrue(math.isclose(concentration.value, 2.0))

        self.assertTrue(math.isclose(calc_enthalpy_flow_rate([1.0, 2.0], [100.0, 200.0]), 500.0))
        self.assertTrue(math.isclose(calc_additive_liquid_volume([1.0, 2.0], [0.018, 0.046], [1000.0, 800.0]), 0.000133))
        self.assertTrue(math.isclose(calc_additive_liquid_volumetric_flow_rate([1.0, 2.0], [0.018, 0.046], [1000.0, 800.0]), 0.000133))

        np.testing.assert_allclose(
            _calc_enthalpy_flow_rate([[1.0, 2.0], [3.0, 4.0]], [[100.0, 200.0], [100.0, 200.0]]),
            np.array([500.0, 1100.0]),
        )

    def test_reaction_heat_generation_sign_convention(self):
        self.assertTrue(math.isclose(calc_reaction_volumetric_heat_source([-100000.0], [0.01]), 1000.0))
        self.assertTrue(math.isclose(calc_reaction_heat_rate([-100000.0], [0.01], 2.0), 2000.0))
        self.assertTrue(math.isclose(_calc_reaction_volumetric_heat_source([-100000.0, 50000.0], [0.01, 0.02]), -0.0))

    def test_invalid_shapes_and_domains_raise(self):
        with self.assertRaises(ValueError):
            _calc_total_heat_capacity([1.0, -1.0], [10.0, 20.0])

        with self.assertRaises(ValueError):
            _calc_reaction_heat_capacity_change([-1.0, 1.0], [10.0])


if __name__ == "__main__":
    unittest.main()
