# import libs
# import packages/modules
import logging
import warnings
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pythermodb_settings.models import Pressure, Temperature, CustomProp, Volume, Component, CustomProperty
# locals
from pythermocalcdb.models import ComponentEnthalpy
from pythermocalcdb.docs.thermo import calc_En, calc_En_IG_ref
# NOTE: for example
# ! model sources
from model_source import model_source


# ====================================================
# SECTION: Components
# ====================================================
# methanol gas
H2_gas = Component(
    name='hydrogen',
    formula='H2',
    state='g',
)

H2_liquid = Component(
    name='hydrogen',
    formula='H2',
    state='l',
)

# ====================================================
# SECTION: Thermodynamic properties
# ====================================================
# temperature
temperature = Temperature(
    value=25.0,
    unit='K',
)
temperature_value = temperature.value
print(f'Temperature: {temperature_value} K')

# enthalpy of vaporization
# Hydrogen EnVap correlation uses t = 1 - T/Tc, so T must be below Tc.
hydrogen_critical_temperature = 33.19  # K
if temperature_value >= hydrogen_critical_temperature:
    raise ValueError(
        f"Hydrogen EnVap correlation is only valid for T < Tc "
        f"(Tc={hydrogen_critical_temperature} K, received {temperature_value} K)."
    )

res_0 = model_source.equation_source['hydrogen-H2']['EnVap'].cal(
    T=temperature_value
)
print(f'Enthalpy of vaporization at {temperature_value}:')
print(res_0)


# enthalpy at temperature
res_1: ComponentEnthalpy | None = calc_En(
    component=H2_gas,
    model_source=model_source,
    temperature=temperature,
)
print(f'Enthalpy of gas at {temperature_value}:')
print(res_1)

# liquid enthalpy at temperature
res_2: ComponentEnthalpy | None = calc_En_IG_ref(
    component=H2_liquid,
    model_source=model_source,
    temperature=temperature,
    component_key="Name-Formula"
)
print(f'Enthalpy of liquid at {temperature_value}:')
print(res_2)

# check difference (should be close to enthalpy of vaporization)
if res_1 is not None and res_2 is not None:
    delta_En = res_1.value - res_2.value
    print(f'Enthalpy of vaporization (gas - liquid): {delta_En} J/mol')
