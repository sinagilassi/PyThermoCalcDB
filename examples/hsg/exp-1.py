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
methanol_gas = Component(
    name='methanol',
    formula='CH3OH',
    state='g',
)

# methanol liquid
methanol_liquid = Component(
    name='methanol',
    formula='CH3OH',
    state='l',
)

# ====================================================
# SECTION: Thermodynamic properties
# ====================================================
# temperature
temperature = Temperature(
    value=398.15,
    unit='K',
)
temperature_value = temperature.value
print(f'Temperature: {temperature_value} K')

# enthalpy of vaporization
res_0 = model_source.equation_source['methanol-CH3OH']['EnVap'].cal(
    T=temperature_value)
print(f'Enthalpy of vaporization of methanol at {temperature_value}:')
print(res_0)


# enthalpy at temperature
res_1: ComponentEnthalpy | None = calc_En(
    component=methanol_gas,
    model_source=model_source,
    temperature=temperature,
)
print(f'Enthalpy of methanol gas at {temperature_value}:')
print(res_1)

# liquid enthalpy at temperature
res_2: ComponentEnthalpy | None = calc_En_IG_ref(
    component=methanol_liquid,
    model_source=model_source,
    temperature=temperature,
    component_key="Name-Formula"
)
print(f'Enthalpy of methanol liquid at {temperature_value}:')
print(res_2)

# check difference (should be close to enthalpy of vaporization)
if res_1 is not None and res_2 is not None:
    delta_En = res_1.value - res_2.value
    print(f'Enthalpy of vaporization (gas - liquid): {delta_En} J/mol')
