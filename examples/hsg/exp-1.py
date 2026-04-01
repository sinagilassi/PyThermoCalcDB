# import libs
# import packages/modules
import logging
import warnings
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pythermodb_settings.models import Pressure, Temperature, CustomProp, Volume, Component, CustomProperty
# locals
from pyThermoCalcDB.docs.thermo import calc_En, calc_En_IG_ref
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
    value=298.15,
    unit='K',
)

# enthalpy at temperature
res_1 = calc_En(
    component=methanol_gas,
    model_source=model_source,
    temperature=temperature,
)
print(f'Enthalpy of methanol gas at {temperature}:')
print(res_1)

# liquid enthalpy at temperature
res_2 = calc_En_IG_ref(
    component=methanol_liquid,
    model_source=model_source,
    temperature=temperature,
    component_key="Name-Formula"
)
print(f'Enthalpy of methanol liquid at {temperature}:')
print(res_2)
