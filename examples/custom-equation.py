# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
# local
from pyThermoCalcDB.docs import (
    mkeq
)

# check version
print(ptdb.__version__)
print(ptdblink.__version__)

# =======================================
# ! LOAD THERMODB
# =======================================
# NOTE: parent directory
parent_dir = os.path.dirname(os.path.abspath(__file__))
print(parent_dir)

# NOTE: thermodb directory
thermodb_dir = os.path.join(parent_dir, 'thermodb')
print(thermodb_dir)

# ! CO2
# thermodb file name
CO2_thermodb_file = os.path.join(thermodb_dir, 'carbon dioxide-g.pkl')

# ! methane
# thermodb file name
methane_thermodb_file = os.path.join(thermodb_dir, 'methane-g.pkl')

# =======================================
# SECTION: COMPONENTS THERMODB SOURCE
# =======================================
# NOTE: carbon dioxide
CO2_comp = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g'
)

CO2_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CO2_comp,
    source=CO2_thermodb_file
)

# NOTE: methane
CH4_Comp = Component(
    name='methane',
    formula='CH4',
    state='g'
)

methane_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CH4_Comp,
    source=methane_thermodb_file
)

# NOTE: component thermodb source
_component_thermodb: list = [
    CO2_component_thermodb,
    methane_component_thermodb,
]
# =======================================
# SECTION: BUILD THERMODB MODEL SOURCE
# =======================================
# update thermodb rule
thermodb_rules: Dict[str, Dict[str, ComponentRule]] = {
    'ALL': {
        'DATA': {
            'critical-pressure': 'Pc',
            'critical-temperature': 'Tc',
            'acentric-factor': 'AcFa',
            'Enthalpy-of-Formation': 'EnFo_IG',
            'Gibbs-Energy-of-Formation': 'GiEnFo_IG',
        },
        'EQUATIONS': {
            'CUSTOM-REF-1::vapor_pressure': 'VaPr',
            'CUSTOM-REF-1::ideal_gas_heat_capacity_at_constant_pressure': 'Cp_IG',
            'CUSTOM-REF-1::liquid_heat_capacity_at_constant_pressure': 'Cp_LIQ',
            'CUSTOM-REF-1::liquid_density': 'rho_LIQ',
        }
    },
}

model_source_: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=_component_thermodb,
    rules=thermodb_rules,
)
print(model_source_)

# get data source and equation source
datasource = model_source_.data_source
equationsource = model_source_.equation_source

# ------------------------------------------------
# ! THERMODYNAMIC PROPERTIES
# ------------------------------------------------
# vapor pressure
VaPr = equationsource['carbon dioxide-g']['VaPr'].cal(T=300.1)
print(VaPr)
# heat capacity
Cp_IG = equationsource['carbon dioxide-g']['Cp_IG'].cal(T=300.1)
print(Cp_IG)

# =======================================
# SECTION: MAKE EQUATION OBJECT
# =======================================
# NOTE: make equation object for vapor pressure
VaPr_eq = mkeq(
    prop_name='VaPr',
    component=CO2_comp,
    model_source=model_source_,
    component_key='Name-State'
)
print(VaPr_eq)
