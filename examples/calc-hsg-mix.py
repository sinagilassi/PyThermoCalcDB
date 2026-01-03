# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
# local
from pyThermoCalcDB.docs.thermo import (
    calc_enthalpy_of_formation_at_temperature,
    calc_gibbs_energy_of_formation_at_temperature,
    calc_enthalpy_of_formation_range,
    calc_gibbs_energy_of_formation_range,
    calc_enthalpy_change,
    calc_entropy_change,
    calc_mixture_enthalpy
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

# ! CO
# FIXME: correct file name
CO_thermodb_file = os.path.join(thermodb_dir, 'carbon dioxide-g.pkl')

# =======================================
# SECTION: COMPONENTS THERMODB SOURCE
# =======================================
# NOTE: carbon dioxide
# ! component
CO2_comp = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g'
)
# ! component thermodb source
CO2_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CO2_comp,
    source=CO2_thermodb_file
)

# ! CO
CO_comp = Component(
    name='carbon monoxide',
    formula='CO',
    state='g'
)
# ! component thermodb source
CO_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CO_comp,
    source=CO_thermodb_file
)

# NOTE: component thermodb source
_component_thermodb: list = [
    CO2_component_thermodb,
    CO_component_thermodb
]

# =======================================
# SECTION: BUILD THERMODB MODEL SOURCE
# =======================================
# update thermodb rule
thermodb_rules: Dict[str, Dict[str, ComponentRule]] = {
    'ALL_HIDDEN': {
        'DATA': {
            'critical-pressure': 'Pc',
            'critical-temperature': 'Tc',
            'acentric-factor': 'AcFa',
            'Enthalpy-of-Formation': 'EnFo_IG',
            'Gibbs-Energy-of-Formation': 'GiEnFo_IG',
        },
        'EQUATIONS': {
            'CUSTOM-REF-1::vapor-pressure': 'VaPr',
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG'
        }
    },
}

model_source_: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=_component_thermodb,
    rules=thermodb_rules,
    original_equation_label=False
)
print(model_source_)

# get data source and equation source
datasource = model_source_.data_source
equationsource = model_source_.equation_source

# =======================================
# SECTION: CALCULATE MIXTURE ENTHALPY
# =======================================
# ! mixture components
mixture_components = [
    CO2_comp,
    CO_comp
]

# ! mole fractions
CO2_comp.mole_fraction = 0.7
CO_comp.mole_fraction = 0.3
# log
print(mixture_components)

# ! temperature and pressure
mixture_temperature = Temperature(value=320.0, unit='K')
mixture_pressure = Pressure(value=5.0, unit='atm')

# ! calculate mixture enthalpy
mixture_enthalpy_res = calc_mixture_enthalpy(
    components=mixture_components,
    model_source=model_source_,
    temperature=mixture_temperature,
    pressure=mixture_pressure,
    phase='IG',
    departure_enthalpy=None,
    excess_enthalpy=None,
    output_unit='kJ/mol',
    mode='log'
)
# log
print(mixture_enthalpy_res)
