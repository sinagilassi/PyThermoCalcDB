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

# NOTE: component thermodb source
_component_thermodb: list = [
    CO2_component_thermodb,
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
# SECTION: CALCULATE ENTHALPY OF FORMATION
# =======================================
# ! temperature
temperature_ref = Temperature(value=298.15, unit='K')
temperature_ = Temperature(value=300.0, unit='K')

# ! pressure
pressure_ref = Pressure(value=1.0, unit='atm')
pressure_ = Pressure(value=1.0, unit='atm')

# ! calculate enthalpy change from reference temperature to temperature
dEn_IG_res = calc_enthalpy_change(
    component=CO2_comp,
    model_source=model_source_,
    temperature_initial=temperature_ref,
    temperature_final=temperature_,
    mode='log'
)
print(dEn_IG_res)

# ! calculate entropy change from reference temperature to temperature
dS_IG_res = calc_entropy_change(
    component=CO2_comp,
    model_source=model_source_,
    temperature_initial=temperature_ref,
    temperature_final=temperature_,
    pressure_initial=pressure_ref,
    pressure_final=pressure_,
    phase='IG',
    mode='log'
)
print(dS_IG_res)

# ! calculate enthalpy of formation at temperature
EnFo_IG_T_res = calc_enthalpy_of_formation_at_temperature(
    component=CO2_comp,
    model_source=model_source_,
    temperature=temperature_,
    mode='attach'
)
print(EnFo_IG_T_res)

EnFo_IG_T_res = calc_enthalpy_of_formation_at_temperature(
    component=CO2_comp,
    model_source=model_source_,
    temperature=temperature_,
    mode='silent'
)
print(EnFo_IG_T_res)

# ! calculate gibbs free energy of formation at temperature
GiEnFo_IG_T_res = calc_gibbs_energy_of_formation_at_temperature(
    component=CO2_comp,
    model_source=model_source_,
    temperature=temperature_,
    mode='log'
)
print(GiEnFo_IG_T_res)

# ! calculate enthalpy of formation range
temperature_range_ = [
    Temperature(value=295.0, unit='K'),
    Temperature(value=300.0, unit='K'),
    Temperature(value=350.0, unit='K'),
    Temperature(value=370.0, unit='K'),
]

EnFo_IG_range_res = calc_enthalpy_of_formation_range(
    component=CO2_comp,
    model_source=model_source_,
    temperatures=temperature_range_,
    mode='attach'
)
print(EnFo_IG_range_res)

# ! calculate gibbs free energy of formation range
GiEnFo_IG_range_res = calc_gibbs_energy_of_formation_range(
    component=CO2_comp,
    model_source=model_source_,
    temperatures=temperature_range_,
    mode='log'
)
print(GiEnFo_IG_range_res)
