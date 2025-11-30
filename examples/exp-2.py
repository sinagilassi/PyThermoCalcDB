# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
# local
from pyThermoCalcDB.docs.thermo import calc_enthalpy_of_formation_at_temperature, calc_gibbs_energy_of_formation_at_temperature

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

# ! acetylene
# thermodb file name
acetylene_thermodb_file = os.path.join(thermodb_dir, 'acetylene-g.pkl')

# ! n-butane
# thermodb file name
n_butane_thermodb_file = os.path.join(thermodb_dir, 'n-butane-g.pkl')

# ! ethanol
# thermodb file name
ethanol_thermodb_file = os.path.join(thermodb_dir, 'ethanol-l.pkl')

# ! methanol
# thermodb file name
methanol_thermodb_file = os.path.join(thermodb_dir, 'methanol-g.pkl')

# ! 1-butene
# thermodb file name
butene_thermodb_file = os.path.join(thermodb_dir, '1-butene-g.pkl')

# ! propane
# thermodb file name
propane_thermodb_file = os.path.join(thermodb_dir, 'propane-g.pkl')

# ! methane
# thermodb file name
methane_thermodb_file = os.path.join(thermodb_dir, 'methane-g.pkl')

# =======================================
# SECTION: COMPONENTS THERMODB SOURCE
# =======================================
# NOTE: carbon dioxide
CO2_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='carbon dioxide',
        formula='CO2',
        state='g'
    ),
    source=CO2_thermodb_file
)

# NOTE: acetylene
acetylene_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='acetylene',
        formula='C2H2',
        state='g'
    ),
    source=acetylene_thermodb_file
)

# NOTE: n-butane
n_butane_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='n-butane',
        formula='C4H10',
        state='g'
    ),
    source=n_butane_thermodb_file
)

# NOTE: ethanol
ethanol_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='ethanol',
        formula='C2H5OH',
        state='l'
    ),
    source=ethanol_thermodb_file
)

# NOTE: methanol
methanol_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='methanol',
        formula='CH3OH',
        state='g'
    ),
    source=methanol_thermodb_file
)

# NOTE: 1-butene
butene_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=Component(
        name='1-butene',
        formula='C4H8',
        state='g'
    ),
    source=butene_thermodb_file
)

# NOTE: propane
C3H8_Comp = Component(
    name='propane',
    formula='C3H8',
    state='g'
)
propane_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=C3H8_Comp,
    source=propane_thermodb_file
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
    acetylene_component_thermodb,
    n_butane_component_thermodb,
    ethanol_component_thermodb,
    methanol_component_thermodb,
    butene_component_thermodb,
    propane_component_thermodb,
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
            'CUSTOM-REF-1::vapor-pressure': 'VaPr',
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG'
        }
    },
    'CH4-g': {
        'DATA': {
            'critical-pressure': 'Pc',
            'critical-temperature': 'Tc',
            'acentric-factor': 'AcFa'
        },
        'EQUATIONS': {
            'CUSTOM-REF-1::vapor-pressure': 'VaPr',
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG'
        }
    },
    'Methane-g': {
        'DATA': {
            'critical-pressure': 'Pc',
            'critical-temperature': 'Tc',
            'acentric-factor': 'AcFa'
        },
        'EQUATIONS': {
            'CUSTOM-REF-1::vapor-pressure': 'VaPr',
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG'
        }
    }
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
# SECTION: CALCULATE ENTHALPY OF FORMATION
# =======================================
# ! component
CO2 = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g'
)

# ! temperature
temperature_ = Temperature(value=300.0, unit='K')

# ! calculate enthalpy of formation at temperature
EnFo_IG_T_res = calc_enthalpy_of_formation_at_temperature(
    component=CO2,
    model_source=model_source_,
    temperature=temperature_
)
print(EnFo_IG_T_res)


# ! calculate gibbs free energy of formation at temperature
GiEnFo_IG_T_res = calc_gibbs_energy_of_formation_at_temperature(
    component=CO2,
    model_source=model_source_,
    temperature=temperature_
)
print(GiEnFo_IG_T_res)
