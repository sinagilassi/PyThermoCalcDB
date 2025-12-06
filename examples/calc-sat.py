# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
from pyThermoDB.core import TableEquation
# local
from pyThermoCalcDB.docs.sat import (
    calc_vapor_pressure_at_temperature,
    calc_enthalpy_of_vaporization_at_temperature,
    calc_saturated_temperature_at_pressure,
    calc_vapor_pressure_sensitivity_at_temperature
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
# ! component
CO2_comp = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g'
)

CO2_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CO2_comp,
    source=CO2_thermodb_file
)

# NOTE: ethanol
# ! component
ethanol_comp = Component(
    name='ethanol',
    formula='C2H5OH',
    state='l'
)

ethanol_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=ethanol_comp,
    source=ethanol_thermodb_file
)

# NOTE: component thermodb source
_component_thermodb: list = [
    CO2_component_thermodb,
    ethanol_component_thermodb
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
VaPr = equationsource['carbon dioxide-g']['VaPr']
if not isinstance(VaPr, TableEquation):
    raise ValueError("VaPr is not an EquationModel")

VaPr_res = VaPr.cal(T=300.1)
print(VaPr_res)
# >> variable range
VaPr_range = VaPr.get_variable_range_values()
print(VaPr_range)

# heat capacity
Cp_IG = equationsource['carbon dioxide-g']['Cp_IG'].cal(T=300.1)
print(Cp_IG)

# =======================================
# SECTION: CALCULATE ENTHALPY OF FORMATION
# =======================================
# ! range of temperature
temperature_range_ = [
    Temperature(value=295.0, unit='K'),
    Temperature(value=300.0, unit='K'),
    Temperature(value=350.0, unit='K'),
    Temperature(value=370.0, unit='K'),
]

# ! calculate vapor pressure at temperatures
for T in temperature_range_:
    Pvap = calc_vapor_pressure_at_temperature(
        component=ethanol_comp,
        temperature=T,
        model_source=model_source_,
        mode='attach'
    )
    if Pvap is not None:
        print(Pvap)

# ! calculate enthalpy of vaporization at temperatures
for T in temperature_range_:
    Hvap = calc_enthalpy_of_vaporization_at_temperature(
        component=ethanol_comp,
        temperature=T,
        model_source=model_source_,
        mode='attach'
    )
    if Hvap is not None:
        print(Hvap)

# ! calculate saturated temperature at pressures
pressure_range_ = [
    Pressure(value=5000.0, unit='kPa'),
    Pressure(value=10000.0, unit='kPa'),
    Pressure(value=15000.0, unit='kPa'),
    Pressure(value=20000.0, unit='kPa'),
]
# {'T': {'min': {'value': 216.58, 'unit': 'K', 'symbol': 'Tmin'}, 'max': {'value': 304.21, 'unit': 'K', 'symbol': 'Tmax'}}}

for P in pressure_range_:
    Tsat = calc_saturated_temperature_at_pressure(
        component=ethanol_comp,
        pressure=P,
        model_source=model_source_,
        temperature_guess=Temperature(value=250.0, unit='K'),
        T_bracket=(
            Temperature(value=216.58, unit='K'),
            Temperature(value=304.21, unit='K')
        ),
        method='least_squares',
        mode='attach'
    )
    if Tsat is not None:
        print(Tsat)

# ! calculate vapor pressure sensitivity at temperatures
for T in temperature_range_:
    dPvap_dT = calc_vapor_pressure_sensitivity_at_temperature(
        component=ethanol_comp,
        temperature=T,
        model_source=model_source_,
        mode='attach'
    )
    if dPvap_dT is not None:
        print(dPvap_dT)
