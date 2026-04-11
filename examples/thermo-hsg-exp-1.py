# import packages/modules
import os
from typing import Dict
from rich import print
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import (
    Component,
    ComponentRule,
    ComponentThermoDBSource,
    Temperature,
    Pressure,
)
# local
from pyThermoCalcDB.docs.thermo import (
    build_hsg_properties,
    calc_dEn,
    calc_dEn_hsg,
    calc_En,
    calc_En_hsg,
    calc_En_IG_ref,
    calc_En_IG_ref_hsg,
    calc_GiFrEn,
    calc_GiFrEn_hsg,
    calc_En_range,
    calc_En_range_hsg,
    calc_GiFrEn_range,
    calc_GiFrEn_range_hsg,
    calc_dEnt,
    calc_dEnt_hsg,
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
CO2_thermodb_file = os.path.join(thermodb_dir, 'carbon dioxide-g.pkl')

# ! CO
CO_thermodb_file = os.path.join(thermodb_dir, 'carbon monoxide-g.pkl')

# =======================================
# SECTION: COMPONENTS THERMODB SOURCE
# =======================================
CO2_comp = Component(
    name='carbon dioxide',
    formula='CO2',
    state='g',
)

CO2_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CO2_comp,
    source=CO2_thermodb_file,
)

CO_comp = Component(
    name='carbon monoxide',
    formula='CO',
    state='g',
)

CO_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CO_comp,
    source=CO_thermodb_file,
)

_component_thermodb: list = [
    CO2_component_thermodb,
    CO_component_thermodb,
]

# =======================================
# SECTION: BUILD THERMODB MODEL SOURCE
# =======================================
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
            'CUSTOM-REF-1::ideal-gas-heat-capacity': 'Cp_IG',
        },
    },
}

model_source_: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=_component_thermodb,
    rules=thermodb_rules,
    original_equation_label=False,
)
print(model_source_)

# =======================================
# SECTION: BUILD HSG PROPERTIES
# =======================================
hsg_props = build_hsg_properties(
    component=CO2_comp,
    model_source=model_source_,
    component_key='Name-Formula',
)
print(hsg_props)

if hsg_props is None:
    raise ValueError('Failed to build HSGProperties instance.')

# =======================================
# SECTION: CALCULATE THERMODYNAMIC PROPERTIES USING *_hsg
# =======================================
temperature_ref = Temperature(value=298.15, unit='K')
temperature_ = Temperature(value=300.0, unit='K')
pressure_ref = Pressure(value=1.0, unit='atm')
pressure_ = Pressure(value=1.0, unit='atm')

temperature_range_ = [
    Temperature(value=295.0, unit='K'),
    Temperature(value=300.0, unit='K'),
    Temperature(value=350.0, unit='K'),
    Temperature(value=370.0, unit='K'),
]

# enthalpy change
dEn_res = calc_dEn_hsg(
    hsg_props=hsg_props,
    temperature_initial=temperature_ref,
    temperature_final=temperature_,
    mode='log',
)
print(dEn_res)

# enthalpy
En_res = calc_En_hsg(
    hsg_props=hsg_props,
    temperature=temperature_,
    mode='log',
)
print(En_res)

# reference enthalpy (phase inferred from hsg_props.component.state)
En_ref_res = calc_En_IG_ref_hsg(
    hsg_props=hsg_props,
    temperature=temperature_,
    mode='log',
)
print(En_ref_res)

# gibbs free energy
GiFrEn_res = calc_GiFrEn_hsg(
    hsg_props=hsg_props,
    temperature=temperature_,
    phase='IG',
    mode='log',
)
print(GiFrEn_res)

# enthalpy range
En_range_res = calc_En_range_hsg(
    hsg_props=hsg_props,
    temperatures=temperature_range_,
    mode='log',
)
print(En_range_res)

# gibbs free energy range
GiFrEn_range_res = calc_GiFrEn_range_hsg(
    hsg_props=hsg_props,
    temperatures=temperature_range_,
    mode='log',
)
print(GiFrEn_range_res)

# entropy change
dEnt_res = calc_dEnt_hsg(
    hsg_props=hsg_props,
    temperature_initial=temperature_ref,
    temperature_final=temperature_,
    pressure_initial=pressure_ref,
    pressure_final=pressure_,
    phase='IG',
    mode='log',
)
print(dEnt_res)

# =======================================
# SECTION: CALCULATE THERMODYNAMIC PROPERTIES USING NORMAL FUNCTIONS
# =======================================
dEn_res_normal = calc_dEn(
    component=CO2_comp,
    model_source=model_source_,
    temperature_initial=temperature_ref,
    temperature_final=temperature_,
    mode='log',
)
print(dEn_res_normal)

En_res_normal = calc_En(
    component=CO2_comp,
    model_source=model_source_,
    temperature=temperature_,
    mode='log',
)
print(En_res_normal)

En_ref_res_normal = calc_En_IG_ref(
    component=CO2_comp,
    model_source=model_source_,
    temperature=temperature_,
    mode='log',
)
print(En_ref_res_normal)

GiFrEn_res_normal = calc_GiFrEn(
    component=CO2_comp,
    model_source=model_source_,
    temperature=temperature_,
    phase='IG',
    mode='log',
)
print(GiFrEn_res_normal)

En_range_res_normal = calc_En_range(
    component=CO2_comp,
    model_source=model_source_,
    temperatures=temperature_range_,
    mode='log',
)
print(En_range_res_normal)

GiFrEn_range_res_normal = calc_GiFrEn_range(
    component=CO2_comp,
    model_source=model_source_,
    temperatures=temperature_range_,
    mode='log',
)
print(GiFrEn_range_res_normal)

dEnt_res_normal = calc_dEnt(
    component=CO2_comp,
    model_source=model_source_,
    temperature_initial=temperature_ref,
    temperature_final=temperature_,
    pressure_initial=pressure_ref,
    pressure_final=pressure_,
    phase='IG',
    mode='log',
)
print(dEnt_res_normal)

# =======================================
# SECTION: COMPARE NORMAL VS *_hsg RESULTS
# =======================================
print('--- compare scalar results (normal - hsg) ---')
if dEn_res is not None and dEn_res_normal is not None:
    print('dEn delta:', dEn_res_normal.value - dEn_res.value)
if En_res is not None and En_res_normal is not None:
    print('En delta:', En_res_normal.value - En_res.value)
if En_ref_res is not None and En_ref_res_normal is not None:
    print('En_ref delta:', En_ref_res_normal.value - En_ref_res.value)
if GiFrEn_res is not None and GiFrEn_res_normal is not None:
    print('GiFrEn delta:', GiFrEn_res_normal.value - GiFrEn_res.value)
if dEnt_res is not None and dEnt_res_normal is not None:
    print('dEnt delta:', dEnt_res_normal.value - dEnt_res.value)

print('--- compare range results (normal - hsg) ---')
if En_range_res is not None and En_range_res_normal is not None:
    En_range_deltas = [
        normal.value - hsg.value
        for normal, hsg in zip(En_range_res_normal, En_range_res)
    ]
    print('En_range deltas:', En_range_deltas)

if GiFrEn_range_res is not None and GiFrEn_range_res_normal is not None:
    GiFrEn_range_deltas = [
        normal.value - hsg.value
        for normal, hsg in zip(GiFrEn_range_res_normal, GiFrEn_range_res)
    ]
    print('GiFrEn_range deltas:', GiFrEn_range_deltas)
