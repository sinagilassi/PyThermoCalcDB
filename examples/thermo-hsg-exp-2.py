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
    CustomProperty
)
# local
from pythermocalcdb.docs.thermo import (
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
    calc_En_IG_ref_hsg_plus,
    calc_En_LIQ_ref_hsg_plus
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

# ! ethanol
CH3OH_thermodb_file = os.path.join(thermodb_dir, 'methanol.pkl')

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

CH3OH = Component(
    name='methanol',
    formula='CH3OH',
    state='g',
)

CH3OH_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CH3OH,
    source=CH3OH_thermodb_file,
)

_component_thermodb: list = [
    CO2_component_thermodb,
    CO_component_thermodb,
    CH3OH_component_thermodb,
]

# =======================================
# SECTION: BUILD THERMODB MODEL SOURCE
# =======================================
# NOTE: build model source from thermodb sources
model_source_: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=_component_thermodb,
    rules=None,  # use default rules
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

# =======================================
# SECTION: Calculate enthalpy at temperature with reference enthalpy at 298.15 K and 1 atm
# =======================================
# NOTE: enthalpy at temperature with reference enthalpy at 298.15 K and 1 atm
EnFo_IG = hsg_props.EnFo_IG
if EnFo_IG is None:
    raise ValueError(
        'Reference enthalpy (EnFo_IG) is not available in HSGProperties.')
print(EnFo_IG)

En_IG_ref_hsg_plus_res = calc_En_IG_ref_hsg_plus(
    hsg_props=hsg_props,
    EnFo_IG=EnFo_IG,
    temperature=temperature_,
    mode='log',
)
print(En_IG_ref_hsg_plus_res)


# NOTE: enthalpy at temperature with reference enthalpy at 298.15 K and 1 atm
EnFo_LIQ = hsg_props.EnFo_LIQ
if EnFo_LIQ is None:
    # set
    EnFo_LIQ = CustomProperty(
        name='Enthalpy of formation of liquid phase at 298 K',
        description='Estimated enthalpy of formation for the liquid phase at 298 K, used as a reference for calculating liquid phase enthalpy at other temperatures.',
        value=EnFo_IG.value + 100,
        unit=EnFo_IG.unit,
        symbol='EnFo_LIQ'
    )
print(EnFo_LIQ)

En_LIQ_ref_hsg_plus_res = calc_En_LIQ_ref_hsg_plus(
    hsg_props=hsg_props,
    EnFo_LIQ=EnFo_LIQ,
    temperature=temperature_,
    mode='log',
)
print(En_LIQ_ref_hsg_plus_res)
