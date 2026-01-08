# import libs
from rich import print
from pyThermoCalcDB.reactions.reactions import Reaction, dEn_rxn_STD
from pythermodb_settings.models import Component, Temperature
import pyThermoDB as ptdb
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource
from typing import Dict
import os


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
CO_thermodb_file = os.path.join(thermodb_dir, 'carbon monoxide-g.pkl')

# ! CH4
CH4_thermodb_file = os.path.join(thermodb_dir, 'methane-g.pkl')

# ! H2O
H2O_thermodb_file = os.path.join(thermodb_dir, 'water-g.pkl')

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

# ! CH4
CH4_comp = Component(
    name='methane',
    formula='CH4',
    state='g'
)
# ! component thermodb source
CH4_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=CH4_comp,
    source=CH4_thermodb_file
)

# ! H2O
H2O_comp = Component(
    name='water',
    formula='H2O',
    state='g'
)
# ! component thermodb source
H2O_component_thermodb: ComponentThermoDBSource = ComponentThermoDBSource(
    component=H2O_comp,
    source=H2O_thermodb_file
)

# NOTE: component thermodb source
_component_thermodb: list = [
    CO2_component_thermodb,
    CO_component_thermodb,
    CH4_component_thermodb,
    H2O_component_thermodb
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

# SECTION: Components
CH4 = Component(
    name="Methane",
    formula="CH4",
    state="g",
)
O2 = Component(
    name="Oxygen",
    formula="O2",
    state="g",
)
CO2 = Component(
    name="Carbon Dioxide",
    formula="CO2",
    state="g",
)
H2O = Component(
    name="Water",
    formula="H2O",
    state="g",
)

components = [CH4, O2, CO2, H2O]

# SECTION: Define reaction
reaction = Reaction(
    name="Combustion of Methane",
    reaction="CH4(g) + 2O2(g) -> CO2(g) + 2H2O(g)",
    components=components,
)

# SECTION: Calculate standard enthalpy change of reaction at 298.15 K
# NOTE: temperature
T = Temperature(value=298.15, unit="K")

# NOTE: calculate standard enthalpy change of reaction
dH_rxn = dEn_rxn_STD(
    reaction=reaction,
    temperature=T,
    model_source=model_source_
)
print(dH_rxn)
