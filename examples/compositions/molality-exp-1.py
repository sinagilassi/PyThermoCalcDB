from pythermocalcdb.compositions.molality import (
    calc_component_molalities_from_props,
    calc_molalities_from_mapping,
    calc_molalities_from_props,
)
from pythermodb_settings.models import Component, CustomProp
from rich import print


components = [
    Component(name="carbon dioxide", formula="CO2", state="g"),
    Component(name="methane", formula="CH4", state="g"),
    Component(name="oxygen", formula="O2", state="g"),
]

# NOTE: Basic dictionary input with numeric solvent mass in kg.
raw_moles = {"A": 2.0, "B": 3.0}
raw_molality = calc_molalities_from_mapping(
    component_moles=raw_moles,
    solvent_mass=10.0,
)
raw_molality_dict = raw_molality.value

print("Raw molality with numeric solvent mass:")
print(raw_molality_dict)

assert raw_molality_dict == {"A": 0.2, "B": 0.3}

# NOTE: Raw dictionary input with CustomProp mole values and CustomProp solvent mass.
custom_moles = {
    "A": CustomProp(value=2.0, unit="mol"),
    "B": CustomProp(value=3.0, unit="mol"),
}
solvent_mass = CustomProp(value=10.0, unit="kg")
custom_mass_molality = calc_molalities_from_props(
    component_moles=custom_moles,
    solvent_mass=solvent_mass,
)
custom_mass_molality_dict = custom_mass_molality.value

print("Raw molality with CustomProp moles and solvent mass:")
print(custom_mass_molality_dict)

assert custom_mass_molality_dict == raw_molality_dict

# NOTE: Input values are intentionally not in component order.
component_moles = {
    "oxygen": CustomProp(value=2.0, unit="mol"),
    "carbon dioxide": CustomProp(value=1.0, unit="mol"),
    "methane": CustomProp(value=1.0, unit="mol"),
}

component_molality = calc_component_molalities_from_props(
    component_moles=component_moles,
    solvent_mass=CustomProp(value=2.0, unit="kg"),
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

component_molality_dict = component_molality.value

print("Component molality keyed by formula-state:")
print(component_molality_dict)

assert component_molality_dict == {
    "CO2-g": 0.5,
    "CH4-g": 0.5,
    "O2-g": 1.0,
}

# NOTE: With component_key=None and components=None, unit-aware component molality uses raw keys.
raw_component_molality = calc_component_molalities_from_props(
    component_moles=custom_moles,
    solvent_mass=solvent_mass,
    components=None,
)

raw_component_molality_dict = raw_component_molality.value

print("Raw component molality with units:")
print(raw_component_molality_dict)

assert raw_component_molality_dict == raw_molality_dict
