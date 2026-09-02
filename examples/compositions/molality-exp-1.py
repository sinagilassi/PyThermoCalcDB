from pythermocalcdb.compositions.molality import molality2, molality3, molality4
from pythermodb_settings.models import Component, CustomProp


components = [
    Component(name="carbon dioxide", formula="CO2", state="g"),
    Component(name="methane", formula="CH4", state="g"),
    Component(name="oxygen", formula="O2", state="g"),
]

# NOTE: Basic dictionary input with numeric solvent mass in kg.
raw_moles = {"A": 2.0, "B": 3.0}
raw_molality_dict, raw_molality_list = molality2(
    component_moles=raw_moles,
    solvent_mass=10.0,
)

print("Raw molality with numeric solvent mass:")
print(raw_molality_dict)
print(raw_molality_list)

assert raw_molality_dict == {"A": 0.2, "B": 0.3}
assert raw_molality_list == [0.2, 0.3]

# NOTE: Raw dictionary input with CustomProp mole values and CustomProp solvent mass.
custom_moles = {
    "A": CustomProp(value=2.0, unit="mol"),
    "B": CustomProp(value=3.0, unit="mol"),
}
solvent_mass = CustomProp(value=10.0, unit="kg")
custom_mass_molality_dict, custom_mass_molality_list = molality3(
    component_moles=custom_moles,
    solvent_mass=solvent_mass,
)

print("Raw molality with CustomProp moles and solvent mass:")
print(custom_mass_molality_dict)
print(custom_mass_molality_list)

assert custom_mass_molality_dict == raw_molality_dict
assert custom_mass_molality_list == raw_molality_list

# NOTE: Numeric mole values still work and are assumed to be in the output mole unit.
numeric_molality3_dict, numeric_molality3_list = molality3(
    component_moles=raw_moles,
    solvent_mass=solvent_mass,
)

assert numeric_molality3_dict == raw_molality_dict
assert numeric_molality3_list == raw_molality_list

# NOTE: Input values are intentionally not in component order.
component_moles = {
    "oxygen": CustomProp(value=2.0, unit="mol"),
    "carbon dioxide": CustomProp(value=1.0, unit="mol"),
    "methane": CustomProp(value=1.0, unit="mol"),
}

component_molality = molality4(
    component_moles=component_moles,
    solvent_mass=CustomProp(value=2.0, unit="kg"),
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

if component_molality is None:
    raise RuntimeError("Failed to calculate component molality.")

component_molality_dict, component_molality_list = component_molality

print("Component molality keyed by formula-state:")
print(component_molality_dict)
print("Component molality in component order:")
print(component_molality_list)

assert component_molality_dict == {
    "CO2-g": 0.5,
    "CH4-g": 0.5,
    "O2-g": 1.0,
}
assert component_molality_list == [0.5, 0.5, 1.0]

# NOTE: With component_key=None and components=None, molality4 uses raw keys.
raw_molality4 = molality4(
    component_moles=raw_moles,
    solvent_mass=solvent_mass,
    components=None,
)

if raw_molality4 is None:
    raise RuntimeError("Failed to calculate raw molality.")

raw_molality4_dict, raw_molality4_list = raw_molality4

print("Raw molality4 values:")
print(raw_molality4_dict)
print(raw_molality4_list)

assert raw_molality4_dict == raw_molality_dict
assert raw_molality4_list == raw_molality_list
