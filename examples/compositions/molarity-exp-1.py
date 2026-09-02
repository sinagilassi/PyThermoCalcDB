from pythermocalcdb.compositions.molarity import molarity2, molarity3, molarity4
from pythermodb_settings.models import Component, CustomProp


components = [
    Component(name="carbon dioxide", formula="CO2", state="g"),
    Component(name="methane", formula="CH4", state="g"),
    Component(name="oxygen", formula="O2", state="g"),
]

# NOTE: Basic dictionary input with numeric solution volume.
raw_moles = {"A": 2.0, "B": 3.0}
raw_molarity_dict, raw_molarity_list = molarity2(
    component_moles=raw_moles,
    solution_volume=10.0,
)

print("Raw molarity with numeric volume:")
print(raw_molarity_dict)
print(raw_molarity_list)

assert raw_molarity_dict == {"A": 0.2, "B": 0.3}
assert raw_molarity_list == [0.2, 0.3]

# NOTE: Raw dictionary input with CustomProp mole values and CustomProp volume.
custom_moles = {
    "A": CustomProp(value=2.0, unit="mol"),
    "B": CustomProp(value=3.0, unit="mol"),
}
solution_volume = CustomProp(value=10.0, unit="L")
custom_volume_molarity_dict, custom_volume_molarity_list = molarity3(
    component_moles=custom_moles,
    solution_volume=solution_volume,
)

print("Raw molarity with CustomProp moles and volume:")
print(custom_volume_molarity_dict)
print(custom_volume_molarity_list)

assert custom_volume_molarity_dict == raw_molarity_dict
assert custom_volume_molarity_list == raw_molarity_list

# NOTE: Numeric mole values still work and are assumed to be in the output mole unit.
numeric_molarity3_dict, numeric_molarity3_list = molarity3(
    component_moles=raw_moles,
    solution_volume=solution_volume,
)

assert numeric_molarity3_dict == raw_molarity_dict
assert numeric_molarity3_list == raw_molarity_list

# NOTE: Input values are intentionally not in component order.
component_moles = {
    "oxygen": CustomProp(value=2.0, unit="mol"),
    "carbon dioxide": CustomProp(value=1.0, unit="mol"),
    "methane": CustomProp(value=1.0, unit="mol"),
}

component_molarity = molarity4(
    component_moles=component_moles,
    solution_volume=CustomProp(value=2.0, unit="L"),
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

if component_molarity is None:
    raise RuntimeError("Failed to calculate component molarity.")

component_molarity_dict, component_molarity_list = component_molarity

print("Component molarity keyed by formula-state:")
print(component_molarity_dict)
print("Component molarity in component order:")
print(component_molarity_list)

assert component_molarity_dict == {
    "CO2-g": 0.5,
    "CH4-g": 0.5,
    "O2-g": 1.0,
}
assert component_molarity_list == [0.5, 0.5, 1.0]

# NOTE: With component_key=None and components=None, molarity4 uses raw keys.
raw_molarity4 = molarity4(
    component_moles=raw_moles,
    solution_volume=solution_volume,
    components=None,
)

if raw_molarity4 is None:
    raise RuntimeError("Failed to calculate raw molarity.")

raw_molarity4_dict, raw_molarity4_list = raw_molarity4

print("Raw molarity4 values:")
print(raw_molarity4_dict)
print(raw_molarity4_list)

assert raw_molarity4_dict == raw_molarity_dict
assert raw_molarity4_list == raw_molarity_list
