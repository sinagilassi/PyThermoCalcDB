# import libs
from pythermocalcdb.compositions.concentration import (
    concentration_amount_volume_2,
    concentration_amount_volume_3,
    concentration_amount_volume_4,
    molar_concentration3,
)
from pythermodb_settings.models import Component, CustomProp
from rich import print


components = [
    Component(name="carbon dioxide", formula="CO2", state="g"),
    Component(name="methane", formula="CH4", state="g"),
    Component(name="oxygen", formula="O2", state="g"),
]

# NOTE: Basic dictionary input with numeric solution volume.
raw_masses = {"A": 2.0, "B": 3.0}
raw_concentration_dict, raw_concentration_list = concentration_amount_volume_2(
    component_amounts=raw_masses,
    solution_volume=10.0,
)

print("Raw mass concentration with numeric volume:")
print(raw_concentration_dict)
print(raw_concentration_list)

assert raw_concentration_dict == {"A": 0.2, "B": 0.3}
assert raw_concentration_list == [0.2, 0.3]

# NOTE: Raw dictionary input with CustomProp mass values and CustomProp volume.
custom_masses = {
    "A": CustomProp(value=2.0, unit="g"),
    "B": CustomProp(value=3.0, unit="g"),
}
solution_volume = CustomProp(value=10.0, unit="L")
custom_volume_concentration_dict, custom_volume_concentration_list = concentration_amount_volume_3(
    component_amounts=custom_masses,
    solution_volume=solution_volume,
    output_unit="g/L",
)

print("Raw mass concentration with CustomProp masses and volume:")
print(custom_volume_concentration_dict)
print(custom_volume_concentration_list)

assert custom_volume_concentration_dict == raw_concentration_dict
assert custom_volume_concentration_list == raw_concentration_list

# NOTE: Numeric amount values still work and are assumed to be in the output numerator unit.
numeric_concentration3_dict, numeric_concentration3_list = concentration_amount_volume_3(
    component_amounts=raw_masses,
    solution_volume=solution_volume,
    output_unit="g/L",
)

assert numeric_concentration3_dict == raw_concentration_dict
assert numeric_concentration3_list == raw_concentration_list

# NOTE: Input values are intentionally not in component order.
component_masses = {
    "oxygen": CustomProp(value=2.0, unit="g"),
    "carbon dioxide": CustomProp(value=1.0, unit="g"),
    "methane": CustomProp(value=1.0, unit="g"),
}

component_concentration = concentration_amount_volume_4(
    component_amounts=component_masses,
    solution_volume=CustomProp(value=2.0, unit="L"),
    output_unit="g/L",
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

if component_concentration is None:
    raise RuntimeError("Failed to calculate component concentration.")

component_concentration_dict, component_concentration_list = component_concentration

print("Component mass concentration keyed by formula-state:")
print(component_concentration_dict)
print("Component mass concentration in component order:")
print(component_concentration_list)

assert component_concentration_dict == {
    "CO2-g": 0.5,
    "CH4-g": 0.5,
    "O2-g": 1.0,
}
assert component_concentration_list == [0.5, 0.5, 1.0]

# NOTE: With component_key=None and components=None, concentration4 uses raw keys.
raw_concentration4 = concentration_amount_volume_4(
    component_amounts=raw_masses,
    solution_volume=solution_volume,
    output_unit="g/L",
    components=None,
)

if raw_concentration4 is None:
    raise RuntimeError("Failed to calculate raw concentration.")

raw_concentration4_dict, raw_concentration4_list = raw_concentration4

print("Raw concentration4 values:")
print(raw_concentration4_dict)
print(raw_concentration4_list)

assert raw_concentration4_dict == raw_concentration_dict
assert raw_concentration4_list == raw_concentration_list

# NOTE: Molar concentration uses the same amount/volume machinery with molar units.
custom_moles = {
    "A": CustomProp(value=2.0, unit="mol"),
    "B": CustomProp(value=3.0, unit="mol"),
}
molar_concentration_dict, molar_concentration_list = molar_concentration3(
    component_moles=custom_moles,
    solution_volume=solution_volume,
    output_unit="mol/L",
)

print("Raw molar concentration with CustomProp moles and volume:")
print(molar_concentration_dict)
print(molar_concentration_list)

assert molar_concentration_dict == raw_concentration_dict
assert molar_concentration_list == raw_concentration_list
