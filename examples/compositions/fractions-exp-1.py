from pythermocalcdb.compositions.fractions import fr2, fr3
from pythermodb_settings.models import Component
from rich import print


components = [
    Component(name="carbon dioxide", formula="CO2", state="g"),
    Component(name="methane", formula="CH4", state="g"),
    Component(name="oxygen", formula="O2", state="g"),
]

# NOTE: Input values are intentionally not in component order.
mole_amounts = {
    "oxygen": 2.0,
    "carbon dioxide": 1.0,
    "methane": 1.0,
}

component_fractions = fr3(
    values=mole_amounts,
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

if component_fractions is None:
    raise RuntimeError("Failed to calculate component fractions.")

component_fractions_dict, component_fractions_list = component_fractions

print("Component fractions keyed by formula-state:")
print(component_fractions_dict)
print("Component fractions in component order:")
print(component_fractions_list)

assert component_fractions_dict == {
    "CO2-g": 0.25,
    "CH4-g": 0.25,
    "O2-g": 0.5,
}
assert component_fractions_list == [0.25, 0.25, 0.5]

# NOTE: With component_key=None and components=None, fr3 normalizes like fr2.
raw_values = {"A": 2.0, "B": 3.0}
raw_fr2 = fr2(raw_values)
raw_fr3 = fr3(
    values=raw_values,
)

if raw_fr3 is None:
    raise RuntimeError("Failed to calculate raw fractions.")

raw_fractions_dict, raw_fractions_list = raw_fr3

print("Raw fr2 fractions:")
print(raw_fr2)
print("Raw fr3 fractions:")
print(raw_fractions_dict)
print(raw_fractions_list)

assert raw_fractions_dict == raw_fr2
assert raw_fractions_list == [0.4, 0.6]
