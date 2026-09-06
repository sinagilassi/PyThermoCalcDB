from typing import List, Dict
from pythermocalcdb.compositions.molarity import (
    _calc_molarities_from_sequence,
    calc_molarities_from_sequence,
    _calc_molarities_from_mapping,
    calc_molarities_from_mapping,
    _calc_molarities_from_props,
    calc_molarities_from_props,
    _calc_component_molarities_from_props,
    calc_component_molarities_from_props,
)
from pythermodb_settings.models import Component, CustomProp, AnnotatedValue
from rich import print

components = [
    Component(name="carbon dioxide", formula="CO2", state="g"),
    Component(name="methane", formula="CH4", state="g"),
    Component(name="oxygen", formula="O2", state="g"),
]

# NOTE: Basic dictionary input with numeric solution volume.
raw_moles = {"A": 2.0, "B": 3.0}
raw_moles_list = list(raw_moles.values())

# ! from mapping
raw_molarity_dict: AnnotatedValue[Dict[str, float]] = calc_molarities_from_mapping(
    component_moles=raw_moles,
    solution_volume=10.0,
)
print(f"[blue]Raw molarity from mapping using AnnotatedValue:[/blue]")
print(raw_molarity_dict)

# ! from mapping normal return
raw_molarity_dict_normal = _calc_molarities_from_mapping(
    component_moles=raw_moles,
    solution_volume=10.0,
)
print(f"[yellow]Raw molarity from mapping normal return:[/yellow]")
print(raw_molarity_dict_normal)

# ! from sequence
raw_molarity_list: AnnotatedValue[List[float]] = calc_molarities_from_sequence(
    component_moles=raw_moles_list,
    solution_volume=10.0,
)
print(f"[blue]Raw molarity from sequence using AnnotatedValue:[/blue]")
print(raw_molarity_list)

# ! from sequence normal return

raw_molarity_list_normal = _calc_molarities_from_sequence(
    component_moles=raw_moles_list,
    solution_volume=10.0,
)
print(f"[yellow]Raw molarity from sequence normal return:[/yellow]")
print(raw_molarity_list_normal)


# NOTE: Raw dictionary input with CustomProp mole values and CustomProp volume.
custom_moles = {
    "A": CustomProp(value=2.0, unit="mol"),
    "B": CustomProp(value=3.0, unit="mol"),
}
solution_volume = CustomProp(value=10.0, unit="L")

custom_volume_molarity: Dict[str, float] = _calc_molarities_from_props(
    component_moles=custom_moles,
    solution_volume=solution_volume,
    output_unit="mol/L",
)
print(f"[green]Raw molarity with CustomProp moles and volume using normal return:[/green]")
print(custom_volume_molarity)

# ! annotated
custom_volume_molarity_dict: AnnotatedValue[Dict[str, float]] = calc_molarities_from_props(
    component_moles=custom_moles,
    solution_volume=solution_volume,
)


print(f"[lime]Raw molarity with CustomProp moles and volume using AnnotatedValue:[/lime]")
print(custom_volume_molarity_dict)

# NOTE: Input values are intentionally not in component order.
component_moles = {
    "oxygen": CustomProp(value=2.0, unit="mol"),
    "carbon dioxide": CustomProp(value=1.0, unit="mol"),
    "methane": CustomProp(value=1.0, unit="mol"),
}

component_molarity = _calc_component_molarities_from_props(
    component_moles=component_moles,
    solution_volume=CustomProp(value=2.0, unit="L"),
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

if component_molarity is None:
    raise RuntimeError("Failed to calculate component molarity.")

print(component_molarity)

# NOTE: With component_key=None and components=None, unit-aware component molarity uses raw keys.
raw_molarity4 = calc_component_molarities_from_props(
    component_moles=custom_moles,
    solution_volume=solution_volume,
    components=None,
)

if raw_molarity4 is None:
    raise RuntimeError("Failed to calculate raw molarity.")

print(raw_molarity4)


# NOTE: Annotated component molarity
raw_molarity5: AnnotatedValue | None = calc_component_molarities_from_props(
    component_moles=component_moles,
    solution_volume=CustomProp(value=2.0, unit="L"),
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

if raw_molarity5 is None:
    raise RuntimeError("Failed to calculate annotated component molarity.")

print("Annotated component molarity with units:")
print(raw_molarity5)

# ! output unit
raw_molarity6: AnnotatedValue | None = calc_component_molarities_from_props(
    component_moles=component_moles,
    solution_volume=CustomProp(value=2.0, unit="L"),
    output_unit="kmol/L",
    components=components,
    component_key="Formula-State",
    case_sensitive=False,
    sort_by_components_order=True,
)

if raw_molarity6 is None:
    raise RuntimeError(
        "Failed to calculate annotated component molarity with output unit.")

print("Annotated component molarity with output unit:")
print(raw_molarity6)
