# NOTE: import libs
from typing import Dict, List, Mapping, Optional
from pythermodb_settings.models import CustomProp
from pycuc import convert_from_to

# NOTE: type aliases
ComponentAmounts = Mapping[str, CustomProp | float | int]
ComponentMoles = Mapping[str, CustomProp | float | int]


# ! ::: Check if units match
def _same_unit(
        unit: str,
        expected_unit: str
) -> bool:
    return unit.strip().lower() == expected_unit.strip().lower()


# ! ::: Check if all CustomProp values have the expected unit
def _all_valid_units(
        values: Dict[str, CustomProp],
        expected_unit: str
) -> bool:
    if not values:
        return True

    all_units = list(set([value.unit for value in values.values()]))

    # NOTE: more than one unique unit
    if len(all_units) > 1:
        return False

    # NOTE: check if the single unit matches the expected unit
    return _same_unit(all_units[0], expected_unit)


# ! ::: Convert a CustomProp/numeric mapping to scalar values
def _to_custom_props_mapping(
        values: Mapping[str, CustomProp | float | int],
        output_unit: Optional[str] = None
) -> Dict[str, float]:
    converted_values: dict[str, float] = {}

    custom_prop_values = {
        key: value
        for key, value in values.items()
        if isinstance(value, CustomProp)
    }

    # SECTION: Check all CustomProp values already use the output unit
    unit_valid = False

    # NOTE: check
    if output_unit is not None and len(custom_prop_values) == len(values):
        unit_valid = _all_valid_units(custom_prop_values, output_unit)

    # NOTE: expected unit
    if unit_valid is True:
        converted_values = {
            key: float(value.value) for key, value in custom_prop_values.items()
        }

        return converted_values

    # SECTION: Convert values to the desired output unit
    for key, value in values.items():
        if isinstance(value, CustomProp):
            # NOTE: set
            val_ = value.value
            unit_ = value.unit

            # ! to output unit
            if output_unit and not _same_unit(unit_, output_unit):
                val_ = convert_from_to(val_, unit_, output_unit)

            # NOTE: set
            converted_values[key] = float(val_)
        else:
            converted_values[key] = float(value)

    return converted_values


# ! ::: Convert component amounts to the requested output unit
def _to_amounts(
        component_amounts: ComponentAmounts,
        output_unit: Optional[str] = None
) -> Dict[str, float]:
    """
    Convert a dictionary of component amounts to float values.

    Parameters
    ----------
    component_amounts : ComponentAmounts
        A dictionary mapping component names to amounts. Numeric values are assumed to already be in output_unit.
    output_unit : str, optional
        The unit to which CustomProp component amounts should be converted. Default is None.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping component names to their respective amounts as floats.
    """
    return _to_custom_props_mapping(component_amounts, output_unit)


# ! ::: Convert component moles to the requested output unit
def _to_moles(
        component_moles: ComponentMoles,
        output_unit: Optional[str] = None
) -> Dict[str, float]:
    """
    Convert a dictionary of component moles to float values.

    Parameters
    ----------
    component_moles : ComponentMoles
        A dictionary mapping component names to mole amounts. Numeric values are assumed to already be in output_unit.
    output_unit : str, optional
        The unit to which CustomProp component moles should be converted. Default is None.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping component names to their respective moles as floats.
    """
    return _to_custom_props_mapping(component_moles, output_unit)


# ! ::: Convert a CustomProp scalar to the requested output unit
def _to_custom_prop_scalar(
        prop: CustomProp,
        output_unit: Optional[str] = None
) -> float:
    val_ = prop.value
    unit_ = prop.unit

    if output_unit and not _same_unit(unit_, output_unit):
        val_ = convert_from_to(val_, unit_, output_unit)

    return float(val_)


# ! ::: Convert to desired volume unit
def _to_volume(
        solution_volume: CustomProp,
        output_unit: Optional[str] = None
) -> float:
    """
    Convert a solution volume defined as a CustomProp object to a float value in the desired unit.

    Parameters
    ----------
    solution_volume : CustomProp
        The solution volume as a CustomProp object.
    output_unit : str, optional
        The unit to which the solution volume should be converted. Default is None.

    Returns
    -------
    float
        The solution volume in the desired unit as a float.
    """
    return _to_custom_prop_scalar(solution_volume, output_unit)


# ! ::: Convert to desired mass unit
def _to_mass(
        solvent_mass: CustomProp,
        output_unit: Optional[str] = None
) -> float:
    """
    Convert a solvent mass defined as a CustomProp object to a float value in the desired unit.

    Parameters
    ----------
    solvent_mass : CustomProp
        The solvent mass as a CustomProp object.
    output_unit : str, optional
        The unit to which the solvent mass should be converted. Default is None.

    Returns
    -------
    float
        The solvent mass in the desired unit as a float.
    """
    return _to_custom_prop_scalar(solvent_mass, output_unit)


# ! ::: Helper function to split unit string into individual units
def _to_units(unit: str) -> List[str]:
    # NOTE: validation
    if not unit or not isinstance(unit, str):
        raise ValueError("Invalid unit string provided.")

    # NOTE: only support units separated by '/'
    if '/' not in unit:
        raise ValueError("Unit string must contain '/' to separate units.")

    # NOTE: more than one '/'
    if unit.count('/') != 1:
        raise ValueError(
            "Unit string must contain exactly one '/' to separate units."
        )

    # NOTE: split the unit string by '/' to get individual units
    return [unit_.strip() for unit_ in unit.strip().split('/')]
