# import libs
import logging
from typing import Any, Dict, List, Mapping, Optional, Tuple
from pythermodb_settings.models import Component, ComponentKey, CustomProp
from pythermodb_settings.utils import config_components_values
from pycuc import convert_from_to
# locals

# NOTE: logger setup
logger = logging.getLogger(__name__)

ComponentAmounts = Mapping[str, CustomProp | float | int]


# ! ::: Check if all component amounts have the expected unit
def _same_unit(
        unit: str,
        expected_unit: str
) -> bool:
    return unit.strip().lower() == expected_unit.strip().lower()


def _all_valid_units(
        values: Dict[str, CustomProp],
        expected_unit: str
) -> bool:
    if not values:
        return True

    all_units = list(set([value.unit for value in values.values()]))

    # more than one unique unit
    if len(all_units) > 1:
        return False

    # check if the single unit matches the expected unit
    return _same_unit(all_units[0], expected_unit)


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
    val_ = solution_volume.value
    unit_ = solution_volume.unit

    if output_unit and not _same_unit(unit_, output_unit):
        val_ = convert_from_to(val_, unit_, output_unit)

    return float(val_)


# ! ::: Convert component amounts to the requested numerator unit
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
    converted_amounts: dict[str, float] = {}

    custom_prop_amounts = {
        key: value
        for key, value in component_amounts.items()
        if isinstance(value, CustomProp)
    }

    # SECTION: Check all CustomProp component amounts already use the output unit
    unit_valid = False

    # >> check
    if output_unit is not None and len(custom_prop_amounts) == len(component_amounts):
        unit_valid = _all_valid_units(custom_prop_amounts, output_unit)

    # NOTE: expected unit
    if unit_valid is True:
        converted_amounts = {
            key: float(value.value) for key, value in custom_prop_amounts.items()
        }

        return converted_amounts

    # SECTION: Convert component amounts to the desired output unit
    for key, value in component_amounts.items():
        if isinstance(value, CustomProp):
            # set
            val_ = value.value
            unit_ = value.unit

            # ! to output unit
            if output_unit and not _same_unit(unit_, output_unit):
                val_ = convert_from_to(val_, unit_, output_unit)

            # set
            converted_amounts[key] = float(val_)
        else:
            converted_amounts[key] = float(value)

    return converted_amounts


# ! ::: Helper function to split unit string into individual units
def _to_units(unit: str) -> List[str]:
    # NOTE: validation
    if not unit or not isinstance(unit, str):
        raise ValueError("Invalid unit string provided.")

    # only support units separated by '/'
    if '/' not in unit:
        raise ValueError("Unit string must contain '/' to separate units.")

    # more than one '/'
    if unit.count('/') != 1:
        raise ValueError(
            "Unit string must contain exactly one '/' to separate units."
        )

    # NOTE: split the unit string by '/' to get individual units
    return [unit_.strip() for unit_ in unit.strip().split('/')]


# ! ::: Concentration [amount/volume]
def concentration_amount_volume_1(
        component_amounts: List[float],
        solution_volume: float,
) -> List[float]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume.

    Parameters
    ----------
    component_amounts : List[float]
        A list of amounts for each component.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    List[float]
        A list of concentration values for each component.
    """
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")

    return [amount / solution_volume for amount in component_amounts]


def concentration_amount_volume_2(
    component_amounts: Dict[str, float | int],
    solution_volume: float,
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume.

    Parameters
    ----------
    component_amounts : Dict[str, float | int]
        A dictionary mapping component names to their respective amounts.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    Tuple[Dict[str, float], List[float]]
        A tuple containing a dictionary of component concentrations and a list of concentration values.
    """
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")

    # ! dict
    component_concentration_dict = {
        key: value / solution_volume for key, value in component_amounts.items()
    }
    # ! list
    component_concentration_list = list(component_concentration_dict.values())

    return component_concentration_dict, component_concentration_list


# ! ::: Concentration [amount/volume] with solution volume as CustomProp
def concentration_amount_volume_3(
    component_amounts: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume as a CustomProp.

    Parameters
    ----------
    component_amounts : ComponentAmounts
        A dictionary mapping component names to their respective amounts. Numeric values are assumed to already be in the numerator unit from output_unit.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output concentration values. Defaults to 'kg/m^3'.

    Returns
    -------
    Tuple[Dict[str, float], List[float]]
        A tuple containing a dictionary of component concentrations and a list of concentration values.
    """
    # SECTION: set default units for amount and volume
    units_ = _to_units(output_unit)
    # >> set
    amount_unit = units_[0]
    volume_unit = units_[1]

    # NOTE: component amounts
    component_amounts_dict: Dict[str, float] = _to_amounts(
        component_amounts=component_amounts,
        output_unit=amount_unit
    )

    # NOTE: solution volume unit should match output unit denominator
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit
    )

    # SECTION: calculate concentration for each component
    component_concentration = concentration_amount_volume_2(
        component_amounts=component_amounts_dict,
        solution_volume=solution_volume_scalar,
    )
    # unpack the result
    component_concentration_dict, component_concentration_list = component_concentration

    return component_concentration_dict, component_concentration_list


# ! ::: Concentration [amount/volume] with component ID mapping and sorting
def concentration_amount_volume_4(
    component_amounts: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume as a CustomProp. The component concentration list and dictionary will be ordered according to the components list if sort_by_components_order is True.

    Parameters
    ----------
    component_amounts : ComponentAmounts
        A dictionary mapping component names to their respective amounts or CustomProp objects representing the amounts.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output concentration values. Defaults to 'kg/m^3'.
    components : Optional[List[Component]], optional
        A list of Component objects to map the component amounts to, by default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component amounts to components, by default None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default True.
    sort_by_components_order : bool, optional
        Whether to sort the component concentrations by the order of components, by default True.

    Returns
    -------
    Optional[Tuple[Dict[str, float], List[float]]]
        A tuple containing a dictionary of component concentrations and a list of concentration values, or None if the calculation could not be performed.
    """
    # SECTION: Unit validation
    units_ = _to_units(output_unit)
    # >> set
    amount_unit = units_[0]
    volume_unit = units_[1]

    # SECTION: convert component amounts if they are CustomProp objects
    # ! component amounts
    component_amounts_dict = _to_amounts(
        component_amounts=component_amounts,
        output_unit=amount_unit
    )

    # ! volume
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit
    )

    # SECTION: get component values
    if component_key is not None:
        if not components:
            logger.warning(
                "Component key is provided but components list is empty.")
            return None

        component_values: Tuple[
            Dict[str, Any],
            List[Any]
        ] | None = config_components_values(
            values=component_amounts_dict,
            components=components,
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order
        )
        # >> check
        if component_values is None:
            logger.warning("Failed to configure component values")
            return None

        # unpack
        component_values_dict, _ = component_values
    else:
        component_values_dict = component_amounts_dict

    # SECTION: calculate concentration for each component
    component_concentration = concentration_amount_volume_2(
        component_amounts=component_values_dict,
        solution_volume=solution_volume_scalar,
    )
    # unpack
    component_concentration_dict, component_concentration_list = component_concentration

    return component_concentration_dict, component_concentration_list


# ! ::: Mass concentration wrappers
def mass_concentration1(
        component_mass: List[float],
        solution_volume: float,
) -> List[float]:
    return concentration_amount_volume_1(
        component_amounts=component_mass,
        solution_volume=solution_volume,
    )


def mass_concentration2(
    component_mass: Dict[str, float | int],
    solution_volume: float,
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_2(
        component_amounts=component_mass,
        solution_volume=solution_volume,
    )


def mass_concentration3(
    component_mass: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_3(
        component_amounts=component_mass,
        solution_volume=solution_volume,
        output_unit=output_unit,
    )


def mass_concentration4(
    component_mass: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    return concentration_amount_volume_4(
        component_amounts=component_mass,
        solution_volume=solution_volume,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )


# ! ::: Molar concentration wrappers
def molar_concentration1(
        component_moles: List[float],
        solution_volume: float,
) -> List[float]:
    return concentration_amount_volume_1(
        component_amounts=component_moles,
        solution_volume=solution_volume,
    )


def molar_concentration2(
    component_moles: Dict[str, float | int],
    solution_volume: float,
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_2(
        component_amounts=component_moles,
        solution_volume=solution_volume,
    )


def molar_concentration3(
    component_moles: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_3(
        component_amounts=component_moles,
        solution_volume=solution_volume,
        output_unit=output_unit,
    )


def molar_concentration4(
    component_moles: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    return concentration_amount_volume_4(
        component_amounts=component_moles,
        solution_volume=solution_volume,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )


# SECTION: Aliases
# ! generic concentration
calculate_concentrations = concentration_amount_volume_1
calculate_keyed_concentrations = concentration_amount_volume_2
calculate_component_concentrations = concentration_amount_volume_4
