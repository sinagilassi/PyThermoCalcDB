# import libs
import logging
from typing import Any, Dict, List, Mapping, Optional, Tuple
from pythermodb_settings.models import Component, ComponentKey, CustomProp
from pythermodb_settings.utils import (
    config_components_values,
)
from pycuc import convert_from_to
# locals

# NOTE: logger setup
logger = logging.getLogger(__name__)

ComponentMoles = Mapping[str, CustomProp | float | int]

# ! ::: Check if all component moles have the expected unit


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

# ! ::: Convert component moles to moles [mol]


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
    converted_moles: dict[str, float] = {}

    custom_prop_moles = {
        key: value
        for key, value in component_moles.items()
        if isinstance(value, CustomProp)
    }

    # SECTION: Check all CustomProp component moles already use the output unit
    unit_valid = False

    # >> check
    if output_unit is not None and len(custom_prop_moles) == len(component_moles):
        unit_valid = _all_valid_units(custom_prop_moles, output_unit)

    # NOTE: expected unit
    if unit_valid is True:
        converted_moles = {
            key: float(value.value) for key, value in custom_prop_moles.items()
        }

        return converted_moles

    # SECTION: Convert component moles to the desired output unit

    # iterate over component moles
    for key, value in component_moles.items():
        if isinstance(value, CustomProp):
            # set
            val_ = value.value
            unit_ = value.unit

            # ! to output unit
            if output_unit and not _same_unit(unit_, output_unit):
                val_ = convert_from_to(val_, unit_, output_unit)

            # set
            converted_moles[key] = float(val_)
        else:
            converted_moles[key] = float(value)

    return converted_moles


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


# ! ::: Molarity [m_i]


def molarity1(
        component_moles: List[float],
        solution_volume: float,
) -> List[float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume.

    Parameters
    ----------
    component_moles : List[float]
        A list of moles for each component.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    List[float]
        A list of molarity values for each component.
    """
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")
    return [moles / solution_volume for moles in component_moles]


def molarity2(
    component_moles: Dict[str, float | int],
    solution_volume: float,
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume.

    Parameters
    ----------
    component_moles : Dict[str, float | int]
        A dictionary mapping component names to their respective moles.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    List[float]
        A list of molarity values for each component.
    """
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")

    # ! dict
    component_molarity_dict = {
        key: value / solution_volume for key, value in component_moles.items()
    }
    # ! list
    component_molarity_list = list(component_molarity_dict.values())

    return component_molarity_dict, component_molarity_list


# ! ::: Molarity [m_i] with solution volume as CustomProp
def molarity3(
    component_moles: ComponentMoles,
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume as a CustomProp. The default
    volume unit is litre (L).

    Parameters
    ----------
    component_moles : ComponentMoles
        A dictionary mapping component names to their respective moles. Numeric values are assumed to already be in the mole unit from output_unit.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output molarity values. Defaults to 'mol/L'.

    Returns
    -------
    Tuple[Dict[str, float], List[float]]
        A tuple containing a dictionary of component molarities and a list of molarity values.

    Notes
    -----
    - The solution volume is expected to be provided as a CustomProp object. If the output_unit is not specified, it defaults to mol/L.
    - Numeric component mole values are assumed to already be in the mole unit from output_unit.
    """
    # SECTION: set default units for moles and volume
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    volume_unit = units_[1]

    # NOTE: component moles
    # ! convert component moles to the specified unit if necessary
    component_moles_dict: Dict[str, float] = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
    )

    # NOTE: solution volume unit should match output unit denominator
    # ! convert solution volume to the specified unit
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit
    )

    # SECTION: calculate molarity for each component
    component_molarity = molarity2(
        component_moles=component_moles_dict,
        solution_volume=solution_volume_scalar,
    )
    # unpack the result
    component_molarity_dict, component_molarity_list = component_molarity

    return component_molarity_dict, component_molarity_list

# ! ::: Molarity [m_i] with component ID mapping and sorting


def molarity4(
    component_moles: ComponentMoles,
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume as a CustomProp. The default
    volume unit is litre (L). The component molarity list and dictionary will be ordered according to the components list if sort_by_components_order is True.

    Parameters
    ----------
    component_moles : ComponentMoles
        A dictionary mapping component names to their respective moles or CustomProp objects representing the moles.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output molarity values. Defaults to 'mol/L'.
    components : Optional[List[Component]], optional
        A list of Component objects to map the component moles to, by default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component moles to components, by default None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default True.
    sort_by_components_order : bool, optional
        Whether to sort the component molarities by the order of components, by default True.

    Returns
    -------
    Optional[Tuple[Dict[str, float], List[float]]]
        A tuple containing a dictionary of component molarities and a list of molarity values, or None if the calculation could not be performed.

    """
    # SECTION: Unit validation
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    volume_unit = units_[1]

    # SECTION: convert component moles to moles if they are CustomProp objects
    # ! component moles
    component_moles_dict = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
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
            values=component_moles_dict,
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
        component_values_dict = component_moles_dict

    # SECTION: calculate molarity for each component
    component_molarity = molarity2(
        component_moles=component_values_dict,
        solution_volume=solution_volume_scalar,
    )
    # unpack
    component_molarity_dict, component_molarity_list = component_molarity

    return component_molarity_dict, component_molarity_list
