# import libs
import logging
from typing import Any, Dict, List, Mapping, Optional, Tuple
from pythermodb_settings.models import Component, ComponentKey, CustomProp
from pythermodb_settings.utils import config_components_values
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
    val_ = solvent_mass.value
    unit_ = solvent_mass.unit

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


# ! ::: Molality [mol/kg]
def molality1(
        component_moles: List[float],
        solvent_mass: float,
) -> List[float]:
    """
    Calculate the molality of each component in a solution given the component moles and the solvent mass in kg.

    Parameters
    ----------
    component_moles : List[float]
        A list of moles for each component.
    solvent_mass : float
        The solvent mass in kg.

    Returns
    -------
    List[float]
        A list of molality values for each component.
    """
    return [moles / solvent_mass for moles in component_moles]


def molality2(
    component_moles: Dict[str, float | int],
    solvent_mass: float,
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the molality of each component in a solution given the component moles and the solvent mass in kg.

    Parameters
    ----------
    component_moles : Dict[str, float | int]
        A dictionary mapping component names to their respective moles.
    solvent_mass : float
        The solvent mass in kg.

    Returns
    -------
    Tuple[Dict[str, float], List[float]]
        A tuple containing a dictionary of component molalities and a list of molality values.
    """
    # ! dict
    component_molality_dict = {
        key: value / solvent_mass for key, value in component_moles.items()
    }
    # ! list
    component_molality_list = list(component_molality_dict.values())

    return component_molality_dict, component_molality_list


# ! ::: Molality [mol/kg] with solvent mass as CustomProp
def molality3(
    component_moles: ComponentMoles,
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the molality of each component in a solution given the component moles and the solvent mass as a CustomProp.

    Parameters
    ----------
    component_moles : ComponentMoles
        A dictionary mapping component names to their respective moles. Numeric values are assumed to already be in the mole unit from output_unit.
    solvent_mass : CustomProp
        The solvent mass as a CustomProp object.
    output_unit : str, optional
        The unit for the output molality values. Defaults to 'mol/kg'.

    Returns
    -------
    Tuple[Dict[str, float], List[float]]
        A tuple containing a dictionary of component molalities and a list of molality values.

    Notes
    -----
    - The solvent mass is expected to be provided as a CustomProp object. If the output_unit is not specified, it defaults to mol/kg.
    - Numeric component mole values are assumed to already be in the mole unit from output_unit.
    """
    # SECTION: set default units for moles and mass
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    mass_unit = units_[1]

    # NOTE: component moles
    # ! convert component moles to the specified unit if necessary
    component_moles_dict: Dict[str, float] = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
    )

    # NOTE: solvent mass unit should match output unit denominator
    # ! convert solvent mass to the specified unit
    solvent_mass_scalar = _to_mass(
        solvent_mass=solvent_mass,
        output_unit=mass_unit
    )

    # SECTION: calculate molality for each component
    component_molality = molality2(
        component_moles=component_moles_dict,
        solvent_mass=solvent_mass_scalar,
    )
    # unpack the result
    component_molality_dict, component_molality_list = component_molality

    return component_molality_dict, component_molality_list


# ! ::: Molality [mol/kg] with component ID mapping and sorting
def molality4(
    component_moles: ComponentMoles,
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    """
    Calculate the molality of each component in a solution given the component moles and the solvent mass as a CustomProp. The component molality list and dictionary will be ordered according to the components list if sort_by_components_order is True.

    Parameters
    ----------
    component_moles : ComponentMoles
        A dictionary mapping component names to their respective moles or CustomProp objects representing the moles.
    solvent_mass : CustomProp
        The solvent mass as a CustomProp object.
    output_unit : str, optional
        The unit for the output molality values. Defaults to 'mol/kg'.
    components : Optional[List[Component]], optional
        A list of Component objects to map the component moles to, by default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component moles to components, by default None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default True.
    sort_by_components_order : bool, optional
        Whether to sort the component molalities by the order of components, by default True.

    Returns
    -------
    Optional[Tuple[Dict[str, float], List[float]]]
        A tuple containing a dictionary of component molalities and a list of molality values, or None if the calculation could not be performed.
    """
    # SECTION: Unit validation
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    mass_unit = units_[1]

    # SECTION: convert component moles to moles if they are CustomProp objects
    # ! component moles
    component_moles_dict = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
    )

    # ! mass
    solvent_mass_scalar = _to_mass(
        solvent_mass=solvent_mass,
        output_unit=mass_unit
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

    # SECTION: calculate molality for each component
    component_molality = molality2(
        component_moles=component_values_dict,
        solvent_mass=solvent_mass_scalar,
    )
    # unpack
    component_molality_dict, component_molality_list = component_molality

    return component_molality_dict, component_molality_list
