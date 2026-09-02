# import libs
import logging
from typing import Any, Dict, List, Optional, Tuple
from pythermodb_settings.models import Component, ComponentKey, CustomProp
from pythermodb_settings.utils import (
    config_components_values,
)
# locals
from .utils import ComponentMoles, _to_moles, _to_units, _to_volume

# NOTE: logger setup
logger = logging.getLogger(__name__)


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
