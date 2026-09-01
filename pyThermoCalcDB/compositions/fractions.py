# import libs
import logging
from typing import List, Dict, Optional, Tuple, Any
from pythermodb_settings.models import Component, ComponentKey
from pythermodb_settings.utils import (
    find_component_by_id,
    set_component_id,
    config_components_values,
    measure_time
)
# locals

# NOTE: logger set
logger = logging.getLogger(__name__)

# ! ::: Calculate Fractions


@measure_time
def fr1(
        values: List[float | int],
        **kwargs
) -> List[float] | None:
    """
    Calculate the fractions of a list of values.

    Parameters
    ----------
    values : List[float | int]
        A list of numerical values.
    kwargs : dict
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    List[float] | None
        A list of fractions corresponding to the input values, or None if the total is zero.
    """
    if any(value < 0 for value in values):
        logger.warning("Negative values found")
        return None

    # total mole
    total = sum(values)

    # >> check
    if total == 0:
        return None

    # calculate fractions
    return [value / total for value in values]


@measure_time
def fr2(
        values: Dict[str, float | int],
        **kwargs
) -> Optional[Dict[str, float]]:
    """
    Calculate the fractions of a dictionary of values.

    Parameters
    ----------
    values : Dict[str, float | int]
        A dictionary of component IDs and their corresponding values.
    kwargs : dict
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    Optional[Dict[str, float]]
        A dictionary of component fractions, or None if the total is zero or negative values are found.
    """
    if any(value < 0 for value in values.values()):
        logger.warning("Negative values found")
        return None

    total = sum(values.values())
    if total == 0:
        return None

    return {comp_id: value / total for comp_id, value in values.items()}


@measure_time
def fr3(
        values: Dict[str, float | int],
        components: Optional[List[Component]] = None,
        component_key: Optional[ComponentKey] = None,
        case_sensitive: bool = True,
        sort_by_components_order: bool = True,
        **kwargs
) -> Optional[tuple[Dict[str, float], List[float]]]:
    """
    Calculate the fractions of components based on their values. When component_key is provided, the component IDs are matched against the provided list of Component objects. When component_key is None, the values are normalized directly, similar to fr2.

    Parameters
    ----------
    values : Dict[str, float | int]
        A dictionary of component IDs and their corresponding values.
    components : Optional[List[Component]], optional
        A list of Component objects. Required when component_key is provided.
    component_key : Optional[ComponentKey], optional
        The key to use for identifying components. Defaults to None.
    case_sensitive : bool, optional
        Whether the component IDs are case-sensitive. Defaults to True.
    kwargs : dict
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    Optional[tuple[Dict[str, float], List[float]]]
        A tuple containing a dictionary of component fractions and a list of component fractions, or None if the input is invalid.

    """
    # SECTION: validate input
    if not values:
        logger.warning("No values provided")
        return None

    # positive values check
    if any(value < 0 for value in values.values()):
        logger.warning("Negative values found")
        return None

    # SECTION: get components values
    # ! configure component values if component_key is provided otherwise
    # ! otherwise use the original values dictionary
    if component_key is not None:
        if not components:
            logger.warning("No components provided")
            return None

        component_values: Tuple[
            Dict[str, Any],
            List[Any]
        ] | None = config_components_values(
            values=values,
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
        component_values_dict = values

    # SECTION: Calculate fractions
    total = sum(component_values_dict.values())
    if total == 0:
        return None

    # component fractions
    # ! dict
    component_fractions: Dict[str, float] = {
        comp_id: value / total for comp_id,
        value in component_values_dict.items()
    }
    # ! list
    component_fractions_list = list(component_fractions.values())

    return component_fractions, component_fractions_list


# SECTION: Aliases
# ! list
calculate_fractions = fr1
# ! dict
calculate_keyed_fractions = fr2
# ! reorder & create keys
calculate_component_fractions = fr3

