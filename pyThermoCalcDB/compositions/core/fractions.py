# import libs
import logging
import numpy as np
from numpy.typing import NDArray
from collections.abc import Mapping, Sequence
from typing import List, Dict, Optional, Tuple, Any, cast
from pythermodb_settings.models import Component, ComponentKey
# locals
from pythermodb_settings.utils import (
    config_components_values,
)

# NOTE: logger set
logger = logging.getLogger(__name__)


# ======================================================================
# *** Helper functions
# ======================================================================
def _validate_values(
        values: Sequence[float | int] | NDArray[np.number]
) -> None:
    # check for negative values
    if any(value < 0 for value in values):
        logger.error("Negative values found")
        raise ValueError("Negative values found")

    total = sum(values)

    # check if total is zero
    if total == 0:
        logger.error("Total of values is zero")
        raise ValueError("Total of values is zero")


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

# ! ::: Fractions from sequence & numpy
def _calc_fractions(
        values: Sequence[float | int] | NDArray[np.number],
) -> NDArray[np.float64]:
    """
    Calculate the fractions of a list of values.

    Parameters
    ----------
    values : List[float | int] | NDArray[np.number]
        A list of numerical values or a numpy array of numerical values.

    Returns
    -------
    NDArray[np.float64]
        A numpy array of fractions corresponding to the input values.
    """
    # SECTION: validation
    _validate_values(values)

    # SECTION: to array
    values = np.array(values, dtype=np.float64)
    total = np.sum(values)

    # calculate fractions
    return cast(NDArray[np.float64], values / total)

# ! ::: Fractions from mapping


def _calc_fractions_from_mapping(
        values: Mapping[str, float | int],
) -> Dict[str, float]:
    """
    Calculate the fractions of a dictionary of values.

    Parameters
    ----------
    values : Mapping[str, float | int]
        A dictionary of component IDs and their corresponding values.

    Returns
    -------
    Optional[Dict[str, float]]
        A dictionary of component fractions, or None if the total is zero or negative values are found.
    """
    res = _calc_fractions(
        values=list(values.values())
    )

    return dict(zip(values.keys(), res.tolist()))


def _calc_component_fractions(
        values: Mapping[str, float | int],
        components: Optional[List[Component]] = None,
        component_key: Optional[ComponentKey] = None,
        case_sensitive: bool = True,
        sort_by_components_order: bool = True,
) -> Dict[str, float]:
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

    Returns
    -------
    Dict[str, float]
        A dictionary of component fractions, or None if the input is invalid.
    """
    # SECTION: get components values
    # ! configure component values if component_key is provided otherwise
    # ! otherwise use the original values dictionary
    if component_key is not None:
        if not components:
            logger.error("No components provided")
            raise ValueError("No components provided")

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
            logger.error("Failed to configure component values")
            raise ValueError("Failed to configure component values")

        # unpack
        component_values_dict, _ = component_values
    else:
        component_values_dict = values

    # SECTION: Calculate fractions
    total = sum(component_values_dict.values())
    if total == 0:
        logger.error("Total of component values is zero")
        raise ValueError("Total of component values is zero")

    # component fractions
    return _calc_fractions_from_mapping(
        values=component_values_dict
    )


# export
__all__ = [
    "_calc_fractions",
    "_calc_fractions_from_mapping",
    "_calc_component_fractions",
]
