# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional, cast, overload, Literal
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp
from pythermodb_settings.utils import (
    config_components_values,
)
# locals
from ...utils.conversions import _to_moles, _to_units, _to_volume

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Helper functions
# ======================================================================
def _validate_moles_and_volume(
    moles: NDArray[np.float64],
    volume: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Validate the shapes and values of moles and volume arrays.

    Parameters
    ----------
    moles : NDArray[np.number]
        Array of component moles.
    volume : NDArray[np.number]
        Array of solution volumes.

    Returns
    -------
    NDArray[np.number]
        Validated volume array. A 1-D per-state volume for 2-D moles is
        returned as ``(n_states, 1)`` so NumPy broadcasts across components.

    Raises
    ------
    ValueError
        If the shapes of moles and volume are incompatible or if any volume is zero.
    """
    if moles.ndim not in (1, 2):
        raise ValueError("component_moles must be a 1-D or 2-D array.")

    if volume.ndim > 2:
        raise ValueError(
            "solution_volume must be a scalar, 1-D array, or 2-D array."
        )

    if moles.ndim == 1:
        if volume.ndim not in (0, 1):
            raise ValueError(
                "For 1-D component_moles, solution_volume must be a scalar "
                "or have the same shape as component_moles."
            )
        if volume.ndim == 1 and volume.shape != moles.shape:
            raise ValueError(
                "For 1-D component_moles, solution_volume must be a scalar "
                "or have the same shape as component_moles."
            )

    elif moles.ndim == 2:
        if volume.ndim == 0:
            pass
        elif volume.ndim == 1:
            if volume.shape[0] != moles.shape[0]:
                raise ValueError(
                    "For 2-D component_moles, 1-D solution_volume must have "
                    "one entry per state."
                )
            volume = volume[:, None]

        elif volume.ndim == 2:
            if volume.shape not in (moles.shape, (moles.shape[0], 1)):
                raise ValueError(
                    "For 2-D component_moles, solution_volume must be a "
                    "scalar, have shape (n_states,), shape (n_states, 1), "
                    "or the same shape as component_moles."
                )

    # NOTE: volume must be finite and greater than zero
    # REVIEW
    if not np.all(np.isfinite(volume)):
        raise ValueError("solution_volume must contain finite values.")

    if np.any(volume <= 0):
        raise ValueError("solution_volume must be greater than zero.")

    return volume

# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

# ! ::: Molarity from array-like inputs


@overload
def _calc_molarities(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solution_volume: float | int | NDArray[np.number],
    *,
    as_list: Literal[False] = False,
) -> NDArray[np.float64]:
    ...


@overload
def _calc_molarities(
    component_moles: Sequence[float | int],
    solution_volume: float | int,
    *,
    as_list: Literal[True],
) -> list[float]:
    ...


def _calc_molarities(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solution_volume: float | int | NDArray[np.number],
    *,
    as_list: bool = False,
) -> NDArray[np.float64] | list[float]:
    """
    Calculate molarities using NumPy vectorization.

    Parameters
    ----------
    component_moles : Sequence[float | int] | NDArray[np.number]
        Component mole amounts. May be a Python sequence or a
        1-D/2-D NumPy array.
    solution_volume : float | int | NDArray[np.number]
        Solution volume. Must be a scalar or have the same shape as
        ``component_moles``.
    as_list : bool, optional
        If True, the result will be returned as a list of floats. Default is False, which returns a NumPy array.

    Returns
    -------
    NDArray[np.floating]
        Molarities with the same shape as ``component_moles``. The default return type is a NumPy array, but if ``as_list`` is True, a list of floats will be returned instead.
    """
    # set
    moles: NDArray[np.float64] = np.asarray(component_moles, dtype=np.float64)
    volume: NDArray[np.float64] = np.asarray(solution_volume, dtype=np.float64)

    # validate
    volume = _validate_moles_and_volume(moles, volume)
    # molarities
    molarities = cast(NDArray[np.float64], moles / volume)

    # check as list
    if as_list:
        return cast(list[float], molarities.tolist())
    return molarities


# ! ::: Molarity from mapping


def _calc_molarities_from_mapping(
    component_moles: Mapping[str, float | int],
    solution_volume: float,
) -> dict[str, float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume.

    Parameters
    ----------
    component_moles : Mapping[str, float | int]
        A mapping of component names to their respective moles.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molarity values.
    """
    # calc
    molarities_ = _calc_molarities(
        component_moles=list(component_moles.values()),
        solution_volume=solution_volume
    )

    # to dict
    component_molarity_dict = dict(
        zip(component_moles.keys(), molarities_.tolist())
    )

    return component_molarity_dict

# ! ::: Molarity with solution volume as CustomProp


def _calc_molarities_from_props(
    component_moles: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume as a CustomProp. The default
    volume unit is litre (L). The component molarity list and dictionary will be ordered according to the components list if sort_by_components_order is True.

    Parameters
    ----------
    component_moles : Mapping[str, CustomProp]
        A mapping of component names to CustomProp objects representing the moles.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output molarity values. Defaults to 'mol/L'.
    components : Optional[Sequence[Component]], optional
        A sequence of Component objects to map the component moles to, by default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component moles to components, by default None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default True.
    sort_by_components_order : bool, optional
        Whether to sort the component molarities by the order of components, by default True.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molarity values, or None if the calculation could not be performed.
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
            logger.error(
                "Component key is provided but components list is empty."
            )
            raise ValueError(
                "Component key is provided but components list is empty."
            )

        component_values: tuple[
            dict[str, Any],
            list[Any]
        ] | None = config_components_values(
            values=component_moles_dict,
            components=list(components),
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order
        )
        # >> check
        if component_values is None:
            logger.error(
                "Failed to configure component values."
            )
            raise ValueError(
                "Failed to configure component values."
            )

        # unpack
        component_values_dict, _ = component_values
    else:
        component_values_dict = component_moles_dict

    # SECTION: calculate molarity for each component
    return _calc_molarities_from_mapping(
        component_moles=component_values_dict,
        solution_volume=solution_volume_scalar,
    )


# all
__all__ = [
    "_calc_molarities",
    "_calc_molarities_from_mapping",
    "_calc_molarities_from_props",
]
