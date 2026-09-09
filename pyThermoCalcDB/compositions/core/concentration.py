# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional, cast, Dict, List, Tuple
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp, UnitConversionFn
from pythermodb_settings.utils import config_components_values
# locals
from ...utils.conversions import _to_units, _to_volume, _to_amounts

# NOTE: logger setup
logger = logging.getLogger(__name__)

# ======================================================================
# *** Helper functions
# ======================================================================


def _validate(
    amount: NDArray[np.float64],
    volume: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Validate the shapes and values of component amount and solution volume arrays.

    Parameters
    ----------
    amount : NDArray[np.number]
        Array of component amounts.
    volume : NDArray[np.number]
        Array of solution volumes.

    Returns
    -------
    NDArray[np.number]
        Validated volume array. A 1-D per-state volume for 2-D amounts is returned
        as ``(n_states, 1)`` so NumPy broadcasts across components.

    Raises
    ------
    ValueError
        If the shapes of amounts and volumes are incompatible or if any volume is
        not finite and greater than zero.
    """
    if amount.ndim not in (1, 2):
        raise ValueError("component_moles must be a 1-D or 2-D array.")

    if volume.ndim > 2:
        raise ValueError(
            "solvent_mass must be a scalar, 1-D array, or 2-D array.")

    if amount.ndim == 1:
        if volume.ndim not in (0, 1):
            raise ValueError(
                "For 1-D component_moles, solvent_mass must be a scalar "
                "or have the same shape as component_moles."
            )
        if volume.ndim == 1 and volume.shape != amount.shape:
            raise ValueError(
                "For 1-D component_moles, solvent_mass must be a scalar "
                "or have the same shape as component_moles."
            )

    elif amount.ndim == 2:
        if volume.ndim == 0:
            pass
        elif volume.ndim == 1:
            if volume.shape[0] != amount.shape[0]:
                raise ValueError(
                    "For 2-D component_moles, 1-D solvent_mass must have "
                    "one entry per state."
                )
            volume = volume[:, None]

        elif volume.ndim == 2:
            if volume.shape not in (amount.shape, (amount.shape[0], 1)):
                raise ValueError(
                    "For 2-D component_moles, solvent_mass must be a scalar, "
                    "have shape (n_states,), shape (n_states, 1), or the "
                    "same shape as component_moles."
                )

    if not np.all(np.isfinite(volume)):
        raise ValueError("solvent_mass must contain finite values.")

    if np.any(volume <= 0):
        raise ValueError("solvent_mass must be greater than zero.")

    return volume

# ======================================================================
# *** Internal deterministic calculations
# ======================================================================
# ! ::: numpy


def _calc_concentrations(
    component_amounts: Sequence[float | int] | NDArray[np.number],
    solution_volume: float | int | NDArray[np.number],
):
    """
    Calculate concentrations using NumPy vectorization.

    Parameters
    ----------
    component_amounts : Sequence[float | int] | NDArray[np.number]
        Component amounts. May be a Python sequence or a 1-D/2-D NumPy array.
    solution_volume : float | int | NDArray[np.number]
        Solution volume. Must be a scalar or have a shape compatible with
        ``component_amounts``.

    Returns
    -------
    NDArray[np.float64]
        Concentrations with the same broadcasted shape as ``component_amounts``.
    """
    # set
    moles: NDArray[np.float64] = np.asarray(
        component_amounts,
        dtype=np.float64
    )
    volume: NDArray[np.float64] = np.asarray(
        solution_volume,
        dtype=np.float64
    )

    # validate
    volume = _validate(moles, volume)

    return cast(NDArray[np.float64], moles / volume)

# ! ::: Concentration [amount/volume] from sequence


def _calc_concentrations_from_sequence(
        component_amounts: Sequence[float | int],
        solution_volume: float | int,
) -> List[float]:
    """
    Calculate the concentration of each component in a solution.

    Parameters
    ----------
    component_amounts : Sequence[float | int]
        A sequence of amounts for each component.
    solution_volume : float | int
        The volume of the solution.

    Returns
    -------
    List[float]
        A list of concentration values for each component.
    """
    return _calc_concentrations(
        component_amounts,
        solution_volume
    ).tolist()


# ! ::: Concentration [amount/volume] from mapping

def _calc_concentrations_from_mapping(
    component_amounts: Mapping[str, float | int],
    solution_volume: float | int,
) -> Dict[str, float]:
    """
    Calculate the concentration of each keyed component in a solution.

    Parameters
    ----------
    component_amounts : Mapping[str, float | int]
        A mapping of component names to their respective amounts.
    solution_volume : float | int
        The volume of the solution.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping component names to their concentration values.
    """
    # calc
    concentrations_ = _calc_concentrations(
        component_amounts=list(component_amounts.values()),
        solution_volume=solution_volume
    )
    return dict(zip(component_amounts.keys(), concentrations_))


# ! ::: Concentration [amount/volume] with solution volume as CustomProp
def _calc_concentrations_from_props(
    component_amounts: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str,
) -> Dict[str, float]:
    """
    Calculate keyed concentrations from unit-aware component amounts and solution volume.

    Parameters
    ----------
    component_amounts : Mapping[str, CustomProp]
        A mapping of component names to their unit-aware amounts.
    solution_volume : CustomProp
        The unit-aware volume of the solution.
    output_unit : str
        The unit for the output concentration values.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping component names to concentration values in
        ``output_unit``.
    """
    # SECTION: set default units for amount and volume
    units_ = _to_units(output_unit)
    # >> set
    amount_unit = units_[0]
    volume_unit = units_[1]

    # NOTE: component amounts
    component_amounts_dict: Dict[str, float] = _to_amounts(
        component_amounts=component_amounts,
        output_unit=amount_unit,
    )

    # NOTE: solution volume unit should match output unit denominator
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit,
    )

    # SECTION: calculate concentration for each component
    return _calc_concentrations_from_mapping(
        component_amounts=component_amounts_dict,
        solution_volume=solution_volume_scalar,
    )


# ! ::: Concentration [amount/volume] with component ID mapping and sorting
def _calc_component_concentrations_from_props(
    component_amounts: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Dict[str, float]:
    """
    Calculate component concentrations with optional component-key remapping and ordering.

    Parameters
    ----------
    component_amounts : Mapping[str, CustomProp]
        A mapping of component identifiers to unit-aware amounts.
    solution_volume : CustomProp
        The unit-aware volume of the solution.
    output_unit : str
        The unit for the output concentration values.
    components : Optional[List[Component]], optional
        Component definitions used to resolve and order component identifiers,
        by default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component amounts to components, by default None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default True.
    sort_by_components_order : bool, optional
        Whether to sort the component concentrations by the order of components, by default True.
    unit_conversion_fn : UnitConversionFn, optional
        Reserved for API consistency. Unit conversion is handled by the
        conversion helpers used in this module.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping resolved component identifiers to concentration
        values in ``output_unit``.
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
        output_unit=amount_unit,
    )

    # ! volume
    solution_volume_scalar: float = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit,
    )

    # SECTION: get component values
    if component_key is not None:
        if not components:
            logger.error(
                "Component key is provided but components list is empty.")
            raise ValueError(
                "Component key is provided but components list is empty.")

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
            logger.error("Failed to configure component values")
            raise ValueError("Failed to configure component values")

        # unpack
        component_values_dict, _ = component_values
    else:
        component_values_dict = component_amounts_dict

    # SECTION: calculate concentration for each component
    return _calc_concentrations_from_mapping(
        component_amounts=component_values_dict,
        solution_volume=solution_volume_scalar,
    )
