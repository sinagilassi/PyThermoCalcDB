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
    Validate the shapes and values of moles and solvent mass arrays.

    Parameters
    ----------
    amount : NDArray[np.number]
        Array of component amount.
    volume : NDArray[np.number]
        Array of solution volume.

    Returns
    -------
    NDArray[np.number]
        Validated mass array. A 1-D per-state mass for 2-D moles is returned
        as ``(n_states, 1)`` so NumPy broadcasts across components.

    Raises
    ------
    ValueError
        If the shapes of moles and mass are incompatible or if any mass is
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
    Calculate the concentration of each component in a solution given component amounts and solution volume as a CustomProp.

    Parameters
    ----------
    component_amounts : ComponentAmounts
        A dictionary mapping component names to their respective amounts. Numeric values are assumed to already be in the numerator unit from output_unit.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output concentration values. Defaults to 'kg/m^3'.
    unit_conversion_fn : UnitConversionFn, optional
        The function to use for unit conversion. Defaults to None. Then it will use the default conversion function `pycuc.convert_from_to`.

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
    Calculate the concentration of each component in a solution given component amounts and solution volume as a CustomProp. The component concentration list and dictionary will be ordered according to the components list if sort_by_components_order is True.

    Parameters
    ----------
    component_amounts : Mapping[str, CustomProp]
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
    unit_conversion_fn : UnitConversionFn, optional
        The function to use for unit conversion. Defaults to None. Then it will use the default conversion function `pycuc.convert_from_to`.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping component names to their calculated concentration values.
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
