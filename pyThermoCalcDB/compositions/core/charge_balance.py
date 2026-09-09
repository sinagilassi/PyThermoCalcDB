"""Charge-balance calculations."""

# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Optional, cast

import numpy as np
from numpy.typing import NDArray

# >> pythermodb-settings
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict
from pythermodb_settings.utils.validators import non_negative, same_shape

# locals
from ...utils.conversions import (
    _configure_component_values,

    _resolve_unit_conversion_fn,
    _validate_same_keys,
)
from ...utils.tools import to_annotated_value

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Helper functions
# ======================================================================
def _validate_charge_balance_inputs(
    concentrations: NDArray[np.float64],
    charges: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """
    Validate and broadcast charge-balance input arrays.

    Parameters
    ----------
    concentrations : NDArray[np.float64]
        Numeric concentration values. Scalars, 1-D arrays, and 2-D arrays are
        supported.
    charges : NDArray[np.float64]
        Numeric charge values aligned with ``concentrations``.

    Returns
    -------
    tuple[NDArray[np.float64], NDArray[np.float64]]
        Broadcast-compatible concentration and charge arrays.
    """
    # SECTION: Validate array dimensions
    if concentrations.ndim > 2:
        raise ValueError(
            "concentrations must be a scalar, 1-D array, or 2-D array."
        )

    if charges.ndim > 2:
        raise ValueError("charges must be a scalar, 1-D array, or 2-D array.")

    if not np.all(np.isfinite(concentrations)):
        raise ValueError("concentrations must contain finite values.")

    if not np.all(np.isfinite(charges)):
        raise ValueError("charges must contain finite values.")

    # ! Concentrations cannot be negative for charge-balance composition input.
    if np.any(concentrations < 0):
        raise ValueError("concentrations must be non-negative.")

    # NOTE: Broadcast inputs to compatible shapes before calculation.
    try:
        concentrations, charges = np.broadcast_arrays(
            concentrations,
            charges,
        )
    except ValueError as exc:
        raise ValueError(
            "concentrations and charges must have broadcast-compatible shapes."
        ) from exc

    return (
        cast(NDArray[np.float64], concentrations),
        cast(NDArray[np.float64], charges),
    )


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================
# ! ::: Numeric array-like core


def _calc_charge_balance(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """
    Calculate charge-balance residual from numeric array-like inputs.

    Parameters
    ----------
    concentrations : float | int | Sequence[float | int] | NDArray[np.number]
        Species concentrations on a common concentration basis.
    charges : float | int | Sequence[float | int] | NDArray[np.number]
        Species charges aligned with ``concentrations``.

    Returns
    -------
    float | NDArray[np.float64]
        Charge-balance residual ``sum_i(z_i*c_i)``. Scalar and 1-D inputs
        return a scalar; 2-D inputs return one residual per row.
    """
    # SECTION: Convert inputs to numeric arrays
    concentration_values: NDArray[np.float64] = np.asarray(
        concentrations,
        dtype=np.float64,
    )
    charge_values: NDArray[np.float64] = np.asarray(
        charges,
        dtype=np.float64,
    )

    # NOTE: validation
    concentration_values, charge_values = _validate_charge_balance_inputs(
        concentration_values,
        charge_values,
    )

    # SECTION: Calculate residual contributions
    values = concentration_values * charge_values

    # NOTE: Scalar and 1-D inputs return a scalar; 2-D inputs return row sums.
    if values.ndim == 0:
        return float(values)

    if values.ndim == 1:
        return float(np.sum(values))

    return cast(NDArray[np.float64], np.sum(values, axis=1))


# ! ::: Sequence adapter

def _calc_charge_balance_from_sequence(
    concentrations: Sequence[float | int],
    charges: Sequence[float | int],
) -> float:
    """
    Calculate charge-balance residual from numeric sequence inputs.

    Parameters
    ----------
    concentrations : Sequence[float | int]
        Species concentrations in sequence order.
    charges : Sequence[float | int]
        Species charges aligned with ``concentrations``.

    Returns
    -------
    float
        Charge-balance residual for the sequence.
    """
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)

    # SECTION: Calculate charge-balance residual
    return cast(
        float,
        _calc_charge_balance(
            concentrations,
            charges,
        ),
    )


# ! ::: Mapping adapter

def _calc_charge_balance_from_mapping(
    concentrations: Mapping[str, float | int],
    charges: Mapping[str, float | int],
) -> float:
    """
    Calculate charge-balance residual from numeric mapping inputs.

    Parameters
    ----------
    concentrations : Mapping[str, float | int]
        Species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Species charges keyed by the same identifiers as ``concentrations``.

    Returns
    -------
    float
        Charge-balance residual for the keyed species set.
    """
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)
    _validate_same_keys(concentrations, charges)

    # SECTION: Calculate charge-balance residual
    return cast(
        float,
        _calc_charge_balance(
            list(concentrations.values()),
            [charges[key] for key in concentrations],
        ),
    )


# ! ::: Unit-aware mapping adapter

def _calc_charge_balance_from_props(
    concentrations: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_concentration_unit: str,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """
    Calculate charge-balance residual from unit-aware concentration inputs.

    Parameters
    ----------
    concentrations : Mapping[str, CustomProp]
        Unit-aware species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Numeric species charges keyed by component identifier.
    output_concentration_unit : str
        Unit used to normalize concentration values before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to the configured conversion helper.
    components : Sequence[Component], optional
        Component metadata used for key remapping and ordering.
    component_key : ComponentKey, optional
        Component identifier format used for mapping keys.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether values should follow the order of ``components``.

    Returns
    -------
    float
        Charge-balance residual on the normalized concentration basis.
    """
    # SECTION: Normalize unit-aware concentrations
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    normalized_concentrations = to_dict(
        concentrations,
        output_concentration_unit,
        unit_conversion_fn=conversion_fn,
    )

    # ? Component metadata can remap keys and enforce component order.
    normalized_concentrations = _configure_component_values(
        normalized_concentrations,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "concentrations",
    )
    normalized_charges = _configure_component_values(
        dict(charges),
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "charges",
    )

    # ! Mapping inputs must pair the same species.
    _validate_same_keys(normalized_concentrations, normalized_charges)

    # SECTION: Calculate charge-balance residual
    return cast(
        float,
        _calc_charge_balance(
            list(normalized_concentrations.values()),
            [normalized_charges[key] for key in normalized_concentrations],
        ),
    )

# ! ::: Numeric array-like check


def _check_electroneutrality(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
    tolerance: float = 1e-12,
) -> bool | NDArray[np.bool_]:
    """
    Check whether numeric charge-balance residuals are within tolerance.

    Parameters
    ----------
    concentrations : float | int | Sequence[float | int] | NDArray[np.number]
        Species concentrations on a common concentration basis.
    charges : float | int | Sequence[float | int] | NDArray[np.number]
        Species charges aligned with ``concentrations``.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.

    Returns
    -------
    bool | NDArray[np.bool_]
        ``True`` where ``abs(sum_i(z_i*c_i)) <= tolerance``.
    """
    # ! Tolerance must be a valid non-negative residual magnitude.
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate charge-balance residual
    result = _calc_charge_balance(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Check electroneutrality
    is_neutral = np.abs(result) <= tolerance

    if isinstance(is_neutral, np.ndarray):
        return cast(NDArray[np.bool_], is_neutral)

    return bool(is_neutral)


# ! ::: Sequence check adapter

def _check_electroneutrality_from_sequence(
    concentrations: Sequence[float | int],
    charges: Sequence[float | int],
    tolerance: float = 1e-12,
) -> bool:
    """
    Check electroneutrality from numeric sequence inputs.

    Parameters
    ----------
    concentrations : Sequence[float | int]
        Species concentrations in sequence order.
    charges : Sequence[float | int]
        Species charges aligned with ``concentrations``.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.

    Returns
    -------
    bool
        ``True`` when the sequence residual is within tolerance.
    """
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)

    # SECTION: Check electroneutrality
    return cast(
        bool,
        _check_electroneutrality(
            concentrations=concentrations,
            charges=charges,
            tolerance=tolerance,
        ),
    )


# ! ::: Mapping check adapter

def _check_electroneutrality_from_mapping(
    concentrations: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    tolerance: float = 1e-12,
) -> bool:
    """
    Check electroneutrality from numeric mapping inputs.

    Parameters
    ----------
    concentrations : Mapping[str, float | int]
        Species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Species charges keyed by the same identifiers as ``concentrations``.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.

    Returns
    -------
    bool
        ``True`` when the keyed residual is within tolerance.
    """
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)
    _validate_same_keys(concentrations, charges)

    # SECTION: Check electroneutrality
    return cast(
        bool,
        _check_electroneutrality(
            concentrations=list(concentrations.values()),
            charges=[charges[key] for key in concentrations],
            tolerance=tolerance,
        ),
    )


# ! ::: Unit-aware mapping check adapter

def _check_electroneutrality_from_props(
    concentrations: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_concentration_unit: str,
    tolerance: float = 1e-12,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> bool:
    """
    Check electroneutrality from unit-aware concentration inputs.

    Parameters
    ----------
    concentrations : Mapping[str, CustomProp]
        Unit-aware species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Numeric species charges keyed by component identifier.
    output_concentration_unit : str
        Unit used to normalize concentration values before checking.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to the configured conversion helper.
    components : Sequence[Component], optional
        Component metadata used for key remapping and ordering.
    component_key : ComponentKey, optional
        Component identifier format used for mapping keys.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether values should follow the order of ``components``.

    Returns
    -------
    bool
        ``True`` when the normalized residual is within tolerance.
    """
    # SECTION: Normalize unit-aware concentrations
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    normalized_concentrations = to_dict(
        concentrations,
        output_concentration_unit,
        unit_conversion_fn=conversion_fn,
    )

    # ? Component metadata can remap keys and enforce component order.
    normalized_concentrations = _configure_component_values(
        normalized_concentrations,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "concentrations",
    )
    normalized_charges = _configure_component_values(
        dict(charges),
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "charges",
    )

    # ! Mapping inputs must pair the same species.
    _validate_same_keys(normalized_concentrations, normalized_charges)

    # SECTION: Check electroneutrality
    return cast(
        bool,
        _check_electroneutrality(
            concentrations=list(normalized_concentrations.values()),
            charges=[
                normalized_charges[key]
                for key in normalized_concentrations
            ],
            tolerance=tolerance,
        ),
    )


# SECTION: Public exports
__all__ = [
    # internal
    "_calc_charge_balance",
    "_calc_charge_balance_from_sequence",
    "_calc_charge_balance_from_mapping",
    "_calc_charge_balance_from_props",
    "_check_electroneutrality",
    "_check_electroneutrality_from_sequence",
    "_check_electroneutrality_from_mapping",
    "_check_electroneutrality_from_props",
]
