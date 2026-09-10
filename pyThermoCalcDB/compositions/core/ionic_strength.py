"""Shared numeric ionic-strength core."""

# import libs
from collections.abc import Mapping, Sequence
from typing import cast, Optional, Literal
from pythermodb_settings.models import Component, ComponentKey, CustomProp, UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.components import extract_components_values
from pythermodb_settings.utils.quantity import to_dict
from pythermodb_settings.utils.validators import non_negative, same_shape

import numpy as np
from numpy.typing import NDArray
from ...utils.conversions import (
    _resolve_unit_conversion_fn,
    _validate_same_keys,
    _configure_component_values
)


# =====================================================================
# *** Helper functions
# =====================================================================

def _validate_concentrations_and_charges(
    concentrations: NDArray[np.float64],
    charges: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Validate and broadcast concentration and charge arrays."""
    # SECTION: Validate array dimensions
    if concentrations.ndim > 2:
        raise ValueError(
            "concentrations must be a scalar, 1-D array, or 2-D array."
        )

    if charges.ndim > 2:
        raise ValueError("charges must be a scalar, 1-D array, or 2-D array.")

    # SECTION: Validate finite numeric values
    if not np.all(np.isfinite(concentrations)):
        raise ValueError("concentrations must contain finite values.")

    if not np.all(np.isfinite(charges)):
        raise ValueError("charges must contain finite values.")

    # ! Concentrations cannot be negative for ionic strength.
    if np.any(concentrations < 0):
        raise ValueError("concentrations must be non-negative.")

    # ? NumPy handles scalar, 1-D, and 2-D broadcasting consistently here.
    try:
        concentrations, charges = np.broadcast_arrays(concentrations, charges)
    except ValueError as exc:
        raise ValueError(
            "concentrations and charges must have broadcast-compatible shapes."
        ) from exc

    return (
        cast(NDArray[np.float64], concentrations),
        cast(NDArray[np.float64], charges),
    )


# ======================================================================
# *** Core calculation
# ======================================================================

def _calc_ionic_strength_core(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate ionic strength: ``I = 0.5 * sum_i(c_i * z_i**2)``."""
    # SECTION: Convert inputs to numeric arrays
    concentration_values: NDArray[np.float64] = np.asarray(
        concentrations,
        dtype=np.float64,
    )
    charge_values: NDArray[np.float64] = np.asarray(
        charges,
        dtype=np.float64,
    )

    # SECTION: Validate inputs
    concentration_values, charge_values = _validate_concentrations_and_charges(
        concentration_values,
        charge_values,
    )

    # SECTION: Calculate component ionic-strength contributions
    values = 0.5 * concentration_values * charge_values**2

    # NOTE: Scalar and 1-D inputs return a scalar; 2-D inputs return row sums.
    if values.ndim == 0:
        return float(values)

    if values.ndim == 1:
        return float(np.sum(values))

    return cast(NDArray[np.float64], np.sum(values, axis=1))


# ======================================================================
# *** Internal deterministic calculations (molality)
# ======================================================================
# ! ::: Numeric array-like core

def _calc_ionic_strength(
    values: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate molality-based ionic strength from numeric array-like inputs."""
    # SECTION: Calculate ionic strength
    return _calc_ionic_strength_core(
        concentrations=values,
        charges=charges,
    )


# ! ::: Mapping adapter

def _calc_ionic_strength_from_mapping(
    values: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    mode: Literal['molarity', 'molality']
) -> float:
    """Calculate molarity- or molality-based ionic strength from numeric mapping inputs."""
    # SECTION: Validate inputs
    non_negative(values, mode)
    same_shape(values, charges)
    _validate_same_keys(values, charges)

    # SECTION: Normalize molalities and charges for consistent ordering
    normalized_molalities = list(dict(sorted(values.items())).values())
    normalized_charges = list(dict(sorted(charges.items())).values())

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength_core(
            normalized_molalities,
            normalized_charges,
        ),
    )

# ! ::: Unit-aware mapping adapter


def _calc_ionic_strength_from_props(
    values: Mapping[str, CustomProp],
    charges: Mapping[str, CustomProp],
    mode: Literal["molarity", "molality"],
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ionic strength from unit-aware molarity inputs."""
    # SECTION: Normalize unit-aware molarities
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    normalized_values = to_dict(
        values,
        output_unit,
        unit_conversion_fn=conversion_fn,
    )
    normalized_charges = to_dict(
        charges,
        None,
        None,
    )

    # ? Component metadata can remap keys and enforce component order.
    normalized_values = _configure_component_values(
        normalized_values,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        mode,
    )
    normalized_charges = _configure_component_values(
        normalized_charges,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "charges",
    )

    # ! Mapping inputs must pair the same species.
    _validate_same_keys(normalized_values, normalized_charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength_core(
            list(normalized_values.values()),
            [normalized_charges[key] for key in normalized_values],
        ),
    )


# SECTION: Public exports
__all__ = [
    "_calc_ionic_strength",
    "_calc_ionic_strength_from_mapping",
    "_calc_ionic_strength_from_props",
]
