"""Core partial-molar property relations."""

# import libs
from collections.abc import Mapping
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...utils.conversions import (
    NumericArrayInput,
    _as_float_array,
    _return_scalar_if_zero_dim,
    _validate_fraction_array,
    _validate_same_array_shape,
    _validate_same_mapping_keys,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Core partial-molar calculations

def _calc_total_property_from_partial_molar_properties(
    amounts: NumericInput,
    partial_molar_properties: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate extensive property ``M = sum_i n_i*Mbar_i``."""
    n = _as_float_array(amounts, "amounts")
    partials = _as_float_array(
        partial_molar_properties,
        "partial_molar_properties",
    )
    _validate_same_array_shape(n, partials, "amounts", "partial_molar_properties")
    if np.any(n < 0.0):
        raise ValueError("amounts must be non-negative.")
    return _return_scalar_if_zero_dim(np.sum(n * partials, axis=-1))


def _calc_molar_property_from_partial_molar_properties(
    mole_fractions: NumericInput,
    partial_molar_properties: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate molar property ``M_m = sum_i x_i*Mbar_i``."""
    x = _as_float_array(mole_fractions, "mole_fractions")
    partials = _as_float_array(
        partial_molar_properties,
        "partial_molar_properties",
    )
    _validate_same_array_shape(x, partials, "mole_fractions", "partial_molar_properties")
    _validate_fraction_array(x, "mole_fractions")
    return _return_scalar_if_zero_dim(np.sum(x * partials, axis=-1))


def _calc_binary_partial_molar_properties(
    molar_property: NumericInput,
    mole_fraction_1: NumericInput,
    dmolar_property_dx1: NumericInput,
) -> tuple[float, float] | tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Calculate binary partial molar properties from ``M`` and ``dM/dx1``."""
    m = np.asarray(molar_property, dtype=np.float64)
    x1 = np.asarray(mole_fraction_1, dtype=np.float64)
    dmdx1 = np.asarray(dmolar_property_dx1, dtype=np.float64)
    for name, arr in (
        ("molar_property", m),
        ("mole_fraction_1", x1),
        ("dmolar_property_dx1", dmdx1),
    ):
        if arr.ndim > 2:
            raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} values must be finite.")
    # ! Binary composition coordinate must remain inside the closed interval.
    if np.any((x1 < 0.0) | (x1 > 1.0)):
        raise ValueError("mole_fraction_1 must be between 0 and 1.")
    mbar1 = np.asarray(m + (1.0 - x1) * dmdx1, dtype=np.float64)
    mbar2 = np.asarray(m - x1 * dmdx1, dtype=np.float64)
    # NOTE: Keep the two tuple entries in the same representation for typing.
    if mbar1.ndim == 0 and mbar2.ndim == 0:
        return float(mbar1), float(mbar2)
    return cast(NDArray[np.float64], mbar1), cast(NDArray[np.float64], mbar2)


# SECTION: Mapping adapters

def _calc_total_property_from_partial_molar_mapping(
    amounts: Mapping[str, float | int],
    partial_molar_properties: Mapping[str, float | int],
) -> float:
    """Calculate total property from aligned component mappings."""
    # ? Component identity is preserved by validating and iterating one key order.
    _validate_same_mapping_keys(
        amounts,
        partial_molar_properties,
        "amounts",
        "partial_molar_properties",
    )
    keys = list(amounts)
    return float(_calc_total_property_from_partial_molar_properties(
        [amounts[key] for key in keys],
        [partial_molar_properties[key] for key in keys],
    ))


def _calc_molar_property_from_partial_molar_mapping(
    mole_fractions: Mapping[str, float | int],
    partial_molar_properties: Mapping[str, float | int],
) -> float:
    """Calculate molar property from aligned component mappings."""
    _validate_same_mapping_keys(
        mole_fractions,
        partial_molar_properties,
        "mole_fractions",
        "partial_molar_properties",
    )
    keys = list(mole_fractions)
    return float(_calc_molar_property_from_partial_molar_properties(
        [mole_fractions[key] for key in keys],
        [partial_molar_properties[key] for key in keys],
    ))


# SECTION: Core exports
__all__ = [
    "_calc_total_property_from_partial_molar_properties",
    "_calc_molar_property_from_partial_molar_properties",
    "_calc_binary_partial_molar_properties",
    "_calc_total_property_from_partial_molar_mapping",
    "_calc_molar_property_from_partial_molar_mapping",
]
