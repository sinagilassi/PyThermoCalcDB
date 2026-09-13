"""Core specific-volume and density conversion calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, ScalarValue, UnitConversionFn

# locals
from pythermocalcdb.utils.conversions import (
    NumericArrayInput,
    _pos,
    _return_scalar_if_zero_dim,
    _validate_custom_prop_scalar,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators

def _as_positive_array(
    values: NumericInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert reciprocal conversion inputs to finite positive float64 arrays."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, 1-D, or 2-D values.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return cast(NDArray[np.float64], arr)


# SECTION: Core numeric calculations

def _calc_density_to_specific_volume(
    density: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reciprocal specific volume ``v = 1/rho``."""
    # SECTION: Normalize and validate
    rho = _as_positive_array(density, "density")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], 1.0 / rho))


def _calc_specific_volume_to_density(
    specific_volume: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reciprocal density ``rho = 1/v``."""
    # SECTION: Normalize and validate
    v = _as_positive_array(specific_volume, "specific_volume")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], 1.0 / v))


# SECTION: Props adapters

def _calc_density_to_specific_volume_from_props(
    density: CustomProp,
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate specific volume from a unit-aware density scalar."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(density, "density")

    # SECTION: Normalize and calculate
    rho = _pos(density, "density", output_density_unit, unit_conversion_fn)
    return float(_calc_density_to_specific_volume(rho))


def _calc_specific_volume_to_density_from_props(
    specific_volume: CustomProp,
    output_specific_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate density from a unit-aware specific-volume scalar."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(specific_volume, "specific_volume")

    # SECTION: Normalize and calculate
    v = _pos(
        specific_volume,
        "specific_volume",
        output_specific_volume_unit,
        unit_conversion_fn,
    )
    return float(_calc_specific_volume_to_density(v))


# SECTION: Scalar adapters

def _calc_density_to_specific_volume_from_scalars(
    density: ScalarValue,
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize a scalar public density input and calculate specific volume."""
    rho = _pos(density, "density", output_density_unit, unit_conversion_fn)
    return float(_calc_density_to_specific_volume(rho))


def _calc_specific_volume_to_density_from_scalars(
    specific_volume: ScalarValue,
    output_specific_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize a scalar public specific-volume input and calculate density."""
    v = _pos(
        specific_volume,
        "specific_volume",
        output_specific_volume_unit,
        unit_conversion_fn,
    )
    return float(_calc_specific_volume_to_density(v))


# SECTION: Core exports
__all__ = [
    "_calc_density_to_specific_volume",
    "_calc_specific_volume_to_density",
    "_calc_density_to_specific_volume_from_props",
    "_calc_specific_volume_to_density_from_props",
    "_calc_density_to_specific_volume_from_scalars",
    "_calc_specific_volume_to_density_from_scalars",
]
