"""Core extensive and intensive property conversion calculations."""

# SECTION: Imports
# import libs
from collections.abc import Sequence
from typing import cast
import numpy as np
from numpy.typing import NDArray


# SECTION: Type aliases
NumericInput = float | int | Sequence[float | int] | NDArray[np.number]


# SECTION: Validation helpers
# ! ::: Normalize numeric inputs
def _as_float_array(value: NumericInput, name: str) -> NDArray[np.float64]:
    """Normalize scalar or array-like numeric input."""
    value_array: NDArray[np.float64] = np.asarray(value, dtype=np.float64)
    # NOTE: Extensive/intensive conversions support scalar, vector, and matrix states.
    if value_array.ndim > 2:
        raise ValueError(f"{name} must be a scalar, 1-D array, or 2-D array.")
    if not np.all(np.isfinite(value_array)):
        raise ValueError(f"{name} must contain finite values.")
    return value_array


# ! ::: Validate NumPy broadcasting semantics
def _validate_broadcastable(
    left: NDArray[np.float64],
    right: NDArray[np.float64],
    left_name: str,
    right_name: str,
) -> None:
    """Validate that two numeric inputs can broadcast intentionally."""
    try:
        np.broadcast_shapes(left.shape, right.shape)
    except ValueError as exc:
        # ? Shape mismatches usually indicate inconsistent state/property inputs.
        raise ValueError(
            f"{left_name} and {right_name} must be broadcast-compatible."
        ) from exc


# ! ::: Preserve scalar return ergonomics
def _return_scalar_if_zero_dim(
    value: NDArray[np.float64],
) -> float | NDArray[np.float64]:
    """Return Python float for scalar calculations, ndarray otherwise."""
    if value.ndim == 0:
        return float(value)
    return value


# SECTION: Molar extensive/intensive conversions
# ! ::: Convert molar property to total property
def _calc_molar_property_to_total(
    moles: NumericInput,
    molar_property: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate a total extensive property from amount and molar property.

    Equation
        Y_total = n * Y_molar

    Scalars, 1-D arrays, and 2-D arrays are supported. Array inputs must be
    broadcast-compatible; 2-D arrays conventionally represent states along
    axis 0 and components/properties along axis 1.
    """
    n = _as_float_array(moles, "moles")
    y_molar = _as_float_array(molar_property, "molar_property")
    # NOTE: Amount is a divisor in the inverse conversion and must define material presence.
    if np.any(n <= 0):
        raise ValueError("moles must be greater than zero.")
    _validate_broadcastable(n, y_molar, "moles", "molar_property")
    result = cast(NDArray[np.float64], n * y_molar)
    return _return_scalar_if_zero_dim(result)


# ! ::: Convert total property to molar property
def _calc_total_to_molar_property(
    total_property: NumericInput,
    moles: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate a molar property from total property and amount.

    Equation
        Y_molar = Y_total / n

    Scalars, 1-D arrays, and 2-D arrays are supported. Array inputs must be
    broadcast-compatible; 2-D arrays conventionally represent states along
    axis 0 and components/properties along axis 1.
    """
    y_total = _as_float_array(total_property, "total_property")
    n = _as_float_array(moles, "moles")
    # NOTE: Positive moles prevent division by zero and negative material amounts.
    if np.any(n <= 0):
        raise ValueError("moles must be greater than zero.")
    _validate_broadcastable(y_total, n, "total_property", "moles")
    result = cast(NDArray[np.float64], y_total / n)
    return _return_scalar_if_zero_dim(result)


# SECTION: Mass extensive/intensive conversions
# ! ::: Convert mass-specific property to total property
def _calc_specific_property_to_total(
    mass: NumericInput,
    specific_property: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate a total extensive property from mass and specific property.

    Equation
        Y_total = m * Y_specific

    Scalars, 1-D arrays, and 2-D arrays are supported. Array inputs must be
    broadcast-compatible; 2-D arrays conventionally represent states along
    axis 0 and components/properties along axis 1.
    """
    mass_value = _as_float_array(mass, "mass")
    y_specific = _as_float_array(specific_property, "specific_property")
    # NOTE: Mass is treated as a strictly positive material basis.
    if np.any(mass_value <= 0):
        raise ValueError("mass must be greater than zero.")
    _validate_broadcastable(mass_value, y_specific, "mass", "specific_property")
    result = cast(NDArray[np.float64], mass_value * y_specific)
    return _return_scalar_if_zero_dim(result)


# ! ::: Convert total property to mass-specific property
def _calc_total_to_specific_property(
    total_property: NumericInput,
    mass: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate a mass-specific property from total property and mass.

    Equation
        Y_specific = Y_total / m

    Scalars, 1-D arrays, and 2-D arrays are supported. Array inputs must be
    broadcast-compatible; 2-D arrays conventionally represent states along
    axis 0 and components/properties along axis 1.
    """
    y_total = _as_float_array(total_property, "total_property")
    mass_value = _as_float_array(mass, "mass")
    # NOTE: Positive mass prevents division by zero and invalid material basis.
    if np.any(mass_value <= 0):
        raise ValueError("mass must be greater than zero.")
    _validate_broadcastable(y_total, mass_value, "total_property", "mass")
    result = cast(NDArray[np.float64], y_total / mass_value)
    return _return_scalar_if_zero_dim(result)


# SECTION: Public exports
__all__ = [
    "_calc_molar_property_to_total",
    "_calc_specific_property_to_total",
    "_calc_total_to_molar_property",
    "_calc_total_to_specific_property",
]
