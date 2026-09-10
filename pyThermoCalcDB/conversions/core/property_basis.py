"""Core property basis conversion calculations."""

# SECTION: Imports
# import libs
from typing import cast
import numpy as np
from numpy.typing import NDArray
# locals
from .extensive_intensive import (
    NumericInput,
    _as_float_array,
    _return_scalar_if_zero_dim,
    _validate_broadcastable,
)


# SECTION: Molar and mass-specific basis conversions
# ! ::: Convert molar property to mass-specific property
def _calc_molar_to_mass_specific(
    molar_property: NumericInput,
    molecular_weight: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate a mass-specific property from a molar property.

    Equation
        Y_mass = Y_molar / M

    Scalars, 1-D arrays, and 2-D arrays are supported. Array inputs must be
    broadcast-compatible; 2-D arrays conventionally represent states along
    axis 0 and components/properties along axis 1.
    """
    y_molar = _as_float_array(molar_property, "molar_property")
    mw = _as_float_array(molecular_weight, "molecular_weight")
    # NOTE: Molecular weight is a material basis and denominator here.
    if np.any(mw <= 0):
        raise ValueError("molecular_weight must be greater than zero.")
    _validate_broadcastable(y_molar, mw, "molar_property", "molecular_weight")
    result = cast(NDArray[np.float64], y_molar / mw)
    return _return_scalar_if_zero_dim(result)


# ! ::: Convert mass-specific property to molar property
def _calc_mass_specific_to_molar(
    mass_specific_property: NumericInput,
    molecular_weight: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate a molar property from a mass-specific property.

    Equation
        Y_molar = Y_mass * M

    Scalars, 1-D arrays, and 2-D arrays are supported. Array inputs must be
    broadcast-compatible; 2-D arrays conventionally represent states along
    axis 0 and components/properties along axis 1.
    """
    y_mass = _as_float_array(mass_specific_property, "mass_specific_property")
    mw = _as_float_array(molecular_weight, "molecular_weight")
    # NOTE: Molecular weight must represent a positive mass per amount basis.
    if np.any(mw <= 0):
        raise ValueError("molecular_weight must be greater than zero.")
    _validate_broadcastable(y_mass, mw, "mass_specific_property", "molecular_weight")
    result = cast(NDArray[np.float64], y_mass * mw)
    return _return_scalar_if_zero_dim(result)


# SECTION: Public exports
__all__ = [
    "_calc_molar_to_mass_specific",
    "_calc_mass_specific_to_molar",
]
