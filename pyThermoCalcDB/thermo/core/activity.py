"""Core activity transformations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim


# SECTION: Numeric helpers

def _as_state_array(values: NumericArrayInput, name: str) -> NDArray[np.float64]:
    """Convert scalar or array-like input to a finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_non_negative(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative values."""
    # ! Fractions and concentrations cannot be negative.
    if np.any(values < 0.0):
        raise ValueError(f"{name} values must be non-negative.")


def _validate_positive(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive values."""
    # ! Coefficients and reference concentrations are denominators or scales.
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


# SECTION: Core activity calculations

def _calc_activity_from_mole_fraction(
    mole_fraction: NumericArrayInput,
    activity_coefficient: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate activity from mole fraction: ``a_i = gamma_i*x_i``."""
    x = _as_state_array(mole_fraction, "mole_fraction")
    gamma = _as_state_array(activity_coefficient, "activity_coefficient")
    _validate_non_negative(x, "mole_fraction")
    _validate_positive(gamma, "activity_coefficient")
    try:
        return _return_scalar_if_zero_dim(x * gamma)
    except ValueError as exc:
        raise ValueError("mole_fraction and activity_coefficient must be broadcast-compatible.") from exc


def _calc_activity_from_concentration(
    concentration: NumericArrayInput,
    activity_coefficient: NumericArrayInput,
    reference_concentration: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate concentration-basis activity: ``a_i = gamma_i*c_i/c_ref``."""
    c = _as_state_array(concentration, "concentration")
    gamma = _as_state_array(activity_coefficient, "activity_coefficient")
    c_ref = _as_state_array(reference_concentration, "reference_concentration")
    _validate_non_negative(c, "concentration")
    _validate_positive(gamma, "activity_coefficient")
    _validate_positive(c_ref, "reference_concentration")
    try:
        return _return_scalar_if_zero_dim(gamma * c / c_ref)
    except ValueError as exc:
        raise ValueError("concentration, activity_coefficient, and reference_concentration must be broadcast-compatible.") from exc


def _calc_effective_concentration(
    activity: NumericArrayInput,
    reference_concentration: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate effective concentration from activity: ``c_eff = a_i*c_ref``."""
    a = _as_state_array(activity, "activity")
    c_ref = _as_state_array(reference_concentration, "reference_concentration")
    _validate_non_negative(a, "activity")
    _validate_positive(c_ref, "reference_concentration")
    try:
        return _return_scalar_if_zero_dim(a * c_ref)
    except ValueError as exc:
        raise ValueError("activity and reference_concentration must be broadcast-compatible.") from exc



def _calc_activity_coefficient_from_fugacity(
    liquid_fugacity: NumericArrayInput,
    mole_fraction: NumericArrayInput,
    standard_state_fugacity: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate activity coefficient from fugacity definition.

    Equation: ``gamma_i = f_i^L/(x_i*f_i^0)``. Fugacities must share a pressure
    unit basis, mole fraction is dimensionless, and all denominator terms must
    be positive.
    """
    f_l = _as_state_array(liquid_fugacity, "liquid_fugacity")
    x = _as_state_array(mole_fraction, "mole_fraction")
    f0 = _as_state_array(standard_state_fugacity, "standard_state_fugacity")
    _validate_positive(f_l, "liquid_fugacity")
    _validate_positive(x, "mole_fraction")
    _validate_positive(f0, "standard_state_fugacity")
    try:
        return _return_scalar_if_zero_dim(f_l / (x * f0))
    except ValueError as exc:
        raise ValueError("liquid_fugacity, mole_fraction, and standard_state_fugacity must be broadcast-compatible.") from exc
# SECTION: Core exports
__all__ = [
    "_calc_activity_from_mole_fraction",
    "_calc_activity_from_concentration",
    "_calc_effective_concentration",
    "_calc_activity_coefficient_from_fugacity",
]


