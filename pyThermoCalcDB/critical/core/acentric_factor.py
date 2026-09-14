"""Core acentric-factor calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_finite_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to finite float64 array form."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _as_positive_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to finite positive float64 array form."""
    arr = _as_finite_array(values, name)
    # ! Vapor pressures and pressure ratios used in logarithms must be positive.
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return arr


# SECTION: Core numeric calculations
def _calc_acentric_factor_from_reduced_vapor_pressure(
    reduced_vapor_pressure: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Pitzer acentric factor from reduced vapor pressure.

    Equation: ``omega = -log10(Pr_sat) - 1`` at ``Tr = 0.7``. The reduced
    vapor pressure is dimensionless and positive. This definitional relation
    returns a dimensionless scalar or ``float64`` array.
    """
    pr_sat = _as_positive_array(reduced_vapor_pressure, "reduced_vapor_pressure")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], -np.log10(pr_sat) - 1.0))


def _calc_acentric_factor_from_vapor_pressure(
    saturation_pressure: NumericInput,
    critical_pressure: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Pitzer acentric factor from vapor and critical pressure.

    Equation: ``omega = -log10(Psat/Pc) - 1`` at ``Tr = 0.7``. Pressures must
    share a unit basis. The relation is definitional and preserves scalar/array
    behavior.
    """
    psat = _as_positive_array(saturation_pressure, "saturation_pressure")
    pc = _as_positive_array(critical_pressure, "critical_pressure")
    try:
        reduced = psat / pc
    except ValueError as exc:
        raise ValueError("saturation_pressure and critical_pressure must be broadcast-compatible.") from exc
    return _calc_acentric_factor_from_reduced_vapor_pressure(cast(NDArray[np.float64], reduced))


def _calc_reduced_vapor_pressure_from_acentric_factor(
    acentric_factor: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reduced vapor pressure from Pitzer acentric factor.

    Equation: ``Pr_sat = 10**(-omega - 1)`` at ``Tr = 0.7``. This exact inverse
    accepts finite dimensionless acentric factors and preserves scalar/array
    behavior.
    """
    omega = _as_finite_array(acentric_factor, "acentric_factor")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], np.power(10.0, -omega - 1.0)))


__all__ = [
    "_calc_acentric_factor_from_reduced_vapor_pressure",
    "_calc_acentric_factor_from_vapor_pressure",
    "_calc_reduced_vapor_pressure_from_acentric_factor",
]
