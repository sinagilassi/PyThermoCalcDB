"""Core viscosity correlations and mixing rules."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_state_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert scalar or array-like input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_positive(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive values."""
    # ! Viscosities and absolute temperatures must be positive.
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


def _validate_fractions(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative fractions that close along the last axis."""
    if np.any(values < 0.0):
        raise ValueError(f"{name} values must be non-negative.")
    if not np.allclose(np.sum(values, axis=-1), 1.0):
        raise ValueError(f"{name} must sum to 1.0 along the component axis.")


# SECTION: Core numeric calculations
def _calc_liquid_mixture_viscosity_log_rule(
    mole_fractions: NumericInput,
    component_viscosities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate liquid-mixture viscosity with the logarithmic mixing rule.

    Equation: ``eta = exp(sum_i x_i*ln(eta_i))``. Mole fractions are
    dimensionless and close along the last axis. Component viscosities must be
    positive and use one shared unit basis, commonly Pa*s. The output uses that
    same viscosity unit basis.
    """
    x = _as_state_array(mole_fractions, "mole_fractions")
    eta = _as_state_array(component_viscosities, "component_viscosities")
    _validate_fractions(x, "mole_fractions")
    _validate_positive(eta, "component_viscosities")
    if x.shape != eta.shape:
        raise ValueError("mole_fractions and component_viscosities must have the same shape.")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], np.exp(np.sum(x * np.log(eta), axis=-1))))


def _calc_viscosity_exponential_correlation(
    temperature: NumericInput,
    A: NumericInput,
    B: NumericInput,
    C: NumericInput,
    D: NumericInput,
    E: NumericInput,
) -> float | NDArray[np.float64]:
    """Evaluate the generic exponential viscosity correlation.

    Equation: ``eta = exp(A + B/T + C*ln(T) + D*T**E)``. Temperature is K.
    Coefficients are supplied externally and determine the output viscosity unit
    basis. This function implements only the equation form, not coefficient data.
    """
    t = _as_state_array(temperature, "temperature")
    a = _as_state_array(A, "A")
    b = _as_state_array(B, "B")
    c = _as_state_array(C, "C")
    d = _as_state_array(D, "D")
    e = _as_state_array(E, "E")
    _validate_positive(t, "temperature")
    try:
        viscosity = np.exp(a + b / t + c * np.log(t) + d * np.power(t, e))
    except ValueError as exc:
        raise ValueError("temperature and viscosity-correlation coefficients must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], viscosity))


__all__ = [
    "_calc_liquid_mixture_viscosity_log_rule",
    "_calc_viscosity_exponential_correlation",
]
