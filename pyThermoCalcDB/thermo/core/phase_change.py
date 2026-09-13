"""Core phase-change thermodynamic identities."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import (
    NumericArrayInput,
    _return_scalar_if_zero_dim,
    _validate_positive_array,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_finite_float_array(
    values: NumericInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert numeric input to a finite scalar/array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(
            f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


# SECTION: Core numeric calculations
def _calc_phase_transition_entropy(
    transition_enthalpy: NumericInput,
    transition_temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``delta_S_tr = delta_H_tr / T_tr``."""
    h_tr = _as_finite_float_array(
        transition_enthalpy, "transition_enthalpy")
    t_tr = _as_finite_float_array(
        transition_temperature, "transition_temperature")
    _validate_positive_array(t_tr, "transition_temperature")
    return _return_scalar_if_zero_dim(h_tr / t_tr)


def _calc_enthalpy_of_sublimation(
    enthalpy_fusion: NumericInput,
    enthalpy_vaporization: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``delta_H_sub = delta_H_fus + delta_H_vap``."""
    h_fus = _as_finite_float_array(enthalpy_fusion, "enthalpy_fusion")
    h_vap = _as_finite_float_array(
        enthalpy_vaporization, "enthalpy_vaporization")
    return _return_scalar_if_zero_dim(h_fus + h_vap)


def _calc_clapeyron_slope(
    transition_enthalpy: NumericInput,
    temperature: NumericInput,
    delta_molar_volume: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``dP/dT = delta_H_tr / (T * delta_V_tr)``."""
    h_tr = _as_finite_float_array(
        transition_enthalpy, "transition_enthalpy")
    t = _as_finite_float_array(temperature, "temperature")
    dv = _as_finite_float_array(delta_molar_volume, "delta_molar_volume")
    _validate_positive_array(t, "temperature")
    if np.any(dv == 0.0):
        raise ValueError("delta_molar_volume must not be zero.")
    return _return_scalar_if_zero_dim(h_tr / (t * dv))


__all__ = [
    "_calc_phase_transition_entropy",
    "_calc_enthalpy_of_sublimation",
    "_calc_clapeyron_slope",
]
