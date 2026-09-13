"""Core derivative-property thermodynamic identities."""

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
def _calc_thermal_expansion_coefficient(
    volume: NumericInput,
    dvolume_dtemperature_at_pressure: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``alpha = (1/V) * (dV/dT)_P``."""
    v = _as_finite_float_array(volume, "volume")
    dvd_t = _as_finite_float_array(
        dvolume_dtemperature_at_pressure,
        "dvolume_dtemperature_at_pressure",
    )
    _validate_positive_array(v, "volume")
    return _return_scalar_if_zero_dim(dvd_t / v)


def _calc_isothermal_compressibility(
    volume: NumericInput,
    dvolume_dpressure_at_temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``kappa_T = -(1/V) * (dV/dP)_T``."""
    v = _as_finite_float_array(volume, "volume")
    dvd_p = _as_finite_float_array(
        dvolume_dpressure_at_temperature,
        "dvolume_dpressure_at_temperature",
    )
    _validate_positive_array(v, "volume")
    return _return_scalar_if_zero_dim(-dvd_p / v)


def _calc_joule_thomson_coefficient(
    temperature: NumericInput,
    molar_volume: NumericInput,
    heat_capacity_cp: NumericInput,
    dvolume_dtemperature_at_pressure: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``mu_JT = (T*(dV/dT)_P - V) / Cp``."""
    t = _as_finite_float_array(temperature, "temperature")
    v = _as_finite_float_array(molar_volume, "molar_volume")
    cp = _as_finite_float_array(heat_capacity_cp, "heat_capacity_cp")
    dvd_t = _as_finite_float_array(
        dvolume_dtemperature_at_pressure,
        "dvolume_dtemperature_at_pressure",
    )
    _validate_positive_array(t, "temperature")
    _validate_positive_array(v, "molar_volume")
    if np.any(cp == 0.0):
        raise ValueError("heat_capacity_cp must not be zero.")
    return _return_scalar_if_zero_dim((t * dvd_t - v) / cp)


__all__ = [
    "_calc_thermal_expansion_coefficient",
    "_calc_isothermal_compressibility",
    "_calc_joule_thomson_coefficient",
]
