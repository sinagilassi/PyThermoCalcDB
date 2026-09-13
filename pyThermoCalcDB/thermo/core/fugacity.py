"""Core fugacity and Poynting-factor transformations."""

# import libs
import math
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim


# SECTION: Numeric helpers

def _as_state_array(values: NumericArrayInput, name: str) -> NDArray[np.float64]:
    """Convert scalar or array-like input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_positive(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive values."""
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


def _validate_non_negative(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative values."""
    if np.any(values < 0.0):
        raise ValueError(f"{name} values must be non-negative.")


# SECTION: Core fugacity calculations

def _calc_poynting_factor_incompressible(
    liquid_molar_volume: NumericArrayInput,
    pressure: NumericArrayInput,
    saturation_pressure: NumericArrayInput,
    temperature: NumericArrayInput,
    gas_constant: NumericArrayInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate incompressible-liquid Poynting factor.

    Equation: ``F = exp(V_L*(P - P_sat)/(R*T))``.
    """
    v_l = _as_state_array(liquid_molar_volume, "liquid_molar_volume")
    p = _as_state_array(pressure, "pressure")
    p_sat = _as_state_array(saturation_pressure, "saturation_pressure")
    t = _as_state_array(temperature, "temperature")
    r = _as_state_array(gas_constant, "gas_constant")
    _validate_non_negative(v_l, "liquid_molar_volume")
    _validate_positive(t, "temperature")
    _validate_positive(r, "gas_constant")
    exponent = v_l * (p - p_sat) / (r * t)
    if np.any(exponent > math.log(np.finfo(np.float64).max)):
        raise OverflowError("Poynting-factor exponent exceeds float64 range.")
    return _return_scalar_if_zero_dim(np.exp(exponent))


def _calc_poynting_factor_from_integral(
    integral_vdp: NumericArrayInput,
    temperature: NumericArrayInput,
    gas_constant: NumericArrayInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate Poynting factor from a supplied ``integral(V dP)`` value."""
    integral = _as_state_array(integral_vdp, "integral_vdp")
    t = _as_state_array(temperature, "temperature")
    r = _as_state_array(gas_constant, "gas_constant")
    _validate_positive(t, "temperature")
    _validate_positive(r, "gas_constant")
    exponent = integral / (r * t)
    if np.any(exponent > math.log(np.finfo(np.float64).max)):
        raise OverflowError("Poynting-factor exponent exceeds float64 range.")
    return _return_scalar_if_zero_dim(np.exp(exponent))


def _calc_liquid_fugacity_coefficient(
    activity_coefficient: NumericArrayInput,
    saturated_fugacity_coefficient: NumericArrayInput,
    saturation_pressure: NumericArrayInput,
    pressure: NumericArrayInput,
    poynting_factor: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate liquid fugacity coefficient from supplied model factors."""
    gamma = _as_state_array(activity_coefficient, "activity_coefficient")
    phi_sat = _as_state_array(saturated_fugacity_coefficient, "saturated_fugacity_coefficient")
    p_sat = _as_state_array(saturation_pressure, "saturation_pressure")
    p = _as_state_array(pressure, "pressure")
    f_poynting = _as_state_array(poynting_factor, "poynting_factor")
    _validate_positive(gamma, "activity_coefficient")
    _validate_positive(phi_sat, "saturated_fugacity_coefficient")
    _validate_positive(p_sat, "saturation_pressure")
    _validate_positive(p, "pressure")
    _validate_positive(f_poynting, "poynting_factor")
    return _return_scalar_if_zero_dim(gamma * phi_sat * (p_sat / p) * f_poynting)


def _calc_liquid_partial_fugacity(
    mole_fraction: NumericArrayInput,
    activity_coefficient: NumericArrayInput,
    saturated_fugacity_coefficient: NumericArrayInput,
    saturation_pressure: NumericArrayInput,
    poynting_factor: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate liquid partial fugacity: ``f_i = x_i*gamma_i*phi_sat*P_sat*F``."""
    x = _as_state_array(mole_fraction, "mole_fraction")
    gamma = _as_state_array(activity_coefficient, "activity_coefficient")
    phi_sat = _as_state_array(saturated_fugacity_coefficient, "saturated_fugacity_coefficient")
    p_sat = _as_state_array(saturation_pressure, "saturation_pressure")
    f_poynting = _as_state_array(poynting_factor, "poynting_factor")
    _validate_non_negative(x, "mole_fraction")
    _validate_positive(gamma, "activity_coefficient")
    _validate_positive(phi_sat, "saturated_fugacity_coefficient")
    _validate_positive(p_sat, "saturation_pressure")
    _validate_positive(f_poynting, "poynting_factor")
    return _return_scalar_if_zero_dim(x * gamma * phi_sat * p_sat * f_poynting)


# SECTION: Core exports
__all__ = [
    "_calc_poynting_factor_incompressible",
    "_calc_poynting_factor_from_integral",
    "_calc_liquid_fugacity_coefficient",
    "_calc_liquid_partial_fugacity",
]
