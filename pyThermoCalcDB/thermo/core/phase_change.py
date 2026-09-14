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


def _calc_transition_enthalpy_from_clapeyron(
    temperature: NumericInput,
    transition_volume_change: NumericInput,
    dpressure_dtemperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``delta_H_tr = T*delta_V_tr*(dP/dT)``."""
    t = _as_finite_float_array(temperature, "temperature")
    dv = _as_finite_float_array(
        transition_volume_change,
        "transition_volume_change",
    )
    dpdt = _as_finite_float_array(dpressure_dtemperature, "dpressure_dtemperature")
    _validate_positive_array(t, "temperature")
    return _return_scalar_if_zero_dim(t * dv * dpdt)


def _calc_enthalpy_vaporization_watson(
    enthalpy_vaporization_reference: NumericInput,
    temperature_reference: NumericInput,
    temperature: NumericInput,
    critical_temperature: NumericInput,
    exponent: NumericInput = 0.38,
) -> float | NDArray[np.float64]:
    """Calculate Watson vaporization-enthalpy temperature correction."""
    h_ref = _as_finite_float_array(
        enthalpy_vaporization_reference,
        "enthalpy_vaporization_reference",
    )
    t_ref = _as_finite_float_array(temperature_reference, "temperature_reference")
    t = _as_finite_float_array(temperature, "temperature")
    tc = _as_finite_float_array(critical_temperature, "critical_temperature")
    n = _as_finite_float_array(exponent, "exponent")
    _validate_positive_array(h_ref, "enthalpy_vaporization_reference")
    _validate_positive_array(t_ref, "temperature_reference")
    _validate_positive_array(t, "temperature")
    _validate_positive_array(tc, "critical_temperature")
    # ! Classical Watson correction is evaluated below the critical point.
    if np.any(t_ref >= tc) or np.any(t >= tc):
        raise ValueError("temperature_reference and temperature must be below critical_temperature.")
    denominator = 1.0 - t_ref / tc
    if np.any(denominator <= 0.0):
        raise ValueError("reference reduced-temperature denominator must be positive.")
    factor = ((1.0 - t / tc) / denominator) ** n
    return _return_scalar_if_zero_dim(h_ref * factor)


def _calc_transition_enthalpy_from_constant_delta_cp(
    transition_enthalpy_reference: NumericInput,
    temperature_reference: NumericInput,
    temperature: NumericInput,
    delta_heat_capacity: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Kirchhoff correction with constant ``delta_Cp``."""
    h_ref = _as_finite_float_array(
        transition_enthalpy_reference,
        "transition_enthalpy_reference",
    )
    t_ref = _as_finite_float_array(temperature_reference, "temperature_reference")
    t = _as_finite_float_array(temperature, "temperature")
    delta_cp = _as_finite_float_array(delta_heat_capacity, "delta_heat_capacity")
    _validate_positive_array(t_ref, "temperature_reference")
    _validate_positive_array(t, "temperature")
    return _return_scalar_if_zero_dim(h_ref + delta_cp * (t - t_ref))


def _calc_transition_enthalpy_from_cp_integral(
    transition_enthalpy_reference: NumericInput,
    delta_cp_integral: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Kirchhoff correction from supplied ``integral(delta_Cp dT)``."""
    h_ref = _as_finite_float_array(
        transition_enthalpy_reference,
        "transition_enthalpy_reference",
    )
    integral = _as_finite_float_array(delta_cp_integral, "delta_cp_integral")
    return _return_scalar_if_zero_dim(h_ref + integral)


def _calc_sublimation_pressure_clapeyron(
    pressure_reference: NumericInput,
    sublimation_enthalpy: NumericInput,
    temperature_reference: NumericInput,
    temperature: NumericInput,
    gas_constant: NumericInput = 8.314462618,
) -> float | NDArray[np.float64]:
    """Calculate sublimation pressure from integrated Clapeyron relation."""
    p_ref = _as_finite_float_array(pressure_reference, "pressure_reference")
    h_sub = _as_finite_float_array(sublimation_enthalpy, "sublimation_enthalpy")
    t_ref = _as_finite_float_array(temperature_reference, "temperature_reference")
    t = _as_finite_float_array(temperature, "temperature")
    r = _as_finite_float_array(gas_constant, "gas_constant")
    _validate_positive_array(p_ref, "pressure_reference")
    _validate_positive_array(h_sub, "sublimation_enthalpy")
    _validate_positive_array(t_ref, "temperature_reference")
    _validate_positive_array(t, "temperature")
    _validate_positive_array(r, "gas_constant")
    ln_ratio = -(h_sub / r) * (1.0 / t - 1.0 / t_ref)
    return _return_scalar_if_zero_dim(p_ref * np.exp(ln_ratio))


__all__ = [
    "_calc_phase_transition_entropy",
    "_calc_enthalpy_of_sublimation",
    "_calc_clapeyron_slope",
    "_calc_transition_enthalpy_from_clapeyron",
    "_calc_enthalpy_vaporization_watson",
    "_calc_transition_enthalpy_from_constant_delta_cp",
    "_calc_transition_enthalpy_from_cp_integral",
    "_calc_sublimation_pressure_clapeyron",
]
