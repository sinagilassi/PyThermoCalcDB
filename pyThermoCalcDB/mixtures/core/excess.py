"""Core excess-property and Gibbs-Duhem relations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import (
    NumericArrayInput,
    _as_float_array,
    _return_scalar_if_zero_dim,
    _validate_fraction_array,
    _validate_positive_array,
    _validate_same_array_shape,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Core excess-property calculations

def _calc_excess_property(
    real_property: NumericInput,
    ideal_property: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate excess property ``M^E = M - M_ideal``."""
    real = np.asarray(real_property, dtype=np.float64)
    ideal = np.asarray(ideal_property, dtype=np.float64)
    for name, arr in (("real_property", real), ("ideal_property", ideal)):
        if arr.ndim > 2:
            raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} values must be finite.")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], real - ideal))


def _calc_excess_gibbs_energy_from_activity_coefficients(
    mole_fractions: NumericInput,
    activity_coefficients: NumericInput,
    temperature: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate molar excess Gibbs energy ``G^E = R*T*sum_i x_i*ln(gamma_i)``."""
    x = _as_float_array(mole_fractions, "mole_fractions")
    gamma = _as_float_array(activity_coefficients, "activity_coefficients")
    t = np.asarray(temperature, dtype=np.float64)
    r = np.asarray(gas_constant, dtype=np.float64)
    _validate_same_array_shape(x, gamma, "mole_fractions", "activity_coefficients")
    _validate_fraction_array(x, "mole_fractions")
    _validate_positive_array(gamma, "activity_coefficients")
    for name, arr in (("temperature", t), ("gas_constant", r)):
        if arr.ndim > 2:
            raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
        if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
            raise ValueError(f"{name} values must be finite and greater than zero.")
    weighted_sum = np.sum(x * np.log(gamma), axis=-1)
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], r * t * weighted_sum))


def _calc_excess_entropy_from_gibbs_enthalpy(
    excess_gibbs_energy: NumericInput,
    excess_enthalpy: NumericInput,
    temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate excess entropy from ``G^E = H^E - T*S^E``."""
    g_e = np.asarray(excess_gibbs_energy, dtype=np.float64)
    h_e = np.asarray(excess_enthalpy, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    for name, arr in (
        ("excess_gibbs_energy", g_e),
        ("excess_enthalpy", h_e),
        ("temperature", t),
    ):
        if arr.ndim > 2:
            raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} values must be finite.")
    if np.any(t <= 0.0):
        raise ValueError("temperature values must be greater than zero.")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], (h_e - g_e) / t))


def _calc_gibbs_duhem_residual(
    mole_fractions: NumericInput,
    dlog_activity_coefficients: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Gibbs-Duhem residual ``sum_i x_i*dln(gamma_i)``."""
    x = _as_float_array(mole_fractions, "mole_fractions")
    dln_gamma = _as_float_array(
        dlog_activity_coefficients,
        "dlog_activity_coefficients",
    )
    _validate_same_array_shape(x, dln_gamma, "mole_fractions", "dlog_activity_coefficients")
    _validate_fraction_array(x, "mole_fractions")
    return _return_scalar_if_zero_dim(np.sum(x * dln_gamma, axis=-1))


def _check_gibbs_duhem_consistency(
    mole_fractions: NumericInput,
    dlog_activity_coefficients: NumericInput,
    tolerance: float = 1.0e-8,
) -> bool | NDArray[np.bool_]:
    """Check whether Gibbs-Duhem residual is close to zero."""
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")
    residual = np.asarray(
        _calc_gibbs_duhem_residual(mole_fractions, dlog_activity_coefficients),
        dtype=np.float64,
    )
    result = np.abs(residual) <= tolerance
    if result.ndim == 0:
        return bool(result)
    return cast(NDArray[np.bool_], result)


# SECTION: Core exports
__all__ = [
    "_calc_excess_property",
    "_calc_excess_gibbs_energy_from_activity_coefficients",
    "_calc_excess_entropy_from_gibbs_enthalpy",
    "_calc_gibbs_duhem_residual",
    "_check_gibbs_duhem_consistency",
]
