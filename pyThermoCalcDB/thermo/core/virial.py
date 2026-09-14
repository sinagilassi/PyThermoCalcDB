"""Core virial equation-of-state transformations."""

# import libs
from collections.abc import Sequence
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim


# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators

def _as_state_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert scalar or array-like input to a finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_positive(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive numeric values."""
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


def _as_coefficients(
    virial_coefficients: Sequence[float | int] | NDArray[np.number],
) -> NDArray[np.float64]:
    """Convert virial coefficients to a finite one-dimensional float64 array."""
    coeffs = np.asarray(virial_coefficients, dtype=np.float64)
    if coeffs.ndim != 1 or coeffs.size == 0:
        raise ValueError("virial_coefficients must be a non-empty one-dimensional sequence.")
    if not np.all(np.isfinite(coeffs)):
        raise ValueError("virial_coefficients values must be finite.")
    return cast(NDArray[np.float64], coeffs)


# SECTION: Second-virial transformations

def _calc_compressibility_from_second_virial(
    second_virial_coefficient: NumericInput,
    pressure: NumericInput,
    temperature: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate ``Z = 1 + B*P/(R*T)`` from a second virial coefficient."""
    b = _as_state_array(second_virial_coefficient, "second_virial_coefficient")
    p = _as_state_array(pressure, "pressure")
    t = _as_state_array(temperature, "temperature")
    r = _as_state_array(gas_constant, "gas_constant")
    _validate_positive(p, "pressure")
    _validate_positive(t, "temperature")
    _validate_positive(r, "gas_constant")
    return _return_scalar_if_zero_dim(1.0 + b * p / (r * t))


def _calc_second_virial_from_compressibility(
    compressibility_factor: NumericInput,
    pressure: NumericInput,
    temperature: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate ``B = (Z - 1)*R*T/P`` from compressibility factor."""
    z = _as_state_array(compressibility_factor, "compressibility_factor")
    p = _as_state_array(pressure, "pressure")
    t = _as_state_array(temperature, "temperature")
    r = _as_state_array(gas_constant, "gas_constant")
    _validate_positive(z, "compressibility_factor")
    _validate_positive(p, "pressure")
    _validate_positive(t, "temperature")
    _validate_positive(r, "gas_constant")
    return _return_scalar_if_zero_dim((z - 1.0) * r * t / p)


# SECTION: General virial forms

def _calc_compressibility_from_virial_density_form(
    molar_density: NumericInput,
    virial_coefficients: Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate density-form virial ``Z = 1 + B*rho + C*rho^2 + ...``."""
    rho = _as_state_array(molar_density, "molar_density")
    coeffs = _as_coefficients(virial_coefficients)
    _validate_positive(rho, "molar_density")
    z = np.ones_like(rho, dtype=np.float64)
    # NOTE: Coefficients are ordered as [B, C, D, ...], powers start at one.
    for power, coefficient in enumerate(coeffs, start=1):
        z = z + coefficient * np.power(rho, power)
    return _return_scalar_if_zero_dim(z)


def _calc_pressure_from_virial_density_form(
    molar_density: NumericInput,
    temperature: NumericInput,
    virial_coefficients: Sequence[float | int] | NDArray[np.number],
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate pressure from density-form virial EOS: ``P = Z*rho*R*T``."""
    rho = _as_state_array(molar_density, "molar_density")
    t = _as_state_array(temperature, "temperature")
    r = _as_state_array(gas_constant, "gas_constant")
    _validate_positive(rho, "molar_density")
    _validate_positive(t, "temperature")
    _validate_positive(r, "gas_constant")
    z = np.asarray(
        _calc_compressibility_from_virial_density_form(rho, virial_coefficients),
        dtype=np.float64,
    )
    return _return_scalar_if_zero_dim(z * rho * r * t)


def _calc_compressibility_from_virial_pressure_form(
    pressure: NumericInput,
    virial_coefficients: Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate pressure-form virial ``Z = 1 + Bp*P + Cp*P^2 + ...``."""
    p = _as_state_array(pressure, "pressure")
    coeffs = _as_coefficients(virial_coefficients)
    _validate_positive(p, "pressure")
    z = np.ones_like(p, dtype=np.float64)
    # NOTE: Coefficients are ordered as [B', C', D', ...], powers start at one.
    for power, coefficient in enumerate(coeffs, start=1):
        z = z + coefficient * np.power(p, power)
    return _return_scalar_if_zero_dim(z)


# SECTION: Core exports
__all__ = [
    "_calc_compressibility_from_second_virial",
    "_calc_second_virial_from_compressibility",
    "_calc_compressibility_from_virial_density_form",
    "_calc_pressure_from_virial_density_form",
    "_calc_compressibility_from_virial_pressure_form",
]
