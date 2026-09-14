"""Core critical-property relation calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_positive_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to a finite positive float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return cast(NDArray[np.float64], arr)


# SECTION: Core numeric calculations
def _calc_critical_compressibility_factor(
    critical_pressure: NumericInput,
    critical_molar_volume: NumericInput,
    critical_temperature: NumericInput,
    gas_constant: NumericInput = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate critical compressibility factor.

    Equation: ``Zc = Pc * Vc / (R * Tc)``. Inputs are expected in consistent
    SI units: Pa, m3/mol, K, and J/(mol.K). The output is dimensionless. This
    exact definition supports scalar, 1-D, and 2-D broadcast-compatible inputs.
    """
    pc = _as_positive_array(critical_pressure, "critical_pressure")
    vc = _as_positive_array(critical_molar_volume, "critical_molar_volume")
    tc = _as_positive_array(critical_temperature, "critical_temperature")
    r = _as_positive_array(gas_constant, "gas_constant")
    try:
        result = pc * vc / (r * tc)
    except ValueError as exc:
        raise ValueError("critical-property inputs must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], result))


def _calc_critical_volume_from_zc(
    critical_compressibility_factor: NumericInput,
    critical_temperature: NumericInput,
    critical_pressure: NumericInput,
    gas_constant: NumericInput = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate critical molar volume from ``Zc = Pc*Vc/(R*Tc)``.

    Inputs use SI-consistent units and positive finite values. The result is
    critical molar volume in m3/mol and preserves scalar/array behavior.
    """
    zc = _as_positive_array(critical_compressibility_factor, "critical_compressibility_factor")
    tc = _as_positive_array(critical_temperature, "critical_temperature")
    pc = _as_positive_array(critical_pressure, "critical_pressure")
    r = _as_positive_array(gas_constant, "gas_constant")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], zc * r * tc / pc))


def _calc_critical_pressure_from_zc(
    critical_compressibility_factor: NumericInput,
    critical_temperature: NumericInput,
    critical_molar_volume: NumericInput,
    gas_constant: NumericInput = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate critical pressure from ``Zc = Pc*Vc/(R*Tc)``.

    Inputs use SI-consistent units and positive finite values. The result is
    critical pressure in Pa and preserves scalar/array behavior.
    """
    zc = _as_positive_array(critical_compressibility_factor, "critical_compressibility_factor")
    tc = _as_positive_array(critical_temperature, "critical_temperature")
    vc = _as_positive_array(critical_molar_volume, "critical_molar_volume")
    r = _as_positive_array(gas_constant, "gas_constant")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], zc * r * tc / vc))


def _calc_critical_temperature_from_zc(
    critical_pressure: NumericInput,
    critical_molar_volume: NumericInput,
    critical_compressibility_factor: NumericInput,
    gas_constant: NumericInput = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate critical temperature from ``Zc = Pc*Vc/(R*Tc)``.

    Inputs use SI-consistent units and positive finite values. The result is
    critical temperature in K and preserves scalar/array behavior.
    """
    pc = _as_positive_array(critical_pressure, "critical_pressure")
    vc = _as_positive_array(critical_molar_volume, "critical_molar_volume")
    zc = _as_positive_array(critical_compressibility_factor, "critical_compressibility_factor")
    r = _as_positive_array(gas_constant, "gas_constant")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], pc * vc / (zc * r)))


__all__ = [
    "_calc_critical_compressibility_factor",
    "_calc_critical_volume_from_zc",
    "_calc_critical_pressure_from_zc",
    "_calc_critical_temperature_from_zc",
]
