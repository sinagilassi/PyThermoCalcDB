"""Core gas flow relations."""

# import libs
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput
from ._common import (
    _as_flow_float_array,
    _return_scalar_if_zero_dim,
    _validate_non_negative,
    _validate_positive,
)

# SECTION: Type aliases
NumericInput: TypeAlias = NumericArrayInput


def _calc_ideal_gas_volumetric_flow_rate(
    molar_flow_rate: NumericInput,
    temperature: NumericInput,
    pressure: NumericInput,
    gas_constant: float = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate ideal-gas volumetric flow rate: Vdot = ndot*R*T/P."""
    f = _as_flow_float_array(molar_flow_rate, "molar_flow_rate")
    t = _as_flow_float_array(temperature, "temperature")
    p = _as_flow_float_array(pressure, "pressure")
    r = _as_flow_float_array(gas_constant, "gas_constant")
    _validate_non_negative(f, "molar_flow_rate")
    _validate_positive(t, "temperature")
    _validate_positive(p, "pressure")
    _validate_positive(r, "gas_constant")
    return _return_scalar_if_zero_dim(f * r * t / p)


def _calc_ideal_gas_molar_flow_rate(
    volumetric_flow_rate: NumericInput,
    temperature: NumericInput,
    pressure: NumericInput,
    gas_constant: float = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate ideal-gas molar flow rate: ndot = P*Vdot/(R*T)."""
    q = _as_flow_float_array(volumetric_flow_rate, "volumetric_flow_rate")
    t = _as_flow_float_array(temperature, "temperature")
    p = _as_flow_float_array(pressure, "pressure")
    r = _as_flow_float_array(gas_constant, "gas_constant")
    _validate_non_negative(q, "volumetric_flow_rate")
    _validate_positive(t, "temperature")
    _validate_positive(p, "pressure")
    _validate_positive(r, "gas_constant")
    return _return_scalar_if_zero_dim(p * q / (r * t))


def _calc_gas_volumetric_flow_rate_from_z(
    molar_flow_rate: NumericInput,
    temperature: NumericInput,
    pressure: NumericInput,
    compressibility_factor: NumericInput,
    gas_constant: float = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate gas volumetric flow rate using supplied Z: Vdot = Z*ndot*R*T/P."""
    z = _as_flow_float_array(compressibility_factor, "compressibility_factor")
    _validate_positive(z, "compressibility_factor")
    return _return_scalar_if_zero_dim(
        z * np.asarray(_calc_ideal_gas_volumetric_flow_rate(molar_flow_rate, temperature, pressure, gas_constant))
    )


def _calc_gas_molar_flow_rate_from_z(
    volumetric_flow_rate: NumericInput,
    temperature: NumericInput,
    pressure: NumericInput,
    compressibility_factor: NumericInput,
    gas_constant: float = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate gas molar flow rate using supplied Z: ndot = P*Vdot/(Z*R*T)."""
    z = _as_flow_float_array(compressibility_factor, "compressibility_factor")
    _validate_positive(z, "compressibility_factor")
    q = _as_flow_float_array(volumetric_flow_rate, "volumetric_flow_rate")
    t = _as_flow_float_array(temperature, "temperature")
    p = _as_flow_float_array(pressure, "pressure")
    r = _as_flow_float_array(gas_constant, "gas_constant")
    _validate_non_negative(q, "volumetric_flow_rate")
    _validate_positive(t, "temperature")
    _validate_positive(p, "pressure")
    _validate_positive(r, "gas_constant")
    return _return_scalar_if_zero_dim(p * q / (z * r * t))


__all__ = [
    "_calc_ideal_gas_volumetric_flow_rate",
    "_calc_ideal_gas_molar_flow_rate",
    "_calc_gas_volumetric_flow_rate_from_z",
    "_calc_gas_molar_flow_rate_from_z",
]
