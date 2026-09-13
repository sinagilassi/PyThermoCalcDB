"""Gas flow relations."""

# import libs
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, Pressure, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn

# locals
from ..utils.conversions import _pos, _resolve_unit_conversion_fn
from .core.gas import (
    _calc_gas_molar_flow_rate_from_z,
    _calc_gas_volumetric_flow_rate_from_z,
    _calc_ideal_gas_molar_flow_rate,
    _calc_ideal_gas_volumetric_flow_rate,
)


def _as_public_scalar(value: float | NDArray[np.float64], name: str) -> float:
    """Return a scalar result for public scalar wrappers."""
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return float(value)
        raise ValueError(f"{name} must be scalar for this public API.")
    return float(value)


def _temperature_k(temperature: Temperature, unit_conversion_fn: UnitConversionFn | None = None) -> float:
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    return float(conversion_fn(float(temperature.value), temperature.unit, "K")) if temperature.unit != "K" else float(temperature.value)


def _pressure_pa(pressure: Pressure, unit_conversion_fn: UnitConversionFn | None = None) -> float:
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    return float(conversion_fn(float(pressure.value), pressure.unit, "Pa")) if pressure.unit != "Pa" else float(pressure.value)


def calc_ideal_gas_volumetric_flow_rate(
    molar_flow_rate: ScalarValue,
    temperature: Temperature,
    pressure: Pressure,
    output_molar_flow_unit: str | None = "mol/s",
    output_unit: str = "m3/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> CustomProp:
    """Calculate ideal-gas volumetric flow rate from molar flow rate."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    f = _pos(molar_flow_rate, "molar_flow_rate", output_molar_flow_unit, conversion_fn)
    q = _calc_ideal_gas_volumetric_flow_rate(f, _temperature_k(temperature, conversion_fn), _pressure_pa(pressure, conversion_fn))
    q_value = _as_public_scalar(q, "volumetric_flow_rate")
    if output_unit != "m3/s":
        q_value = conversion_fn(q_value, "m3/s", output_unit)
    return CustomProp(value=q_value, unit=output_unit)


def calc_ideal_gas_molar_flow_rate(
    volumetric_flow_rate: ScalarValue,
    temperature: Temperature,
    pressure: Pressure,
    output_volumetric_flow_unit: str | None = "m3/s",
    output_unit: str = "mol/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> CustomProp:
    """Calculate ideal-gas molar flow rate from volumetric flow rate."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    q = _pos(volumetric_flow_rate, "volumetric_flow_rate", output_volumetric_flow_unit, conversion_fn)
    f = _calc_ideal_gas_molar_flow_rate(q, _temperature_k(temperature, conversion_fn), _pressure_pa(pressure, conversion_fn))
    f_value = _as_public_scalar(f, "molar_flow_rate")
    if output_unit != "mol/s":
        f_value = conversion_fn(f_value, "mol/s", output_unit)
    return CustomProp(value=f_value, unit=output_unit)


def calc_gas_volumetric_flow_rate_from_z(
    molar_flow_rate: ScalarValue,
    temperature: Temperature,
    pressure: Pressure,
    compressibility_factor: ScalarValue,
    output_molar_flow_unit: str | None = "mol/s",
    output_unit: str = "m3/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> CustomProp:
    """Calculate gas volumetric flow rate from a supplied compressibility factor."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    f = _pos(molar_flow_rate, "molar_flow_rate", output_molar_flow_unit, conversion_fn)
    z = _pos(compressibility_factor, "compressibility_factor")
    q = _calc_gas_volumetric_flow_rate_from_z(f, _temperature_k(temperature, conversion_fn), _pressure_pa(pressure, conversion_fn), z)
    q_value = _as_public_scalar(q, "volumetric_flow_rate")
    if output_unit != "m3/s":
        q_value = conversion_fn(q_value, "m3/s", output_unit)
    return CustomProp(value=q_value, unit=output_unit)


def calc_gas_molar_flow_rate_from_z(
    volumetric_flow_rate: ScalarValue,
    temperature: Temperature,
    pressure: Pressure,
    compressibility_factor: ScalarValue,
    output_volumetric_flow_unit: str | None = "m3/s",
    output_unit: str = "mol/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> CustomProp:
    """Calculate gas molar flow rate from volumetric flow rate and supplied Z."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    q = _pos(volumetric_flow_rate, "volumetric_flow_rate", output_volumetric_flow_unit, conversion_fn)
    z = _pos(compressibility_factor, "compressibility_factor")
    f = _calc_gas_molar_flow_rate_from_z(q, _temperature_k(temperature, conversion_fn), _pressure_pa(pressure, conversion_fn), z)
    f_value = _as_public_scalar(f, "molar_flow_rate")
    if output_unit != "mol/s":
        f_value = conversion_fn(f_value, "mol/s", output_unit)
    return CustomProp(value=f_value, unit=output_unit)


__all__ = [
    "calc_ideal_gas_volumetric_flow_rate",
    "calc_ideal_gas_molar_flow_rate",
    "calc_gas_volumetric_flow_rate_from_z",
    "calc_gas_molar_flow_rate_from_z",
]
