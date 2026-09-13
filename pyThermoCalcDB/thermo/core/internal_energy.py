"""Core internal-energy identity calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import (
    CustomProp,
    ScalarValue,
    Temperature,
    UnitConversionFn,
)

# locals
from pythermocalcdb.utils.conversions import (
    NumericArrayInput,
    _generic_temperature,
    _pos,
    _resolve_unit_conversion_fn,
    _return_scalar_if_zero_dim,
    _scalar,
    _validate_custom_prop_scalar,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators

def _as_finite_array(
    values: NumericInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert numeric inputs to finite float64 arrays."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, 1-D, or 2-D values.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _as_positive_array(
    values: NumericInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert numeric inputs to finite positive float64 arrays."""
    arr = _as_finite_array(values, name)
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return arr


# SECTION: Core numeric calculations

def _calc_internal_energy(
    enthalpy: NumericInput,
    pressure: NumericInput,
    volume: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``U = H - P*V`` using finite numeric inputs.

    Inputs may be scalars, 1-D arrays, or 2-D arrays when their shapes are
    broadcast-compatible. For 2-D inputs, axis 0 is states and axis 1 is
    component/property values.
    """
    # SECTION: Normalize and validate
    h = _as_finite_array(enthalpy, "enthalpy")
    p = _as_positive_array(pressure, "pressure")
    v = _as_positive_array(volume, "volume")

    try:
        internal_energy = h - p * v
    except ValueError as exc:
        raise ValueError(
            "enthalpy, pressure, and volume must be broadcast-compatible."
        ) from exc

    return _return_scalar_if_zero_dim(
        cast(NDArray[np.float64], internal_energy)
    )


def _calc_ideal_gas_internal_energy(
    molar_enthalpy: NumericInput,
    temperature: NumericInput,
    universal_gas_constant: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ideal-gas ``U_molar = H_molar - R*T``."""
    # SECTION: Normalize and validate
    h = _as_finite_array(molar_enthalpy, "molar_enthalpy")
    t = _as_positive_array(temperature, "temperature")
    r = _as_positive_array(universal_gas_constant, "universal_gas_constant")

    try:
        internal_energy = h - r * t
    except ValueError as exc:
        raise ValueError(
            "molar_enthalpy, temperature, and universal_gas_constant must be "
            "broadcast-compatible."
        ) from exc

    return _return_scalar_if_zero_dim(
        cast(NDArray[np.float64], internal_energy)
    )


# SECTION: Props adapters

def _calc_internal_energy_from_props(
    enthalpy: CustomProp,
    pressure: CustomProp,
    volume: CustomProp,
    output_enthalpy_unit: str | None = None,
    output_pressure_unit: str | None = None,
    output_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate internal energy from unit-aware scalar inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(enthalpy, "enthalpy")
    _validate_custom_prop_scalar(pressure, "pressure")
    _validate_custom_prop_scalar(volume, "volume")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    h = _scalar(enthalpy, "enthalpy", output_enthalpy_unit, conversion_fn)
    p = _pos(pressure, "pressure", output_pressure_unit, conversion_fn)
    v = _pos(volume, "volume", output_volume_unit, conversion_fn)

    # SECTION: Calculate internal energy
    return float(_calc_internal_energy(h, p, v))


def _calc_ideal_gas_internal_energy_from_props(
    molar_enthalpy: CustomProp,
    temperature: Temperature,
    output_molar_enthalpy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    universal_gas_constant: ScalarValue = 8.31446261815324,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal-gas internal energy from unit-aware scalar inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(molar_enthalpy, "molar_enthalpy")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    h = _scalar(
        molar_enthalpy,
        "molar_enthalpy",
        output_molar_enthalpy_unit,
        conversion_fn,
    )
    t = _generic_temperature(
        temperature,
        output_temperature_unit,
        conversion_fn,
    )
    r = _pos(universal_gas_constant, "universal_gas_constant")

    # SECTION: Calculate ideal-gas internal energy
    return float(_calc_ideal_gas_internal_energy(h, t, r))


# SECTION: Scalar adapters

def _calc_internal_energy_from_scalars(
    enthalpy: ScalarValue,
    pressure: ScalarValue,
    volume: ScalarValue,
    output_enthalpy_unit: str | None = None,
    output_pressure_unit: str | None = None,
    output_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate internal energy."""
    # SECTION: Normalize scalar inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    h = _scalar(enthalpy, "enthalpy", output_enthalpy_unit, conversion_fn)
    p = _pos(pressure, "pressure", output_pressure_unit, conversion_fn)
    v = _pos(volume, "volume", output_volume_unit, conversion_fn)

    # SECTION: Calculate internal energy
    return float(_calc_internal_energy(h, p, v))


def _calc_ideal_gas_internal_energy_from_scalars(
    molar_enthalpy: ScalarValue,
    temperature: Temperature,
    output_molar_enthalpy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    universal_gas_constant: ScalarValue = 8.31446261815324,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate ideal-gas internal energy."""
    # SECTION: Normalize scalar inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    h = _scalar(
        molar_enthalpy,
        "molar_enthalpy",
        output_molar_enthalpy_unit,
        conversion_fn,
    )
    t = _generic_temperature(
        temperature,
        output_temperature_unit,
        conversion_fn,
    )
    r = _pos(universal_gas_constant, "universal_gas_constant")

    # SECTION: Calculate ideal-gas internal energy
    return float(_calc_ideal_gas_internal_energy(h, t, r))


# SECTION: Core exports
__all__ = [
    "_calc_internal_energy",
    "_calc_ideal_gas_internal_energy",
    "_calc_internal_energy_from_props",
    "_calc_ideal_gas_internal_energy_from_props",
    "_calc_internal_energy_from_scalars",
    "_calc_ideal_gas_internal_energy_from_scalars",
]
