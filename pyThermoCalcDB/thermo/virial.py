"""Virial equation-of-state transformations."""

# import libs
from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, Pressure, Temperature
from pythermodb_settings.models.units import UnitConversionFn
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import _pos, _resolve_unit_conversion_fn, _scalar, _to_kelvin
from .core.virial import (
    _calc_compressibility_from_second_virial,
    _calc_compressibility_from_virial_density_form,
    _calc_compressibility_from_virial_pressure_form,
    _calc_pressure_from_virial_density_form,
    _calc_second_virial_from_compressibility,
)


# SECTION: Unit helpers

def _pressure_to_pa(
    pressure,
    name: str,
    unit_conversion_fn,
) -> float:
    """Normalize pressure-like scalar input to Pa."""
    if isinstance(pressure, Pressure):
        value = float(pressure.value)
        if pressure.unit != "Pa":
            value = float(unit_conversion_fn(value, pressure.unit, "Pa"))
        if value <= 0.0:
            raise ValueError(f"{name} must be greater than zero.")
        return value
    return _pos(pressure, name, "Pa", unit_conversion_fn) if isinstance(pressure, CustomProp) else _pos(pressure, name)


# SECTION: Public second-virial transformations

def calc_compressibility_from_second_virial(
    second_virial_coefficient,
    pressure,
    temperature,
    gas_constant: float = R_J_molK,
    output_second_virial_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate compressibility factor from a supplied second virial coefficient.

    Numeric ``B`` is assumed to be in ``m3/mol`` unless a unit-aware
    ``CustomProp`` and ``output_second_virial_unit`` are supplied.
    """
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    b_unit = output_second_virial_unit or "m3/mol"
    b = _scalar(
        second_virial_coefficient,
        "second_virial_coefficient",
        b_unit if isinstance(second_virial_coefficient, CustomProp) else None,
        conversion_fn,
    )
    if b_unit != "m3/mol":
        b = float(conversion_fn(b, b_unit, "m3/mol"))
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    r = _pos(gas_constant, "gas_constant")
    return _calc_compressibility_from_second_virial(b, p, t, r)


def calc_second_virial_from_compressibility(
    compressibility_factor,
    pressure,
    temperature,
    gas_constant: float = R_J_molK,
    output_unit: str = "m3/mol",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate second virial coefficient from supplied compressibility factor."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    z = _pos(compressibility_factor, "compressibility_factor")
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    r = _pos(gas_constant, "gas_constant")
    value = float(_calc_second_virial_from_compressibility(z, p, t, r))
    if output_unit != "m3/mol":
        value = float(conversion_fn(value, "m3/mol", output_unit))
    return value


# SECTION: Public general virial forms

def calc_compressibility_from_virial_density_form(
    molar_density,
    virial_coefficients: Sequence[float | int] | NDArray[np.number],
    output_molar_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate compressibility factor from density-form virial coefficients."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    rho_unit = output_molar_density_unit or "mol/m3"
    rho = _pos(
        molar_density,
        "molar_density",
        rho_unit if isinstance(molar_density, CustomProp) else None,
        conversion_fn,
    )
    if rho_unit != "mol/m3":
        rho = float(conversion_fn(rho, rho_unit, "mol/m3"))
    return _calc_compressibility_from_virial_density_form(rho, virial_coefficients)


def calc_pressure_from_virial_density_form(
    molar_density,
    temperature,
    virial_coefficients: Sequence[float | int] | NDArray[np.number],
    gas_constant: float = R_J_molK,
    output_pressure_unit: str = "Pa",
    output_molar_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate pressure from density-form virial EOS."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    rho_unit = output_molar_density_unit or "mol/m3"
    rho = _pos(
        molar_density,
        "molar_density",
        rho_unit if isinstance(molar_density, CustomProp) else None,
        conversion_fn,
    )
    if rho_unit != "mol/m3":
        rho = float(conversion_fn(rho, rho_unit, "mol/m3"))
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    r = _pos(gas_constant, "gas_constant")
    value = float(_calc_pressure_from_virial_density_form(rho, t, virial_coefficients, r))
    if output_pressure_unit != "Pa":
        value = float(conversion_fn(value, "Pa", output_pressure_unit))
    return value


def calc_compressibility_from_virial_pressure_form(
    pressure,
    virial_coefficients: Sequence[float | int] | NDArray[np.number],
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate compressibility factor from pressure-form virial coefficients."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    return _calc_compressibility_from_virial_pressure_form(p, virial_coefficients)


# SECTION: Public exports
__all__ = [
    "calc_compressibility_from_second_virial",
    "calc_second_virial_from_compressibility",
    "calc_compressibility_from_virial_density_form",
    "calc_pressure_from_virial_density_form",
    "calc_compressibility_from_virial_pressure_form",
]
