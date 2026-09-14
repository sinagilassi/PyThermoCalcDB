"""Model-neutral fugacity transformations."""

# import libs
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, Pressure, Temperature
from pythermodb_settings.models.units import UnitConversionFn
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import _pos, _resolve_unit_conversion_fn, _scalar, _to_kelvin
from .core.fugacity import (
    _calc_liquid_fugacity_coefficient,
    _calc_liquid_partial_fugacity,
    _calc_poynting_factor_from_integral,
    _calc_poynting_factor_incompressible,
    _calc_fugacity_coefficient,
    _calc_fugacity_from_coefficient,
    _calc_phase_equilibrium_fugacity_residual,
)


# SECTION: Public Poynting helpers

def _pressure_to_pa(
    pressure,
    name: str,
    unit_conversion_fn,
) -> float:
    """Normalize pressure-like scalar input to Pa."""
    # SECTION: Pressure model normalization
    if isinstance(pressure, Pressure):
        value = float(pressure.value)
        if pressure.unit != "Pa":
            value = float(unit_conversion_fn(value, pressure.unit, "Pa"))
        if value <= 0.0:
            raise ValueError(f"{name} must be greater than zero.")
        return value
    return _pos(pressure, name, "Pa", unit_conversion_fn) if isinstance(pressure, CustomProp) else _pos(pressure, name)


def calc_poynting_factor_incompressible(
    liquid_molar_volume,
    pressure,
    saturation_pressure,
    temperature,
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate the incompressible-liquid Poynting factor.

    Unit-aware scalar inputs are normalized to ``m3/mol``, ``Pa``, and ``K``.
    Numeric values are assumed to already be on those bases.
    """
    # SECTION: Normalize unit-bearing scalar inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    v_l = _pos(liquid_molar_volume, "liquid_molar_volume",
               "m3/mol", conversion_fn)
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    p_sat = _pressure_to_pa(saturation_pressure,
                            "saturation_pressure", conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature, "temperature", "K", conversion_fn)
    r = _pos(gas_constant, "gas_constant")
    return _calc_poynting_factor_incompressible(v_l, p, p_sat, t, r)


def calc_poynting_factor_from_integral(
    integral_vdp,
    temperature,
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate Poynting factor from a supplied ``integral(V dP)`` value."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    integral = _scalar(integral_vdp, "integral_vdp", "J/mol", conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature, "temperature", "K", conversion_fn)
    r = _pos(gas_constant, "gas_constant")
    return _calc_poynting_factor_from_integral(integral, t, r)


# SECTION: Public fugacity transformations

def calc_liquid_fugacity_coefficient(
    activity_coefficient,
    saturated_fugacity_coefficient,
    saturation_pressure,
    pressure,
    poynting_factor,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate liquid fugacity coefficient from supplied model-neutral factors."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    gamma = _pos(activity_coefficient, "activity_coefficient")
    phi_sat = _pos(saturated_fugacity_coefficient,
                   "saturated_fugacity_coefficient")
    p_sat = _pressure_to_pa(saturation_pressure,
                            "saturation_pressure", conversion_fn)
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    f_poynting = _pos(poynting_factor, "poynting_factor")
    return _calc_liquid_fugacity_coefficient(gamma, phi_sat, p_sat, p, f_poynting)


def calc_liquid_partial_fugacity(
    mole_fraction,
    activity_coefficient,
    saturated_fugacity_coefficient,
    saturation_pressure,
    poynting_factor,
    output_pressure_unit: str = "Pa",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | CustomProp | NDArray[np.float64]:
    """Calculate liquid partial fugacity from supplied activity and fugacity factors."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = _scalar(mole_fraction, "mole_fraction")
    gamma = _pos(activity_coefficient, "activity_coefficient")
    phi_sat = _pos(saturated_fugacity_coefficient,
                   "saturated_fugacity_coefficient")
    p_sat = _pressure_to_pa(saturation_pressure,
                            "saturation_pressure", conversion_fn)
    f_poynting = _pos(poynting_factor, "poynting_factor")
    value = _calc_liquid_partial_fugacity(x, gamma, phi_sat, p_sat, f_poynting)
    if output_pressure_unit != "Pa":
        value = conversion_fn(value=float(
            value), from_unit="Pa", to_unit=output_pressure_unit)
    if isinstance(saturation_pressure, (CustomProp, Pressure)) or output_pressure_unit != "Pa":
        return CustomProp(value=float(value), unit=output_pressure_unit)
    return value



def calc_fugacity_coefficient(
    fugacity,
    mole_fraction,
    pressure,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate fugacity coefficient ``phi_i = f_i/(y_i*P)``."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    f = _pressure_to_pa(fugacity, "fugacity", conversion_fn)
    y = _pos(mole_fraction, "mole_fraction")
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    return _calc_fugacity_coefficient(f, y, p)


def calc_fugacity_from_coefficient(
    fugacity_coefficient,
    mole_fraction,
    pressure,
    output_pressure_unit: str = "Pa",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | CustomProp | NDArray[np.float64]:
    """Calculate fugacity from coefficient ``f_i = phi_i*y_i*P``."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    phi = _pos(fugacity_coefficient, "fugacity_coefficient")
    y = _scalar(mole_fraction, "mole_fraction")
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    value = _calc_fugacity_from_coefficient(phi, y, p)
    if output_pressure_unit != "Pa":
        value = conversion_fn(value=float(value), from_unit="Pa", to_unit=output_pressure_unit)
    if isinstance(pressure, (CustomProp, Pressure)) or output_pressure_unit != "Pa":
        return CustomProp(value=float(value), unit=output_pressure_unit)
    return value


def calc_phase_equilibrium_fugacity_residual(
    fugacity_phase_1,
    fugacity_phase_2,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate phase-equilibrium fugacity residual ``f1 - f2``."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    f1 = _pressure_to_pa(fugacity_phase_1, "fugacity_phase_1", conversion_fn)
    f2 = _pressure_to_pa(fugacity_phase_2, "fugacity_phase_2", conversion_fn)
    return _calc_phase_equilibrium_fugacity_residual(f1, f2)
# SECTION: Public exports
__all__ = [
    "calc_poynting_factor_incompressible",
    "calc_poynting_factor_from_integral",
    "calc_liquid_fugacity_coefficient",
    "calc_liquid_partial_fugacity",
    "calc_fugacity_coefficient",
    "calc_fugacity_from_coefficient",
    "calc_phase_equilibrium_fugacity_residual",
]


