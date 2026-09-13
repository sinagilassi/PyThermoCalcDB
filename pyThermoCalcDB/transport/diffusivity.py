"""Public gas and liquid diffusivity correlations."""

# import libs
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, Pressure, Temperature
from pythermodb_settings.models.units import UnitConversionFn
# locals
from ..utils.conversions import _pos, _resolve_unit_conversion_fn, _to_kelvin
from .core.diffusivity import (
    _calc_fuller_schettler_giddings_diffusivity,
    _calc_wilke_chang_diffusivity,
    _calc_wilke_lee_diffusivity,
)


# SECTION: Public gas diffusivity

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


def calc_fuller_schettler_giddings_diffusivity(
    temperature,
    pressure,
    molecular_weight_i,
    molecular_weight_j,
    diffusion_volume_i,
    diffusion_volume_j,
    output_unit: str = "m2/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | CustomProp | NDArray[np.float64]:
    """Calculate Fuller-Schettler-Giddings dilute binary gas diffusivity.

    This empirical correlation assumes dilute, low-pressure gas behavior and
    Fuller diffusion volumes supplied by the caller.
    """
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature, "temperature", "K", conversion_fn)
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    mw_i = _pos(molecular_weight_i, "molecular_weight_i",
                "kg/mol", conversion_fn)
    mw_j = _pos(molecular_weight_j, "molecular_weight_j",
                "kg/mol", conversion_fn)
    v_i = _pos(diffusion_volume_i, "diffusion_volume_i")
    v_j = _pos(diffusion_volume_j, "diffusion_volume_j")
    value = _calc_fuller_schettler_giddings_diffusivity(
        t, p, mw_i, mw_j, v_i, v_j)
    if output_unit != "m2/s":
        value = conversion_fn(value=float(
            value), from_unit="m2/s", to_unit=output_unit)
    if any(isinstance(item, CustomProp) for item in (molecular_weight_i, molecular_weight_j)) or output_unit != "m2/s":
        return CustomProp(value=float(value), unit=output_unit)
    return value


def calc_wilke_lee_diffusivity(
    temperature,
    pressure,
    molecular_weight_i,
    molecular_weight_j,
    sigma_i,
    sigma_j,
    epsilon_over_k_i,
    epsilon_over_k_j,
    output_unit: str = "m2/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | CustomProp | NDArray[np.float64]:
    """Calculate Wilke-Lee dilute binary gas diffusivity."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature, "temperature", "K", conversion_fn)
    p = _pressure_to_pa(pressure, "pressure", conversion_fn)
    mw_i = _pos(molecular_weight_i, "molecular_weight_i",
                "kg/mol", conversion_fn)
    mw_j = _pos(molecular_weight_j, "molecular_weight_j",
                "kg/mol", conversion_fn)
    value = _calc_wilke_lee_diffusivity(
        t, p, mw_i, mw_j, sigma_i, sigma_j, epsilon_over_k_i, epsilon_over_k_j)
    if output_unit != "m2/s":
        value = conversion_fn(value=float(
            value), from_unit="m2/s", to_unit=output_unit)
    if any(isinstance(item, CustomProp) for item in (molecular_weight_i, molecular_weight_j)) or output_unit != "m2/s":
        return CustomProp(value=float(value), unit=output_unit)
    return value


# SECTION: Public liquid diffusivity

def calc_wilke_chang_diffusivity(
    temperature,
    solvent_viscosity,
    solvent_molecular_weight,
    solute_molar_volume_at_normal_boiling_point,
    solvent_association_factor: float = 1.0,
    output_unit: str = "m2/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | CustomProp | NDArray[np.float64]:
    """Calculate Wilke-Chang infinite-dilution liquid diffusivity."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature, "temperature", "K", conversion_fn)
    mu = _pos(solvent_viscosity, "solvent_viscosity", "Pa*s", conversion_fn)
    mw = _pos(solvent_molecular_weight,
              "solvent_molecular_weight", "kg/mol", conversion_fn)
    v_a = _pos(
        solute_molar_volume_at_normal_boiling_point,
        "solute_molar_volume_at_normal_boiling_point",
        "m3/mol",
        conversion_fn,
    )
    phi = _pos(solvent_association_factor, "solvent_association_factor")
    value = _calc_wilke_chang_diffusivity(t, mu, mw, v_a, phi)
    if output_unit != "m2/s":
        value = conversion_fn(value=float(
            value), from_unit="m2/s", to_unit=output_unit)
    if any(isinstance(item, CustomProp) for item in (solvent_viscosity, solvent_molecular_weight, solute_molar_volume_at_normal_boiling_point)) or output_unit != "m2/s":
        return CustomProp(value=float(value), unit=output_unit)
    return value


# SECTION: Public exports
__all__ = [
    "calc_fuller_schettler_giddings_diffusivity",
    "calc_wilke_lee_diffusivity",
    "calc_wilke_chang_diffusivity",
]
