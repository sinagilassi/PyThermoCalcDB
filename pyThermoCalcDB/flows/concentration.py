"""Concentration and flow identities."""

# import libs
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, ScalarValue
from pythermodb_settings.models.units import UnitConversionFn

# locals
from ..utils.conversions import _pos, _resolve_unit_conversion_fn
from .core.concentration import (
    _calc_concentration_from_molar_flow_rate,
    _calc_molar_flow_rate_from_concentration,
)


def _as_public_scalar(value: float | NDArray[np.float64], name: str) -> float:
    """Return a scalar result for public scalar wrappers."""
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return float(value)
        raise ValueError(f"{name} must be scalar for this public API.")
    return float(value)


def calc_molar_flow_rate_from_concentration(
    concentration: ScalarValue,
    volumetric_flow_rate: ScalarValue,
    output_concentration_unit: str | None = "mol/m3",
    output_volumetric_flow_unit: str | None = "m3/s",
    output_unit: str = "mol/s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> CustomProp:
    """Calculate molar flow rate from concentration and volumetric flow rate."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    c = _pos(concentration, "concentration",
             output_concentration_unit, conversion_fn)
    q = _pos(volumetric_flow_rate, "volumetric_flow_rate",
             output_volumetric_flow_unit, conversion_fn)
    f = _calc_molar_flow_rate_from_concentration(c, q)
    f_value = _as_public_scalar(f, "molar_flow_rate")
    if output_unit != "mol/s":
        f_value = conversion_fn(
            value=f_value,
            from_unit="mol/s",
            to_unit=output_unit
        )
    return CustomProp(value=f_value, unit=output_unit)


def calc_concentration_from_molar_flow_rate(
    molar_flow_rate: ScalarValue,
    volumetric_flow_rate: ScalarValue,
    output_molar_flow_unit: str | None = "mol/s",
    output_volumetric_flow_unit: str | None = "m3/s",
    output_unit: str = "mol/m3",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> CustomProp:
    """Calculate concentration from molar flow rate and volumetric flow rate."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    f = _pos(molar_flow_rate, "molar_flow_rate",
             output_molar_flow_unit, conversion_fn)
    q = _pos(volumetric_flow_rate, "volumetric_flow_rate",
             output_volumetric_flow_unit, conversion_fn)
    c = _calc_concentration_from_molar_flow_rate(f, q)
    c_value = _as_public_scalar(c, "concentration")
    if output_unit != "mol/m3":
        c_value = conversion_fn(
            value=c_value,
            from_unit="mol/m3",
            to_unit=output_unit
        )
    return CustomProp(value=c_value, unit=output_unit)


__all__ = [
    "calc_molar_flow_rate_from_concentration",
    "calc_concentration_from_molar_flow_rate",
]
