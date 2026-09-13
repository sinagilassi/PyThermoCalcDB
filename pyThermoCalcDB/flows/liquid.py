"""Liquid-flow approximations."""

# import libs
from collections.abc import Mapping, Sequence

from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import positive, same_shape

# locals
from ..utils.conversions import _resolve_unit_conversion_fn
from .core.liquid import (
    _calc_additive_liquid_volumetric_flow_rate,
    _calc_additive_liquid_volumetric_flow_rate_from_mapping,
)


def calc_additive_liquid_volumetric_flow_rate(
    molar_flow_rates: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molecular_weights: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    component_densities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molar_flow_unit: str | None = "mol/s",
    output_molecular_weight_unit: str | None = "kg/mol",
    output_density_unit: str | None = "kg/m3",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate additive liquid volumetric flow approximation.

    Equation: ``Vdot ~= sum_i(F_i*M_i/rho_i)``.
    """
    positive(molecular_weights, "molecular_weights")
    positive(component_densities, "component_densities")
    same_shape(molar_flow_rates, molecular_weights)
    same_shape(molar_flow_rates, component_densities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    if (
        isinstance(molar_flow_rates, Mapping)
        and isinstance(molecular_weights, Mapping)
        and isinstance(component_densities, Mapping)
    ):
        f = to_dict(molar_flow_rates, output_molar_flow_unit, unit_conversion_fn=conversion_fn)
        mw = to_dict(molecular_weights, output_molecular_weight_unit, unit_conversion_fn=conversion_fn)
        rho = to_dict(component_densities, output_density_unit, unit_conversion_fn=conversion_fn)
        return _calc_additive_liquid_volumetric_flow_rate_from_mapping(f, mw, rho)

    if (
        isinstance(molar_flow_rates, Mapping)
        or isinstance(molecular_weights, Mapping)
        or isinstance(component_densities, Mapping)
    ):
        raise TypeError("All component inputs must be mappings or all sequences.")

    f = to_list(molar_flow_rates, output_molar_flow_unit, unit_conversion_fn=conversion_fn)
    mw = to_list(molecular_weights, output_molecular_weight_unit, unit_conversion_fn=conversion_fn)
    rho = to_list(component_densities, output_density_unit, unit_conversion_fn=conversion_fn)
    return float(_calc_additive_liquid_volumetric_flow_rate(f, mw, rho))


__all__ = ["calc_additive_liquid_volumetric_flow_rate"]
