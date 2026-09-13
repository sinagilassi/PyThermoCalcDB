"""Mixture volume approximations."""

# import libs
from collections.abc import Mapping, Sequence

from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import positive, same_shape

# locals
from ..utils.conversions import _resolve_unit_conversion_fn
from .core.volume import (
    _calc_additive_liquid_volume,
    _calc_additive_liquid_volume_from_mapping,
)


# SECTION: Additive liquid volume approximation

def calc_additive_liquid_volume(
    component_moles: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molecular_weights: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    component_densities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_moles_unit: str | None = "mol",
    output_molecular_weight_unit: str | None = "kg/mol",
    output_density_unit: str | None = "kg/m3",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate additive liquid mixture volume.

    This is an additive-volume approximation:
    ``V_mix ~= sum_i(n_i*M_i/rho_i)``. It is not an exact nonideal-mixture
    relation; rigorous mixture volumes require partial molar volumes.
    """
    # SECTION: Validate inputs
    positive(molecular_weights, "molecular_weights")
    positive(component_densities, "component_densities")
    same_shape(component_moles, molecular_weights)
    same_shape(component_moles, component_densities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if (
        isinstance(component_moles, Mapping)
        and isinstance(molecular_weights, Mapping)
        and isinstance(component_densities, Mapping)
    ):
        n = to_dict(component_moles, output_moles_unit, unit_conversion_fn=conversion_fn)
        mw = to_dict(molecular_weights, output_molecular_weight_unit, unit_conversion_fn=conversion_fn)
        rho = to_dict(component_densities, output_density_unit, unit_conversion_fn=conversion_fn)
        return _calc_additive_liquid_volume_from_mapping(n, mw, rho)

    if (
        isinstance(component_moles, Mapping)
        or isinstance(molecular_weights, Mapping)
        or isinstance(component_densities, Mapping)
    ):
        raise TypeError("All component inputs must be mappings or all sequences.")

    # SECTION: Sequence implementation
    n = to_list(component_moles, output_moles_unit, unit_conversion_fn=conversion_fn)
    mw = to_list(molecular_weights, output_molecular_weight_unit, unit_conversion_fn=conversion_fn)
    rho = to_list(component_densities, output_density_unit, unit_conversion_fn=conversion_fn)
    return float(_calc_additive_liquid_volume(n, mw, rho))


# SECTION: Public exports
__all__ = ["calc_additive_liquid_volume"]
