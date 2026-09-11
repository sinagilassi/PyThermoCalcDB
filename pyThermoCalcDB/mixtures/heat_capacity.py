"""Ideal mixture heat-capacity rules."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import fractions, positive, same_shape
# locals
from ..utils.conversions import (
    _all_custom_props,
    _configure_component_values,
    _resolve_unit_conversion_fn,
)
from .core.heat_capacity import (
    _calc_ideal_mixture_heat_capacity,
    _calc_ideal_mixture_heat_capacity_from_mapping,
    _calc_ideal_mixture_heat_capacity_from_props,
)


# SECTION: Ideal heat-capacity mixing rule

def calc_ideal_mixture_heat_capacity(
    fraction_values: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    heat_capacities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal mixture heat capacity from component heat capacities.

    Parameters
    ----------
    fraction_values : mapping or sequence of float | int | CustomProp
        Mole fractions for molar heat capacities, or mass fractions for
        mass-specific heat capacities.
    heat_capacities : mapping or sequence of float | int | CustomProp
        Component heat capacities on the same basis, for example J/mol/K or
        J/kg/K.
    output_heat_capacity_unit : str, optional
        Unit used to normalize ``heat_capacities`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.
    components : list[Component], optional
        Components used to remap and order mapping-input values when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format for remapping mapping-input keys.
    case_sensitive : bool, optional
        Whether component ID matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether mapping values should follow the order of ``components``.

    Returns
    -------
    float
        Ideal mixture heat capacity on the same basis as component heat capacities.

    Notes
    -----
    Equation
        `Cp_mix = sum_i(f_i*Cp_i), where f_i is x_i or w_i.`
    """
    # SECTION: Validate inputs
    fractions(fraction_values, "fractions")
    positive(heat_capacities, "heat_capacities")
    same_shape(fraction_values, heat_capacities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(fraction_values, Mapping) and isinstance(heat_capacities, Mapping):
        if _all_custom_props(fraction_values) and _all_custom_props(heat_capacities):
            return _calc_ideal_mixture_heat_capacity_from_props(
                fraction_values,
                heat_capacities,
                output_heat_capacity_unit,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )

        # SECTION: Normalize mixed/numeric mapping inputs
        fr = to_dict(fraction_values, unit_conversion_fn=conversion_fn)
        cp = to_dict(
            heat_capacities,
            output_heat_capacity_unit,
            unit_conversion_fn=conversion_fn,
        )
        fr = _configure_component_values(
            fr, components, component_key, case_sensitive, sort_by_components_order, "fractions")
        cp = _configure_component_values(
            cp, components, component_key, case_sensitive, sort_by_components_order, "heat_capacities")
        return _calc_ideal_mixture_heat_capacity_from_mapping(fr, cp)

    if isinstance(fraction_values, Mapping) or isinstance(heat_capacities, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    fr = to_list(fraction_values, unit_conversion_fn=conversion_fn)
    cp = to_list(
        heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn=conversion_fn
    )
    return float(_calc_ideal_mixture_heat_capacity(fr, cp))


# SECTION: Public exports
__all__ = ["calc_ideal_mixture_heat_capacity"]
