"""Ideal mixture heat-capacity rules."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List
from pycuc import convert_from_to
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import fractions, positive, same_shape


# SECTION: Internal helpers
def _resolve_unit_conversion_fn(
    unit_conversion_fn: UnitConversionFn | None,
) -> UnitConversionFn:
    """Return the provided converter or the module default converter."""
    return convert_from_to if unit_conversion_fn is None else unit_conversion_fn


def _configure_component_values(
    values: Mapping[str, float],
    components: Optional[List[Component]],
    component_key: Optional[ComponentKey],
    case_sensitive: bool,
    sort_by_components_order: bool,
    name: str,
) -> dict[str, float]:
    """Remap and order mapping values using component metadata when requested."""
    # ! If no component key is provided, preserve original keys and order.
    if component_key is None:
        return dict(values)
    if not components:
        raise ValueError(
            f"component_key is provided but components is empty for {name}.")
    component_values = config_components_values(
        values=dict(values),
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    if component_values is None:
        raise ValueError(f"Failed to configure {name} component values.")
    component_values_dict, _ = component_values
    return component_values_dict


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
        fr = to_dict(fraction_values, unit_conversion_fn=conversion_fn)
        cp = to_dict(
            heat_capacities,
            output_heat_capacity_unit,
            unit_conversion_fn=conversion_fn
        )
        fr = _configure_component_values(
            fr, components, component_key, case_sensitive, sort_by_components_order, "fractions")
        cp = _configure_component_values(
            cp, components, component_key, case_sensitive, sort_by_components_order, "heat_capacities")
        return sum(fr[key] * cp[key] for key in fr)

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
    return sum(fr_i * cp_i for fr_i, cp_i in zip(fr, cp))


# SECTION: Public exports
__all__ = ["calc_ideal_mixture_heat_capacity"]
