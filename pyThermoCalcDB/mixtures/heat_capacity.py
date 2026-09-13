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
    _get_all_custom_props,
    _configure_component_values,
    _resolve_unit_conversion_fn,
)
from .core.heat_capacity import (
    _calc_ideal_mixture_heat_capacity,
    _calc_ideal_mixture_heat_capacity_from_mapping,
    _calc_ideal_mixture_heat_capacity_from_props,
    _calc_total_heat_capacity,
    _calc_total_heat_capacity_from_mapping,
    _calc_total_heat_capacity_from_props,
)


# SECTION: Ideal heat-capacity mixing rule

def calc_ideal_mixture_heat_capacity_from_alls(
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
            fraction_custom_props: Mapping[str, CustomProp] = dict(
                zip(fraction_values.keys(), _get_all_custom_props(fraction_values, return_type="list"))
            )
            heat_capacity_custom_props: Mapping[str, CustomProp] = dict(
                zip(heat_capacities.keys(), _get_all_custom_props(heat_capacities, return_type="list"))
            )
            return _calc_ideal_mixture_heat_capacity_from_props(
                fraction_custom_props,
                heat_capacity_custom_props,
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

# ! ::: Ideal heat-capacity mixing rule from sequences

def calc_ideal_mixture_heat_capacity_from_sequence(
    fraction_values: Sequence[float | int | CustomProp],
    heat_capacities: Sequence[float | int | CustomProp],
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal mixture heat capacity from sequence inputs."""
    if isinstance(fraction_values, Mapping) or isinstance(heat_capacities, Mapping):
        raise TypeError("Both component inputs must be sequences.")

    return calc_ideal_mixture_heat_capacity_from_alls(
        fraction_values,
        heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn,
    )


# ! ::: Ideal heat-capacity mixing rule from mappings

def calc_ideal_mixture_heat_capacity_from_mapping(
    fraction_values: Mapping[str, float | int],
    heat_capacities: Mapping[str, float | int],
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal mixture heat capacity from numeric mapping inputs."""
    if not isinstance(fraction_values, Mapping) or not isinstance(heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings.")
    if _all_custom_props(fraction_values) or _all_custom_props(heat_capacities):
        raise TypeError("CustomProp mappings must use calc_ideal_mixture_heat_capacity_from_props.")

    return calc_ideal_mixture_heat_capacity_from_alls(
        fraction_values,
        heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


# ! ::: Ideal heat-capacity mixing rule from props

def calc_ideal_mixture_heat_capacity_from_props(
    fraction_values: Mapping[str, CustomProp],
    heat_capacities: Mapping[str, CustomProp],
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal mixture heat capacity from unit-aware mapping inputs."""
    if not isinstance(fraction_values, Mapping) or not isinstance(heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings of CustomProp instances.")
    if not _all_custom_props(fraction_values) or not _all_custom_props(heat_capacities):
        raise TypeError("Both component mappings must contain only CustomProp instances.")

    return calc_ideal_mixture_heat_capacity_from_alls(
        fraction_values,
        heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


# SECTION: Backwards-compatible aliases
calc_ideal_mixture_heat_capacity_from_all = calc_ideal_mixture_heat_capacity_from_alls
calc_ideal_mixture_heat_capacity = calc_ideal_mixture_heat_capacity_from_all


# SECTION: Extensive heat-capacity calculation

def calc_total_heat_capacity_from_alls(
    component_moles: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molar_heat_capacities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_moles_unit: str | None = "mol",
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total heat capacity for a finite component inventory.

    Equation: ``Cp_total = sum_i(n_i*Cp_i)``. With moles in ``mol`` and molar
    heat capacities in ``J/mol.K``, the result is in ``J/K``.
    """
    # SECTION: Validate inputs
    positive(molar_heat_capacities, "molar_heat_capacities")
    same_shape(component_moles, molar_heat_capacities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(component_moles, Mapping) and isinstance(molar_heat_capacities, Mapping):
        if _all_custom_props(component_moles) and _all_custom_props(molar_heat_capacities):
            moles_custom_props: Mapping[str, CustomProp] = dict(
                zip(component_moles.keys(), _get_all_custom_props(component_moles, return_type="list"))
            )
            heat_capacity_custom_props: Mapping[str, CustomProp] = dict(
                zip(molar_heat_capacities.keys(), _get_all_custom_props(molar_heat_capacities, return_type="list"))
            )
            return _calc_total_heat_capacity_from_props(
                moles_custom_props,
                heat_capacity_custom_props,
                output_moles_unit,
                output_heat_capacity_unit,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )

        # SECTION: Normalize mixed/numeric mapping inputs
        n = to_dict(
            component_moles,
            output_moles_unit,
            unit_conversion_fn=conversion_fn,
        )
        cp = to_dict(
            molar_heat_capacities,
            output_heat_capacity_unit,
            unit_conversion_fn=conversion_fn,
        )
        n = _configure_component_values(
            n, components, component_key, case_sensitive, sort_by_components_order, "component_moles")
        cp = _configure_component_values(
            cp, components, component_key, case_sensitive, sort_by_components_order, "molar_heat_capacities")
        return _calc_total_heat_capacity_from_mapping(n, cp)

    if isinstance(component_moles, Mapping) or isinstance(molar_heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    n = to_list(
        component_moles,
        output_moles_unit,
        unit_conversion_fn=conversion_fn,
    )
    cp = to_list(
        molar_heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn=conversion_fn,
    )
    return float(_calc_total_heat_capacity(n, cp))


def calc_total_heat_capacity_from_sequence(
    component_moles: Sequence[float | int | CustomProp],
    molar_heat_capacities: Sequence[float | int | CustomProp],
    output_moles_unit: str | None = "mol",
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate total heat capacity from sequence inputs."""
    if isinstance(component_moles, Mapping) or isinstance(molar_heat_capacities, Mapping):
        raise TypeError("Both component inputs must be sequences.")

    return calc_total_heat_capacity_from_alls(
        component_moles,
        molar_heat_capacities,
        output_moles_unit,
        output_heat_capacity_unit,
        unit_conversion_fn,
    )


def calc_total_heat_capacity_from_mapping(
    component_moles: Mapping[str, float | int],
    molar_heat_capacities: Mapping[str, float | int],
    output_moles_unit: str | None = "mol",
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total heat capacity from numeric mapping inputs."""
    if not isinstance(component_moles, Mapping) or not isinstance(molar_heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings.")
    if _all_custom_props(component_moles) or _all_custom_props(molar_heat_capacities):
        raise TypeError("CustomProp mappings must use calc_total_heat_capacity_from_props.")

    return calc_total_heat_capacity_from_alls(
        component_moles,
        molar_heat_capacities,
        output_moles_unit,
        output_heat_capacity_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


def calc_total_heat_capacity_from_props(
    component_moles: Mapping[str, CustomProp],
    molar_heat_capacities: Mapping[str, CustomProp],
    output_moles_unit: str | None = "mol",
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total heat capacity from unit-aware mapping inputs."""
    if not isinstance(component_moles, Mapping) or not isinstance(molar_heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings of CustomProp instances.")
    if not _all_custom_props(component_moles) or not _all_custom_props(molar_heat_capacities):
        raise TypeError("Both component mappings must contain only CustomProp instances.")

    return calc_total_heat_capacity_from_alls(
        component_moles,
        molar_heat_capacities,
        output_moles_unit,
        output_heat_capacity_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


calc_total_heat_capacity_from_all = calc_total_heat_capacity_from_alls
calc_total_heat_capacity = calc_total_heat_capacity_from_all


# SECTION: Public exports
__all__ = [
    "calc_ideal_mixture_heat_capacity_from_alls",
    "calc_ideal_mixture_heat_capacity_from_all",
    "calc_ideal_mixture_heat_capacity_from_sequence",
    "calc_ideal_mixture_heat_capacity_from_mapping",
    "calc_ideal_mixture_heat_capacity_from_props",
    "calc_ideal_mixture_heat_capacity",
    "calc_total_heat_capacity_from_alls",
    "calc_total_heat_capacity_from_all",
    "calc_total_heat_capacity_from_sequence",
    "calc_total_heat_capacity_from_mapping",
    "calc_total_heat_capacity_from_props",
    "calc_total_heat_capacity",
]
