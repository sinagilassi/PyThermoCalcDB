"""Ideal mixture density rules."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List
from pycuc import convert_from_to
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import fractions, positive, same_shape
# locals
from ..utils.conversions import _resolve_unit_conversion_fn


# SECTION: Internal helpers

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


# SECTION: Ideal density mixing rule

def calc_ideal_mixture_density(
    mass_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    densities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal mixture density from mass fractions and pure densities.

    Parameters
    ----------
    mass_fractions : mapping or sequence of float | int | CustomProp
        Component mass fractions.
    densities : mapping or sequence of float | int | CustomProp
        Pure-component densities at the same temperature and pressure.
    output_density_unit : str, optional
        Unit used to normalize ``densities`` before calculation.
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
        Ideal mixture density in the normalized density unit.

    Notes
    -----
    Assumption
        Component volumes are additive and densities are evaluated at the same T,P.

    Equation
        `rho_mix = 1 / sum_i(w_i/rho_i)`
    """
    # SECTION: Validate inputs
    fractions(mass_fractions, "mass_fractions")
    positive(densities, "densities")
    same_shape(mass_fractions, densities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(mass_fractions, Mapping) and isinstance(densities, Mapping):
        w = to_dict(mass_fractions, unit_conversion_fn=conversion_fn)
        rho = to_dict(
            densities,
            output_density_unit,
            unit_conversion_fn=conversion_fn
        )
        w = _configure_component_values(
            w, components, component_key, case_sensitive, sort_by_components_order, "mass_fractions")
        rho = _configure_component_values(
            rho, components, component_key, case_sensitive, sort_by_components_order, "densities")
        denominator = sum(w[key] / rho[key] for key in w)
        if denominator <= 0.0:
            raise ValueError(
                "The ideal-volume density denominator must be positive.")
        return 1.0 / denominator

    if isinstance(mass_fractions, Mapping) or isinstance(densities, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    w = to_list(mass_fractions, unit_conversion_fn=conversion_fn)
    rho = to_list(
        densities,
        output_density_unit,
        unit_conversion_fn=conversion_fn
    )
    denominator = sum(w_i / rho_i for w_i, rho_i in zip(w, rho))
    if denominator <= 0.0:
        raise ValueError(
            "The ideal-volume density denominator must be positive.")
    return 1.0 / denominator


# SECTION: Public exports
__all__ = ["calc_ideal_mixture_density"]
