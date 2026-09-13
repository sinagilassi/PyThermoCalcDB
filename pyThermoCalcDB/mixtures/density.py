"""Ideal mixture density rules."""

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
from .core.density import (
    _calc_ideal_mixture_density,
    _calc_ideal_mixture_density_from_mapping,
    _calc_ideal_mixture_density_from_props,
)


# ! ::: Ideal density mixing rule from sequences

def calc_ideal_mixture_density_from_sequence(
    mass_fractions: Sequence[float | int | CustomProp],
    densities: Sequence[float | int | CustomProp],
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """
    Calculate ideal mixture density from mass fractions and pure densities using sequence.

    Parameters
    ----------
    mass_fractions : sequence of float | int | CustomProp
        Component mass fractions.
    densities : sequence of float | int | CustomProp
        Pure-component densities at the same temperature and pressure.
    output_density_unit : str, optional
        Unit used to normalize ``densities`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

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

    # SECTION: Sequence implementation
    w = to_list(
        mass_fractions,
        unit_conversion_fn=conversion_fn
    )
    rho = to_list(
        densities,
        output_density_unit,
        unit_conversion_fn=conversion_fn
    )
    return float(_calc_ideal_mixture_density(w, rho))

# ! ::: Ideal density mixing rule from mappings


def calc_ideal_mixture_density_from_mapping(
    mass_fractions: Mapping[str, float | int],
    densities: Mapping[str, float | int],
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """
    Calculate ideal mixture density from mass fractions and pure densities using mappings.

    Parameters
    ----------
    mass_fractions : mapping or sequence of float | int
        Component mass fractions.
    densities : mapping or sequence of float | int
        Pure-component densities at the same temperature and pressure.
    output_density_unit : str, optional
        Unit used to normalize ``densities`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

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

    # SECTION: Normalize mixed/numeric mapping inputs
    w = to_dict(mass_fractions)
    rho = to_dict(
        densities,
        output_density_unit,
        unit_conversion_fn=conversion_fn,
    )

    return _calc_ideal_mixture_density_from_mapping(w, rho)

# !::: Ideal mixture density from props


def calc_ideal_mixture_density_from_props(
    mass_fractions: Mapping[str, CustomProp],
    densities: Mapping[str, CustomProp],
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """
    Calculate ideal mixture density from mass fractions and pure densities using mappings.

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
        if _all_custom_props(mass_fractions) and _all_custom_props(densities):
            return _calc_ideal_mixture_density_from_props(
                mass_fractions,
                densities,
                output_density_unit,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )
        else:
            raise ValueError(
                "Unsupported input types for mass_fractions and densities.")
    else:
        raise ValueError(
            "Unsupported input types for mass_fractions and densities.")


# SECTION: Public exports
__all__ = [
    "calc_ideal_mixture_density_from_sequence",
    "calc_ideal_mixture_density_from_mapping",
    "calc_ideal_mixture_density_from_props",
]
