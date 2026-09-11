"""Volume fraction and ideal-volume-additivity conversion helpers."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import fractions, non_negative, positive, same_shape
# locals
from ..utils.conversions import (
    _all_custom_props,
    _configure_component_values,
    _resolve_unit_conversion_fn,
)
from .core.volume_fraction import (
    _calc_mass_fraction_to_volume_fraction,
    _calc_mass_fraction_to_volume_fraction_from_mapping,
    _calc_mass_fraction_to_volume_fraction_from_props,
    _calc_volume_fractions,
    _calc_volume_fractions_from_mapping,
    _calc_volume_fractions_from_props,
)


# SECTION: Direct volume fractions

def calc_volume_fractions(
    volumes: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float] | list[float]:
    """Calculate volume fractions from component volumes.

    Parameters
    ----------
    volumes : mapping or sequence of float | int | CustomProp
        Component volumes. Values must be non-negative and at least one value
        must be positive.
    output_volume_unit : str, optional
        Unit used to normalize ``volumes`` before calculation.
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
    dict[str, float] or list[float]
        Component volume fractions.

    Notes
    -----
    Equation
        `phi_i = V_i / sum_j(V_j)`
    """
    # SECTION: Validate inputs
    non_negative(volumes, "volumes")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(volumes, Mapping):
        if _all_custom_props(volumes):
            return _calc_volume_fractions_from_props(
                volumes,
                output_volume_unit,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )

        # SECTION: Normalize mixed/numeric mapping inputs
        volume_values = to_dict(
            volumes,
            output_volume_unit,
            unit_conversion_fn=conversion_fn,
        )
        volume_values = _configure_component_values(
            volume_values,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "volumes",
        )
        return _calc_volume_fractions_from_mapping(volume_values)

    # SECTION: Sequence implementation
    volume_values = to_list(volumes, output_volume_unit,
                            unit_conversion_fn=conversion_fn)
    return _calc_volume_fractions(volume_values).tolist()


# SECTION: Ideal volume-additivity conversion

def mass_fraction_to_volume_fraction(
    mass_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    densities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float] | list[float]:
    """Convert mass fractions to volume fractions under ideal volume additivity.

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
    dict[str, float] or list[float]
        Component volume fractions.

    Notes
    -----
    Assumption
        Component volumes are additive and densities are evaluated at the same T,P.

    Equation
        `phi_i = (w_i/rho_i) / sum_j(w_j/rho_j)`
    """
    # SECTION: Validate inputs
    fractions(mass_fractions, "mass_fractions")
    positive(densities, "densities")
    same_shape(mass_fractions, densities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(mass_fractions, Mapping) and isinstance(densities, Mapping):
        if _all_custom_props(mass_fractions) and _all_custom_props(densities):
            return _calc_mass_fraction_to_volume_fraction_from_props(
                mass_fractions,
                densities,
                output_density_unit,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )

        # SECTION: Normalize mixed/numeric mapping inputs
        w = to_dict(mass_fractions, unit_conversion_fn=conversion_fn)
        rho = to_dict(
            densities,
            output_density_unit,
            unit_conversion_fn=conversion_fn,
        )
        w = _configure_component_values(
            w, components, component_key, case_sensitive, sort_by_components_order, "mass_fractions")
        rho = _configure_component_values(
            rho, components, component_key, case_sensitive, sort_by_components_order, "densities")
        return _calc_mass_fraction_to_volume_fraction_from_mapping(w, rho)

    if isinstance(mass_fractions, Mapping) or isinstance(densities, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    w = to_list(mass_fractions, unit_conversion_fn=conversion_fn)
    rho = to_list(
        densities, output_density_unit,
        unit_conversion_fn=conversion_fn
    )
    return _calc_mass_fraction_to_volume_fraction(w, rho).tolist()


# SECTION: Public exports
__all__ = ["calc_volume_fractions", "mass_fraction_to_volume_fraction"]
