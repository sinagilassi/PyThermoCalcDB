"""Volume fraction and ideal-volume-additivity conversion helpers."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List
from pycuc import convert_from_to
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import fractions, non_negative, positive, same_shape


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
        volume_values = to_dict(
            volumes, output_volume_unit, unit_conversion_fn=conversion_fn)
        volume_values = _configure_component_values(
            volume_values,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "volumes"
        )
        total = sum(volume_values.values())
        if total <= 0.0:
            raise ValueError("Total volume must be positive.")
        return {key: value / total for key, value in volume_values.items()}

    # SECTION: Sequence implementation
    volume_values = to_list(volumes, output_volume_unit,
                            unit_conversion_fn=conversion_fn)
    total = sum(volume_values)
    if total <= 0.0:
        raise ValueError("Total volume must be positive.")
    return [value / total for value in volume_values]


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
        partial_volumes = {key: w[key] / rho[key] for key in w}
        total_volume = sum(partial_volumes.values())
        if total_volume <= 0.0:
            raise ValueError("Total ideal partial volume must be positive.")
        return {key: value / total_volume for key, value in partial_volumes.items()}

    if isinstance(mass_fractions, Mapping) or isinstance(densities, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    w = to_list(mass_fractions, unit_conversion_fn=conversion_fn)
    rho = to_list(
        densities, output_density_unit,
        unit_conversion_fn=conversion_fn
    )
    partial_volumes = [w_i / rho_i for w_i, rho_i in zip(w, rho)]
    total_volume = sum(partial_volumes)
    if total_volume <= 0.0:
        raise ValueError("Total ideal partial volume must be positive.")
    return [value / total_volume for value in partial_volumes]


# SECTION: Public exports
__all__ = ["calc_volume_fractions", "mass_fraction_to_volume_fraction"]
