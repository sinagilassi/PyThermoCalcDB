"""Core volume-fraction calculations."""

# import libs
from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp, UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict
# locals
from ...utils.conversions import (
    NumericArrayInput,
    _as_float_array,
    _configure_component_values,
    _resolve_unit_conversion_fn,
    _validate_custom_prop_mapping,
    _validate_fraction_array,
    _validate_non_negative_array,
    _validate_positive_array,
    _validate_same_array_shape,
    _validate_same_mapping_keys,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Direct volume fractions

def _calc_volume_fractions(
    volumes: NumericInput,
) -> NDArray[np.float64]:
    """Calculate phi_i = V_i / sum_j(V_j).

    For 2-D inputs, axis 0 is states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    volume_values = _as_float_array(volumes, "volumes")
    _validate_non_negative_array(volume_values, "volumes")

    # SECTION: Calculate normalized volume fractions
    total = np.sum(volume_values, axis=-1, keepdims=True)
    if np.any(total <= 0.0):
        raise ValueError("Total volume must be positive.")
    return volume_values / total


# SECTION: Direct mapping adapter

def _calc_volume_fractions_from_mapping(
    volumes: Mapping[str, float | int],
) -> dict[str, float]:
    """Calculate keyed volume fractions from keyed volumes."""
    # NOTE: Single-input mapping order is preserved in the result.
    fractions = _calc_volume_fractions(list(volumes.values())).tolist()
    return dict(zip(volumes.keys(), fractions))


# SECTION: Ideal volume-additivity conversion

def _calc_mass_fraction_to_volume_fraction(
    mass_fractions: NumericInput,
    densities: NumericInput,
) -> NDArray[np.float64]:
    """Calculate phi_i = (w_i / rho_i) / sum_j(w_j / rho_j).

    For 2-D inputs, axis 0 is states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    w = _as_float_array(mass_fractions, "mass_fractions")
    rho = _as_float_array(densities, "densities")
    _validate_same_array_shape(w, rho, "mass_fractions", "densities")
    _validate_fraction_array(w, "mass_fractions")
    _validate_positive_array(rho, "densities")

    # SECTION: Calculate ideal partial volumes and normalize
    partial_volumes = w / rho
    total = np.sum(partial_volumes, axis=-1, keepdims=True)
    if np.any(total <= 0.0):
        raise ValueError("Total ideal partial volume must be positive.")
    return partial_volumes / total


# SECTION: Ideal conversion mapping adapter

def _calc_mass_fraction_to_volume_fraction_from_mapping(
    mass_fractions: Mapping[str, float | int],
    densities: Mapping[str, float | int],
) -> dict[str, float]:
    """Calculate keyed volume fractions from keyed mass fractions and densities."""
    # SECTION: Align by caller mass-fraction order after key validation
    _validate_same_mapping_keys(mass_fractions, densities, "mass_fractions", "densities")
    keys = list(mass_fractions)
    fractions = _calc_mass_fraction_to_volume_fraction(
        [mass_fractions[key] for key in keys],
        [densities[key] for key in keys],
    ).tolist()
    return dict(zip(keys, fractions))


# SECTION: Props adapters

def _calc_volume_fractions_from_props(
    volumes: Mapping[str, CustomProp],
    output_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Calculate keyed volume fractions from unit-aware keyed volumes."""
    # SECTION: Validate props input contract
    _validate_custom_prop_mapping(volumes, "volumes")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    volume_values = to_dict(
        volumes,
        output_volume_unit,
        unit_conversion_fn=conversion_fn,
    )

    # SECTION: Remap component keys
    volume_values = _configure_component_values(
        volume_values,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "volumes",
    )

    # SECTION: Calculate from mapping
    return _calc_volume_fractions_from_mapping(volume_values)


def _calc_mass_fraction_to_volume_fraction_from_props(
    mass_fractions: Mapping[str, CustomProp],
    densities: Mapping[str, CustomProp],
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Calculate keyed volume fractions from unit-aware mass fractions and densities."""
    # SECTION: Validate props input contract
    _validate_custom_prop_mapping(mass_fractions, "mass_fractions")
    _validate_custom_prop_mapping(densities, "densities")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    w = to_dict(mass_fractions, unit_conversion_fn=conversion_fn)
    rho = to_dict(
        densities,
        output_density_unit,
        unit_conversion_fn=conversion_fn,
    )

    # SECTION: Remap component keys
    w = _configure_component_values(
        w,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "mass_fractions",
    )
    rho = _configure_component_values(
        rho,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "densities",
    )

    # SECTION: Calculate from aligned mapping
    return _calc_mass_fraction_to_volume_fraction_from_mapping(w, rho)


# SECTION: Core exports
__all__ = [
    "_calc_volume_fractions",
    "_calc_volume_fractions_from_mapping",
    "_calc_mass_fraction_to_volume_fraction",
    "_calc_mass_fraction_to_volume_fraction_from_mapping",
    "_calc_volume_fractions_from_props",
    "_calc_mass_fraction_to_volume_fraction_from_props",
]
