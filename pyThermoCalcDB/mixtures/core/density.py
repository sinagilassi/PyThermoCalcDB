"""Core ideal mixture density calculations."""

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
    _return_scalar_if_zero_dim,
    _validate_custom_prop_mapping,
    _validate_fraction_array,
    _validate_positive_array,
    _validate_same_array_shape,
    _validate_same_mapping_keys,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Ideal density mixing rule

def _calc_ideal_mixture_density(
    mass_fractions: NumericInput,
    densities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate rho_mix = 1 / sum_i(w_i / rho_i).

    For 2-D inputs, axis 0 is states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    w = _as_float_array(mass_fractions, "mass_fractions")
    rho = _as_float_array(densities, "densities")
    _validate_same_array_shape(w, rho, "mass_fractions", "densities")
    _validate_fraction_array(w, "mass_fractions")
    _validate_positive_array(rho, "densities")

    # SECTION: Calculate reciprocal ideal-volume density
    denominator = np.sum(w / rho, axis=-1)
    if np.any(denominator <= 0.0):
        raise ValueError("The ideal-volume density denominator must be positive.")
    return _return_scalar_if_zero_dim(1.0 / denominator)


# SECTION: Mapping adapter

def _calc_ideal_mixture_density_from_mapping(
    mass_fractions: Mapping[str, float | int],
    densities: Mapping[str, float | int],
) -> float:
    """Calculate ideal mixture density from aligned keyed inputs."""
    # SECTION: Align by caller mass-fraction order after key validation
    _validate_same_mapping_keys(mass_fractions, densities, "mass_fractions", "densities")
    keys = list(mass_fractions)
    return float(
        _calc_ideal_mixture_density(
            [mass_fractions[key] for key in keys],
            [densities[key] for key in keys],
        )
    )


# SECTION: Props adapter

def _calc_ideal_mixture_density_from_props(
    mass_fractions: Mapping[str, CustomProp],
    densities: Mapping[str, CustomProp],
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal mixture density from unit-aware keyed inputs."""
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
    return _calc_ideal_mixture_density_from_mapping(w, rho)


# SECTION: Core exports
__all__ = [
    "_calc_ideal_mixture_density",
    "_calc_ideal_mixture_density_from_mapping",
    "_calc_ideal_mixture_density_from_props",
]
