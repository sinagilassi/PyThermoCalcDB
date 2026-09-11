"""Core mixture molecular-weight calculations."""

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


# SECTION: Mole-fraction basis

def _calc_mixture_molecular_weight_from_mole_fractions(
    mole_fractions: NumericInput,
    molecular_weights: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate M_mix = sum_i(x_i * M_i).

    For 2-D inputs, axis 0 is states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    x = _as_float_array(mole_fractions, "mole_fractions")
    mw = _as_float_array(molecular_weights, "molecular_weights")
    _validate_same_array_shape(x, mw, "mole_fractions", "molecular_weights")
    _validate_fraction_array(x, "mole_fractions")
    _validate_positive_array(mw, "molecular_weights")

    # SECTION: Calculate mole-fraction weighted molecular weight
    return _return_scalar_if_zero_dim(np.sum(x * mw, axis=-1))


# SECTION: Mass-fraction basis

def _calc_mixture_molecular_weight_from_mass_fractions(
    mass_fractions: NumericInput,
    molecular_weights: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate M_mix = 1 / sum_i(w_i / M_i).

    For 2-D inputs, axis 0 is states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    w = _as_float_array(mass_fractions, "mass_fractions")
    mw = _as_float_array(molecular_weights, "molecular_weights")
    _validate_same_array_shape(w, mw, "mass_fractions", "molecular_weights")
    _validate_fraction_array(w, "mass_fractions")
    _validate_positive_array(mw, "molecular_weights")

    # SECTION: Calculate reciprocal mass-fraction weighted molecular weight
    denominator = np.sum(w / mw, axis=-1)
    if np.any(denominator <= 0.0):
        raise ValueError("The reciprocal molecular-weight sum must be positive.")
    return _return_scalar_if_zero_dim(1.0 / denominator)


# SECTION: Mapping adapters

def _calc_mixture_molecular_weight_from_mole_fraction_mapping(
    mole_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
) -> float:
    """Calculate mole-fraction-basis mixture molecular weight from keyed inputs."""
    # SECTION: Align by caller mole-fraction order after key validation
    _validate_same_mapping_keys(mole_fractions, molecular_weights, "mole_fractions", "molecular_weights")
    keys = list(mole_fractions)
    return float(
        _calc_mixture_molecular_weight_from_mole_fractions(
            [mole_fractions[key] for key in keys],
            [molecular_weights[key] for key in keys],
        )
    )


def _calc_mixture_molecular_weight_from_mass_fraction_mapping(
    mass_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
) -> float:
    """Calculate mass-fraction-basis mixture molecular weight from keyed inputs."""
    # SECTION: Align by caller mass-fraction order after key validation
    _validate_same_mapping_keys(mass_fractions, molecular_weights, "mass_fractions", "molecular_weights")
    keys = list(mass_fractions)
    return float(
        _calc_mixture_molecular_weight_from_mass_fractions(
            [mass_fractions[key] for key in keys],
            [molecular_weights[key] for key in keys],
        )
    )


# SECTION: Props adapters

def _calc_mixture_molecular_weight_from_mole_fractions_from_props(
    mole_fractions: Mapping[str, CustomProp],
    molecular_weights: Mapping[str, CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate mole-fraction-basis mixture molecular weight from unit-aware keyed inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_mapping(mole_fractions, "mole_fractions")
    _validate_custom_prop_mapping(molecular_weights, "molecular_weights")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)
    mw = to_dict(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn=conversion_fn,
    )

    # SECTION: Remap component keys
    x = _configure_component_values(
        x,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "mole_fractions",
    )
    mw = _configure_component_values(
        mw,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molecular_weights",
    )

    # SECTION: Calculate from aligned mapping
    return _calc_mixture_molecular_weight_from_mole_fraction_mapping(x, mw)


def _calc_mixture_molecular_weight_from_mass_fractions_from_props(
    mass_fractions: Mapping[str, CustomProp],
    molecular_weights: Mapping[str, CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate mass-fraction-basis mixture molecular weight from unit-aware keyed inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_mapping(mass_fractions, "mass_fractions")
    _validate_custom_prop_mapping(molecular_weights, "molecular_weights")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    w = to_dict(mass_fractions, unit_conversion_fn=conversion_fn)
    mw = to_dict(
        molecular_weights,
        output_molecular_weight_unit,
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
    mw = _configure_component_values(
        mw,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molecular_weights",
    )

    # SECTION: Calculate from aligned mapping
    return _calc_mixture_molecular_weight_from_mass_fraction_mapping(w, mw)


# SECTION: Compatibility aliases

# NOTE: Preserve the initially-added basis-specific adapter names.
_calc_mixture_molecular_weight_from_mole_fraction_props = (
    _calc_mixture_molecular_weight_from_mole_fractions_from_props
)
_calc_mixture_molecular_weight_from_mass_fraction_props = (
    _calc_mixture_molecular_weight_from_mass_fractions_from_props
)


# SECTION: Core exports
__all__ = [
    "_calc_mixture_molecular_weight_from_mole_fractions",
    "_calc_mixture_molecular_weight_from_mass_fractions",
    "_calc_mixture_molecular_weight_from_mole_fraction_mapping",
    "_calc_mixture_molecular_weight_from_mass_fraction_mapping",
    "_calc_mixture_molecular_weight_from_mole_fractions_from_props",
    "_calc_mixture_molecular_weight_from_mass_fractions_from_props",
    "_calc_mixture_molecular_weight_from_mole_fraction_props",
    "_calc_mixture_molecular_weight_from_mass_fraction_props",
]
