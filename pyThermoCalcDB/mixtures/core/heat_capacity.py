"""Core ideal mixture heat-capacity calculations."""

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


# SECTION: Ideal heat-capacity mixing rule

def _calc_ideal_mixture_heat_capacity(
    fraction_values: NumericInput,
    heat_capacities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Cp_mix = sum_i(f_i * Cp_i).

    For 2-D inputs, axis 0 is states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    fr = _as_float_array(fraction_values, "fraction_values")
    cp = _as_float_array(heat_capacities, "heat_capacities")
    _validate_same_array_shape(fr, cp, "fraction_values", "heat_capacities")
    _validate_fraction_array(fr, "fraction_values")
    _validate_positive_array(cp, "heat_capacities")

    # SECTION: Calculate weighted-average heat capacity
    return _return_scalar_if_zero_dim(np.sum(fr * cp, axis=-1))


# SECTION: Mapping adapter

def _calc_ideal_mixture_heat_capacity_from_mapping(
    fraction_values: Mapping[str, float | int],
    heat_capacities: Mapping[str, float | int],
) -> float:
    """Calculate ideal mixture heat capacity from aligned keyed inputs."""
    # SECTION: Align by caller fraction order after key validation
    _validate_same_mapping_keys(fraction_values, heat_capacities, "fraction_values", "heat_capacities")
    keys = list(fraction_values)
    return float(
        _calc_ideal_mixture_heat_capacity(
            [fraction_values[key] for key in keys],
            [heat_capacities[key] for key in keys],
        )
    )


# SECTION: Props adapter

def _calc_ideal_mixture_heat_capacity_from_props(
    fraction_values: Mapping[str, CustomProp],
    heat_capacities: Mapping[str, CustomProp],
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal mixture heat capacity from unit-aware keyed inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_mapping(fraction_values, "fraction_values")
    _validate_custom_prop_mapping(heat_capacities, "heat_capacities")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    fr = to_dict(fraction_values, unit_conversion_fn=conversion_fn)
    cp = to_dict(
        heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn=conversion_fn,
    )

    # SECTION: Remap component keys
    fr = _configure_component_values(
        fr,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "fractions",
    )
    cp = _configure_component_values(
        cp,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "heat_capacities",
    )

    # SECTION: Calculate from aligned mapping
    return _calc_ideal_mixture_heat_capacity_from_mapping(fr, cp)


# SECTION: Core exports
__all__ = [
    "_calc_ideal_mixture_heat_capacity",
    "_calc_ideal_mixture_heat_capacity_from_mapping",
    "_calc_ideal_mixture_heat_capacity_from_props",
]
