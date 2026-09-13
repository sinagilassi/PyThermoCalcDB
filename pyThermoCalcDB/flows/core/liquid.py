"""Core liquid-flow approximations."""

# import libs
from collections.abc import Mapping
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _validate_same_mapping_keys
from ._common import (
    _as_flow_float_array,
    _return_scalar_if_zero_dim,
    _validate_non_negative,
    _validate_positive,
    _validate_same_shape,
)

# SECTION: Type aliases
NumericInput: TypeAlias = NumericArrayInput


def _calc_additive_liquid_volumetric_flow_rate(
    molar_flow_rates: NumericInput,
    molecular_weights: NumericInput,
    component_densities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate additive liquid volumetric flow approximation.

    Equation: ``Vdot ~= sum_i(F_i*M_i/rho_i)``.
    """
    f = _as_flow_float_array(molar_flow_rates, "molar_flow_rates")
    mw = _as_flow_float_array(molecular_weights, "molecular_weights")
    rho = _as_flow_float_array(component_densities, "component_densities")
    _validate_same_shape(f, mw, "molar_flow_rates", "molecular_weights")
    _validate_same_shape(f, rho, "molar_flow_rates", "component_densities")
    _validate_non_negative(f, "molar_flow_rates")
    _validate_positive(mw, "molecular_weights")
    _validate_positive(rho, "component_densities")
    return _return_scalar_if_zero_dim(np.sum(f * mw / rho, axis=-1))


def _calc_additive_liquid_volumetric_flow_rate_from_mapping(
    molar_flow_rates: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    component_densities: Mapping[str, float | int],
) -> float:
    """Calculate additive liquid volumetric flow from aligned keyed inputs."""
    _validate_same_mapping_keys(molar_flow_rates, molecular_weights, "molar_flow_rates", "molecular_weights")
    _validate_same_mapping_keys(molar_flow_rates, component_densities, "molar_flow_rates", "component_densities")
    keys = list(molar_flow_rates)
    return float(
        _calc_additive_liquid_volumetric_flow_rate(
            [molar_flow_rates[key] for key in keys],
            [molecular_weights[key] for key in keys],
            [component_densities[key] for key in keys],
        )
    )


__all__ = [
    "_calc_additive_liquid_volumetric_flow_rate",
    "_calc_additive_liquid_volumetric_flow_rate_from_mapping",
]
