"""Core mixture volume calculations."""

# import libs
from collections.abc import Mapping

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import (
    NumericArrayInput,
    _as_float_array,
    _return_scalar_if_zero_dim,
    _validate_non_negative_array,
    _validate_positive_array,
    _validate_same_array_shape,
    _validate_same_mapping_keys,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Additive liquid volume approximation

def _calc_additive_liquid_volume(
    component_moles: NumericInput,
    molecular_weights: NumericInput,
    component_densities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate additive liquid volume approximation.

    Equation: ``V_mix ~= sum_i(n_i*M_i/rho_i)``. For 2-D inputs, axis 0 is
    states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    n = _as_float_array(component_moles, "component_moles")
    mw = _as_float_array(molecular_weights, "molecular_weights")
    rho = _as_float_array(component_densities, "component_densities")
    _validate_same_array_shape(n, mw, "component_moles", "molecular_weights")
    _validate_same_array_shape(n, rho, "component_moles", "component_densities")
    _validate_non_negative_array(n, "component_moles")
    _validate_positive_array(mw, "molecular_weights")
    _validate_positive_array(rho, "component_densities")

    # SECTION: Calculate additive volume
    return _return_scalar_if_zero_dim(np.sum(n * mw / rho, axis=-1))


def _calc_additive_liquid_volume_from_mapping(
    component_moles: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    component_densities: Mapping[str, float | int],
) -> float:
    """Calculate additive liquid volume from aligned keyed inputs."""
    # SECTION: Align by caller mole order after key validation
    _validate_same_mapping_keys(component_moles, molecular_weights, "component_moles", "molecular_weights")
    _validate_same_mapping_keys(component_moles, component_densities, "component_moles", "component_densities")
    keys = list(component_moles)
    return float(
        _calc_additive_liquid_volume(
            [component_moles[key] for key in keys],
            [molecular_weights[key] for key in keys],
            [component_densities[key] for key in keys],
        )
    )


# SECTION: Core exports
__all__ = [
    "_calc_additive_liquid_volume",
    "_calc_additive_liquid_volume_from_mapping",
]
