"""Core vapor-liquid surface-tension mixing rules."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim


# SECTION: Numeric helpers

def _as_state_array(values: NumericArrayInput, name: str) -> NDArray[np.float64]:
    """Convert scalar or array-like input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_positive(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive values."""
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


def _validate_fractions(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative fractions that close along the last axis."""
    if np.any(values < 0.0):
        raise ValueError(f"{name} values must be non-negative.")
    if not np.allclose(np.sum(values, axis=-1), 1.0):
        raise ValueError(f"{name} must sum to 1.0 along the component axis.")


# SECTION: Surface-tension mixing rules

def _calc_ideal_vapor_liquid_surface_tension(
    liquid_mole_fractions: NumericArrayInput,
    pure_surface_tensions: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate ideal vapor-liquid surface tension: ``sigma = sum_i x_i*sigma_i``."""
    x = _as_state_array(liquid_mole_fractions, "liquid_mole_fractions")
    sigma = _as_state_array(pure_surface_tensions, "pure_surface_tensions")
    _validate_fractions(x, "liquid_mole_fractions")
    _validate_positive(sigma, "pure_surface_tensions")
    if x.shape != sigma.shape:
        raise ValueError("liquid_mole_fractions and pure_surface_tensions must have the same shape.")
    return _return_scalar_if_zero_dim(np.sum(x * sigma, axis=-1))


def _calc_winterfeld_vapor_liquid_surface_tension(
    liquid_mole_fractions: NumericArrayInput,
    pure_surface_tensions: NumericArrayInput,
    liquid_molar_densities: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate Winterfeld-Scriven-Davis vapor-liquid surface tension.

    ``liquid_molar_densities`` are pure-component liquid molar densities in a
    consistent unit, such as mol/m3. The rule is normalized by the mixture
    molar volume squared.
    """
    x = _as_state_array(liquid_mole_fractions, "liquid_mole_fractions")
    sigma = _as_state_array(pure_surface_tensions, "pure_surface_tensions")
    rho_m = _as_state_array(liquid_molar_densities, "liquid_molar_densities")
    _validate_fractions(x, "liquid_mole_fractions")
    _validate_positive(sigma, "pure_surface_tensions")
    _validate_positive(rho_m, "liquid_molar_densities")
    if x.shape != sigma.shape or x.shape != rho_m.shape:
        raise ValueError("liquid_mole_fractions, pure_surface_tensions, and liquid_molar_densities must have the same shape.")

    # NOTE: V_i = 1/rho_i and V_mix = sum_i x_i*V_i.
    molar_volumes = 1.0 / rho_m
    weighted_volumes = x * molar_volumes
    v_mix = np.sum(weighted_volumes, axis=-1)
    if np.any(v_mix <= 0.0):
        raise ValueError("mixture molar volume must be greater than zero.")
    pair_factor = np.sqrt(np.expand_dims(sigma, -1) * np.expand_dims(sigma, -2))
    weighted_pairs = np.expand_dims(weighted_volumes, -1) * np.expand_dims(weighted_volumes, -2)
    numerator = np.sum(weighted_pairs * pair_factor, axis=(-2, -1))
    return _return_scalar_if_zero_dim(numerator / np.power(v_mix, 2.0))


# SECTION: Core exports
__all__ = [
    "_calc_ideal_vapor_liquid_surface_tension",
    "_calc_winterfeld_vapor_liquid_surface_tension",
]
