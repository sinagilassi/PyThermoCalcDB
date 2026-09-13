"""Core collision-property helpers for transport correlations."""

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
    """Validate strictly positive collision-property values."""
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


# SECTION: Collision helpers

def _calc_lennard_jones_pair_diameter(
    sigma_i: NumericArrayInput,
    sigma_j: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate Lennard-Jones pair diameter: ``sigma_ij = (sigma_i + sigma_j)/2``."""
    sigma_i_arr = _as_state_array(sigma_i, "sigma_i")
    sigma_j_arr = _as_state_array(sigma_j, "sigma_j")
    _validate_positive(sigma_i_arr, "sigma_i")
    _validate_positive(sigma_j_arr, "sigma_j")
    return _return_scalar_if_zero_dim(0.5 * (sigma_i_arr + sigma_j_arr))


def _calc_lennard_jones_pair_energy(
    epsilon_over_k_i: NumericArrayInput,
    epsilon_over_k_j: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate Lennard-Jones pair energy by geometric mean."""
    eps_i = _as_state_array(epsilon_over_k_i, "epsilon_over_k_i")
    eps_j = _as_state_array(epsilon_over_k_j, "epsilon_over_k_j")
    _validate_positive(eps_i, "epsilon_over_k_i")
    _validate_positive(eps_j, "epsilon_over_k_j")
    return _return_scalar_if_zero_dim(np.sqrt(eps_i * eps_j))


def _calc_reduced_collision_temperature(
    temperature: NumericArrayInput,
    epsilon_over_k_ij: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate reduced collision temperature: ``T* = T/(epsilon_ij/k_B)``."""
    t = _as_state_array(temperature, "temperature")
    eps_ij = _as_state_array(epsilon_over_k_ij, "epsilon_over_k_ij")
    _validate_positive(t, "temperature")
    _validate_positive(eps_ij, "epsilon_over_k_ij")
    return _return_scalar_if_zero_dim(t / eps_ij)


def _calc_neufeld_diffusion_collision_integral(
    reduced_temperature: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate Neufeld diffusion collision integral for Lennard-Jones gases."""
    t_star = _as_state_array(reduced_temperature, "reduced_temperature")
    _validate_positive(t_star, "reduced_temperature")
    omega = (
        1.06036 / np.power(t_star, 0.15610)
        + 0.19300 * np.exp(-0.47635 * t_star)
        + 1.03587 * np.exp(-1.52996 * t_star)
        + 1.76474 * np.exp(-3.89411 * t_star)
    )
    return _return_scalar_if_zero_dim(omega)


# SECTION: Core exports
__all__ = [
    "_calc_lennard_jones_pair_diameter",
    "_calc_lennard_jones_pair_energy",
    "_calc_reduced_collision_temperature",
    "_calc_neufeld_diffusion_collision_integral",
]
