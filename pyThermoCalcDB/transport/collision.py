"""Public collision-property helpers for transport correlations."""

# import libs
import numpy as np
from numpy.typing import NDArray
# locals
from .core.collision import (
    _calc_lennard_jones_pair_diameter,
    _calc_lennard_jones_pair_energy,
    _calc_neufeld_diffusion_collision_integral,
    _calc_reduced_collision_temperature,
)


# SECTION: Public collision helpers

def calc_lennard_jones_pair_diameter(
    sigma_i,
    sigma_j,
) -> float | NDArray[np.float64]:
    """Calculate the Lennard-Jones pair diameter from two pure-species diameters."""
    return _calc_lennard_jones_pair_diameter(sigma_i, sigma_j)


def calc_lennard_jones_pair_energy(
    epsilon_over_k_i,
    epsilon_over_k_j,
) -> float | NDArray[np.float64]:
    """Calculate pair ``epsilon/k`` from pure-species values by geometric mean."""
    return _calc_lennard_jones_pair_energy(epsilon_over_k_i, epsilon_over_k_j)


def calc_reduced_collision_temperature(
    temperature,
    epsilon_over_k_ij,
) -> float | NDArray[np.float64]:
    """Calculate reduced collision temperature for a pair interaction."""
    return _calc_reduced_collision_temperature(temperature, epsilon_over_k_ij)


def calc_neufeld_diffusion_collision_integral(
    reduced_temperature,
) -> float | NDArray[np.float64]:
    """Calculate the Neufeld diffusion collision integral."""
    return _calc_neufeld_diffusion_collision_integral(reduced_temperature)


# SECTION: Public exports
__all__ = [
    "calc_lennard_jones_pair_diameter",
    "calc_lennard_jones_pair_energy",
    "calc_reduced_collision_temperature",
    "calc_neufeld_diffusion_collision_integral",
]
