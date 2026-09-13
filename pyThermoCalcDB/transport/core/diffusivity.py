"""Core gas and liquid diffusivity correlations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim
from .collision import (
    _calc_lennard_jones_pair_diameter,
    _calc_lennard_jones_pair_energy,
    _calc_neufeld_diffusion_collision_integral,
    _calc_reduced_collision_temperature,
)


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
    """Validate strictly positive diffusivity-correlation inputs."""
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


# SECTION: Gas diffusivity correlations

def _calc_fuller_schettler_giddings_diffusivity(
    temperature: NumericArrayInput,
    pressure: NumericArrayInput,
    molecular_weight_i: NumericArrayInput,
    molecular_weight_j: NumericArrayInput,
    diffusion_volume_i: NumericArrayInput,
    diffusion_volume_j: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate Fuller-Schettler-Giddings binary gas diffusivity in ``m2/s``.

    Inputs use SI for state variables and molecular weights: ``T`` in K,
    ``P`` in Pa, molecular weights in kg/mol. Fuller diffusion volumes are the
    tabulated dimensionless group increments used by the correlation.
    """
    t = _as_state_array(temperature, "temperature")
    p = _as_state_array(pressure, "pressure")
    mw_i = _as_state_array(molecular_weight_i, "molecular_weight_i")
    mw_j = _as_state_array(molecular_weight_j, "molecular_weight_j")
    v_i = _as_state_array(diffusion_volume_i, "diffusion_volume_i")
    v_j = _as_state_array(diffusion_volume_j, "diffusion_volume_j")
    for name, values in (
        ("temperature", t),
        ("pressure", p),
        ("molecular_weight_i", mw_i),
        ("molecular_weight_j", mw_j),
        ("diffusion_volume_i", v_i),
        ("diffusion_volume_j", v_j),
    ):
        _validate_positive(values, name)
    # NOTE: The empirical coefficient is in cm2/s, atm, g/mol convention.
    p_atm = p / 101325.0
    mw_i_g_mol = mw_i * 1000.0
    mw_j_g_mol = mw_j * 1000.0
    denominator = p_atm * np.power(np.cbrt(v_i) + np.cbrt(v_j), 2.0)
    d_cm2_s = 0.00143 * np.power(t, 1.75) * np.sqrt(1.0 / mw_i_g_mol + 1.0 / mw_j_g_mol) / denominator
    return _return_scalar_if_zero_dim(d_cm2_s * 1.0e-4)


def _calc_wilke_lee_diffusivity(
    temperature: NumericArrayInput,
    pressure: NumericArrayInput,
    molecular_weight_i: NumericArrayInput,
    molecular_weight_j: NumericArrayInput,
    sigma_i: NumericArrayInput,
    sigma_j: NumericArrayInput,
    epsilon_over_k_i: NumericArrayInput,
    epsilon_over_k_j: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate Wilke-Lee gas binary diffusivity in ``m2/s``.

    State variables are SI, molecular weights are kg/mol, Lennard-Jones
    diameters are Angstrom, and epsilon values are expressed as ``epsilon/k`` K.
    """
    t = _as_state_array(temperature, "temperature")
    p = _as_state_array(pressure, "pressure")
    mw_i = _as_state_array(molecular_weight_i, "molecular_weight_i")
    mw_j = _as_state_array(molecular_weight_j, "molecular_weight_j")
    _validate_positive(t, "temperature")
    _validate_positive(p, "pressure")
    _validate_positive(mw_i, "molecular_weight_i")
    _validate_positive(mw_j, "molecular_weight_j")
    sigma_ij = np.asarray(_calc_lennard_jones_pair_diameter(sigma_i, sigma_j), dtype=np.float64)
    eps_ij = np.asarray(_calc_lennard_jones_pair_energy(epsilon_over_k_i, epsilon_over_k_j), dtype=np.float64)
    t_star = np.asarray(_calc_reduced_collision_temperature(t, eps_ij), dtype=np.float64)
    omega = np.asarray(_calc_neufeld_diffusion_collision_integral(t_star), dtype=np.float64)
    p_atm = p / 101325.0
    mw_i_g_mol = mw_i * 1000.0
    mw_j_g_mol = mw_j * 1000.0
    m_ab = 2.0 / (1.0 / mw_i_g_mol + 1.0 / mw_j_g_mol)
    coefficient = (3.03 - 0.98 / np.sqrt(m_ab)) * 1.0e-3
    d_cm2_s = coefficient * np.power(t, 1.5) / (p_atm * np.power(sigma_ij, 2.0) * omega * np.sqrt(m_ab))
    return _return_scalar_if_zero_dim(d_cm2_s * 1.0e-4)


# SECTION: Liquid diffusivity correlations

def _calc_wilke_chang_diffusivity(
    temperature: NumericArrayInput,
    solvent_viscosity: NumericArrayInput,
    solvent_molecular_weight: NumericArrayInput,
    solute_molar_volume_at_normal_boiling_point: NumericArrayInput,
    solvent_association_factor: NumericArrayInput = 1.0,
) -> float | NDArray[np.float64]:
    """Calculate Wilke-Chang liquid diffusivity in ``m2/s``.

    Inputs are SI: K, Pa*s, kg/mol, and m3/mol. The empirical correlation is
    evaluated in cP, g/mol, and cm3/mol internally.
    """
    t = _as_state_array(temperature, "temperature")
    mu = _as_state_array(solvent_viscosity, "solvent_viscosity")
    mw = _as_state_array(solvent_molecular_weight, "solvent_molecular_weight")
    v_a = _as_state_array(solute_molar_volume_at_normal_boiling_point, "solute_molar_volume_at_normal_boiling_point")
    phi = _as_state_array(solvent_association_factor, "solvent_association_factor")
    for name, values in (
        ("temperature", t),
        ("solvent_viscosity", mu),
        ("solvent_molecular_weight", mw),
        ("solute_molar_volume_at_normal_boiling_point", v_a),
        ("solvent_association_factor", phi),
    ):
        _validate_positive(values, name)
    mu_cp = mu * 1000.0
    mw_g_mol = mw * 1000.0
    v_cm3_mol = v_a * 1.0e6
    d_cm2_s = 7.4e-8 * np.sqrt(phi * mw_g_mol) * t / (mu_cp * np.power(v_cm3_mol, 0.6))
    return _return_scalar_if_zero_dim(d_cm2_s * 1.0e-4)


# SECTION: Core exports
__all__ = [
    "_calc_fuller_schettler_giddings_diffusivity",
    "_calc_wilke_lee_diffusivity",
    "_calc_wilke_chang_diffusivity",
]
