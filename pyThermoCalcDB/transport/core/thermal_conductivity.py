"""Core thermal-conductivity correlations and mixing rules."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
# locals
from ...configs.constants import R_J_molK
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


# SECTION: Thermal conductivity calculations

def _calc_stiel_thodos_gas_thermal_conductivity(
    viscosity: NumericArrayInput,
    molecular_weight: NumericArrayInput,
    molar_heat_capacity: NumericArrayInput,
    gas_constant: NumericArrayInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate Stiel-Thodos gas thermal conductivity.

    SI inputs produce W/(m*K): viscosity in Pa*s, molecular weight in kg/mol,
    heat capacity and gas constant in J/(mol*K).
    """
    mu = _as_state_array(viscosity, "viscosity")
    mw = _as_state_array(molecular_weight, "molecular_weight")
    cp = _as_state_array(molar_heat_capacity, "molar_heat_capacity")
    r = _as_state_array(gas_constant, "gas_constant")
    _validate_positive(mu, "viscosity")
    _validate_positive(mw, "molecular_weight")
    _validate_positive(cp, "molar_heat_capacity")
    _validate_positive(r, "gas_constant")
    return _return_scalar_if_zero_dim((mu / mw) * (1.15 * cp + 0.88 * r))


def _calc_ideal_liquid_thermal_conductivity(
    mole_fractions: NumericArrayInput,
    thermal_conductivities: NumericArrayInput,
) -> float | NDArray[np.float64]:
    """Calculate ideal liquid thermal-conductivity mixing rule."""
    x = _as_state_array(mole_fractions, "mole_fractions")
    k = _as_state_array(thermal_conductivities, "thermal_conductivities")
    _validate_fractions(x, "mole_fractions")
    _validate_positive(k, "thermal_conductivities")
    if x.shape != k.shape:
        raise ValueError("mole_fractions and thermal_conductivities must have the same shape.")
    return _return_scalar_if_zero_dim(np.sum(x * k, axis=-1))


# SECTION: Core exports
__all__ = [
    "_calc_stiel_thodos_gas_thermal_conductivity",
    "_calc_ideal_liquid_thermal_conductivity",
]
