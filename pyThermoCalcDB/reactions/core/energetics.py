"""Core reaction energetic identity calculations."""

# import libs
from collections.abc import Mapping
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import (
    NumericArrayInput,
    _return_scalar_if_zero_dim,
    _validate_positive_array,
    _validate_same_mapping_keys,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators

def _as_finite_float_array(
    values: NumericInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert numeric input to a finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_component_arrays(
    stoichiometric_coefficients: NDArray[np.float64],
    standard_entropies: NDArray[np.float64],
) -> None:
    """Validate pairwise component arrays for reaction entropy calculations."""
    if stoichiometric_coefficients.ndim not in (1, 2):
        raise ValueError("stoichiometric_coefficients must be a 1-D or 2-D array.")
    if standard_entropies.ndim not in (1, 2):
        raise ValueError("standard_entropies must be a 1-D or 2-D array.")
    if stoichiometric_coefficients.shape != standard_entropies.shape:
        raise ValueError("stoichiometric_coefficients and standard_entropies must have the same shape.")


# SECTION: Standard reaction entropy

def _calc_reaction_entropy_std(
    stoichiometric_coefficients: NumericInput,
    standard_entropies: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate delta_S_rxn_std = sum_i(nu_i*S_i_std).

    For 2-D inputs, axis 0 is states and axis 1 is species/components.
    """
    # SECTION: Normalize and validate
    nu = _as_finite_float_array(stoichiometric_coefficients, "stoichiometric_coefficients")
    entropy = _as_finite_float_array(standard_entropies, "standard_entropies")
    _validate_component_arrays(nu, entropy)

    # SECTION: Calculate standard reaction entropy
    return _return_scalar_if_zero_dim(np.sum(nu * entropy, axis=-1))


def _calc_reaction_entropy_std_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int],
    standard_entropies: Mapping[str, float | int],
) -> float:
    """Calculate standard reaction entropy from keyed species values."""
    _validate_same_mapping_keys(
        stoichiometric_coefficients,
        standard_entropies,
        "stoichiometric_coefficients",
        "standard_entropies",
    )
    return float(
        _calc_reaction_entropy_std(
            [stoichiometric_coefficients[key] for key in stoichiometric_coefficients],
            [standard_entropies[key] for key in stoichiometric_coefficients],
        )
    )


# SECTION: Entropy from enthalpy and Gibbs energy

def _calc_reaction_entropy_std_from_enthalpy_gibbs(
    delta_h_reaction_std: NumericInput,
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate delta_S_rxn_std = (delta_H_rxn_std - delta_G_rxn_std)/T."""
    # SECTION: Normalize and validate
    dh = _as_finite_float_array(delta_h_reaction_std, "delta_h_reaction_std")
    dg = _as_finite_float_array(delta_g_reaction_std, "delta_g_reaction_std")
    t = _as_finite_float_array(temperature, "temperature")
    _validate_positive_array(t, "temperature")

    # SECTION: Calculate reaction entropy
    return _return_scalar_if_zero_dim((dh - dg) / t)


# SECTION: Core exports
__all__ = [
    "_calc_reaction_entropy_std",
    "_calc_reaction_entropy_std_from_mapping",
    "_calc_reaction_entropy_std_from_enthalpy_gibbs",
]
