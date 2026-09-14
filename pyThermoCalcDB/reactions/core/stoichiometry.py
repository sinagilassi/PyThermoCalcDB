"""Core stoichiometric reaction-extent calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_finite_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _as_non_negative_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to finite non-negative float64 array."""
    arr = _as_finite_array(values, name)
    # ! Component mole amounts cannot be negative before applying reaction extents.
    if np.any(arr < 0.0):
        raise ValueError(f"{name} values must be non-negative.")
    return arr


# SECTION: Core numeric calculations
def _calc_component_moles_from_reaction_extent(
    initial_moles: NumericInput,
    stoichiometric_coefficients: NumericInput,
    reaction_extent: NumericInput,
) -> NDArray[np.float64]:
    """Calculate component moles for one reaction extent.

    Equation: ``N_j = N_j0 + nu_j*xi``. Initial moles are mol,
    stoichiometric coefficients are signed dimensionless coefficients, and
    reaction extent is mol. The output is component moles, mol. This exact
    stoichiometric relation preserves component order along the last axis.
    """
    n0 = _as_non_negative_array(initial_moles, "initial_moles")
    nu = _as_finite_array(stoichiometric_coefficients, "stoichiometric_coefficients")
    xi = _as_finite_array(reaction_extent, "reaction_extent")
    if n0.shape != nu.shape:
        raise ValueError("initial_moles and stoichiometric_coefficients must have the same shape.")
    result = n0 + nu * xi
    # ! Negative product means the extent over-consumed at least one component.
    if np.any(result < 0.0):
        raise ValueError("calculated component moles must be non-negative.")
    return cast(NDArray[np.float64], result)


def _calc_component_moles_from_reaction_extents(
    initial_moles: NumericInput,
    stoichiometric_matrix: NumericInput,
    reaction_extents: NumericInput,
) -> NDArray[np.float64]:
    """Calculate component moles for multiple reaction extents.

    Equation: ``N = N0 + xi @ nu`` where rows of ``stoichiometric_matrix`` are
    reactions and columns are components. Initial moles and output are mol;
    reaction extents are mol; stoichiometric coefficients are dimensionless.
    """
    n0 = _as_non_negative_array(initial_moles, "initial_moles")
    nu = _as_finite_array(stoichiometric_matrix, "stoichiometric_matrix")
    xi = _as_finite_array(reaction_extents, "reaction_extents")
    if n0.ndim != 1:
        raise ValueError("initial_moles must be one-dimensional for multiple reactions.")
    if nu.ndim != 2:
        raise ValueError("stoichiometric_matrix must be two-dimensional with reactions on axis 0.")
    if xi.ndim != 1:
        raise ValueError("reaction_extents must be one-dimensional.")
    # ? Matrix columns must align to the caller's component order.
    if nu.shape[1] != n0.shape[0] or nu.shape[0] != xi.shape[0]:
        raise ValueError("stoichiometric_matrix shape must be (n_reactions, n_components).")
    result = n0 + xi @ nu
    if np.any(result < 0.0):
        raise ValueError("calculated component moles must be non-negative.")
    return cast(NDArray[np.float64], result)


__all__ = [
    "_calc_component_moles_from_reaction_extent",
    "_calc_component_moles_from_reaction_extents",
]
