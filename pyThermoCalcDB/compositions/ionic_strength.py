"""Shared numeric ionic-strength core."""

# import libs
from collections.abc import Sequence
from typing import cast

import numpy as np
from numpy.typing import NDArray


# =====================================================================
# *** Helper functions
# =====================================================================

def _validate_concentrations_and_charges(
    concentrations: NDArray[np.float64],
    charges: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Validate and broadcast concentration and charge arrays."""
    # SECTION: Validate array dimensions
    if concentrations.ndim > 2:
        raise ValueError(
            "concentrations must be a scalar, 1-D array, or 2-D array."
        )

    if charges.ndim > 2:
        raise ValueError("charges must be a scalar, 1-D array, or 2-D array.")

    # SECTION: Validate finite numeric values
    if not np.all(np.isfinite(concentrations)):
        raise ValueError("concentrations must contain finite values.")

    if not np.all(np.isfinite(charges)):
        raise ValueError("charges must contain finite values.")

    # ! Concentrations cannot be negative for ionic strength.
    if np.any(concentrations < 0):
        raise ValueError("concentrations must be non-negative.")

    # ? NumPy handles scalar, 1-D, and 2-D broadcasting consistently here.
    try:
        concentrations, charges = np.broadcast_arrays(concentrations, charges)
    except ValueError as exc:
        raise ValueError(
            "concentrations and charges must have broadcast-compatible shapes."
        ) from exc

    return (
        cast(NDArray[np.float64], concentrations),
        cast(NDArray[np.float64], charges),
    )


# ======================================================================
# *** Core calculation
# ======================================================================

def _calc_ionic_strength(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate ionic strength: ``I = 0.5 * sum_i(c_i * z_i**2)``."""
    # SECTION: Convert inputs to numeric arrays
    concentration_values: NDArray[np.float64] = np.asarray(
        concentrations,
        dtype=np.float64,
    )
    charge_values: NDArray[np.float64] = np.asarray(
        charges,
        dtype=np.float64,
    )

    # SECTION: Validate inputs
    concentration_values, charge_values = _validate_concentrations_and_charges(
        concentration_values,
        charge_values,
    )

    # SECTION: Calculate component ionic-strength contributions
    values = 0.5 * concentration_values * charge_values**2

    # NOTE: Scalar and 1-D inputs return a scalar; 2-D inputs return row sums.
    if values.ndim == 0:
        return float(values)

    if values.ndim == 1:
        return float(np.sum(values))

    return cast(NDArray[np.float64], np.sum(values, axis=1))


# SECTION: Public exports
__all__ = ["_calc_ionic_strength"]
