"""Core ideal-mixture enthalpy helpers."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim


def _calc_ideal_enthalpy_of_mixing(
    mole_fractions: NumericArrayInput | None = None,
) -> float | NDArray[np.float64]:
    """Return zero for ideal enthalpy of mixing."""
    if mole_fractions is None:
        return 0.0
    arr = np.asarray(mole_fractions, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError("mole_fractions must be one- or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError("mole_fractions values must be finite.")
    if np.any(arr < 0.0):
        raise ValueError("mole_fractions must be non-negative.")
    if arr.ndim == 0:
        return 0.0
    shape = arr.shape[:-1]
    return _return_scalar_if_zero_dim(np.zeros(shape, dtype=np.float64))


__all__ = ["_calc_ideal_enthalpy_of_mixing"]
