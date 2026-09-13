"""Shared core helpers for flow calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim


def _as_flow_float_array(
    values: NumericArrayInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert scalar or array-like numeric input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_positive(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive array values."""
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


def _validate_non_negative(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative array values."""
    if np.any(values < 0.0):
        raise ValueError(f"{name} values must be non-negative.")


def _validate_same_shape(left: NDArray[np.float64], right: NDArray[np.float64], left_name: str, right_name: str) -> None:
    """Validate identical array shapes."""
    if left.shape != right.shape:
        raise ValueError(f"{left_name} and {right_name} must have the same shape.")


__all__ = [
    "_as_flow_float_array",
    "_return_scalar_if_zero_dim",
    "_validate_non_negative",
    "_validate_positive",
    "_validate_same_shape",
]
