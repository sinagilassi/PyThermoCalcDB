"""Core concentration and flow identities."""

# import libs
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput
from ._common import (
    _as_flow_float_array,
    _return_scalar_if_zero_dim,
    _validate_non_negative,
    _validate_positive,
)

# SECTION: Type aliases
NumericInput: TypeAlias = NumericArrayInput


def _calc_molar_flow_rate_from_concentration(
    concentration: NumericInput,
    volumetric_flow_rate: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate molar flow rate from concentration: ndot_i = C_i*Vdot."""
    c = _as_flow_float_array(concentration, "concentration")
    q = _as_flow_float_array(volumetric_flow_rate, "volumetric_flow_rate")
    _validate_non_negative(c, "concentration")
    _validate_non_negative(q, "volumetric_flow_rate")
    return _return_scalar_if_zero_dim(c * q)


def _calc_concentration_from_molar_flow_rate(
    molar_flow_rate: NumericInput,
    volumetric_flow_rate: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate concentration from molar flow rate: C_i = ndot_i/Vdot."""
    f = _as_flow_float_array(molar_flow_rate, "molar_flow_rate")
    q = _as_flow_float_array(volumetric_flow_rate, "volumetric_flow_rate")
    _validate_non_negative(f, "molar_flow_rate")
    _validate_positive(q, "volumetric_flow_rate")
    return _return_scalar_if_zero_dim(f / q)


__all__ = [
    "_calc_molar_flow_rate_from_concentration",
    "_calc_concentration_from_molar_flow_rate",
]
