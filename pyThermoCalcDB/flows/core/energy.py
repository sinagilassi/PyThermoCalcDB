"""Core stream energy flow calculations."""

# import libs
from collections.abc import Mapping
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _validate_same_mapping_keys
from ._common import (
    _as_flow_float_array,
    _return_scalar_if_zero_dim,
    _validate_non_negative,
    _validate_positive,
    _validate_same_shape,
)

# SECTION: Type aliases
NumericInput: TypeAlias = NumericArrayInput


def _calc_flowing_heat_capacity(
    molar_flow_rates: NumericInput,
    molar_heat_capacities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate flowing heat capacity: Cpdot = sum_i(F_i*Cp_i)."""
    f = _as_flow_float_array(molar_flow_rates, "molar_flow_rates")
    cp = _as_flow_float_array(molar_heat_capacities, "molar_heat_capacities")
    _validate_same_shape(f, cp, "molar_flow_rates", "molar_heat_capacities")
    _validate_non_negative(f, "molar_flow_rates")
    _validate_positive(cp, "molar_heat_capacities")
    return _return_scalar_if_zero_dim(np.sum(f * cp, axis=-1))


def _calc_enthalpy_flow_rate(
    molar_flow_rates: NumericInput,
    molar_enthalpies: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate enthalpy flow rate: Hdot = sum_i(F_i*h_i)."""
    f = _as_flow_float_array(molar_flow_rates, "molar_flow_rates")
    h = _as_flow_float_array(molar_enthalpies, "molar_enthalpies")
    _validate_same_shape(f, h, "molar_flow_rates", "molar_enthalpies")
    _validate_non_negative(f, "molar_flow_rates")
    return _return_scalar_if_zero_dim(np.sum(f * h, axis=-1))


def _calc_flowing_heat_capacity_from_mapping(
    molar_flow_rates: Mapping[str, float | int],
    molar_heat_capacities: Mapping[str, float | int],
) -> float:
    """Calculate flowing heat capacity from aligned keyed inputs."""
    _validate_same_mapping_keys(molar_flow_rates, molar_heat_capacities, "molar_flow_rates", "molar_heat_capacities")
    keys = list(molar_flow_rates)
    return float(
        _calc_flowing_heat_capacity(
            [molar_flow_rates[key] for key in keys],
            [molar_heat_capacities[key] for key in keys],
        )
    )


def _calc_enthalpy_flow_rate_from_mapping(
    molar_flow_rates: Mapping[str, float | int],
    molar_enthalpies: Mapping[str, float | int],
) -> float:
    """Calculate enthalpy flow rate from aligned keyed inputs."""
    _validate_same_mapping_keys(molar_flow_rates, molar_enthalpies, "molar_flow_rates", "molar_enthalpies")
    keys = list(molar_flow_rates)
    return float(
        _calc_enthalpy_flow_rate(
            [molar_flow_rates[key] for key in keys],
            [molar_enthalpies[key] for key in keys],
        )
    )


__all__ = [
    "_calc_flowing_heat_capacity",
    "_calc_enthalpy_flow_rate",
    "_calc_flowing_heat_capacity_from_mapping",
    "_calc_enthalpy_flow_rate_from_mapping",
]
