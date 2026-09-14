"""Core thermodynamic departure-property relations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_finite_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric scalar, sequence, or array input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _as_positive_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to finite positive float64 array."""
    arr = _as_finite_array(values, name)
    # ! Temperatures and gas constants are positive thermodynamic scales.
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return arr


# SECTION: Core numeric calculations
def _calc_enthalpy_departure(
    enthalpy: NumericInput,
    ideal_gas_enthalpy: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate enthalpy departure ``H^R = H - H^id``.

    Enthalpy and ideal-gas enthalpy must share an energy basis, commonly J/mol.
    The output uses the same basis. This is a fundamental departure definition.
    """
    h = _as_finite_array(enthalpy, "enthalpy")
    h_id = _as_finite_array(ideal_gas_enthalpy, "ideal_gas_enthalpy")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], h - h_id))


def _calc_entropy_departure(
    entropy: NumericInput,
    ideal_gas_entropy: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate entropy departure ``S^R = S - S^id``.

    Entropy and ideal-gas entropy must share a basis, commonly J/(mol.K). The
    output uses the same basis. This is a fundamental departure definition.
    """
    s = _as_finite_array(entropy, "entropy")
    s_id = _as_finite_array(ideal_gas_entropy, "ideal_gas_entropy")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], s - s_id))


def _calc_dimensionless_enthalpy_departure(
    enthalpy_departure: NumericInput,
    temperature: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate dimensionless enthalpy departure ``H^R/(R*T)``."""
    h_dep = _as_finite_array(enthalpy_departure, "enthalpy_departure")
    t = _as_positive_array(temperature, "temperature")
    r = _as_positive_array(gas_constant, "gas_constant")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], h_dep / (r * t)))


def _calc_dimensionless_entropy_departure(
    entropy_departure: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate dimensionless entropy departure ``S^R/R``."""
    s_dep = _as_finite_array(entropy_departure, "entropy_departure")
    r = _as_positive_array(gas_constant, "gas_constant")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], s_dep / r))


def _calc_enthalpy_from_ideal_and_departure(
    ideal_gas_enthalpy: NumericInput,
    enthalpy_departure: NumericInput,
) -> float | NDArray[np.float64]:
    """Recover real-fluid enthalpy from ideal-gas and departure terms."""
    h_id = _as_finite_array(ideal_gas_enthalpy, "ideal_gas_enthalpy")
    h_dep = _as_finite_array(enthalpy_departure, "enthalpy_departure")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], h_id + h_dep))


def _calc_entropy_from_ideal_and_departure(
    ideal_gas_entropy: NumericInput,
    entropy_departure: NumericInput,
) -> float | NDArray[np.float64]:
    """Recover real-fluid entropy from ideal-gas and departure terms."""
    s_id = _as_finite_array(ideal_gas_entropy, "ideal_gas_entropy")
    s_dep = _as_finite_array(entropy_departure, "entropy_departure")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], s_id + s_dep))


def _calc_cp_departure_from_eos_derivatives(
    temperature: NumericInput,
    integral_d2p_dt2_dv: NumericInput,
    dpressure_dtemperature_at_volume: NumericInput,
    dpressure_dvolume_at_temperature: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate generic heat-capacity departure from EOS derivative data.

    Equation: ``Cp^R = T*I - T*(dP/dT)_V**2/(dP/dV)_T - R`` where ``I`` is
    the supplied volume integral of ``(d2P/dT2)_V``. EOS-specific derivatives
    and integrals are supplied by the caller/model layer.
    """
    t = _as_positive_array(temperature, "temperature")
    integral = _as_finite_array(integral_d2p_dt2_dv, "integral_d2p_dt2_dv")
    dpdt = _as_finite_array(dpressure_dtemperature_at_volume, "dpressure_dtemperature_at_volume")
    dpdv = _as_finite_array(dpressure_dvolume_at_temperature, "dpressure_dvolume_at_temperature")
    r = _as_positive_array(gas_constant, "gas_constant")
    if np.any(dpdv == 0.0):
        raise ValueError("dpressure_dvolume_at_temperature must not be zero.")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], t * integral - t * dpdt**2 / dpdv - r))


__all__ = [
    "_calc_enthalpy_departure",
    "_calc_entropy_departure",
    "_calc_dimensionless_enthalpy_departure",
    "_calc_dimensionless_entropy_departure",
    "_calc_enthalpy_from_ideal_and_departure",
    "_calc_entropy_from_ideal_and_departure",
    "_calc_cp_departure_from_eos_derivatives",
]
