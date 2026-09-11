"""Core reaction equilibrium calculations."""

# import libs
from collections.abc import Mapping
from typing import cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Temperature
from pythermodb_settings.models.units import UnitConversionFn

# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import (
    NumericArrayInput,
    _resolve_unit_conversion_fn,
    _return_scalar_if_zero_dim,
    _validate_positive_array,
    _validate_positive_scalar,
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
        raise ValueError(
            f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_component_arrays(
    stoichiometric_coefficients: NDArray[np.float64],
    activities: NDArray[np.float64],
) -> None:
    """Validate pairwise component arrays for reaction quotient calculations."""
    if stoichiometric_coefficients.ndim not in (1, 2):
        raise ValueError(
            "stoichiometric_coefficients must be a 1-D or 2-D array.")
    if activities.ndim not in (1, 2):
        raise ValueError("activities must be a 1-D or 2-D array.")
    if stoichiometric_coefficients.shape != activities.shape:
        raise ValueError(
            "stoichiometric_coefficients and activities must have the same shape.")
    _validate_positive_array(activities, "activities")


# SECTION: Temperature adapter

def _temperature_k(
    temperature: Temperature,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Return absolute temperature in K."""
    # SECTION: Normalize temperature
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    value = float(temperature.value)
    unit = temperature.unit.strip()
    if unit != "K":
        value = float(conversion_fn(
            value=value, from_unit=unit, to_unit="K"))

    # ! Log/equilibrium thermodynamic identities require T > 0 K.
    return _validate_positive_scalar(value, "temperature")


# SECTION: Equilibrium constant

def _calc_log_equilibrium_constant(
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate ln(K) = -delta_G_rxn_std/(R*T)."""
    # SECTION: Normalize and validate
    dg = _as_finite_float_array(delta_g_reaction_std, "delta_g_reaction_std")
    t = _as_finite_float_array(temperature, "temperature")
    _validate_positive_array(t, "temperature")
    r = _validate_positive_scalar(gas_constant, "gas_constant")

    # SECTION: Calculate logarithmic equilibrium constant
    return _return_scalar_if_zero_dim(-dg / (r * t))


def _calc_equilibrium_constant(
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate K = exp(-delta_G_rxn_std/(R*T))."""
    return _return_scalar_if_zero_dim(
        np.exp(
            np.asarray(
                _calc_log_equilibrium_constant(
                    delta_g_reaction_std,
                    temperature,
                    gas_constant,
                ),
                dtype=np.float64,
            )
        )
    )


# SECTION: Reaction quotient

def _calc_log_reaction_quotient(
    stoichiometric_coefficients: NumericInput,
    activities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ln(Q) = sum_i(nu_i*ln(a_i)).

    For 2-D inputs, axis 0 is states and axis 1 is species/components.
    """
    # SECTION: Normalize and validate
    nu = _as_finite_float_array(
        stoichiometric_coefficients, "stoichiometric_coefficients")
    a = _as_finite_float_array(activities, "activities")
    _validate_component_arrays(nu, a)

    # SECTION: Calculate logarithmic reaction quotient
    return _return_scalar_if_zero_dim(np.sum(nu * np.log(a), axis=-1))


def _calc_reaction_quotient(
    stoichiometric_coefficients: NumericInput,
    activities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Q = product_i(a_i**nu_i)."""
    return _return_scalar_if_zero_dim(
        np.exp(
            np.asarray(
                _calc_log_reaction_quotient(
                    stoichiometric_coefficients,
                    activities,
                ),
                dtype=np.float64,
            )
        )
    )


def _calc_log_reaction_quotient_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int],
    activities: Mapping[str, float | int],
) -> float:
    """Calculate ln(Q) from keyed stoichiometric coefficients and activities."""
    _validate_same_mapping_keys(
        stoichiometric_coefficients,
        activities,
        "stoichiometric_coefficients",
        "activities",
    )
    return float(
        _calc_log_reaction_quotient(
            [stoichiometric_coefficients[key]
                for key in stoichiometric_coefficients],
            [activities[key] for key in stoichiometric_coefficients],
        )
    )


def _calc_reaction_quotient_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int],
    activities: Mapping[str, float | int],
) -> float:
    """Calculate Q from keyed stoichiometric coefficients and activities."""
    return float(
        np.exp(
            _calc_log_reaction_quotient_from_mapping(
                stoichiometric_coefficients,
                activities,
            )
        )
    )


# SECTION: Actual reaction Gibbs energy

def _calc_reaction_gibbs_energy(
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
    log_reaction_quotient: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate delta_G_rxn = delta_G_rxn_std + R*T*ln(Q)."""
    # SECTION: Normalize and validate
    dg_std = _as_finite_float_array(
        delta_g_reaction_std, "delta_g_reaction_std")
    t = _as_finite_float_array(temperature, "temperature")
    ln_q = _as_finite_float_array(
        log_reaction_quotient, "log_reaction_quotient")
    _validate_positive_array(t, "temperature")
    r = _validate_positive_scalar(gas_constant, "gas_constant")

    # SECTION: Calculate actual reaction Gibbs energy
    return _return_scalar_if_zero_dim(dg_std + r * t * ln_q)


# SECTION: van't Hoff relations

def _calc_dlnK_dT(
    delta_h_reaction_std: NumericInput,
    temperature: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate dlnK/dT = delta_H_rxn_std/(R*T**2)."""
    # SECTION: Normalize and validate
    dh = _as_finite_float_array(delta_h_reaction_std, "delta_h_reaction_std")
    t = _as_finite_float_array(temperature, "temperature")
    _validate_positive_array(t, "temperature")
    r = _validate_positive_scalar(gas_constant, "gas_constant")

    # SECTION: Calculate derivative
    return _return_scalar_if_zero_dim(dh / (r * t ** 2))


def _calc_log_equilibrium_constant_at_temperature(
    equilibrium_constant_initial: float | int,
    delta_h_reaction_std: NumericInput,
    temperature_initial: NumericInput,
    temperature_final: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate ln(K2) from the integrated van't Hoff relation."""
    # SECTION: Normalize and validate
    k_initial = _validate_positive_scalar(
        equilibrium_constant_initial,
        "equilibrium_constant_initial",
    )
    dh = _as_finite_float_array(delta_h_reaction_std, "delta_h_reaction_std")
    t_initial = _as_finite_float_array(
        temperature_initial, "temperature_initial")
    t_final = _as_finite_float_array(temperature_final, "temperature_final")
    _validate_positive_array(t_initial, "temperature_initial")
    _validate_positive_array(t_final, "temperature_final")
    r = _validate_positive_scalar(gas_constant, "gas_constant")

    # SECTION: Calculate final logarithmic equilibrium constant
    return _return_scalar_if_zero_dim(
        np.log(k_initial) - (dh / r) * (1.0 / t_final - 1.0 / t_initial)
    )


def _calc_equilibrium_constant_at_temperature(
    equilibrium_constant_initial: float | int,
    delta_h_reaction_std: NumericInput,
    temperature_initial: NumericInput,
    temperature_final: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate K2 from the integrated van't Hoff relation."""
    return _return_scalar_if_zero_dim(
        np.exp(
            np.asarray(
                _calc_log_equilibrium_constant_at_temperature(
                    equilibrium_constant_initial,
                    delta_h_reaction_std,
                    temperature_initial,
                    temperature_final,
                    gas_constant,
                ),
                dtype=np.float64,
            )
        )
    )


# SECTION: Core exports
__all__ = [
    "_temperature_k",
    "_calc_log_equilibrium_constant",
    "_calc_equilibrium_constant",
    "_calc_log_reaction_quotient",
    "_calc_reaction_quotient",
    "_calc_log_reaction_quotient_from_mapping",
    "_calc_reaction_quotient_from_mapping",
    "_calc_reaction_gibbs_energy",
    "_calc_dlnK_dT",
    "_calc_log_equilibrium_constant_at_temperature",
    "_calc_equilibrium_constant_at_temperature",
]
