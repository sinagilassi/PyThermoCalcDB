"""Core reaction equilibrium calculations."""

# import libs
from collections.abc import Mapping
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import (
    NumericArrayInput,
    _return_scalar_if_zero_dim,
    _validate_positive_array,
    _validate_positive_scalar,
    _validate_same_mapping_keys,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
# ? Check if input is finite and has appropriate dimensions

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

# ? Validate component arrays for reaction quotient calculations


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


# ! ::: Equilibrium constant

def _calc_log_equilibrium_constant(
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """
    Calculate the natural logarithm of the thermodynamic equilibrium constant
    from the standard reaction Gibbs energy.

    The equilibrium constant is related to the standard reaction Gibbs energy by

        ln(K) = -ΔG°_rxn(T) / (R T)

    where K is the dimensionless thermodynamic equilibrium constant,
    ΔG°_rxn(T) is the standard reaction Gibbs energy at temperature T,
    R is the gas constant, and T is the absolute temperature.

    Parameters
    ----------
    delta_g_reaction_std : NumericInput
        Standard reaction Gibbs energy, ΔG°_rxn(T).

        The energy unit must be consistent with `gas_constant`. For example,
        when `gas_constant` is expressed in J/(mol·K), ΔG°_rxn must be
        expressed in J/mol.

    temperature : NumericInput
        Absolute temperature, T, in kelvin. All values must be greater
        than zero.

    gas_constant : float | int, optional
        Universal gas constant, R. The default is `R_J_molK`, expressed
        in J/(mol·K).

    Returns
    -------
    float | NDArray[np.float64]
        Natural logarithm of the dimensionless equilibrium constant, ln(K).

        A scalar is returned for scalar inputs; otherwise a NumPy array
        is returned according to NumPy broadcasting rules.

    Raises
    ------
    ValueError
        If any input contains non-finite values, if temperature is not
        positive, or if `gas_constant` is not positive.

    Notes
    -----
    The ``_std`` suffix denotes a standard-state thermodynamic property
    and does not imply a temperature of 298.15 K.

    The relation is valid at any temperature T provided that
    ΔG°_rxn(T) corresponds to that same temperature.

    The thermodynamic equilibrium constant K is dimensionless because
    it is defined in terms of dimensionless activities.
    """
    # SECTION: Normalize and validate
    dg = _as_finite_float_array(delta_g_reaction_std, "delta_g_reaction_std")
    t = _as_finite_float_array(temperature, "temperature")
    _validate_positive_array(t, "temperature")
    r = _validate_positive_scalar(gas_constant, "gas_constant")

    # SECTION: Calculate logarithmic equilibrium constant
    return _return_scalar_if_zero_dim(-dg / (r * t))

# ! ::: Equilibrium constant (wrapper)


def _calc_equilibrium_constant(
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """
    Calculate the thermodynamic equilibrium constant from the standard
    reaction Gibbs energy.

    The equilibrium constant is calculated as

        K = exp[-ΔG°_rxn(T) / (R T)]

    where K is the dimensionless thermodynamic equilibrium constant,
    ΔG°_rxn(T) is the standard reaction Gibbs energy at temperature T,
    R is the gas constant, and T is the absolute temperature.

    Parameters
    ----------
    delta_g_reaction_std : NumericInput
        Standard reaction Gibbs energy, ΔG°_rxn(T). Its energy unit must
        be consistent with `gas_constant`.

    temperature : NumericInput
        Absolute temperature, T, in kelvin. Must be greater than zero.

    gas_constant : float | int, optional
        Universal gas constant, R. Default is `R_J_molK`.

    Returns
    -------
    float | NDArray[np.float64]
        Dimensionless thermodynamic equilibrium constant, K.

    Notes
    -----
    This function evaluates `_calc_log_equilibrium_constant` and
    exponentiates the result.

    The ``_std`` suffix refers to the standard state and does not imply
    a temperature of 298.15 K.
    """
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


# ! ::: Reaction quotient

def _calc_log_reaction_quotient(
    stoichiometric_coefficients: NumericInput,
    activities: NumericInput,
) -> float | NDArray[np.float64]:
    """
    Calculate the natural logarithm of the reaction quotient from species
    activities and signed stoichiometric coefficients.

    The logarithmic reaction quotient is defined as

        ln(Q) = Σ_i ν_i ln(a_i)

    which is equivalent to

        Q = Π_i a_i**ν_i

    where ν_i is the signed stoichiometric coefficient and a_i is the
    dimensionless activity of species i.

    The signed stoichiometric convention is assumed:

        ν_i < 0  for reactants
        ν_i > 0  for products

    Parameters
    ----------
    stoichiometric_coefficients : NumericInput
        Signed stoichiometric coefficients, ν_i.

        Supported shapes are:

        - 1-D: ``(n_species,)``
        - 2-D: ``(n_states, n_species)``

        For 2-D inputs, each row represents one reaction state and each
        column represents one species.

    activities : NumericInput
        Dimensionless species activities, a_i, corresponding element-wise
        to `stoichiometric_coefficients`.

        Activities must be finite and strictly greater than zero because
        their natural logarithms are evaluated.

    Returns
    -------
    float | NDArray[np.float64]
        Natural logarithm of the dimensionless reaction quotient, ln(Q).

        For 1-D inputs, a scalar is returned.

        For 2-D inputs with shape ``(n_states, n_species)``, a 1-D array
        with shape ``(n_states,)`` is returned.

    Raises
    ------
    ValueError
        If the inputs contain non-finite values, are not 1-D or 2-D,
        do not have matching shapes, or if any activity is less than
        or equal to zero.

    Notes
    -----
    Activities must be dimensionless thermodynamic activities.

    Depending on the phase and thermodynamic model, an activity may be
    constructed from quantities such as fugacity, partial pressure,
    concentration, mole fraction, or activity coefficients relative to
    an appropriate standard state.
    """
    # SECTION: Normalize and validate
    nu = _as_finite_float_array(
        stoichiometric_coefficients, "stoichiometric_coefficients")
    a = _as_finite_float_array(activities, "activities")
    _validate_component_arrays(nu, a)

    # SECTION: Calculate logarithmic reaction quotient
    return _return_scalar_if_zero_dim(np.sum(nu * np.log(a), axis=-1))

# ! ::: Reaction quotient (wrapper)


def _calc_reaction_quotient(
    stoichiometric_coefficients: NumericInput,
    activities: NumericInput,
) -> float | NDArray[np.float64]:
    """
    Calculate the thermodynamic reaction quotient from species activities.

    The reaction quotient is defined as

        Q = Π_i a_i**ν_i

    where ν_i is the signed stoichiometric coefficient of species i and
    a_i is its dimensionless activity.

    The signed stoichiometric convention is assumed:

        ν_i < 0  for reactants
        ν_i > 0  for products

    Parameters
    ----------
    stoichiometric_coefficients : NumericInput
        Signed stoichiometric coefficients of the reaction.

    activities : NumericInput
        Dimensionless species activities corresponding element-wise to
        `stoichiometric_coefficients`. All activities must be strictly
        greater than zero.

    Returns
    -------
    float | NDArray[np.float64]
        Dimensionless reaction quotient, Q.

    Notes
    -----
    This function evaluates `_calc_log_reaction_quotient` and
    exponentiates the resulting ln(Q).

    At thermodynamic equilibrium,

        Q = K

    while Q < K or Q > K indicates a thermodynamic driving force away
    from the current composition.
    """
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

# ! ::: Reaction quotient from mapping


def _calc_log_reaction_quotient_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int],
    activities: Mapping[str, float | int],
) -> float:
    """
    Calculate the natural logarithm of the reaction quotient from
    species-keyed mappings.

    The logarithmic reaction quotient is

        ln(Q) = Σ_i ν_i ln(a_i)

    where ν_i is the signed stoichiometric coefficient and a_i is the
    dimensionless activity of species i.

    Parameters
    ----------
    stoichiometric_coefficients : Mapping[str, float | int]
        Mapping of species identifiers to signed stoichiometric
        coefficients.

        Reactants should have negative coefficients and products positive
        coefficients.

    activities : Mapping[str, float | int]
        Mapping of species identifiers to dimensionless activities.

        The species keys must exactly match those in
        `stoichiometric_coefficients`, and all activities must be
        strictly greater than zero.

    Returns
    -------
    float
        Natural logarithm of the dimensionless reaction quotient, ln(Q).

    Raises
    ------
    ValueError
        If the two mappings do not contain the same species keys, or if
        the underlying activity validation fails.

    Notes
    -----
    This function validates the species mappings and delegates the
    numerical calculation to `_calc_log_reaction_quotient`.
    """
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

# ! ::: Reaction quotient from mapping


def _calc_reaction_quotient_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int],
    activities: Mapping[str, float | int],
) -> float:
    """
    Calculate the thermodynamic reaction quotient from species-keyed
    stoichiometric coefficients and activities.

    The reaction quotient is

        Q = Π_i a_i**ν_i

    where ν_i is the signed stoichiometric coefficient and a_i is the
    dimensionless activity of species i.

    Parameters
    ----------
    stoichiometric_coefficients : Mapping[str, float | int]
        Mapping of species identifiers to signed stoichiometric
        coefficients.

    activities : Mapping[str, float | int]
        Mapping of the same species identifiers to dimensionless
        activities. All activities must be strictly greater than zero.

    Returns
    -------
    float
        Dimensionless reaction quotient, Q.

    Notes
    -----
    This function evaluates `_calc_log_reaction_quotient_from_mapping`
    and exponentiates the resulting ln(Q).
    """
    return float(
        np.exp(
            _calc_log_reaction_quotient_from_mapping(
                stoichiometric_coefficients,
                activities,
            )
        )
    )


# ! ::: Actual reaction Gibbs energy

def _calc_reaction_gibbs_energy(
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
    log_reaction_quotient: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """
    Calculate the reaction Gibbs energy at the actual system composition.

    The reaction Gibbs energy is related to the standard reaction Gibbs
    energy and reaction quotient by

        ΔG_rxn(T) = ΔG°_rxn(T) + R T ln(Q)

    where ΔG°_rxn(T) is the standard reaction Gibbs energy and Q is the
    dimensionless reaction quotient describing the current composition.

    Parameters
    ----------
    delta_g_reaction_std : NumericInput
        Standard reaction Gibbs energy, ΔG°_rxn(T).

        Its energy unit must be consistent with `gas_constant`.

    temperature : NumericInput
        Absolute temperature, T, in kelvin. Must be greater than zero.

    log_reaction_quotient : NumericInput
        Natural logarithm of the dimensionless reaction quotient, ln(Q).

    gas_constant : float | int, optional
        Universal gas constant, R. Default is `R_J_molK`.

    Returns
    -------
    float | NDArray[np.float64]
        Reaction Gibbs energy, ΔG_rxn(T), under the specified composition
        or activity state.

    Raises
    ------
    ValueError
        If any input contains non-finite values, temperature is not
        positive, or `gas_constant` is not positive.

    Notes
    -----
    Unlike ΔG°_rxn, the returned ΔG_rxn represents the thermodynamic
    driving force at the actual system composition.

    At equilibrium,

        Q = K
        ΔG_rxn = 0

    For the reaction direction as written:

        ΔG_rxn < 0  indicates forward thermodynamic favorability
        ΔG_rxn > 0  indicates reverse thermodynamic favorability
        ΔG_rxn = 0  indicates equilibrium
    """
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


# ! ::: van't Hoff relations

def _calc_dlnK_dT(
    delta_h_reaction_std: NumericInput,
    temperature: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """
    Calculate the temperature derivative of the logarithmic equilibrium
    constant using the differential van't Hoff equation.

    The differential van't Hoff relation is

        d ln(K) / dT = ΔH°_rxn(T) / (R T²)

    where ΔH°_rxn(T) is the standard reaction enthalpy at temperature T.

    Parameters
    ----------
    delta_h_reaction_std : NumericInput
        Standard reaction enthalpy, ΔH°_rxn(T).

        Its energy unit must be consistent with `gas_constant`.

    temperature : NumericInput
        Absolute temperature, T, in kelvin. Must be greater than zero.

    gas_constant : float | int, optional
        Universal gas constant, R. Default is `R_J_molK`.

    Returns
    -------
    float | NDArray[np.float64]
        Temperature derivative d ln(K)/dT, typically in K⁻¹.

    Notes
    -----
    This is the differential form of the van't Hoff equation and permits
    ΔH°_rxn to depend on temperature.

    For an endothermic reaction, ΔH°_rxn > 0 and K tends to increase
    with temperature.

    For an exothermic reaction, ΔH°_rxn < 0 and K tends to decrease
    with temperature.
    """
    # SECTION: Normalize and validate
    dh = _as_finite_float_array(delta_h_reaction_std, "delta_h_reaction_std")
    t = _as_finite_float_array(temperature, "temperature")
    _validate_positive_array(t, "temperature")
    r = _validate_positive_scalar(gas_constant, "gas_constant")

    # SECTION: Calculate derivative
    return _return_scalar_if_zero_dim(dh / (r * t ** 2))

# ! ::: Logarithmic equilibrium constant at temperature


def _calc_log_equilibrium_constant_at_temperature(
    equilibrium_constant_initial: float | int,
    delta_h_reaction_std: NumericInput,
    temperature_initial: NumericInput,
    temperature_final: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """
    Calculate ln(K) at a new temperature using the integrated van't Hoff
    equation.

    Assuming the standard reaction enthalpy is constant over the
    temperature interval,

        ln(K2 / K1) = -(ΔH°_rxn / R) (1/T2 - 1/T1)

    or equivalently,

        ln(K2) = ln(K1) - (ΔH°_rxn / R) (1/T2 - 1/T1)

    Parameters
    ----------
    equilibrium_constant_initial : float | int
        Dimensionless equilibrium constant K1 at `temperature_initial`.
        Must be strictly greater than zero.

    delta_h_reaction_std : NumericInput
        Standard reaction enthalpy, ΔH°_rxn, assumed constant over the
        temperature interval from T1 to T2.

        Its energy unit must be consistent with `gas_constant`.

    temperature_initial : NumericInput
        Initial absolute temperature, T1, in kelvin. Must be greater
        than zero.

    temperature_final : NumericInput
        Final absolute temperature, T2, in kelvin. Must be greater
        than zero.

    gas_constant : float | int, optional
        Universal gas constant, R. Default is `R_J_molK`.

    Returns
    -------
    float | NDArray[np.float64]
        Natural logarithm of the equilibrium constant at the final
        temperature, ln(K2).

    Raises
    ------
    ValueError
        If `equilibrium_constant_initial` is not positive, either
        temperature is not positive, any array input contains non-finite
        values, or `gas_constant` is not positive.

    Notes
    -----
    This integrated form assumes that ΔH°_rxn is approximately constant
    between `temperature_initial` and `temperature_final`.

    For larger temperature intervals, or when reaction heat-capacity
    effects are significant, ΔH°_rxn(T) should be treated as
    temperature-dependent and the van't Hoff equation should be
    integrated accordingly.
    """
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

# ! ::: Equilibrium constant at temperature


def _calc_equilibrium_constant_at_temperature(
    equilibrium_constant_initial: float | int,
    delta_h_reaction_std: NumericInput,
    temperature_initial: NumericInput,
    temperature_final: NumericInput,
    gas_constant: float | int = R_J_molK,
) -> float | NDArray[np.float64]:
    """
    Calculate the equilibrium constant at a new temperature using the
    integrated van't Hoff equation.

    Assuming the standard reaction enthalpy remains approximately constant
    over the temperature interval,

        K2 = exp[ln(K1) - (ΔH°_rxn / R) (1/T2 - 1/T1)]

    Parameters
    ----------
    equilibrium_constant_initial : float | int
        Dimensionless equilibrium constant, K1, at the initial
        temperature. Must be strictly greater than zero.

    delta_h_reaction_std : NumericInput
        Standard reaction enthalpy, ΔH°_rxn, assumed constant between
        the initial and final temperatures.

    temperature_initial : NumericInput
        Initial absolute temperature, T1, in kelvin.

    temperature_final : NumericInput
        Final absolute temperature, T2, in kelvin.

    gas_constant : float | int, optional
        Universal gas constant, R. Default is `R_J_molK`.

    Returns
    -------
    float | NDArray[np.float64]
        Dimensionless equilibrium constant, K2, at the final temperature.

    Notes
    -----
    This function evaluates
    `_calc_log_equilibrium_constant_at_temperature` and exponentiates
    the resulting ln(K2).

    The calculation assumes that ΔH°_rxn is constant over the temperature
    interval. For temperature-dependent reaction enthalpy, a more general
    integration of the van't Hoff equation is required.
    """
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
