"""Reaction equilibrium primitives."""

# import libs
import math
from collections.abc import Mapping, Sequence

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos, to_dict, to_list, to_scalar
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import _resolve_unit_conversion_fn


# SECTION: Internal scalar helpers

def _scalar(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert scalar input to float, optionally normalizing units."""
    return to_scalar(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


def _pos(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert scalar input to a positive float, optionally normalizing units."""
    return pos(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


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
        value = float(conversion_fn(value, unit, "K"))

    # ! Log/equilibrium thermodynamic identities require T > 0 K.
    if value <= 0.0:
        raise ValueError("temperature must be greater than zero K.")
    return value


def _same_keys(left: Mapping[str, float], right: Mapping[str, float]) -> None:
    """Validate matching mapping keys."""
    # ? Mismatched keys usually indicate a missing species activity or coefficient.
    if set(left) != set(right):
        raise ValueError("mapping inputs must have the same component keys.")


# SECTION: Equilibrium constant

def calc_log_equilibrium_constant(
    delta_g_reaction_std: ScalarValue,
    temperature: Temperature,
    output_delta_g_unit: str | None = "J/mol",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate the natural logarithm of the equilibrium constant.

    Parameters
    ----------
    delta_g_reaction_std : float | int | CustomProp
        Standard Gibbs energy of reaction. Converted to
        ``output_delta_g_unit`` when supplied as ``CustomProp``.
    temperature : Temperature
        Reaction temperature. Converted to K before calculation.
    output_delta_g_unit : str, optional
        Unit used for ``delta_g_reaction_std``. Defaults to ``J/mol``.
    gas_constant : float, optional
        Gas constant consistent with ``output_delta_g_unit`` per mol per K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        ``ln(K)`` for a dimensionless thermodynamic equilibrium constant.

    Notes
    -----
    Equation: ``ln(K) = -delta_G_rxn_std/(R*T)``. Exposing ``ln(K)`` avoids
    unnecessary overflow/underflow when the exponential form is not needed.

    Raises
    ------
    ValueError
        If temperature or gas constant is not positive.
    """
    # SECTION: Normalize inputs
    dg = _scalar(
        delta_g_reaction_std,
        "delta_g_reaction_std",
        output_delta_g_unit,
        unit_conversion_fn,
    )
    temperature_k = _temperature_k(temperature, unit_conversion_fn)
    r = _pos(gas_constant, "gas_constant")

    # NOTE: ln(K) is exposed to avoid unnecessary exp overflow/underflow.
    return -dg / (r * temperature_k)


def calc_equilibrium_constant(
    delta_g_reaction_std: ScalarValue,
    temperature: Temperature,
    output_delta_g_unit: str | None = "J/mol",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate the dimensionless thermodynamic equilibrium constant.

    Parameters
    ----------
    delta_g_reaction_std : float | int | CustomProp
        Standard Gibbs energy of reaction.
    temperature : Temperature
        Reaction temperature. Converted to K before calculation.
    output_delta_g_unit : str, optional
        Unit used for ``delta_g_reaction_std``. Defaults to ``J/mol``.
    gas_constant : float, optional
        Gas constant in units consistent with energy and K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Dimensionless equilibrium constant ``K``.

    Notes
    -----
    Equation: ``K = exp(-delta_G_rxn_std/(R*T))``. Activities defining ``K``
    are assumed dimensionless relative to their standard states.
    """
    # SECTION: Calculate from logarithmic form
    return math.exp(
        calc_log_equilibrium_constant(
            delta_g_reaction_std,
            temperature,
            output_delta_g_unit,
            gas_constant,
            unit_conversion_fn,
        )
    )


# SECTION: Reaction quotient

def calc_log_reaction_quotient(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    activities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate the natural logarithm of the reaction quotient.

    Parameters
    ----------
    stoichiometric_coefficients : mapping or sequence of float | int | CustomProp
        Stoichiometric coefficients, positive for products and negative for
        reactants.
    activities : mapping or sequence of float | int | CustomProp
        Dimensionless species activities. Values must be greater than zero.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function used for ``CustomProp`` normalization.

    Returns
    -------
    float
        ``ln(Q)``.

    Notes
    -----
    Equation: ``ln(Q) = sum_i(nu_i*ln(a_i))``. This function does not calculate
    activities or activity coefficients.

    Raises
    ------
    TypeError
        If one component input is a mapping and the other is a sequence.
    ValueError
        If activities are non-positive or input shapes/keys differ.
    """
    # SECTION: Resolve conversion function
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(stoichiometric_coefficients, Mapping) and isinstance(activities, Mapping):
        nu = to_dict(stoichiometric_coefficients, unit_conversion_fn=conversion_fn)
        a = to_dict(activities, unit_conversion_fn=conversion_fn)
        _same_keys(nu, a)

        # ! Activities must already be dimensionless and strictly positive.
        if any(value <= 0.0 for value in a.values()):
            raise ValueError("activities must be greater than zero.")
        return sum(nu[key] * math.log(a[key]) for key in nu)

    # ! Mixed mapping/sequence input is ambiguous.
    if isinstance(stoichiometric_coefficients, Mapping) or isinstance(activities, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    nu = to_list(stoichiometric_coefficients, unit_conversion_fn=conversion_fn)
    a = to_list(activities, unit_conversion_fn=conversion_fn)
    if len(nu) != len(a):
        raise ValueError("stoichiometric_coefficients and activities must have the same length.")
    if any(value <= 0.0 for value in a):
        raise ValueError("activities must be greater than zero.")
    return sum(nu_i * math.log(a_i) for nu_i, a_i in zip(nu, a))


def calc_reaction_quotient(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    activities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate the dimensionless reaction quotient from activities.

    Parameters
    ----------
    stoichiometric_coefficients : mapping or sequence of float | int | CustomProp
        Stoichiometric coefficients, positive for products and negative for
        reactants.
    activities : mapping or sequence of float | int | CustomProp
        Dimensionless species activities.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Dimensionless reaction quotient ``Q``.

    Notes
    -----
    Equation: ``Q = product_i(a_i**nu_i)``. The logarithmic form is used
    internally.
    """
    # SECTION: Calculate from logarithmic form
    return math.exp(
        calc_log_reaction_quotient(
            stoichiometric_coefficients,
            activities,
            unit_conversion_fn,
        )
    )


# SECTION: Actual reaction Gibbs energy

def calc_reaction_gibbs_energy(
    delta_g_reaction_std: ScalarValue,
    temperature: Temperature,
    reaction_quotient: ScalarValue | None = None,
    log_reaction_quotient: ScalarValue | None = None,
    output_delta_g_unit: str | None = "J/mol",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate actual reaction Gibbs energy from standard state and ``Q``.

    Parameters
    ----------
    delta_g_reaction_std : float | int | CustomProp
        Standard Gibbs energy of reaction.
    temperature : Temperature
        Reaction temperature. Converted to K before calculation.
    reaction_quotient : float | int | CustomProp, optional
        Dimensionless reaction quotient ``Q``. Must be positive. Provide this
        or ``log_reaction_quotient``, not both.
    log_reaction_quotient : float | int | CustomProp, optional
        Natural logarithm of ``Q``. Preferred when already available.
    output_delta_g_unit : str, optional
        Unit used for ``delta_g_reaction_std``. Defaults to ``J/mol``.
    gas_constant : float, optional
        Gas constant in units consistent with energy and K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Actual reaction Gibbs energy in the normalized energy unit.

    Notes
    -----
    Equation: ``delta_G_rxn = delta_G_rxn_std + R*T*ln(Q)``.
    """
    # SECTION: Validate quotient source
    if reaction_quotient is None and log_reaction_quotient is None:
        raise ValueError("reaction_quotient or log_reaction_quotient must be provided.")
    if reaction_quotient is not None and log_reaction_quotient is not None:
        raise ValueError("provide only one of reaction_quotient or log_reaction_quotient.")

    # SECTION: Normalize thermodynamic inputs
    dg_std = _scalar(
        delta_g_reaction_std,
        "delta_g_reaction_std",
        output_delta_g_unit,
        unit_conversion_fn,
    )
    temperature_k = _temperature_k(temperature, unit_conversion_fn)
    r = _pos(gas_constant, "gas_constant")

    # NOTE: Prefer caller-provided ln(Q) when available for numerical stability.
    if log_reaction_quotient is None:
        q = _pos(reaction_quotient, "reaction_quotient")
        ln_q = math.log(q)
    else:
        ln_q = _scalar(log_reaction_quotient, "log_reaction_quotient")

    # SECTION: Calculate actual reaction Gibbs energy
    return dg_std + r * temperature_k * ln_q


# SECTION: van't Hoff relations

def calc_dlnK_dT(
    delta_h_reaction_std: ScalarValue,
    temperature: Temperature,
    output_delta_h_unit: str | None = "J/mol",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate the van't Hoff derivative ``dlnK/dT``.

    Parameters
    ----------
    delta_h_reaction_std : float | int | CustomProp
        Standard reaction enthalpy. Converted to ``output_delta_h_unit`` when
        supplied as ``CustomProp``.
    temperature : Temperature
        Temperature at which the derivative is evaluated. Converted to K.
    output_delta_h_unit : str, optional
        Unit used for ``delta_h_reaction_std``. Defaults to ``J/mol``.
    gas_constant : float, optional
        Gas constant consistent with enthalpy per mol per K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Temperature derivative ``dlnK/dT`` in ``1/K``.

    Notes
    -----
    Equation: ``dlnK/dT = delta_H_rxn_std/(R*T**2)``.
    """
    # SECTION: Normalize inputs
    dh = _scalar(
        delta_h_reaction_std,
        "delta_h_reaction_std",
        output_delta_h_unit,
        unit_conversion_fn,
    )
    temperature_k = _temperature_k(temperature, unit_conversion_fn)
    r = _pos(gas_constant, "gas_constant")

    # SECTION: Calculate derivative
    return dh / (r * temperature_k ** 2)


def calc_equilibrium_constant_at_temperature(
    equilibrium_constant_initial: ScalarValue,
    delta_h_reaction_std: ScalarValue,
    temperature_initial: Temperature,
    temperature_final: Temperature,
    output_delta_h_unit: str | None = "J/mol",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate equilibrium constant at a new temperature by integrated van't Hoff.

    Parameters
    ----------
    equilibrium_constant_initial : float | int | CustomProp
        Initial dimensionless equilibrium constant ``K1``.
    delta_h_reaction_std : float | int | CustomProp
        Standard reaction enthalpy, assumed constant over the temperature
        interval.
    temperature_initial : Temperature
        Initial temperature ``T1``. Converted to K.
    temperature_final : Temperature
        Final temperature ``T2``. Converted to K.
    output_delta_h_unit : str, optional
        Unit used for ``delta_h_reaction_std``. Defaults to ``J/mol``.
    gas_constant : float, optional
        Gas constant consistent with enthalpy per mol per K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Estimated dimensionless equilibrium constant ``K2``.

    Notes
    -----
    Equation: ``ln(K2/K1) = -delta_H_rxn_std/R * (1/T2 - 1/T1)``. This is the
    integrated van't Hoff relation for approximately constant reaction enthalpy.
    """
    # SECTION: Normalize inputs
    k_initial = _pos(equilibrium_constant_initial, "equilibrium_constant_initial")
    dh = _scalar(
        delta_h_reaction_std,
        "delta_h_reaction_std",
        output_delta_h_unit,
        unit_conversion_fn,
    )
    t_initial = _temperature_k(temperature_initial, unit_conversion_fn)
    t_final = _temperature_k(temperature_final, unit_conversion_fn)
    r = _pos(gas_constant, "gas_constant")

    # NOTE: The logarithmic form is used internally for numerical stability.
    ln_k_final = math.log(k_initial) - (dh / r) * (1.0 / t_final - 1.0 / t_initial)

    # SECTION: Calculate final equilibrium constant
    return math.exp(ln_k_final)


# SECTION: Public exports
__all__ = [
    "calc_log_equilibrium_constant",
    "calc_equilibrium_constant",
    "calc_log_reaction_quotient",
    "calc_reaction_quotient",
    "calc_reaction_gibbs_energy",
    "calc_dlnK_dT",
    "calc_equilibrium_constant_at_temperature",
]
