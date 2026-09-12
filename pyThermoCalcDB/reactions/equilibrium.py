"""Reaction equilibrium primitives."""

# import libs
from collections.abc import Mapping, Sequence
import math

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos, to_dict, to_list, to_scalar
from pycuc.canonical import to_J_per_mol, to_K
from pycuc import convert_from_to
# locals
from ..utils.conversions import _resolve_unit_conversion_fn, _to_kelvin, _to_scalar, _pos
# ! core
from .core.equilibrium import (
    _calc_dlnK_dT,
    _calc_equilibrium_constant,
    _calc_equilibrium_constant_at_temperature,
    _calc_log_equilibrium_constant_at_temperature,
    _calc_log_equilibrium_constant,
    _calc_log_reaction_quotient,
    _calc_log_reaction_quotient_from_mapping,
    _calc_reaction_gibbs_energy,
    _calc_reaction_quotient,
    _calc_reaction_quotient_from_mapping,
)


# SECTION: Equilibrium constant

def calc_log_equilibrium_constant(
    delta_g_reaction_std: CustomProp,
    temperature: Temperature,
) -> float:
    """Calculate the natural logarithm of the equilibrium constant.

    Parameters
    ----------
    delta_g_reaction_std : CustomProp
        Standard Gibbs energy of reaction, provided as a ``CustomProp``. Then converted to J/mol before calculation.
    temperature : Temperature
        Reaction temperature. Converted to K before calculation.

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
    # ! to J/mol
    dg = to_J_per_mol(
        value=delta_g_reaction_std.value,
        from_unit=delta_g_reaction_std.unit,
    )
    # ! to K
    temperature_k = _to_kelvin(temperature)

    # NOTE: ln(K) is exposed to avoid unnecessary exp overflow/underflow.
    return float(_calc_log_equilibrium_constant(dg, temperature_k))

# ! ::: Calculate equilibrium constant


def calc_equilibrium_constant(
    delta_g_reaction_std: CustomProp,
    temperature: Temperature,
) -> float:
    """Calculate the dimensionless thermodynamic equilibrium constant.

    Parameters
    ----------
    delta_g_reaction_std : CustomProp
        Standard Gibbs energy of reaction, provided as a ``CustomProp``. Then converted to J/mol before calculation.
    temperature : Temperature
        Reaction temperature. Converted to K before calculation.

    Returns
    -------
    float
        Dimensionless equilibrium constant ``K``.

    Notes
    -----
    Equation: ``K = exp(-delta_G_rxn_std/(R*T))``. Activities defining ``K``
    are assumed dimensionless relative to their standard states.
    """
    # SECTION: Normalize inputs
    # ! to J/mol
    dg = to_J_per_mol(
        value=delta_g_reaction_std.value,
        from_unit=delta_g_reaction_std.unit,
    )
    # ! to K
    temperature_k = _to_kelvin(temperature)

    # SECTION: Calculate equilibrium constant
    return float(_calc_equilibrium_constant(dg, temperature_k))


# ! ::: Reaction quotient

def calc_log_reaction_quotient(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    activities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
) -> float:
    """Calculate the natural logarithm of the reaction quotient.

    Parameters
    ----------
    stoichiometric_coefficients : mapping or sequence of float | int | CustomProp
        Stoichiometric coefficients, positive for products and negative for
        reactants.
    activities : mapping or sequence of float | int | CustomProp
        Dimensionless species activities. Values must be greater than zero.

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
    # SECTION: Mapping implementation
    if isinstance(stoichiometric_coefficients, Mapping) and isinstance(activities, Mapping):
        nu = to_dict(stoichiometric_coefficients)
        a = to_dict(activities)
        return _calc_log_reaction_quotient_from_mapping(nu, a)

    # ! Mixed mapping/sequence input is ambiguous.
    if isinstance(stoichiometric_coefficients, Mapping) or isinstance(activities, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    nu = to_list(stoichiometric_coefficients)
    a = to_list(activities)
    return float(_calc_log_reaction_quotient(nu, a))

# ! ::: Reaction quotient


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
    # SECTION: Resolve conversion function
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(stoichiometric_coefficients, Mapping) and isinstance(activities, Mapping):
        nu = to_dict(
            stoichiometric_coefficients,
            unit_conversion_fn=conversion_fn
        )
        a = to_dict(
            activities,
            unit_conversion_fn=conversion_fn
        )
        return _calc_reaction_quotient_from_mapping(nu, a)

    # ! Mixed mapping/sequence input is ambiguous.
    if isinstance(stoichiometric_coefficients, Mapping) or isinstance(activities, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    nu = to_list(stoichiometric_coefficients, unit_conversion_fn=conversion_fn)
    a = to_list(activities, unit_conversion_fn=conversion_fn)
    return float(_calc_reaction_quotient(nu, a))


# ! ::: Actual reaction Gibbs energy

def calc_reaction_gibbs_energy(
    delta_g_reaction_std: CustomProp,
    temperature: Temperature,
    reaction_quotient: ScalarValue | None = None,
    log_reaction_quotient: ScalarValue | None = None,
    output_unit: str = "J/mol"
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
        raise ValueError(
            "reaction_quotient or log_reaction_quotient must be provided.")
    if reaction_quotient is not None and log_reaction_quotient is not None:
        raise ValueError(
            "provide only one of reaction_quotient or log_reaction_quotient.")

    # SECTION: Normalize thermodynamic inputs
    # ! to J/mol
    dg_std = to_J_per_mol(
        delta_g_reaction_std.value,
        from_unit=delta_g_reaction_std.unit,
    )
    # ! to K
    temperature_k = to_K(
        temperature.value,
        from_unit=temperature.unit,
    )

    # NOTE: Prefer caller-provided ln(Q) when available for numerical stability.
    if log_reaction_quotient is not None:
        ln_q = _to_scalar(log_reaction_quotient, "log_reaction_quotient")
    else:
        if reaction_quotient is None:
            raise ValueError(
                "reaction_quotient or log_reaction_quotient must be provided.")
        q = _pos(reaction_quotient, "reaction_quotient")
        ln_q = math.log(q)

    # SECTION: Calculate actual reaction Gibbs energy
    res = float(_calc_reaction_gibbs_energy(dg_std, temperature_k, ln_q))

    # NOTE: unit conversion if needed
    if output_unit != "J/mol":
        res = convert_from_to(res, from_unit="J/mol", to_unit=output_unit)
    return res


# ! ::: van't Hoff relations

def calc_dlnK_dT(
    delta_h_reaction_std: CustomProp,
    temperature: Temperature,
) -> float:
    """
    Calculate the temperature derivative of the logarithmic equilibrium
    constant using the differential van't Hoff equation.

    The derivative is calculated as

        d ln(K) / dT = ΔH°_rxn(T) / (R T²)

    Parameters
    ----------
    delta_h_reaction_std : CustomProp
        Standard reaction enthalpy, ΔH°_rxn(T). The value is internally
        converted to J/mol.

    temperature : Temperature
        Absolute temperature at which the derivative is evaluated.
        The value is internally converted to kelvin.

    Returns
    -------
    float
        Temperature derivative d ln(K)/dT, in K⁻¹.

    Notes
    -----
    The ``_std`` suffix denotes a standard-state reaction property and
    does not imply a temperature of 298.15 K.

    For an endothermic reaction, ΔH°_rxn > 0 and K increases with
    temperature locally. For an exothermic reaction, ΔH°_rxn < 0 and K
    decreases with temperature locally.
    """
    # SECTION: Normalize inputs
    # ! to J/mol
    dh = to_J_per_mol(
        delta_h_reaction_std.value,
        from_unit=delta_h_reaction_std.unit,
    )
    # ! to K
    temperature_k = _to_kelvin(temperature)

    # SECTION: Calculate derivative
    return float(_calc_dlnK_dT(dh, temperature_k))

# ! ::: Integrated van't Hoff


def calc_equilibrium_constant_at_temperature(
    equilibrium_constant_initial: float,
    delta_h_reaction_std: CustomProp,
    temperature_initial: Temperature,
    temperature_final: Temperature,
) -> float:
    """
    Calculate the equilibrium constant at a new temperature using the
    integrated van't Hoff equation.

    Assuming the standard reaction enthalpy remains approximately constant
    between the two temperatures,

        ln(K2 / K1) = -(ΔH°_rxn / R) (1/T2 - 1/T1)

    Parameters
    ----------
    equilibrium_constant_initial : float
        Dimensionless equilibrium constant, K1, at `temperature_initial`.
        Must be greater than zero.

    delta_h_reaction_std : CustomProp
        Standard reaction enthalpy, ΔH°_rxn, assumed approximately constant
        over the temperature interval. The value is internally converted
        to J/mol.

    temperature_initial : Temperature
        Initial absolute temperature, T1. Internally converted to kelvin.

    temperature_final : Temperature
        Final absolute temperature, T2. Internally converted to kelvin.

    Returns
    -------
    float
        Dimensionless equilibrium constant, K2, at `temperature_final`.

    Notes
    -----
    The integrated van't Hoff equation used here assumes that ΔH°_rxn is
    approximately constant between T1 and T2.

    For a significant temperature interval, especially when reaction heat
    capacity effects are important, ΔH°_rxn(T) should be treated as
    temperature-dependent and the differential van't Hoff relation should
    be integrated accordingly.
    """
    # SECTION: Normalize inputs
    # ! check that equilibrium_constant_initial is positive
    k_initial = _pos(
        equilibrium_constant_initial,
        "equilibrium_constant_initial"
    )
    # ! to J/mol
    dh = to_J_per_mol(
        delta_h_reaction_std.value,
        from_unit=delta_h_reaction_std.unit,
    )
    # ! to K
    t_initial = _to_kelvin(temperature_initial)
    t_final = _to_kelvin(temperature_final)

    # SECTION: Calculate final equilibrium constant
    return float(
        _calc_equilibrium_constant_at_temperature(
            k_initial,
            dh,
            t_initial,
            t_final,
        )
    )

# ! :::


def calc_log_equilibrium_constant_at_temperature(
    equilibrium_constant_initial: float,
    delta_h_reaction_std: CustomProp,
    temperature_initial: Temperature,
    temperature_final: Temperature,
):
    """
    Calculate the natural logarithm of the equilibrium constant at a new
    temperature using the integrated van't Hoff equation.

    Assuming the standard reaction enthalpy remains approximately constant
    over the temperature interval,

        ln(K2 / K1) = -(ΔH°_rxn / R) (1/T2 - 1/T1)

    or equivalently,

        ln(K2) = ln(K1) - (ΔH°_rxn / R) (1/T2 - 1/T1)

    Parameters
    ----------
    equilibrium_constant_initial : float
        Dimensionless equilibrium constant, K1, at `temperature_initial`.
        Must be strictly greater than zero.

    delta_h_reaction_std : CustomProp
        Standard reaction enthalpy, ΔH°_rxn, assumed approximately constant
        between `temperature_initial` and `temperature_final`.

        The value is internally converted to J/mol before calculation.

    temperature_initial : Temperature
        Initial absolute temperature, T1. The value is internally converted
        to kelvin.

    temperature_final : Temperature
        Final absolute temperature, T2. The value is internally converted
        to kelvin.

    Returns
    -------
    float
        Natural logarithm of the dimensionless equilibrium constant at the
        final temperature, ln(K2).

    Raises
    ------
    ValueError
        If `equilibrium_constant_initial` is less than or equal to zero,
        or if either temperature is invalid or non-positive.

    Notes
    -----
    The ``_std`` suffix denotes a standard-state reaction property and
    does not imply a temperature of 298.15 K.

    This integrated van't Hoff relation assumes that ΔH°_rxn is
    approximately constant over the temperature interval from T1 to T2.

    For large temperature intervals, or when reaction heat-capacity effects
    are significant, ΔH°_rxn(T) should be treated as temperature-dependent
    and the differential van't Hoff equation should be integrated instead.

    Returning ln(K2) can be preferable to returning K2 directly when the
    equilibrium constant is very large or very small, because it avoids
    unnecessary exponential overflow or underflow.
    """
    return float(
        _calc_log_equilibrium_constant_at_temperature(
            _pos(
                equilibrium_constant_initial,
                "equilibrium_constant_initial"
            ),
            to_J_per_mol(
                delta_h_reaction_std.value,
                from_unit=delta_h_reaction_std.unit,
            ),
            _to_kelvin(temperature_initial),
            _to_kelvin(temperature_final),
        )
    )


# SECTION: Public exports
__all__ = [
    "calc_log_equilibrium_constant",
    "calc_equilibrium_constant",
    "calc_log_reaction_quotient",
    "calc_reaction_quotient",
    "calc_reaction_gibbs_energy",
    "calc_dlnK_dT",
    "calc_equilibrium_constant_at_temperature",
    "calc_log_equilibrium_constant_at_temperature",
]
