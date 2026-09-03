"""Reaction energetic identity helpers."""

# import libs
from collections.abc import Mapping, Sequence

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list, to_scalar
# locals
from ..utils.conversions import _resolve_unit_conversion_fn
from .equilibrium import _temperature_k


# SECTION: Internal helpers

def _same_keys(left: Mapping[str, float], right: Mapping[str, float]) -> None:
    """Validate matching mapping keys."""
    # ? Mismatched keys usually indicate a missing species entropy value.
    if set(left) != set(right):
        raise ValueError("mapping inputs must have the same component keys.")


# SECTION: Standard reaction entropy

def calc_reaction_entropy_std(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    standard_entropies: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_entropy_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate standard reaction entropy from species standard entropies.

    Parameters
    ----------
    stoichiometric_coefficients : mapping or sequence of float | int | CustomProp
        Stoichiometric coefficients, positive for products and negative for
        reactants.
    standard_entropies : mapping or sequence of float | int | CustomProp
        Species standard molar entropies. ``CustomProp`` values are converted to
        ``output_entropy_unit`` when provided.
    output_entropy_unit : str, optional
        Unit used to normalize species entropies. Defaults to ``J/mol.K``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Standard reaction entropy in the normalized entropy unit.

    Notes
    -----
    Equation: ``delta_S_rxn_std = sum_i(nu_i*S_i_std)``. Mapping inputs must
    provide the same species keys for coefficients and entropies.
    """
    # SECTION: Resolve conversion function
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(stoichiometric_coefficients, Mapping) and isinstance(standard_entropies, Mapping):
        nu = to_dict(stoichiometric_coefficients, unit_conversion_fn=conversion_fn)
        entropy = to_dict(
            standard_entropies,
            output_entropy_unit,
            unit_conversion_fn=conversion_fn,
        )
        _same_keys(nu, entropy)
        return sum(nu[key] * entropy[key] for key in nu)

    # ! Mixed mapping/sequence input is ambiguous.
    if isinstance(stoichiometric_coefficients, Mapping) or isinstance(standard_entropies, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    nu = to_list(stoichiometric_coefficients, unit_conversion_fn=conversion_fn)
    entropy = to_list(
        standard_entropies,
        output_entropy_unit,
        unit_conversion_fn=conversion_fn,
    )
    if len(nu) != len(entropy):
        raise ValueError("stoichiometric_coefficients and standard_entropies must have the same length.")
    return sum(nu_i * s_i for nu_i, s_i in zip(nu, entropy))


# SECTION: Entropy from enthalpy and Gibbs energy

def calc_reaction_entropy_std_from_enthalpy_gibbs(
    delta_h_reaction_std: ScalarValue,
    delta_g_reaction_std: ScalarValue,
    temperature: Temperature,
    output_delta_h_unit: str | None = "J/mol",
    output_delta_g_unit: str | None = "J/mol",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate standard reaction entropy from standard enthalpy and Gibbs energy.

    Parameters
    ----------
    delta_h_reaction_std : float | int | CustomProp
        Standard reaction enthalpy.
    delta_g_reaction_std : float | int | CustomProp
        Standard reaction Gibbs energy.
    temperature : Temperature
        Temperature at which the relation is evaluated. Converted to K.
    output_delta_h_unit : str, optional
        Unit used to normalize ``delta_h_reaction_std``. Defaults to ``J/mol``.
    output_delta_g_unit : str, optional
        Unit used to normalize ``delta_g_reaction_std``. Defaults to ``J/mol``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Standard reaction entropy, typically J/mol/K.

    Notes
    -----
    Equation: ``delta_S_rxn_std = (delta_H_rxn_std - delta_G_rxn_std)/T``.
    Enthalpy and Gibbs energy must use compatible reaction bases.
    """
    # SECTION: Normalize inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    dh = to_scalar(
        delta_h_reaction_std,
        "delta_h_reaction_std",
        output_delta_h_unit,
        unit_conversion_fn=conversion_fn,
    )
    dg = to_scalar(
        delta_g_reaction_std,
        "delta_g_reaction_std",
        output_delta_g_unit,
        unit_conversion_fn=conversion_fn,
    )
    temperature_k = _temperature_k(temperature, unit_conversion_fn)

    # SECTION: Calculate reaction entropy
    return (dh - dg) / temperature_k


# SECTION: Public exports
__all__ = [
    "calc_reaction_entropy_std",
    "calc_reaction_entropy_std_from_enthalpy_gibbs",
]

