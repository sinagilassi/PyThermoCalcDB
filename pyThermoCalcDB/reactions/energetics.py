"""Reaction energetic identity helpers."""

# import libs
from collections.abc import Mapping, Sequence

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list, to_scalar
from pycuc.canonical import to_K
# locals
from ..utils.conversions import _resolve_unit_conversion_fn
from .core.energetics import (
    _calc_reaction_entropy_std,
    _calc_reaction_entropy_std_from_enthalpy_gibbs,
    _calc_reaction_entropy_std_from_mapping,
)


# ! ::: Standard reaction entropy

def calc_reaction_entropy_std_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp],
    standard_entropies: Mapping[str, float | int | CustomProp],
    output_entropy_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate standard reaction entropy from species standard entropies.

    Parameters
    ----------
    stoichiometric_coefficients : mapping of float | int | CustomProp
        Stoichiometric coefficients, positive for products and negative for
        reactants.
    standard_entropies : mapping of float | int | CustomProp
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
    nu = to_dict(stoichiometric_coefficients,)
    entropy = to_dict(
        standard_entropies,
        output_entropy_unit,
        unit_conversion_fn=conversion_fn,
    )
    return _calc_reaction_entropy_std_from_mapping(nu, entropy)

# ! ::: Standard reaction entropy from sequence


def calc_reaction_entropy_std_from_sequence(
    stoichiometric_coefficients: Sequence[float | int | CustomProp],
    standard_entropies: Sequence[float | int | CustomProp],
    output_entropy_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate standard reaction entropy from species standard entropies.

    Parameters
    ----------
    stoichiometric_coefficients : sequence of float | int | CustomProp
        Stoichiometric coefficients, positive for products and negative for
        reactants.
    standard_entropies : sequence of float | int | CustomProp
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

    # SECTION: Sequence implementation
    nu = to_list(stoichiometric_coefficients)
    entropy = to_list(
        standard_entropies,
        output_entropy_unit,
        unit_conversion_fn=conversion_fn,
    )
    return float(_calc_reaction_entropy_std(nu, entropy))

# ! ::: Entropy from enthalpy and Gibbs energy


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
    temperature_k = to_K(
        value=temperature.value,
        from_unit=temperature.unit,
    )

    # SECTION: Calculate reaction entropy
    return float(
        _calc_reaction_entropy_std_from_enthalpy_gibbs(
            dh,
            dg,
            temperature_k,
        )
    )


# SECTION: Public exports
__all__ = [
    "calc_reaction_entropy_std_from_mapping",
    "calc_reaction_entropy_std_from_sequence",
    "calc_reaction_entropy_std_from_enthalpy_gibbs",
]
