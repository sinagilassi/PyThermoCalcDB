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
    _calc_reaction_enthalpy_from_constant_delta_cp,
    _calc_reaction_entropy_std,
    _calc_reaction_entropy_std_from_enthalpy_gibbs,
    _calc_reaction_entropy_std_from_mapping,
    _calc_reaction_heat_capacity_change,
    _calc_reaction_heat_capacity_change_from_mapping,
    _calc_reaction_heat_rate,
    _calc_reaction_volumetric_heat_source,
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


def calc_reaction_entropy_std(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    standard_entropies: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_entropy_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate standard reaction entropy from species standard entropies."""
    if isinstance(stoichiometric_coefficients, Mapping) and isinstance(standard_entropies, Mapping):
        return calc_reaction_entropy_std_from_mapping(
            stoichiometric_coefficients,
            standard_entropies,
            output_entropy_unit,
            unit_conversion_fn,
        )
    if isinstance(stoichiometric_coefficients, Mapping) or isinstance(standard_entropies, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")
    return calc_reaction_entropy_std_from_sequence(
        stoichiometric_coefficients,
        standard_entropies,
        output_entropy_unit,
        unit_conversion_fn,
    )

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


# SECTION: Reaction heat-capacity change

def calc_reaction_heat_capacity_change_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp],
    component_heat_capacities: Mapping[str, float | int | CustomProp],
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate reaction heat-capacity change from component heat capacities.

    Equation: ``delta_Cp_r = sum_i(nu_i*Cp_i)``. Stoichiometric coefficients
    are signed: negative for reactants and positive for products.
    """
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    nu = to_dict(stoichiometric_coefficients)
    cp = to_dict(
        component_heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn=conversion_fn,
    )
    return _calc_reaction_heat_capacity_change_from_mapping(nu, cp)


def calc_reaction_heat_capacity_change_from_sequence(
    stoichiometric_coefficients: Sequence[float | int | CustomProp],
    component_heat_capacities: Sequence[float | int | CustomProp],
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate reaction heat-capacity change from sequence inputs."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    nu = to_list(stoichiometric_coefficients)
    cp = to_list(
        component_heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn=conversion_fn,
    )
    return float(_calc_reaction_heat_capacity_change(nu, cp))


def calc_reaction_heat_capacity_change(
    stoichiometric_coefficients: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    component_heat_capacities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate reaction heat-capacity change from species heat capacities."""
    if isinstance(stoichiometric_coefficients, Mapping) and isinstance(component_heat_capacities, Mapping):
        return calc_reaction_heat_capacity_change_from_mapping(
            stoichiometric_coefficients,
            component_heat_capacities,
            output_heat_capacity_unit,
            unit_conversion_fn,
        )
    if isinstance(stoichiometric_coefficients, Mapping) or isinstance(component_heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")
    return calc_reaction_heat_capacity_change_from_sequence(
        stoichiometric_coefficients,
        component_heat_capacities,
        output_heat_capacity_unit,
        unit_conversion_fn,
    )


calc_reaction_heat_capacity_change_from_props = calc_reaction_heat_capacity_change_from_mapping


# SECTION: Kirchhoff enthalpy correction

def calc_reaction_enthalpy_from_constant_delta_cp(
    delta_h_reaction_ref: ScalarValue,
    delta_cp_reaction: ScalarValue,
    temperature: Temperature,
    reference_temperature: Temperature,
    output_delta_h_unit: str | None = "J/mol",
    output_delta_cp_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate reaction enthalpy from a constant reaction heat-capacity change."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    dh_ref = to_scalar(
        delta_h_reaction_ref,
        "delta_h_reaction_ref",
        output_delta_h_unit,
        unit_conversion_fn=conversion_fn,
    )
    dcp = to_scalar(
        delta_cp_reaction,
        "delta_cp_reaction",
        output_delta_cp_unit,
        unit_conversion_fn=conversion_fn,
    )
    temperature_k = to_K(value=temperature.value, from_unit=temperature.unit)
    reference_temperature_k = to_K(
        value=reference_temperature.value,
        from_unit=reference_temperature.unit,
    )
    return float(
        _calc_reaction_enthalpy_from_constant_delta_cp(
            dh_ref,
            dcp,
            temperature_k,
            reference_temperature_k,
        )
    )


# SECTION: Reaction heat source/rate

def calc_reaction_volumetric_heat_source(
    reaction_enthalpies: Sequence[float | int | CustomProp],
    reaction_rates: Sequence[float | int | CustomProp],
    output_reaction_enthalpy_unit: str | None = "J/mol",
    output_reaction_rate_unit: str | None = "mol/m3.s",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate generated volumetric reaction heat source.

    Equation: ``q''' = -sum_j(delta_H_r,j*r_j)``. Exothermic reactions have
    ``delta_H_r < 0`` and produce positive generated heat for positive rates.
    """
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    dh = to_list(
        reaction_enthalpies,
        output_reaction_enthalpy_unit,
        unit_conversion_fn=conversion_fn,
    )
    rates = to_list(
        reaction_rates,
        output_reaction_rate_unit,
        unit_conversion_fn=conversion_fn,
    )
    return float(_calc_reaction_volumetric_heat_source(dh, rates))


def calc_reaction_heat_rate(
    reaction_enthalpies: Sequence[float | int | CustomProp],
    reaction_rates: Sequence[float | int | CustomProp],
    volume: ScalarValue,
    output_reaction_enthalpy_unit: str | None = "J/mol",
    output_reaction_rate_unit: str | None = "mol/m3.s",
    output_volume_unit: str | None = "m3",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate total generated reaction heat rate.

    Equation: ``Qdot = -V*sum_j(delta_H_r,j*r_j)``.
    """
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    dh = to_list(
        reaction_enthalpies,
        output_reaction_enthalpy_unit,
        unit_conversion_fn=conversion_fn,
    )
    rates = to_list(
        reaction_rates,
        output_reaction_rate_unit,
        unit_conversion_fn=conversion_fn,
    )
    v = to_scalar(
        volume,
        "volume",
        output_volume_unit,
        unit_conversion_fn=conversion_fn,
    )
    return float(_calc_reaction_heat_rate(dh, rates, v))


# SECTION: Public exports
__all__ = [
    "calc_reaction_entropy_std",
    "calc_reaction_entropy_std_from_mapping",
    "calc_reaction_entropy_std_from_sequence",
    "calc_reaction_entropy_std_from_enthalpy_gibbs",
    "calc_reaction_heat_capacity_change",
    "calc_reaction_heat_capacity_change_from_mapping",
    "calc_reaction_heat_capacity_change_from_sequence",
    "calc_reaction_heat_capacity_change_from_props",
    "calc_reaction_enthalpy_from_constant_delta_cp",
    "calc_reaction_volumetric_heat_source",
    "calc_reaction_heat_rate",
]
