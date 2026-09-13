"""Core reaction energetic identity calculations."""

# import libs
from collections.abc import Mapping, Sequence
from typing import cast, TypeAlias, Optional

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, Component, ComponentKey

# locals
from ...utils.conversions import (
    NumericArrayInput,
    _return_scalar_if_zero_dim,
    _validate_positive_array,
    _validate_same_mapping_keys,
    _configure_component_values,
    _to_values,
)

# SECTION: Type aliases
NumericInput: TypeAlias = NumericArrayInput


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
    component_values: NDArray[np.float64],
) -> None:
    """Validate pairwise component arrays for reaction entropy calculations."""
    if stoichiometric_coefficients.ndim not in (1, 2):
        raise ValueError(
            "stoichiometric_coefficients must be a 1-D or 2-D array.")
    if component_values.ndim not in (1, 2):
        raise ValueError("component_values must be a 1-D or 2-D array.")
    if stoichiometric_coefficients.shape != component_values.shape:
        raise ValueError(
            "stoichiometric_coefficients and component_values must have the same shape.")


# ! ::: Standard reaction entropy

def _calc_reaction_entropy_std(
    stoichiometric_coefficients: NumericInput,
    standard_entropies: NumericInput,
) -> float | NDArray[np.float64]:
    """
    Calculate the standard reaction entropy from stoichiometric coefficients
    and species standard molar entropies.

    The standard reaction entropy is defined as

        ΔS°_rxn(T) = Σ_i ν_i S°_i(T)

    where ν_i is the signed stoichiometric coefficient of species i and
    S°_i(T) is the standard molar entropy of that species at temperature T.

    The signed stoichiometric convention is assumed:

        ν_i < 0  for reactants
        ν_i > 0  for products

    Therefore, the expression is equivalent to

        ΔS°_rxn(T) = Σ_products ν_i S°_i(T) - Σ_reactants |ν_i| S°_i(T)

    Parameters
    ----------
    stoichiometric_coefficients : NumericInput
        Signed stoichiometric coefficients, ν_i.

        Supported input shapes are:

        - scalar
        - 1-D array with shape ``(n_species,)``
        - 2-D array with shape ``(n_states, n_species)``

        For 2-D inputs, each row represents one thermodynamic state or
        reaction evaluation and each column represents a species.

    standard_entropies : NumericInput
        Standard molar entropies, S°_i(T), corresponding element-wise to
        `stoichiometric_coefficients`.

        The input must have the same shape as
        `stoichiometric_coefficients`. All entropy values within a given
        reaction evaluation must correspond to the same temperature and
        standard-state convention.

        A typical unit is J/(mol·K).

    Returns
    -------
    float | NDArray[np.float64]
        Standard reaction entropy, ΔS°_rxn(T).

        For 1-D inputs, a scalar is returned.

        For 2-D inputs with shape ``(n_states, n_species)``, a 1-D array
        with shape ``(n_states,)`` is returned, with one standard reaction
        entropy value for each row.

        The output unit is inherited from `standard_entropies`, typically
        J/(mol·K).

    Raises
    ------
    ValueError
        If either input contains non-finite values, has more than two
        dimensions, or if the two inputs do not have matching shapes.

    Notes
    -----
    The ``_std`` suffix denotes a standard-state thermodynamic quantity.
    It does not imply a temperature of 298.15 K.

    This calculation is valid at any temperature T provided that all
    species standard molar entropies, S°_i(T), refer to that same
    temperature.

    If the supplied values are specifically standard molar entropies at
    298.15 K, then the result is

        ΔS°_rxn(298.15 K) = Σ_i ν_i S°_i(298.15 K)

    commonly written as ΔS°_rxn,298.15.
    """
    # SECTION: Normalize and validate
    # ? Normalize stoichiometric coefficients to a finite float array.
    nu = _as_finite_float_array(
        stoichiometric_coefficients,
        "stoichiometric_coefficients"
    )

    # ? Normalize standard entropies to a finite float array.
    entropy = _as_finite_float_array(standard_entropies, "standard_entropies")

    # >> validate
    _validate_component_arrays(nu, entropy)

    # SECTION: Calculate standard reaction entropy
    return _return_scalar_if_zero_dim(np.sum(nu * entropy, axis=-1))


# ! ::: Calculate reaction entropy using mappings

def _calc_reaction_entropy_std_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int],
    standard_entropies: Mapping[str, float | int],
) -> float:
    """
    Calculate the standard reaction entropy from species-keyed mappings.

    The standard reaction entropy is calculated as

        ΔS°_rxn(T) = Σ_i ν_i S°_i(T)

    where ν_i is the signed stoichiometric coefficient of species i and
    S°_i(T) is its standard molar entropy at temperature T.

    The signed stoichiometric convention is assumed:

        ν_i < 0  for reactants
        ν_i > 0  for products

    Parameters
    ----------
    stoichiometric_coefficients : Mapping[str, float | int]
        Mapping of species identifiers to signed stoichiometric
        coefficients, ν_i.

        The mapping keys identify the species participating in the reaction
        and must exactly match the keys in `standard_entropies`.

    standard_entropies : Mapping[str, float | int]
        Mapping of species identifiers to standard molar entropies,
        S°_i(T).

        All entropy values must correspond to the same temperature and use
        consistent units, typically J/(mol·K).

    Returns
    -------
    float
        Standard reaction entropy, ΔS°_rxn(T), in the same entropy unit
        used by `standard_entropies`, typically J/(mol·K).

    Raises
    ------
    ValueError
        If `stoichiometric_coefficients` and `standard_entropies` do not
        contain the same species keys.

    Notes
    -----
    The ``_std`` suffix denotes a standard-state thermodynamic quantity
    and does not imply a temperature of 298.15 K.

    This calculation is valid at any temperature T provided that all
    species standard molar entropies, S°_i(T), correspond to that same
    temperature.

    If the supplied entropy values are specifically evaluated at
    298.15 K, the result is the standard reaction entropy at 298.15 K:

        ΔS°_rxn,298.15 = Σ_i ν_i S°_i,298.15

    This function validates the species keys and delegates the numerical
    calculation to `_calc_reaction_entropy_std`.
    """
    # >> validate
    _validate_same_mapping_keys(
        stoichiometric_coefficients,
        standard_entropies,
        "stoichiometric_coefficients",
        "standard_entropies",
    )

    # >> calculate reaction entropy from mappings
    return float(
        _calc_reaction_entropy_std(
            [stoichiometric_coefficients[key]
                for key in stoichiometric_coefficients],
            [standard_entropies[key] for key in stoichiometric_coefficients],
        )
    )

# ! ::: Calculate reaction entropy with props (CustomProp)


def _calc_reaction_entropy_std_from_props(
    stoichiometric_coefficients: Mapping[str, CustomProp],
    standard_entropies: Mapping[str, CustomProp],
    output_unit: str = 'J/mol.K',
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """
    Calculate the standard reaction entropy from species standard molar
    entropies and stoichiometric coefficients.

    The standard reaction entropy is calculated as

        ΔS°_rxn(T) = Σ_i ν_i S°_i(T)

    where ν_i is the signed stoichiometric coefficient of species i and
    S°_i(T) is its standard molar entropy at temperature T.

    With the signed stoichiometric convention:

        ν_i < 0  for reactants
        ν_i > 0  for products

    Parameters
    ----------
    stoichiometric_coefficients : Mapping[str, CustomProp]
        Stoichiometric coefficients of the reaction, keyed by species.
    standard_entropies : Mapping[str, CustomProp]
        Standard entropies of the species/components, keyed by species.
    output_unit : str, optional
        The unit to which the standard entropies should be converted. Default is 'J/mol.K'.
    components : Optional[Sequence[Component]], optional
        The list of component metadata. Default is None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping components. Default is None.
    case_sensitive : bool, optional
        Whether the component key matching should be case-sensitive. Default is True.
    sort_by_components_order : bool, optional
        Whether to sort the output by the order of components. Default is True.

    Returns
    -------
    float
        Standard reaction entropy calculated from the CustomProp mappings.

    Notes
    -----
    The `_std` suffix denotes standard-state thermodynamic properties and
    does not imply a temperature of 298.15 K.

    This calculation is valid at any temperature T provided that all
    species standard molar entropies, S°_i(T), are evaluated at the same
    temperature.

    If the supplied entropy values are specifically tabulated at
    298.15 K, the result is the standard reaction entropy at 298.15 K:

        ΔS°_rxn,298.15 = Σ_i ν_i S°_i,298.15
    """
    # SECTION: Convert CustomProp mappings to values
    stoichiometric_coefficients_dict = _to_values(
        data=stoichiometric_coefficients,
        name="stoichiometric_coefficients",
    )
    # ? convert to output unit
    standard_entropies_dict = _to_values(
        data=standard_entropies,
        name="standard_entropies",
        output_unit=output_unit,
    )

    # SECTION: Assign converted values to the dictionaries used for further processing
    if (
        stoichiometric_coefficients_dict is not None and
        components is not None and
        component_key is not None
    ):
        stoichiometric_coefficients_dict = _configure_component_values(
            values=stoichiometric_coefficients,
            components=list(components),
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order,
            name="stoichiometric_coefficients",
        )
        standard_entropies_dict = _configure_component_values(
            values=standard_entropies,
            components=list(components),
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order,
            name="standard_entropies",
        )

    return _calc_reaction_entropy_std_from_mapping(
        stoichiometric_coefficients=stoichiometric_coefficients_dict,
        standard_entropies=standard_entropies_dict,
    )

# ! ::: Entropy from enthalpy and Gibbs energy


def _calc_reaction_entropy_std_from_enthalpy_gibbs(
    delta_h_reaction_std: NumericInput,
    delta_g_reaction_std: NumericInput,
    temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """
    Calculate the standard reaction entropy at temperature T.

    ΔS°_rxn(T) = [ΔH°_rxn(T) - ΔG°_rxn(T)] / T

    Parameters
    ----------
    delta_h_reaction_std : NumericInput
        Standard reaction enthalpy at the specified temperature.
    delta_g_reaction_std : NumericInput
        Standard reaction Gibbs energy at the specified temperature.
    temperature : NumericInput
        Temperature at which the standard reaction entropy is calculated.

    Returns
    -------
    float | NDArray[np.float64]
        Standard reaction entropy, ΔS°_rxn(T). The resulting entropy unit
        is determined by the energy unit of the input properties divided
        by kelvin. For example, if ΔH° and ΔG° are given in J/mol, the
        result is in J/(mol·K).

    Notes
    -----
    The superscript degree symbol (°), represented here by the `_std`
    suffix, denotes a standard-state thermodynamic property and does not
    imply a temperature of 298.15 K.

    Therefore, this relation is valid at any temperature T provided that
    ΔH°_rxn(T) and ΔG°_rxn(T) are evaluated at that same temperature.

    If the supplied properties specifically correspond to 298.15 K, the
    calculated quantity is ΔS°_rxn(298.15 K), which may be denoted
    ΔS°_rxn,298.15.

    """
    # SECTION: Normalize and validate
    dh = _as_finite_float_array(delta_h_reaction_std, "delta_h_reaction_std")
    dg = _as_finite_float_array(delta_g_reaction_std, "delta_g_reaction_std")
    t = _as_finite_float_array(temperature, "temperature")
    _validate_positive_array(t, "temperature")

    # SECTION: Calculate reaction entropy
    return _return_scalar_if_zero_dim((dh - dg) / t)


# SECTION: Reaction heat-capacity change

def _calc_reaction_heat_capacity_change(
    stoichiometric_coefficients: NumericInput,
    component_heat_capacities: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reaction heat-capacity change: delta_Cp_r = sum_i(nu_i*Cp_i)."""
    # SECTION: Normalize and validate
    nu = _as_finite_float_array(stoichiometric_coefficients, "stoichiometric_coefficients")
    cp = _as_finite_float_array(component_heat_capacities, "component_heat_capacities")
    _validate_component_arrays(nu, cp)
    _validate_positive_array(cp, "component_heat_capacities")

    # SECTION: Calculate reaction heat-capacity change
    return _return_scalar_if_zero_dim(np.sum(nu * cp, axis=-1))


def _calc_reaction_heat_capacity_change_from_mapping(
    stoichiometric_coefficients: Mapping[str, float | int],
    component_heat_capacities: Mapping[str, float | int],
) -> float:
    """Calculate reaction heat-capacity change from aligned keyed inputs."""
    _validate_same_mapping_keys(
        stoichiometric_coefficients,
        component_heat_capacities,
        "stoichiometric_coefficients",
        "component_heat_capacities",
    )
    keys = list(stoichiometric_coefficients)
    return float(
        _calc_reaction_heat_capacity_change(
            [stoichiometric_coefficients[key] for key in keys],
            [component_heat_capacities[key] for key in keys],
        )
    )


def _calc_reaction_heat_capacity_change_from_props(
    stoichiometric_coefficients: Mapping[str, CustomProp],
    component_heat_capacities: Mapping[str, CustomProp],
    output_unit: str = "J/mol.K",
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate reaction heat-capacity change from unit-aware keyed inputs."""
    nu = _to_values(
        data=stoichiometric_coefficients,
        name="stoichiometric_coefficients",
    )
    cp = _to_values(
        data=component_heat_capacities,
        name="component_heat_capacities",
        output_unit=output_unit,
    )

    if components is not None and component_key is not None:
        nu = _configure_component_values(
            values=nu,
            components=list(components),
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order,
            name="stoichiometric_coefficients",
        )
        cp = _configure_component_values(
            values=cp,
            components=list(components),
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order,
            name="component_heat_capacities",
        )

    return _calc_reaction_heat_capacity_change_from_mapping(nu, cp)


# SECTION: Kirchhoff enthalpy correction

def _calc_reaction_enthalpy_from_constant_delta_cp(
    delta_h_reaction_ref: NumericInput,
    delta_cp_reaction: NumericInput,
    temperature: NumericInput,
    reference_temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reaction enthalpy using constant reaction heat-capacity change."""
    # SECTION: Normalize and validate
    dh_ref = _as_finite_float_array(delta_h_reaction_ref, "delta_h_reaction_ref")
    dcp = _as_finite_float_array(delta_cp_reaction, "delta_cp_reaction")
    t = _as_finite_float_array(temperature, "temperature")
    t_ref = _as_finite_float_array(reference_temperature, "reference_temperature")
    _validate_positive_array(t, "temperature")
    _validate_positive_array(t_ref, "reference_temperature")

    # SECTION: Calculate Kirchhoff correction
    return _return_scalar_if_zero_dim(dh_ref + dcp * (t - t_ref))


# SECTION: Reaction heat source/rate

def _calc_reaction_volumetric_heat_source(
    reaction_enthalpies: NumericInput,
    reaction_rates: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate generated reaction heat source: q''' = -sum_j(delta_H_j*r_j)."""
    dh = _as_finite_float_array(reaction_enthalpies, "reaction_enthalpies")
    rates = _as_finite_float_array(reaction_rates, "reaction_rates")
    if dh.shape != rates.shape:
        raise ValueError("reaction_enthalpies and reaction_rates must have the same shape.")
    return _return_scalar_if_zero_dim(-np.sum(dh * rates, axis=-1))


def _calc_reaction_heat_rate(
    reaction_enthalpies: NumericInput,
    reaction_rates: NumericInput,
    volume: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate total generated reaction heat rate: Qdot = V*q'''."""
    q_source = np.asarray(
        _calc_reaction_volumetric_heat_source(reaction_enthalpies, reaction_rates),
        dtype=np.float64,
    )
    v = _as_finite_float_array(volume, "volume")
    _validate_positive_array(v, "volume")
    return _return_scalar_if_zero_dim(v * q_source)


# SECTION: Core exports
__all__ = [
    "_calc_reaction_entropy_std",
    "_calc_reaction_entropy_std_from_mapping",
    "_calc_reaction_entropy_std_from_props",
    "_calc_reaction_entropy_std_from_enthalpy_gibbs",
    "_calc_reaction_heat_capacity_change",
    "_calc_reaction_heat_capacity_change_from_mapping",
    "_calc_reaction_heat_capacity_change_from_props",
    "_calc_reaction_enthalpy_from_constant_delta_cp",
    "_calc_reaction_volumetric_heat_source",
    "_calc_reaction_heat_rate",
]
