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
    standard_entropies: NDArray[np.float64],
) -> None:
    """Validate pairwise component arrays for reaction entropy calculations."""
    if stoichiometric_coefficients.ndim not in (1, 2):
        raise ValueError(
            "stoichiometric_coefficients must be a 1-D or 2-D array.")
    if standard_entropies.ndim not in (1, 2):
        raise ValueError("standard_entropies must be a 1-D or 2-D array.")
    if stoichiometric_coefficients.shape != standard_entropies.shape:
        raise ValueError(
            "stoichiometric_coefficients and standard_entropies must have the same shape.")


# ! ::: Standard reaction entropy

def _calc_reaction_entropy_std(
    stoichiometric_coefficients: NumericInput,
    standard_entropies: NumericInput,
) -> float | NDArray[np.float64]:
    """
    Calculate the entropy of a reaction at standard conditions.
        Ent_RXN_STD = sum_i(nu_i*S_i_std).

    For 2-D inputs, axis 0 is states and axis 1 is species/components. Both must obey the same order as the input mappings.

    Parameters
    ----------
    stoichiometric_coefficients : NumericInput
        Stoichiometric coefficients of the reaction.
    standard_entropies : NumericInput
        Standard entropies of the species/components.

    Returns
    -------
    float | NDArray[np.float64]
        Standard reaction entropy.
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
    Calculate standard reaction entropy from keyed species values.

    Parameters
    ----------
    stoichiometric_coefficients : Mapping[str, float | int]
        Stoichiometric coefficients of the reaction, keyed by species.
    standard_entropies : Mapping[str, float | int]
        Standard entropies of the species/components, keyed by species.

    Returns
    -------
    float
        Standard reaction entropy calculated from the keyed species values.
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

# ! ::: Calculate reaction entropy using CustomProp


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


# SECTION: Core exports
__all__ = [
    "_calc_reaction_entropy_std",
    "_calc_reaction_entropy_std_from_mapping",
    "_calc_reaction_entropy_std_from_enthalpy_gibbs",
]
