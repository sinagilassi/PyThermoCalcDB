"""Ideal Gibbs-energy-of-mixing helpers."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, Component, ComponentKey, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos, to_dict, to_list
from pythermodb_settings.utils.validators import fractions
from pycuc.canonical import to_K
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import (
    _all_custom_props,
    _get_all_custom_props,
    _configure_component_values,
    _resolve_unit_conversion_fn,
)
from .core.gibbs import (
    _calc_ideal_gibbs_energy_of_mixing,
    _calc_ideal_gibbs_energy_of_mixing_from_props,
    _calc_ideal_molar_gibbs_energy_of_mixing,
    _calc_ideal_molar_gibbs_energy_of_mixing_from_mapping,
    _calc_ideal_molar_gibbs_energy_of_mixing_from_props,
)


# SECTION: Ideal molar Gibbs energy of mixing
# ! ::: from mapping or sequence
def calc_ideal_molar_gibbs_energy_of_mixing_from_sequence(
    mole_fractions: Sequence[float | int],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
) -> float:
    """Calculate ideal molar Gibbs energy of mixing.

    Parameters
    ----------
    mole_fractions : mapping or sequence of float | int | CustomProp
        Component mole fractions. Zero fractions are allowed and contribute
        zero through the ``x*ln(x)`` limiting behavior.
    temperature : Temperature
        Mixture temperature. Converted to K before calculation.
    gas_constant : float, optional
        Gas constant in energy units per mol per K. Defaults to
        ``8.314462618`` J/mol/K.

    Returns
    -------
    float
        Ideal molar Gibbs energy of mixing, typically J/mol.

    Notes
    -----
    Equation: ``delta_G_mix = R*T*sum_i(x_i*ln(x_i))``. The mixture is assumed
    ideal, so no excess Gibbs energy term is included.
    """
    # SECTION: Validate inputs
    fractions(mole_fractions, "mole_fractions")
    r = pos(gas_constant, "gas_constant")

    # SECTION: Calculate ideal molar Gibbs energy of mixing
    temperature_k = to_K(temperature.value, temperature.unit)
    x = to_list(mole_fractions)
    return float(_calc_ideal_molar_gibbs_energy_of_mixing(x, temperature_k, r))

# ! ::: from mapping


def calc_ideal_molar_gibbs_energy_of_mixing_from_mapping(
    mole_fractions: Mapping[str, float | int],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
) -> float:
    """Calculate ideal molar Gibbs energy of mixing.

    Parameters
    ----------
    mole_fractions : mapping of float | int
        Component mole fractions. Zero fractions are allowed and contribute
        zero through the ``x*ln(x)`` limiting behavior.
    temperature : Temperature
        Mixture temperature. Converted to K before calculation.
    gas_constant : float, optional
        Gas constant in energy units per mol per K. Defaults to
        ``8.314462618`` J/mol/K.

    Returns
    -------
    float
        Ideal molar Gibbs energy of mixing, typically J/mol.

    Notes
    -----
    Equation: ``delta_G_mix = R*T*sum_i(x_i*ln(x_i))``. The mixture is assumed
    ideal, so no excess Gibbs energy term is included.
    """
    # SECTION: Validate inputs
    fractions(mole_fractions, "mole_fractions")
    r = pos(gas_constant, "gas_constant")

    # SECTION: Normalize mixed/numeric mapping inputs
    temperature_k = to_K(temperature.value, temperature.unit)
    x = to_dict(mole_fractions)
    return _calc_ideal_molar_gibbs_energy_of_mixing_from_mapping(x, temperature_k, r)


# ! ::: from props


def calc_ideal_molar_gibbs_energy_of_mixing_from_props(
    mole_fractions: Mapping[str, CustomProp],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal molar Gibbs energy of mixing.

    Parameters
    ----------
    mole_fractions : mapping of float | int | CustomProp
        Component mole fractions. Zero fractions are allowed and contribute
        zero through the ``x*ln(x)`` limiting behavior.
    temperature : Temperature
        Mixture temperature. Converted to K before calculation.
    gas_constant : float, optional
        Gas constant in energy units per mol per K. Defaults to
        ``8.314462618`` J/mol/K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function used by shared quantity helpers.
    components : list[Component], optional
        Component metadata used for mapping-key remapping when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format used for mapping inputs.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether mapping values should follow ``components`` order.

    Returns
    -------
    float
        Ideal molar Gibbs energy of mixing, typically J/mol.

    Notes
    -----
    Equation: ``delta_G_mix = R*T*sum_i(x_i*ln(x_i))``. The mixture is assumed
    ideal, so no excess Gibbs energy term is included.
    """
    # SECTION: Validate inputs
    fractions(mole_fractions, "mole_fractions")
    r = pos(gas_constant, "gas_constant")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    return _calc_ideal_molar_gibbs_energy_of_mixing_from_props(
        mole_fractions,
        temperature,
        r,
        conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


# SECTION: Total ideal Gibbs energy of mixing

# ! form sequence

def calc_ideal_gibbs_energy_of_mixing_from_alls(
    total_moles: ScalarValue,
    mole_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
    output_total_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total ideal Gibbs energy of mixing for a specified amount.

    Parameters
    ----------
    total_moles : float | int | CustomProp
        Total amount of mixture. Must be positive.
    mole_fractions : mapping or sequence of float | int | CustomProp
        Component mole fractions.
    temperature : Temperature
        Mixture temperature. Converted to K before calculation.
    gas_constant : float, optional
        Gas constant in energy units per mol per K.
    output_total_moles_unit : str, optional
        Unit used to normalize ``total_moles`` when supplied as ``CustomProp``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.
    components : list[Component], optional
        Component metadata used for mapping-key remapping.
    component_key : ComponentKey, optional
        Component identifier format used for mapping inputs.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether mapping values should follow ``components`` order.

    Returns
    -------
    float
        Total ideal Gibbs energy of mixing, typically J.

    Notes
    -----
    Equation: ``delta_G_mix,total = n_total*delta_G_mix,molar``.
    """
    # SECTION: Normalize total amount
    r = pos(gas_constant, "gas_constant")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    if isinstance(mole_fractions, Mapping):
        if isinstance(total_moles, CustomProp) and _all_custom_props(mole_fractions):
            # get all CustomProp instances from the mapping
            custom_props_all = _get_all_custom_props(
                mole_fractions, return_type="mapping"
            )
            return _calc_ideal_gibbs_energy_of_mixing_from_props(
                total_moles,
                custom_props_all,
                temperature,
                r,
                output_total_moles_unit,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )

        # SECTION: Normalize mixed/numeric mapping inputs
        n_total = pos(
            total_moles,
            "total_moles",
            output_total_moles_unit,
            unit_conversion_fn=conversion_fn,
        )
        temperature_k = to_K(temperature.value, temperature.unit)
        x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)
        x = _configure_component_values(
            x, components, component_key, case_sensitive, sort_by_components_order, "mole_fractions"
        )
        return float(_calc_ideal_gibbs_energy_of_mixing(n_total, list(x.values()), temperature_k, r))

    # SECTION: Scale molar Gibbs energy by total moles
    n_total = pos(
        total_moles,
        "total_moles",
        output_total_moles_unit,
        unit_conversion_fn=conversion_fn,
    )
    temperature_k = to_K(temperature.value, temperature.unit)
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    return float(_calc_ideal_gibbs_energy_of_mixing(n_total, x, temperature_k, r))

# ! ::: from sequence


def calc_ideal_gibbs_energy_of_mixing_from_sequence(
    total_moles: ScalarValue,
    mole_fractions: Sequence[float | int | CustomProp],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
    output_total_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate total ideal Gibbs energy of mixing from sequence inputs."""
    return calc_ideal_gibbs_energy_of_mixing_from_alls(
        total_moles,
        mole_fractions,
        temperature,
        gas_constant,
        output_total_moles_unit,
        unit_conversion_fn,
    )

# ! ::: from mapping


def calc_ideal_gibbs_energy_of_mixing_from_mapping(
    total_moles: float | int,
    mole_fractions: Mapping[str, float | int],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
) -> float:
    """Calculate total ideal Gibbs energy of mixing from numeric mapping inputs."""
    return calc_ideal_gibbs_energy_of_mixing_from_alls(
        total_moles,
        mole_fractions,
        temperature,
        gas_constant,
    )

# ! ::: from props


def calc_ideal_gibbs_energy_of_mixing_from_props(
    total_moles: CustomProp,
    mole_fractions: Mapping[str, CustomProp],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
    output_total_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total ideal Gibbs energy of mixing from unit-aware mapping inputs."""
    return calc_ideal_gibbs_energy_of_mixing_from_alls(
        total_moles,
        mole_fractions,
        temperature,
        gas_constant,
        output_total_moles_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


# SECTION: Backwards-compatible aliases
calc_ideal_molar_gibbs_energy_of_mixing = calc_ideal_molar_gibbs_energy_of_mixing_from_sequence
calc_ideal_gibbs_energy_of_mixing_from_all = calc_ideal_gibbs_energy_of_mixing_from_alls
calc_ideal_gibbs_energy_of_mixing = calc_ideal_gibbs_energy_of_mixing_from_all


# SECTION: Public exports
__all__ = [
    "calc_ideal_molar_gibbs_energy_of_mixing_from_sequence",
    "calc_ideal_molar_gibbs_energy_of_mixing_from_mapping",
    "calc_ideal_molar_gibbs_energy_of_mixing_from_props",
    "calc_ideal_molar_gibbs_energy_of_mixing",
    "calc_ideal_gibbs_energy_of_mixing_from_alls",
    "calc_ideal_gibbs_energy_of_mixing_from_all",
    "calc_ideal_gibbs_energy_of_mixing_from_sequence",
    "calc_ideal_gibbs_energy_of_mixing_from_mapping",
    "calc_ideal_gibbs_energy_of_mixing_from_props",
    "calc_ideal_gibbs_energy_of_mixing",
]
