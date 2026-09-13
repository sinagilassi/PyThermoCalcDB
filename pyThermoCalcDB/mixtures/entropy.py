"""Ideal entropy-of-mixing helpers."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, Component, ComponentKey, ScalarValue
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos, to_dict, to_list
from pythermodb_settings.utils.validators import fractions
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import (
    _all_custom_props,
    _configure_component_values,
    _resolve_unit_conversion_fn,
)
from .core.entropy import (
    _calc_ideal_entropy_of_mixing,
    _calc_ideal_entropy_of_mixing_from_props,
    _calc_ideal_molar_entropy_of_mixing,
    _calc_ideal_molar_entropy_of_mixing_from_mapping,
    _calc_ideal_molar_entropy_of_mixing_from_props,
)


# ! ::: Ideal molar entropy of mixing from mapping
def calc_ideal_molar_entropy_of_mixing_from_mapping(
    mole_fractions: Mapping[str, float | int],
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal molar entropy of mixing.

    Parameters
    ----------
    mole_fractions : mapping or sequence of float | int
        Component mole fractions. Zero fractions are allowed and contribute
        zero through the ``x*ln(x)`` limiting behavior.
    gas_constant : float, optional
        Gas constant in entropy units per mol per K. Defaults to
        ``8.314462618`` J/mol/K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function used by shared quantity helpers.

    Returns
    -------
    float
        Ideal molar entropy of mixing, typically J/mol/K.

    Notes
    -----
    Equation: ``delta_S_mix = -R*sum_i(x_i*ln(x_i))``. The mixture is assumed
    ideal and no excess entropy term is included.
    """
    # SECTION: Validate inputs
    fractions(mole_fractions, "mole_fractions")
    r = pos(gas_constant, "gas_constant")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Normalize mixed/numeric mapping inputs
    x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)

    return _calc_ideal_molar_entropy_of_mixing_from_mapping(x, r)


# ! ::: Ideal molar entropy of mixing from mapping

def calc_ideal_molar_entropy_of_mixing_from_props(
    mole_fractions: Mapping[str, CustomProp],
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal molar entropy of mixing.

    Parameters
    ----------
    mole_fractions : mapping of CustomProp
        Component mole fractions. Zero fractions are allowed and contribute
        zero through the ``x*ln(x)`` limiting behavior.
    gas_constant : float, optional
        Gas constant in entropy units per mol per K. Defaults to
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
        Ideal molar entropy of mixing, typically J/mol/K.

    Notes
    -----
    Equation: ``delta_S_mix = -R*sum_i(x_i*ln(x_i))``. The mixture is assumed
    ideal and no excess entropy term is included.
    """
    # SECTION: Validate inputs
    fractions(mole_fractions, "mole_fractions")
    r = pos(gas_constant, "gas_constant")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Normalize composition input
    if isinstance(mole_fractions, Mapping):
        if _all_custom_props(mole_fractions):
            return _calc_ideal_molar_entropy_of_mixing_from_props(
                mole_fractions,
                r,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )
        else:
            raise ValueError(
                "Mole fractions must be all CustomProp instances when provided as a mapping.")
    else:
        raise ValueError("Mole fractions must be provided as a mapping.")


# SECTION: Total ideal entropy of mixing

def calc_ideal_entropy_of_mixing(
    total_moles: ScalarValue,
    mole_fractions: Mapping[str, CustomProp] | Sequence[float | int | CustomProp],
    gas_constant: float = R_J_molK,
    output_total_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total ideal entropy of mixing for a specified amount.

    Parameters
    ----------
    total_moles : float | int | CustomProp
        Total amount of mixture. Must be positive.
    mole_fractions : mapping or sequence of float | int | CustomProp
        Component mole fractions.
    gas_constant : float, optional
        Gas constant in entropy units per mol per K.
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
        Total ideal entropy of mixing, typically J/K.

    Notes
    -----
    Equation: ``delta_S_mix,total = n_total*delta_S_mix,molar``.
    """
    # SECTION: Normalize total amount
    r = pos(gas_constant, "gas_constant")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    if isinstance(mole_fractions, Mapping):
        if isinstance(total_moles, CustomProp) and _all_custom_props(mole_fractions):
            return _calc_ideal_entropy_of_mixing_from_props(
                total_moles,
                mole_fractions,
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
        x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)
        x = _configure_component_values(
            x, components, component_key, case_sensitive, sort_by_components_order, "mole_fractions"
        )
        return float(_calc_ideal_entropy_of_mixing(n_total, list(x.values()), r))

    # SECTION: Scale molar entropy by total moles
    n_total = pos(
        total_moles,
        "total_moles",
        output_total_moles_unit,
        unit_conversion_fn=conversion_fn,
    )
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    return float(_calc_ideal_entropy_of_mixing(n_total, x, r))

# ! ::: ideal entropy of mixing using sequence


def calc_ideal_entropy_of_mixing_from_sequence(
    total_moles: float | int | CustomProp,
    mole_fractions: Sequence[float | int | CustomProp],
    gas_constant: float = R_J_molK,
    output_total_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate total ideal entropy of mixing for a specified amount.

    Parameters
    ----------
    total_moles : float | int | CustomProp
        Total amount of mixture. Must be positive.
    mole_fractions : mapping or sequence of float | int | CustomProp
        Component mole fractions.
    gas_constant : float, optional
        Gas constant in entropy units per mol per K.
    output_total_moles_unit : str, optional
        Unit used to normalize ``total_moles`` when supplied as ``CustomProp``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Total ideal entropy of mixing, typically J/K.

    Notes
    -----
    - Equation: ``delta_S_mix,total = n_total*delta_S_mix,molar``.
    - Unit conversion is applied to ``total_moles`` when supplied as ``CustomProp`` using the specified ``unit_conversion_fn``.
    """
    # SECTION: Normalize total amount
    r = pos(gas_constant, "gas_constant")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Scale molar entropy by total moles
    n_total = pos(
        total_moles,
        "total_moles",
        output_total_moles_unit,
        unit_conversion_fn=conversion_fn,
    )
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    return float(_calc_ideal_entropy_of_mixing(n_total, x, r))

# ! ::: ideal entropy of mixing using mapping


def calc_ideal_entropy_of_mixing_from_mapping(
    total_moles: float | int,
    mole_fractions: Mapping[str, float | int],
    gas_constant: float = R_J_molK,
) -> float:
    """Calculate total ideal entropy of mixing for a specified amount.

    Parameters
    ----------
    total_moles : float | int
        Total amount of mixture. Must be positive.
    mole_fractions : mapping of float | int
        Component mole fractions.
    gas_constant : float, optional
        Gas constant in entropy units per mol per K.

    Returns
    -------
    float
        Total ideal entropy of mixing, typically J/K.

    Notes
    -----
    Equation: ``delta_S_mix,total = n_total*delta_S_mix,molar``.
    """
    # SECTION: Normalize total amount
    r = pos(gas_constant, "gas_constant")

    # SECTION: Normalize mixed/numeric mapping inputs
    n_total = pos(
        total_moles,
        "total_moles",
    )
    x = to_dict(mole_fractions)

    return float(_calc_ideal_entropy_of_mixing(n_total, list(x.values()), r))


# ! ::: ideal entropy of mixing using props


def calc_ideal_entropy_of_mixing_from_props(
    total_moles: CustomProp,
    mole_fractions: Mapping[str, CustomProp],
    gas_constant: float = R_J_molK,
    output_total_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total ideal entropy of mixing for a specified amount.

    Parameters
    ----------
    total_moles : float | int | CustomProp
        Total amount of mixture. Must be positive.
    mole_fractions : mapping of float | int | CustomProp
        Component mole fractions.
    gas_constant : float, optional
        Gas constant in entropy units per mol per K.
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
        Total ideal entropy of mixing, typically J/K.

    Notes
    -----
    Equation: ``delta_S_mix,total = n_total*delta_S_mix,molar``.
    """
    # SECTION: Normalize total amount
    r = pos(gas_constant, "gas_constant")
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    if isinstance(mole_fractions, Mapping):
        if isinstance(total_moles, CustomProp) and _all_custom_props(mole_fractions):
            return _calc_ideal_entropy_of_mixing_from_props(
                total_moles,
                mole_fractions,
                r,
                output_total_moles_unit,
                conversion_fn,
                components,
                component_key,
                case_sensitive,
                sort_by_components_order,
            )
        else:
            raise TypeError(
                "mole_fractions must be a mapping of CustomProp instances.")
    else:
        raise TypeError(
            "mole_fractions must be a mapping of CustomProp instances.")


# SECTION: Public exports
__all__ = [
    "calc_ideal_entropy_of_mixing_from_mapping",
    "calc_ideal_entropy_of_mixing_from_props",
    "calc_ideal_entropy_of_mixing",
    "calc_ideal_entropy_of_mixing_from_sequence",
    "calc_ideal_entropy_of_mixing_from_mapping",
    "calc_ideal_entropy_of_mixing_from_props",
]
