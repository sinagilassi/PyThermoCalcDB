"""Ideal entropy-of-mixing helpers."""

# import libs
import math
from collections.abc import Mapping, Sequence
from typing import Optional, List

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, Component, ComponentKey, ScalarValue
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.quantity import pos, to_dict, to_list
from pythermodb_settings.utils.validators import fractions
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import _resolve_unit_conversion_fn


# SECTION: Internal helpers

def _configure_component_values(
    values: Mapping[str, float],
    components: Optional[List[Component]],
    component_key: Optional[ComponentKey],
    case_sensitive: bool,
    sort_by_components_order: bool,
    name: str,
) -> dict[str, float]:
    """Remap and order mapping values using component metadata when requested."""
    # ! When no component key is requested, preserve caller mapping keys/order.
    if component_key is None:
        return dict(values)

    # NOTE: Component metadata is required only for key remapping.
    if not components:
        raise ValueError(f"component_key is provided but components is empty for {name}.")

    # SECTION: Remap values through pythermodb-settings utilities
    component_values = config_components_values(
        values=dict(values),
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    if component_values is None:
        raise ValueError(f"Failed to configure {name} component values.")
    component_values_dict, _ = component_values
    return component_values_dict


def _x_log_x_sum(values: Mapping[str, float] | Sequence[float]) -> float:
    """Return sum(x_i*ln(x_i)), using the x*ln(x) -> 0 limit at x = 0."""
    # NOTE: Skip zero fractions to avoid log(0).
    iterable = values.values() if isinstance(values, Mapping) else values
    return sum(value * math.log(value) for value in iterable if value > 0.0)


# SECTION: Ideal molar entropy of mixing

def calc_ideal_molar_entropy_of_mixing(
    mole_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
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
    mole_fractions : mapping or sequence of float | int | CustomProp
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
        x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)
        x = _configure_component_values(
            x, components, component_key, case_sensitive, sort_by_components_order, "mole_fractions"
        )
    else:
        x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)

    # SECTION: Calculate ideal molar entropy of mixing
    return -r * _x_log_x_sum(x)


# SECTION: Total ideal entropy of mixing

def calc_ideal_entropy_of_mixing(
    total_moles: ScalarValue,
    mole_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
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
    n_total = pos(
        total_moles,
        "total_moles",
        output_total_moles_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )

    # SECTION: Scale molar entropy by total moles
    return n_total * calc_ideal_molar_entropy_of_mixing(
        mole_fractions,
        gas_constant,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


# SECTION: Public exports
__all__ = [
    "calc_ideal_molar_entropy_of_mixing",
    "calc_ideal_entropy_of_mixing",
]

