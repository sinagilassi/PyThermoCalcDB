"""Partial-molar property relations."""

# import libs
from collections.abc import Mapping, Sequence

from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import fractions, non_negative, same_shape
# locals
from ..utils.conversions import _resolve_unit_conversion_fn
from .core.partial_molar import (
    _calc_binary_partial_molar_properties,
    _calc_molar_property_from_partial_molar_mapping,
    _calc_molar_property_from_partial_molar_properties,
    _calc_total_property_from_partial_molar_mapping,
    _calc_total_property_from_partial_molar_properties,
)


# SECTION: Public sequence APIs

def calc_total_property_from_partial_molar_properties(
    amounts: Sequence[float | int | CustomProp],
    partial_molar_properties: Sequence[float | int | CustomProp],
    output_amount_unit: str | None = None,
    output_partial_property_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate an extensive property from amounts and partial molar properties."""
    non_negative(amounts, "amounts")
    same_shape(amounts, partial_molar_properties)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    n = to_list(amounts, output_amount_unit, unit_conversion_fn=conversion_fn)
    partials = to_list(
        partial_molar_properties,
        output_partial_property_unit,
        unit_conversion_fn=conversion_fn,
    )
    return float(_calc_total_property_from_partial_molar_properties(n, partials))


def calc_molar_property_from_partial_molar_properties(
    mole_fractions: Sequence[float | int | CustomProp],
    partial_molar_properties: Sequence[float | int | CustomProp],
    output_partial_property_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate a molar mixture property from mole fractions and partial properties."""
    fractions(mole_fractions, "mole_fractions")
    same_shape(mole_fractions, partial_molar_properties)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    partials = to_list(
        partial_molar_properties,
        output_partial_property_unit,
        unit_conversion_fn=conversion_fn,
    )
    return float(_calc_molar_property_from_partial_molar_properties(x, partials))


def calc_binary_partial_molar_properties(
    molar_property: float | int,
    mole_fraction_1: float | int,
    dmolar_property_dx1: float | int,
) -> tuple[float, float]:
    """Calculate binary partial molar properties from ``M`` and ``dM/dx1``."""
    left, right = _calc_binary_partial_molar_properties(
        molar_property,
        mole_fraction_1,
        dmolar_property_dx1,
    )
    return float(left), float(right)


# SECTION: Public mapping APIs

def calc_total_property_from_partial_molar_mapping(
    amounts: Mapping[str, float | int | CustomProp],
    partial_molar_properties: Mapping[str, float | int | CustomProp],
    output_amount_unit: str | None = None,
    output_partial_property_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate an extensive property from aligned component mappings."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    n = to_dict(amounts, output_amount_unit, unit_conversion_fn=conversion_fn)
    partials = to_dict(
        partial_molar_properties,
        output_partial_property_unit,
        unit_conversion_fn=conversion_fn,
    )
    return _calc_total_property_from_partial_molar_mapping(n, partials)


def calc_molar_property_from_partial_molar_mapping(
    mole_fractions: Mapping[str, float | int | CustomProp],
    partial_molar_properties: Mapping[str, float | int | CustomProp],
    output_partial_property_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate a molar property from aligned component mappings."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)
    partials = to_dict(
        partial_molar_properties,
        output_partial_property_unit,
        unit_conversion_fn=conversion_fn,
    )
    return _calc_molar_property_from_partial_molar_mapping(x, partials)


# SECTION: Public exports
__all__ = [
    "calc_total_property_from_partial_molar_properties",
    "calc_molar_property_from_partial_molar_properties",
    "calc_binary_partial_molar_properties",
    "calc_total_property_from_partial_molar_mapping",
    "calc_molar_property_from_partial_molar_mapping",
]
