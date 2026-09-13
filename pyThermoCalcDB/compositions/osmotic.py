"""Osmotic composition primitives."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp, Temperature
from pythermodb_settings.utils import to_annotated_value
from pythermodb_settings.utils.quantity import to_dict, to_list

# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import (
    _configure_component_values,
    _pos,
    _resolve_unit_conversion_fn,
    _to_kelvin,
    _to_values,
    _validate_custom_prop_mapping,
)
from .core.osmotic import (
    _calc_ideal_osmotic_pressure,
    _calc_osmolality,
    _calc_osmolality_from_mapping,
    _calc_osmolarity,
    _calc_osmolarity_from_mapping,
)


# SECTION: Osmolarity
@calculation_info(
    name="osmolarity",
    description="Calculate osmolarity by summing supplied dissolved particle molarities.",
    equation="C_osm = sum_i C_i",
    inputs={"species_molarities": "Molarities of explicitly supplied dissolved particle species."},
    outputs={"osmolarity": "Total osmolarity."},
    aliases=("osmotic concentration", "osmolar concentration"),
    notes=("No dissociation or speciation is inferred.",),
    tags=("osmolarity", "composition", "low_level"),
)
def calc_osmolarity(
    species_molarities: Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "osmolarity",
    description: str = "Calculate osmolarity from supplied particle molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float | NDArray[np.float64]]:
    """Calculate osmolarity from explicitly supplied particle molarities."""
    return to_annotated_value(
        _calc_osmolarity(species_molarities),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_osmolarity",
    )


def calc_osmolarity_from_sequence(
    species_molarities: Sequence[float | int],
    *,
    name: str = "osmolarity",
    description: str = "Calculate osmolarity from supplied particle molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate osmolarity from a sequence of supplied particle molarities."""
    return to_annotated_value(
        float(_calc_osmolarity(to_list(species_molarities))),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_osmolarity",
    )


def calc_osmolarity_from_mapping(
    species_molarities: Mapping[str, float | int | CustomProp],
    *,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    name: str = "osmolarity",
    description: str = "Calculate osmolarity from supplied particle molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate osmolarity from a species-keyed molarity mapping."""
    values = to_dict(species_molarities)
    values = _configure_component_values(
        values,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "species_molarities",
    )
    return to_annotated_value(
        _calc_osmolarity_from_mapping(values),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_osmolarity_from_mapping",
    )


def calc_osmolarity_from_props(
    species_molarities: Mapping[str, CustomProp],
    output_unit: str = "mol/L",
    *,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn=None,
    name: str = "osmolarity",
    description: str = "Calculate osmolarity from supplied particle molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate unit-aware osmolarity from supplied particle molarities."""
    _validate_custom_prop_mapping(species_molarities, "species_molarities")
    if unit is None:
        unit = output_unit
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})")
    values = _to_values(
        species_molarities,
        "species_molarities",
        output_unit,
        unit_conversion_fn,
    )
    values = _configure_component_values(
        values,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "species_molarities",
    )
    return to_annotated_value(
        _calc_osmolarity_from_mapping(values),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_osmolarity_from_mapping",
    )


# SECTION: Osmolality
@calculation_info(
    name="osmolality",
    description="Calculate osmolality by summing supplied dissolved particle molalities.",
    equation="b_osm = sum_i b_i",
    inputs={"species_molalities": "Molalities of explicitly supplied dissolved particle species."},
    outputs={"osmolality": "Total osmolality."},
    aliases=("osmotic molality",),
    notes=("No dissociation or speciation is inferred.",),
    tags=("osmolality", "composition", "low_level"),
)
def calc_osmolality(
    species_molalities: Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "osmolality",
    description: str = "Calculate osmolality from supplied particle molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float | NDArray[np.float64]]:
    """Calculate osmolality from explicitly supplied particle molalities."""
    return to_annotated_value(
        _calc_osmolality(species_molalities),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_osmolality",
    )


def calc_osmolality_from_sequence(
    species_molalities: Sequence[float | int],
    *,
    name: str = "osmolality",
    description: str = "Calculate osmolality from supplied particle molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate osmolality from a sequence of supplied particle molalities."""
    return to_annotated_value(
        float(_calc_osmolality(to_list(species_molalities))),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_osmolality",
    )


def calc_osmolality_from_mapping(
    species_molalities: Mapping[str, float | int | CustomProp],
    *,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    name: str = "osmolality",
    description: str = "Calculate osmolality from supplied particle molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate osmolality from a species-keyed molality mapping."""
    values = to_dict(species_molalities)
    values = _configure_component_values(
        values,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "species_molalities",
    )
    return to_annotated_value(
        _calc_osmolality_from_mapping(values),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_osmolality_from_mapping",
    )


def calc_osmolality_from_props(
    species_molalities: Mapping[str, CustomProp],
    output_unit: str = "mol/kg",
    *,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn=None,
    name: str = "osmolality",
    description: str = "Calculate osmolality from supplied particle molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate unit-aware osmolality from supplied particle molalities."""
    _validate_custom_prop_mapping(species_molalities, "species_molalities")
    if unit is None:
        unit = output_unit
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})")
    values = _to_values(
        species_molalities,
        "species_molalities",
        output_unit,
        unit_conversion_fn,
    )
    values = _configure_component_values(
        values,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "species_molalities",
    )
    return to_annotated_value(
        _calc_osmolality_from_mapping(values),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_osmolality_from_mapping",
    )


@calculation_info(
    name="ideal_osmotic_pressure",
    description="Calculate dilute ideal osmotic pressure from molar concentration.",
    equation="Pi = i*c*R*T",
    inputs={
        "molar_concentration": "Osmotically active molar concentration.",
        "temperature": "Absolute temperature.",
        "vant_hoff_factor": "Optional van 't Hoff factor.",
    },
    outputs={"osmotic_pressure": "Ideal osmotic pressure."},
    aliases=("van't Hoff osmotic pressure",),
    notes=("No osmotic coefficient or detailed electrolyte model is inferred.",),
    tags=("osmotic_pressure", "composition", "electrolyte", "low_level"),
)
def calc_ideal_osmotic_pressure(
    molar_concentration,
    temperature,
    vant_hoff_factor: float = 1.0,
    gas_constant: float = R_J_molK,
    output_concentration_unit: str | None = None,
    output_pressure_unit: str = "Pa",
    unit_conversion_fn=None,
    *,
    name: str = "ideal_osmotic_pressure",
    description: str = "Calculate dilute ideal osmotic pressure.",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate ideal osmotic pressure ``Pi = i*c*R*T``.

    Numeric concentrations are assumed to be in ``mol/m3`` unless the caller
    supplies a unit-aware ``CustomProp`` and ``output_concentration_unit``.
    """
    # SECTION: Normalize scalar inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    concentration_unit = output_concentration_unit or "mol/m3"
    c = _pos(
        molar_concentration,
        "molar_concentration",
        concentration_unit if isinstance(molar_concentration, CustomProp) else None,
        conversion_fn,
    )
    if concentration_unit != "mol/m3":
        c = float(conversion_fn(c, concentration_unit, "mol/m3"))
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    pressure = float(_calc_ideal_osmotic_pressure(c, t, vant_hoff_factor, gas_constant))
    if output_pressure_unit != "Pa":
        pressure = float(conversion_fn(pressure, "Pa", output_pressure_unit))
    return to_annotated_value(
        pressure,
        name=name,
        description=description,
        unit=output_pressure_unit,
        symbol=symbol,
        implementation="_calc_ideal_osmotic_pressure",
    )


__all__ = [
    "calc_osmolarity",
    "calc_osmolarity_from_sequence",
    "calc_osmolarity_from_mapping",
    "calc_osmolarity_from_props",
    "calc_osmolality",
    "calc_osmolality_from_sequence",
    "calc_osmolality_from_mapping",
    "calc_osmolality_from_props",
    "calc_ideal_osmotic_pressure",
]
