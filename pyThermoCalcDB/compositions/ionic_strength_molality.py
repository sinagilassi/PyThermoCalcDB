"""Molality-based ionic-strength calculations."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, cast

import numpy as np
from numpy.typing import NDArray

# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.components import extract_components_values
from pythermodb_settings.utils.quantity import to_dict
from pythermodb_settings.utils.validators import non_negative, same_shape

# locals
from ..utils.conversions import (
    _resolve_result_unit,
    _resolve_unit_conversion_fn,
    _validate_same_keys,
)
from ..utils.tools import to_annotated_value
from .ionic_strength import _calc_ionic_strength


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================
# ! ::: Numeric array-like core

def _calc_ionic_strength_molality(
    molalities: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate molality-based ionic strength from numeric array-like inputs."""
    # SECTION: Calculate ionic strength
    return _calc_ionic_strength(
        concentrations=molalities,
        charges=charges,
    )


# ! ::: Sequence adapter

def _calc_ionic_strength_molality_from_sequence(
    molalities: Sequence[float | int],
    charges: Sequence[float | int],
) -> float:
    """Calculate molality-based ionic strength from numeric sequence inputs."""
    # SECTION: Validate inputs
    non_negative(molalities, "molalities")
    same_shape(molalities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            molalities,
            charges,
        ),
    )


# ! ::: Mapping adapter

def _calc_ionic_strength_molality_from_mapping(
    molalities: Mapping[str, float | int],
    charges: Mapping[str, float | int],
) -> float:
    """Calculate molality-based ionic strength from numeric mapping inputs."""
    # SECTION: Validate inputs
    non_negative(molalities, "molalities")
    same_shape(molalities, charges)
    _validate_same_keys(molalities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            list(molalities.values()),
            [charges[key] for key in molalities],
        ),
    )


# ! ::: Unit-aware mapping adapter

def _calc_ionic_strength_molality_from_props(
    molalities: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_molality_unit: str = "mol/kg",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate molality-based ionic strength from unit-aware molality inputs."""
    # SECTION: Normalize unit-aware molalities
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    normalized_molalities = to_dict(
        molalities,
        output_molality_unit,
        unit_conversion_fn=conversion_fn,
    )

    # ! Mapping inputs must pair the same species.
    _validate_same_keys(normalized_molalities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            list(normalized_molalities.values()),
            [charges[key] for key in normalized_molalities],
        ),
    )


# ! ::: Mapping input with component metadata

def _calc_ionic_strength_molality_with_components_from_mapping(
    molalities: Mapping[str, float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
) -> float:
    """Calculate molality ionic strength from mapping input and components."""
    # SECTION: Normalize molality values by component metadata
    configured = config_components_values(
        values=dict(molalities),
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=True,
    )

    # ! Component configuration must succeed before charges are extracted.
    if configured is None:
        raise ValueError("Failed to normalize molalities component values.")

    _, molality_values = configured

    # SECTION: Extract component charges
    extracted = extract_components_values(
        attribute_name="net_charge",
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    if extracted is None:
        raise ValueError("Failed to extract component charges.")

    _, charge_values = extracted

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            molality_values,
            charge_values,
        ),
    )


# ! ::: Sequence input with component metadata

def _calc_ionic_strength_molality_with_components_from_sequence(
    molalities: Sequence[float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
) -> float:
    """Calculate molality ionic strength from sequence input and components."""
    # SECTION: Extract component charges
    extracted = extract_components_values(
        attribute_name="net_charge",
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    if extracted is None:
        raise ValueError("Failed to extract component charges.")

    _, charge_values = extracted

    # SECTION: Calculate ionic strength
    return _calc_ionic_strength_molality_from_sequence(
        molalities=molalities,
        charges=charge_values,
    )


# =====================================================================
# *** Public annotated API
# =====================================================================
# ::: annotated for numpy array

@calculation_info(
    name="ionic_strength_molality",
    description="Calculate molality-based ionic strength.",
    equation="I_m = 0.5 * sum_i(m_i * z_i**2)",
    inputs={
        "molalities": "Species molalities.",
        "charges": "Species charges.",
    },
    outputs={
        "ionic_strength": "Molality-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molality",
        "array_like",
        "numpy",
        "numeric",
    ),
)
def calc_ionic_strength_molality(
    molalities: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float | NDArray[np.float64]]:
    """Return annotated molality-based ionic strength from numeric inputs."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molality(
        molalities=molalities,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molality",
    )


# ::: annotated for sequence

@calculation_info(
    name="ionic_strength_molality",
    description="Calculate molality-based ionic strength from sequence inputs.",
    equation="I_m = 0.5 * sum_i(m_i * z_i**2)",
    inputs={
        "molalities": "Species molalities.",
        "charges": "Species charges.",
    },
    outputs={
        "ionic_strength": "Molality-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molality",
        "sequence",
        "numeric",
    ),
)
def calc_ionic_strength_molality_from_sequence(
    molalities: Sequence[float | int],
    charges: Sequence[float | int],
    *,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality-based ionic strength from numeric sequences."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molality_from_sequence(
        molalities=molalities,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molality_from_sequence",
    )


# ::: annotated for mapping

@calculation_info(
    name="ionic_strength_molality",
    description="Calculate molality-based ionic strength from mapping inputs.",
    equation="I_m = 0.5 * sum_i(m_i * z_i**2)",
    inputs={
        "molalities": "Keyed species molalities.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "ionic_strength": "Molality-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molality",
        "mapping",
        "numeric",
    ),
)
def calc_ionic_strength_molality_from_mapping(
    molalities: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    *,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality-based ionic strength from numeric mappings."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molality_from_mapping(
        molalities=molalities,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molality_from_mapping",
    )


# ::: annotated for CustomProp mapping

@calculation_info(
    name="ionic_strength_molality",
    description="Calculate molality-based ionic strength from unit-aware mapping inputs.",
    equation="I_m = 0.5 * sum_i(m_i * z_i**2)",
    inputs={
        "molalities": "Keyed unit-aware species molalities.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "ionic_strength": "Molality-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molality",
        "mapping",
        "unit_aware",
    ),
)
def calc_ionic_strength_molality_from_props(
    molalities: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_molality_unit: str = "mol/kg",
    unit_conversion_fn: UnitConversionFn | None = None,
    *,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength.",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality-based ionic strength from unit-aware mappings."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molality_from_props(
        molalities=molalities,
        charges=charges,
        output_molality_unit=output_molality_unit,
        unit_conversion_fn=unit_conversion_fn,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=_resolve_result_unit(
            "molalities",
            molalities,
            output_molality_unit,
        ),
        symbol=symbol,
        implementation="_calc_ionic_strength_molality_from_props",
    )


# ::: annotated for mapping and component metadata

def calc_ionic_strength_molality_with_components_from_mapping(
    molalities: Mapping[str, float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    *,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength from component charges.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality ionic strength from mapping component metadata."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molality_with_components_from_mapping(
        molalities=molalities,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molality_with_components_from_mapping",
    )


# ::: annotated for sequence and component metadata

def calc_ionic_strength_molality_with_components_from_sequence(
    molalities: Sequence[float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    *,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength from component charges.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality ionic strength from sequence component metadata."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molality_with_components_from_sequence(
        molalities=molalities,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molality_with_components_from_sequence",
    )


# =====================================================================
# *** Aliases
# =====================================================================

calc_mapping_ionic_strength_molality = calc_ionic_strength_molality_from_mapping
calc_sequence_ionic_strength_molality = calc_ionic_strength_molality_from_sequence
calc_mapping_ionic_strength_molality_from_components = (
    calc_ionic_strength_molality_with_components_from_mapping
)
calc_sequence_ionic_strength_molality_from_components = (
    calc_ionic_strength_molality_with_components_from_sequence
)


# SECTION: Public exports
__all__ = [
    # internal
    "_calc_ionic_strength_molality",
    "_calc_ionic_strength_molality_from_sequence",
    "_calc_ionic_strength_molality_from_mapping",
    "_calc_ionic_strength_molality_from_props",
    "_calc_ionic_strength_molality_with_components_from_sequence",
    "_calc_ionic_strength_molality_with_components_from_mapping",

    # public
    "calc_ionic_strength_molality",
    "calc_ionic_strength_molality_from_sequence",
    "calc_ionic_strength_molality_from_mapping",
    "calc_ionic_strength_molality_from_props",
    "calc_ionic_strength_molality_with_components_from_sequence",
    "calc_ionic_strength_molality_with_components_from_mapping",

    # compatibility
    "calc_mapping_ionic_strength_molality",
    "calc_sequence_ionic_strength_molality",
    "calc_mapping_ionic_strength_molality_from_components",
    "calc_sequence_ionic_strength_molality_from_components",
]
