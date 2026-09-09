"""Molarity-based ionic-strength calculations."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional

import numpy as np
from numpy.typing import NDArray

# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp
from pythermodb_settings.models.units import UnitConversionFn

# locals
from ..utils.conversions import (
    _resolve_result_unit,
)
from ..utils.tools import to_annotated_value
from .core.ionic_strength import (
    _calc_ionic_strength_molarity,
    _calc_ionic_strength_molarity_from_sequence,
    _calc_ionic_strength_molarity_from_mapping,
    _calc_ionic_strength_molarity_from_props,
    _calc_ionic_strength_molarity_with_components_from_sequence,
    _calc_ionic_strength_molarity_with_components_from_mapping,
)


# =====================================================================
# *** Public annotated API
# =====================================================================
# ::: annotated for numpy array

@calculation_info(
    name="ionic_strength_molarity",
    description="Calculate molarity-based ionic strength.",
    equation="I_c = 0.5 * sum_i(c_i * z_i**2)",
    inputs={
        "molarities": "Species molarities.",
        "charges": "Species charges.",
    },
    outputs={
        "ionic_strength": "Molarity-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molarity",
        "array_like",
        "numpy",
        "numeric",
    ),
)
def calc_ionic_strength_molarity(
    molarities: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "ionic_strength",
    description: str = "Molarity-based ionic strength.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float | NDArray[np.float64]]:
    """Return annotated molarity-based ionic strength from numeric inputs."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molarity(
        molarities=molarities,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molarity",
    )


# ::: annotated for sequence

@calculation_info(
    name="ionic_strength_molarity",
    description="Calculate molarity-based ionic strength from sequence inputs.",
    equation="I_c = 0.5 * sum_i(c_i * z_i**2)",
    inputs={
        "molarities": "Species molarities.",
        "charges": "Species charges.",
    },
    outputs={
        "ionic_strength": "Molarity-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molarity",
        "sequence",
        "numeric",
    ),
)
def calc_ionic_strength_molarity_from_sequence(
    molarities: Sequence[float | int],
    charges: Sequence[float | int],
    *,
    name: str = "ionic_strength",
    description: str = "Molarity-based ionic strength.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity-based ionic strength from numeric sequences."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molarity_from_sequence(
        molarities=molarities,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molarity_from_sequence",
    )


# ::: annotated for mapping

@calculation_info(
    name="ionic_strength_molarity",
    description="Calculate molarity-based ionic strength from mapping inputs.",
    equation="I_c = 0.5 * sum_i(c_i * z_i**2)",
    inputs={
        "molarities": "Keyed species molarities.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "ionic_strength": "Molarity-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molarity",
        "mapping",
        "numeric",
    ),
)
def calc_ionic_strength_molarity_from_mapping(
    molarities: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    *,
    name: str = "ionic_strength",
    description: str = "Molarity-based ionic strength.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity-based ionic strength from numeric mappings."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molarity_from_mapping(
        molarities=molarities,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_molarity_from_mapping",
    )


# ::: annotated for CustomProp mapping

@calculation_info(
    name="ionic_strength_molarity",
    description="Calculate molarity-based ionic strength from unit-aware mapping inputs.",
    equation="I_c = 0.5 * sum_i(c_i * z_i**2)",
    inputs={
        "molarities": "Keyed unit-aware species molarities.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "ionic_strength": "Molarity-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molarity",
        "mapping",
        "unit_aware",
    ),
)
def calc_ionic_strength_molarity_from_props(
    molarities: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_molarity_unit: str = "mol/L",
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "ionic_strength",
    description: str = "Molarity-based ionic strength.",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity-based ionic strength from unit-aware mappings."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molarity_from_props(
        molarities=molarities,
        charges=charges,
        output_molarity_unit=output_molarity_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=_resolve_result_unit(
            "molarities",
            molarities,
            output_molarity_unit,
        ),
        symbol=symbol,
        implementation="_calc_ionic_strength_molarity_from_props",
    )


# ::: annotated for mapping and component metadata

@calculation_info(
    name="ionic_strength_molarity",
    description="Calculate molarity-based ionic strength from mapping and component metadata.",
    equation="I_c = 0.5 * sum_i(c_i * z_i**2)",
    inputs={
        "molarities": "Keyed species molarities.",
        "components": "Species component metadata used to obtain charges.",
    },
    outputs={
        "ionic_strength": "Molarity-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molarity",
        "mapping",
        "components",
    ),
)
def calc_ionic_strength_molarity_with_components_from_mapping(
    molarities: Mapping[str, float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    *,
    name: str = "ionic_strength",
    description: str = "Molarity-based ionic strength from component charges.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity ionic strength from mapping component metadata."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molarity_with_components_from_mapping(
        molarities=molarities,
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
        implementation="_calc_ionic_strength_molarity_with_components_from_mapping",
    )


# ::: annotated for sequence and component metadata

@calculation_info(
    name="ionic_strength_molarity",
    description="Calculate molarity-based ionic strength from sequence and component metadata.",
    equation="I_c = 0.5 * sum_i(c_i * z_i**2)",
    inputs={
        "molarities": "Species molarities.",
        "components": "Species component metadata used to obtain charges.",
    },
    outputs={
        "ionic_strength": "Molarity-based ionic strength."
    },
    tags=(
        "ionic_strength",
        "molarity",
        "sequence",
        "components",
    ),
)
def calc_ionic_strength_molarity_with_components_from_sequence(
    molarities: Sequence[float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    *,
    name: str = "ionic_strength",
    description: str = "Molarity-based ionic strength from component charges.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity ionic strength from sequence component metadata."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_molarity_with_components_from_sequence(
        molarities=molarities,
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
        implementation="_calc_ionic_strength_molarity_with_components_from_sequence",
    )


# SECTION: Public exports
__all__ = [
    # public
    "calc_ionic_strength_molarity",
    "calc_ionic_strength_molarity_from_sequence",
    "calc_ionic_strength_molarity_from_mapping",
    "calc_ionic_strength_molarity_from_props",
    "calc_ionic_strength_molarity_with_components_from_sequence",
    "calc_ionic_strength_molarity_with_components_from_mapping",
]
