"""Molality-based ionic-strength calculations."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, Literal, cast

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
# ! core
from .core.ionic_strength import (
    _calc_ionic_strength,
    _calc_ionic_strength_from_mapping,
    _calc_ionic_strength_from_props,
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
    value = _calc_ionic_strength(
        values=molalities,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength",
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
    value = float(_calc_ionic_strength(
        values=molalities,
        charges=charges,
    ))

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength",
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
    value = _calc_ionic_strength_from_mapping(
        values=molalities,
        charges=charges,
        mode='molality'
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_ionic_strength_from_mapping",
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
    charges: Mapping[str, CustomProp],
    output_molality_unit: str = "mol/kg",
    unit_conversion_fn: UnitConversionFn | None = None,
    *,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength.",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality-based ionic strength from unit-aware mappings."""
    # SECTION: Calculate value
    value = _calc_ionic_strength_from_props(
        values=molalities,
        charges=charges,
        mode=cast(Literal['molality'], 'molality'),
        output_unit=output_molality_unit,
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
        implementation="_calc_ionic_strength_from_props",
    )


# SECTION: Public exports
__all__ = [
    # public
    "calc_ionic_strength_molality",
    "calc_ionic_strength_molality_from_sequence",
    "calc_ionic_strength_molality_from_mapping",
    "calc_ionic_strength_molality_from_props",
]
