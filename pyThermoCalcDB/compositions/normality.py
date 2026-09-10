"""Normality and equivalent-concentration helpers."""

from collections.abc import Sequence
from typing import Literal, overload, cast

import numpy as np
from numpy.typing import NDArray
# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import to_annotated_value

# locals
from .core.normality import (
    _calc_normality,
    _calc_normality_from_props
)


# ======================================================================
# *** Public annotated API
# ======================================================================

@calculation_info(
    name="normality",
    description="Calculate normality from numeric molarity and an equivalence factor.",
    equation="normality = molarity * equivalence_factor",
    inputs={
        "molarity": "Molar concentration of the solute on the desired basis.",
        "equivalence_factor": "Reaction-context equivalence factor.",
    },
    outputs={
        "normality": "Normality on the same volume basis as molarity."
    },
    aliases=(
        "equivalent concentration",
        "normal concentration",
    ),
    notes=(
        "The equivalence factor must be supplied by the caller; it is not inferred from formula, valence, charge, or stoichiometry.",
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when molarity and equivalence_factor are already expressed on that normality basis.",
    ),
    tags=(
        "normality",
        "equivalent_concentration",
        "numeric",
        "scalar",
    )
)
def calc_normality(
    molarity: NDArray[np.number],
    equivalence_factor: float | int,
    *,
    name: str = "normality",
    description: str = "Calculate normality from molarity and an equivalence factor.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.float64]]:
    """Calculate annotated normality from numeric inputs."""
    value = _calc_normality(
        molarity=molarity,
        equivalence_factor=equivalence_factor,
        as_list=False,
    )

    return cast(
        AnnotatedValue[NDArray[np.float64]],
        to_annotated_value(
            value,
            name=name,
            description=description,
            unit=unit,
            symbol=symbol,
            implementation="_calc_normality",
        ),
    )


@calculation_info(
    name="normality",
    description="Calculate normality values from numeric molarities and an equivalence factor.",
    equation="normality = molarity * equivalence_factor",
    inputs={
        "molarity": "Molar concentrations of the solute on the desired basis.",
        "equivalence_factor": "Reaction-context equivalence factor.",
    },
    outputs={
        "normality": "Normalities on the same volume basis as molarity."
    },
    aliases=(
        "equivalent concentration",
        "normal concentration",
    ),
    notes=(
        "The equivalence factor must be supplied by the caller; it is not inferred from formula, valence, charge, or stoichiometry.",
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when molarity and equivalence_factor are already expressed on that normality basis.",
    ),
    tags=(
        "normality",
        "equivalent_concentration",
        "sequence",
        "numeric",
    )
)
def calc_normality_from_sequence(
    molarity: Sequence[float | int],
    equivalence_factor: float | int | Sequence[float | int],
    *,
    name: str = "normality",
    description: str = "Calculate normality from molarity and an equivalence factor.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Calculate annotated normality values from numeric sequences."""
    return to_annotated_value(
        _calc_normality(
            molarity=molarity,
            equivalence_factor=equivalence_factor,
            as_list=True,
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_normality",
    )


@calculation_info(
    name="normality",
    description="Calculate normality from unit-aware molarity and an equivalence factor.",
    equation="normality = molarity * equivalence_factor",
    inputs={
        "molarity": "Unit-aware molarity of the solute.",
        "equivalence_factor": "Reaction-context equivalence factor.",
        "output_unit": "Normality unit used to normalize molarity and annotate the result.",
    },
    outputs={
        "normality": "Normality in output_unit."
    },
    aliases=(
        "equivalent concentration",
        "normal concentration",
    ),
    notes=(
        "The output_unit must be a ratio such as eq/L.",
        "The denominator of output_unit defines the molarity conversion basis, for example eq/L uses mol/L.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "normality",
        "equivalent_concentration",
        "unit_aware",
        "unit_conversion",
        "scalar",
    )
)
def calc_normality_from_props(
    molarity: CustomProp,
    equivalence_factor: float | int,
    output_unit: str = "eq/L",
    unit_conversion_fn: UnitConversionFn | None = None,
    *,
    name: str = "normality",
    description: str = "Calculate normality from molarity and an equivalence factor.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate annotated normality from unit-aware molarity."""
    # SECTION: set default unit for output if not provided
    if unit is None:
        unit = output_unit

    # check unit & output unit consistency
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    return to_annotated_value(
        _calc_normality_from_props(
            molarity=molarity,
            equivalence_factor=equivalence_factor,
            output_unit=output_unit,
            unit_conversion_fn=unit_conversion_fn,
        ),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_normality_from_props",
    )


# SECTION: Public exports
__all__ = [
    "calc_normality",
    "calc_normality_from_sequence",
    "calc_normality_from_props",
]
