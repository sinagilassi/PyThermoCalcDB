# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional, cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import (
    AnnotatedValue,
    Component,
    ComponentKey,
    CustomProp,
)
from pythermodb_settings.utils import to_annotated_value

# locals
from .core.molality import (
    _calc_molalities,
    _calc_molalities_from_sequence,
    _calc_molalities_from_mapping,
    _calc_molalities_from_props,
    _calc_component_molalities_from_props
)

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Public annotated API
# ======================================================================

# ::: annotated for numpy array
@calculation_info(
    name="molality",
    description="Calculate the molality of each component in a solution.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solvent_mass": "Mass of the solvent."
    },
    outputs={
        "molality": "Molality of each component in the solution."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solvent_mass are already expressed on that molality basis.",
    ),
    tags=(
        "molality",
        "array_like",
        "numpy",
        "component_wise",
        "numeric",
    )
)
def _molality_annotated(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solvent_mass: float | int | NDArray[np.number],
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.floating]]:
    """Calculate annotated molality values from array-like inputs."""
    return to_annotated_value(
        _calc_molalities(
            component_moles=component_moles,
            solvent_mass=solvent_mass
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molalities",
    )


# ::: annotated for sequence
@calculation_info(
    name="molality",
    description="Calculate the molality of each component in a solution.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solvent_mass": "Mass of the solvent."
    },
    outputs={
        "molality": "Molality of each component in the solution."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solvent_mass are already expressed on that molality basis.",
    ),
    tags=(
        "molality",
        "sequence",
        "component_wise",
        "numeric",
    )
)
def _molality_1_annotated(
    component_moles: Sequence[float],
    solvent_mass: float,
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None
) -> AnnotatedValue[list[float]]:
    """Calculate annotated molality values from a sequence of component moles."""
    return to_annotated_value(
        _calc_molalities_from_sequence(
            component_moles=component_moles,
            solvent_mass=solvent_mass
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molalities_from_sequence"
    )


# ::: annotated for mapping
@calculation_info(
    name="molality",
    description="Calculate the molality of each keyed component in a solution.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Mapping of component identifiers to component moles in the solution.",
        "solvent_mass": "Mass of the solvent."
    },
    outputs={
        "molality": "Mapping of component identifiers to molality values."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solvent_mass are already expressed on that molality basis.",
    ),
    tags=(
        "numeric",
        "molality",
        "mapping",
        "component_aware",
        "keyed",
    )
)
def _molality_2_annotated(
    component_moles: Mapping[str, float | int],
    solvent_mass: float,
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molality values from a component-mole mapping."""
    return to_annotated_value(
        _calc_molalities_from_mapping(
            component_moles=component_moles,
            solvent_mass=solvent_mass
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molalities_from_mapping"
    )


# ::: annotated for mapping with custom properties
@calculation_info(
    name="molality",
    description="Calculate keyed molality values from unit-aware component moles and solvent mass.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Mapping of component identifiers to unit-aware component mole amounts.",
        "solvent_mass": "Unit-aware mass of the solvent.",
        "output_unit": "Molality unit used to normalize component moles and solvent mass."
    },
    outputs={
        "molality": "Mapping of component identifiers to molality values in output_unit."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "The output_unit must be a ratio such as mol/kg with amount in the numerator and mass in the denominator.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "molality",
        "mapping",
        "unit_aware",
        "component_aware",
        "unit_conversion",
    )
)
def _molality_3_annotated(
    component_moles: Mapping[str, CustomProp],
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molality values from unit-aware component moles."""
    # SECTION: set default unit for output if not provided
    if unit is None:
        unit = output_unit

    # check unit & output unit consistency
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    return to_annotated_value(
        _calc_molalities_from_props(
            component_moles=component_moles,
            solvent_mass=solvent_mass,
            output_unit=output_unit
        ),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_molalities_from_props"
    )


# ::: annotated for component mapping with custom properties
@calculation_info(
    name="component_molality",
    description="Calculate unit-aware component molalities with optional component-key remapping and ordering.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Mapping of component identifiers to unit-aware component mole amounts.",
        "solvent_mass": "Unit-aware mass of the solvent.",
        "output_unit": "Molality unit used to normalize component moles and solvent mass.",
        "components": "Optional component definitions used to resolve and order component identifiers.",
        "component_key": "Optional component key used for identifier matching."
    },
    outputs={
        "component_molality": "Mapping of resolved component identifiers to molality values in output_unit."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "The output_unit must be a ratio such as mol/kg with amount in the numerator and mass in the denominator.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "component_molality",
        "molality",
        "mapping",
        "component_wise",
        "unit_aware",
        "unit_conversion",
        "component_key",
        "component_ordering",
    )
)
def _molality_4_annotated(
    component_moles: Mapping[str, CustomProp],
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molality values for mapped components with units."""
    # SECTION: set default unit for output if not provided
    if unit is None:
        unit = output_unit

    # check unit & output unit consistency
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    res = _calc_component_molalities_from_props(
        component_moles=component_moles,
        solvent_mass=solvent_mass,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order
    )

    # convert the result to an annotated value
    return to_annotated_value(
        value=res,
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_component_molalities_from_props"
    )


# ======================================================================
# *** Aliases
# ======================================================================
# >> molalities
calc_molalities = _molality_annotated

# >> molalities from sequence
calc_molalities_from_sequence = _molality_1_annotated

# >> molalities from mapping
calc_molalities_from_mapping = _molality_2_annotated

# >> molalities with units
calc_molalities_from_props = _molality_3_annotated

# >> component molalities
calc_component_molalities_from_props = _molality_4_annotated

# all
__all__ = [
    "calc_molalities",
    "calc_molalities_from_sequence",
    "calc_molalities_from_mapping",
    "calc_molalities_from_props",
    "calc_component_molalities_from_props",
]
