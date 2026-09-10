# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Optional
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import (
    AnnotatedValue,
    Component,
    ComponentKey,
    CustomProp,
    UnitConversionFn,
)
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.utils import (
    to_annotated_value,
)
# locals
from .core.concentration import (
    _calc_concentrations,
    _calc_concentrations_from_mapping,
    _calc_concentrations_from_props,
)

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

# ! ::: Mass concentration
@calculation_info(
    name="mass_concentration",
    description="Calculate the mass concentration of each component in a solution.",
    equation="mass_concentration = component_mass / solution_volume",
    inputs={
        "component_amounts": "Component masses in the solution.",
        "solution_volume": "Volume of the solution.",
    },
    outputs={
        "mass_concentration": "Mass concentration of each component in the solution.",
    },
    aliases=(
        "mass concentration",
        "density concentration",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_amounts and solution_volume are already expressed on that concentration basis.",
    ),
    tags=(
        "mass_concentration",
        "array_like",
        "numpy",
        "component_wise",
        "numeric",
    ),
)
def calc_mass_concentrations(
        component_amounts: Sequence[float | int] | NDArray[np.number],
        solution_volume: float | int | NDArray[np.number],
        *,
        name: str = "mass_concentration",
        description: str = "Calculate the mass concentration of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.float64]]:
    """Calculate annotated mass concentrations from array-like inputs."""
    return to_annotated_value(
        _calc_concentrations(
            component_amounts=component_amounts,
            solution_volume=solution_volume,
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_concentrations",
    )


@calculation_info(
    name="mass_concentration",
    description="Calculate the mass concentration of each component in a solution.",
    equation="mass_concentration = component_mass / solution_volume",
    inputs={
        "component_mass": "Component masses in the solution.",
        "solution_volume": "Volume of the solution.",
    },
    outputs={
        "mass_concentration": "Mass concentration of each component in the solution.",
    },
    aliases=(
        "mass concentration",
        "density concentration",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_mass and solution_volume are already expressed on that concentration basis.",
    ),
    tags=(
        "mass_concentration",
        "sequence",
        "component_wise",
        "numeric",
    ),
)
def calc_mass_concentrations_from_sequence(
        component_mass: Sequence[float | int],
        solution_volume: float | int,
        *,
        name: str = "mass_concentration",
        description: str = "Calculate the mass concentration of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Calculate annotated mass concentrations from a sequence of masses."""
    return to_annotated_value(
        _calc_concentrations(
            component_amounts=component_mass,
            solution_volume=solution_volume,
        ).tolist(),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_concentrations",
    )


@calculation_info(
    name="mass_concentration",
    description="Calculate the mass concentration of each keyed component in a solution.",
    equation="mass_concentration = component_mass / solution_volume",
    inputs={
        "component_mass": "Mapping of component identifiers to component masses in the solution.",
        "solution_volume": "Volume of the solution.",
    },
    outputs={
        "mass_concentration": "Mapping of component identifiers to mass concentration values.",
    },
    aliases=(
        "mass concentration",
        "density concentration",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_mass and solution_volume are already expressed on that concentration basis.",
    ),
    tags=(
        "mass_concentration",
        "mapping",
        "component_aware",
        "keyed",
        "numeric",
    ),
)
def calc_mass_concentrations_from_mapping(
    component_mass: Mapping[str, float | int],
    solution_volume: float,
    *,
    name: str = "mass_concentration",
    description: str = "Calculate the mass concentration of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated mass concentrations from a keyed mass mapping."""
    return to_annotated_value(
        _calc_concentrations_from_mapping(
            component_amounts=component_mass,
            solution_volume=solution_volume,
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_concentrations_from_mapping",
    )


@calculation_info(
    name="component_mass_concentration",
    description="Calculate unit-aware component mass concentrations with optional component-key remapping and ordering.",
    equation="mass_concentration = component_mass / solution_volume",
    inputs={
        "component_mass": "Mapping of component identifiers to unit-aware component masses.",
        "solution_volume": "Unit-aware volume of the solution.",
        "output_unit": "Mass concentration unit used to normalize component masses and solution volume.",
        "components": "Optional component definitions used to resolve and order component identifiers.",
        "component_key": "Optional component key used for identifier matching.",
    },
    outputs={
        "component_mass_concentration": "Mapping of resolved component identifiers to mass concentration values in output_unit.",
    },
    aliases=(
        "mass concentration",
        "density concentration",
    ),
    notes=(
        "The output_unit must be a ratio such as kg/m^3 with mass in the numerator and volume in the denominator.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "component_mass_concentration",
        "mass_concentration",
        "mapping",
        "component_wise",
        "unit_aware",
        "unit_conversion",
        "component_key",
        "component_ordering",
    ),
)
def calc_mass_concentrations_from_props(
    component_mass: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
    *,
    name: str = "component_mass_concentration",
    description: str = "Calculate the mass concentration of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated component-aware unit-normalized mass concentrations."""
    if unit is None:
        unit = output_unit

    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    return to_annotated_value(
        _calc_concentrations_from_props(
            component_amounts=component_mass,
            solution_volume=solution_volume,
            output_unit=output_unit,
            components=components,
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order,
            unit_conversion_fn=unit_conversion_fn,
        ),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_concentrations_from_props",
    )


# ! ::: Molar concentration
@calculation_info(
    name="molar_concentration",
    description="Calculate the molar concentration of each component in a solution.",
    equation="molar_concentration = component_moles / solution_volume",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solution_volume": "Volume of the solution.",
    },
    outputs={
        "molar_concentration": "Molar concentration of each component in the solution.",
    },
    aliases=(
        "molar concentration",
        "amount concentration",
        "molarity",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solution_volume are already expressed on that concentration basis.",
    ),
    tags=(
        "molar_concentration",
        "array_like",
        "numpy",
        "component_wise",
        "numeric",
    ),
)
def calc_molar_concentrations(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solution_volume: float | int | NDArray[np.number],
    *,
    name: str = "molar_concentration",
    description: str = "Calculate the molar concentration of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.float64]]:
    """Calculate annotated molar concentrations from array-like inputs."""
    return to_annotated_value(
        _calc_concentrations(
            component_amounts=component_moles,
            solution_volume=solution_volume,
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_concentrations",
    )


@calculation_info(
    name="molar_concentration",
    description="Calculate the molar concentration of each component in a solution.",
    equation="molar_concentration = component_moles / solution_volume",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solution_volume": "Volume of the solution.",
    },
    outputs={
        "molar_concentration": "Molar concentration of each component in the solution.",
    },
    aliases=(
        "molar concentration",
        "amount concentration",
        "molarity",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solution_volume are already expressed on that concentration basis.",
    ),
    tags=(
        "molar_concentration",
        "sequence",
        "component_wise",
        "numeric",
    ),
)
def calc_molar_concentrations_from_sequence(
        component_moles: Sequence[float],
        solution_volume: float,
        *,
        name: str = "molar_concentration",
        description: str = "Calculate the molar concentration of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Calculate annotated molar concentrations from a sequence of moles."""
    return to_annotated_value(
        _calc_concentrations(
            component_amounts=component_moles,
            solution_volume=solution_volume,
        ).tolist(),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_concentrations_from_sequence",
    )


@calculation_info(
    name="molar_concentration",
    description="Calculate the molar concentration of each keyed component in a solution.",
    equation="molar_concentration = component_moles / solution_volume",
    inputs={
        "component_moles": "Mapping of component identifiers to component moles in the solution.",
        "solution_volume": "Volume of the solution.",
    },
    outputs={
        "molar_concentration": "Mapping of component identifiers to molar concentration values.",
    },
    aliases=(
        "molar concentration",
        "amount concentration",
        "molarity",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solution_volume are already expressed on that concentration basis.",
    ),
    tags=(
        "molar_concentration",
        "mapping",
        "component_aware",
        "keyed",
        "numeric",
    ),
)
def calc_molar_concentrations_from_mapping(
    component_moles: Mapping[str, float | int],
    solution_volume: float,
    *,
    name: str = "molar_concentration",
    description: str = "Calculate the molar concentration of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molar concentrations from a keyed mole mapping."""
    return to_annotated_value(
        _calc_concentrations_from_mapping(
            component_amounts=component_moles,
            solution_volume=solution_volume,
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_concentrations_from_mapping",
    )


@calculation_info(
    name="component_molar_concentration",
    description="Calculate unit-aware component molar concentrations with optional component-key remapping and ordering.",
    equation="molar_concentration = component_moles / solution_volume",
    inputs={
        "component_moles": "Mapping of component identifiers to unit-aware component mole amounts.",
        "solution_volume": "Unit-aware volume of the solution.",
        "output_unit": "Molar concentration unit used to normalize component moles and solution volume.",
        "components": "Optional component definitions used to resolve and order component identifiers.",
        "component_key": "Optional component key used for identifier matching.",
    },
    outputs={
        "component_molar_concentration": "Mapping of resolved component identifiers to molar concentration values in output_unit.",
    },
    aliases=(
        "molar concentration",
        "amount concentration",
        "molarity",
    ),
    notes=(
        "The output_unit must be a ratio such as mol/L with amount in the numerator and volume in the denominator.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "component_molar_concentration",
        "molar_concentration",
        "mapping",
        "component_wise",
        "unit_aware",
        "unit_conversion",
        "component_key",
        "component_ordering",
    ),
)
def calc_molar_concentrations_from_props(
    component_moles: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
    *,
    name: str = "component_molar_concentration",
    description: str = "Calculate the molar concentration of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated component-aware unit-normalized molar concentrations."""
    if unit is None:
        unit = output_unit

    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    return to_annotated_value(
        _calc_concentrations_from_props(
            component_amounts=component_moles,
            solution_volume=solution_volume,
            output_unit=output_unit,
            components=components,
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order,
            unit_conversion_fn=unit_conversion_fn,
        ),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_component_concentrations_from_props",
    )


# export functions
__all__ = [
    # mass
    "calc_mass_concentrations",
    "calc_mass_concentrations_from_sequence",
    "calc_mass_concentrations_from_mapping",
    "calc_mass_concentrations_from_props",
    # molar
    "calc_molar_concentrations",
    "calc_molar_concentrations_from_sequence",
    "calc_molar_concentrations_from_mapping",
    "calc_molar_concentrations_from_props",
]
