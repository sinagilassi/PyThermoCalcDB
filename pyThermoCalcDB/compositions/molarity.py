# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Optional
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp, AnnotatedValue
from pythermodb_settings.utils import (
    to_annotated_value,
)
from pythermodb_settings.decorators import calculation_info
# locals
from .core.molarity import (
    _calc_molarities,
    _calc_molarities_from_mapping,
    _calc_molarities_from_props,
)

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Public annotated API
# ======================================================================

# ::: annotated for numpy array


@calculation_info(
    name="molarity",
    description="Calculate the molarity of each component in a solution.",
    equation="molarity = component_moles / solution_volume",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solution_volume": "Volume of the solution."
    },
    outputs={
        "molarity": "Molarity of each component in the solution."
    },
    aliases=(
        "molar concentration",
        "amount concentration",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solution_volume are already expressed on that molarity basis.",
    ),
    tags=(
        "molarity",
        "array_like",
        "numpy",
        "component_wise",
        "numeric",
    )
)
def calc_molarities(
        component_moles: Sequence[float | int] | NDArray[np.number],
        solution_volume: float | int | NDArray[np.number],
        *,
        name: str = "molarity",
        description: str = "Calculate the molarity of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.floating]]:
    """Calculate annotated molarity values from array-like inputs.

    Parameters
    ----------
    component_moles : Sequence[float | int] | NDArray[np.number]
        Component mole amounts. May be a Python sequence or a
        1-D/2-D NumPy array.
    solution_volume : float | int | NDArray[np.number]
        Solution volume.

    Returns
    -------
    AnnotatedValue[NDArray[np.floating]]
        Molarities with the same shape as ``component_moles``.
    """
    return to_annotated_value(
        _calc_molarities(
            component_moles=component_moles,
            solution_volume=solution_volume
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molarities",
    )


# ::: annotated for sequence


@calculation_info(
    name="molarity",
    description="Calculate the molarity of each component in a solution.",
    equation="molarity = component_moles / solution_volume",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solution_volume": "Volume of the solution."
    },
    outputs={
        "molarity": "Molarity of each component in the solution."
    },
    aliases=(
        "molar concentration",
        "amount concentration",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solution_volume are already expressed on that molarity basis.",
    ),
    tags=(
        "molarity",
        "sequence",
        "component_wise",
        "numeric",
    )
)
def calc_molarities_from_sequence(
        component_moles: Sequence[float],
        solution_volume: float,
        *,
        name: str = "molarity",
        description: str = "Calculate the molarity of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None
) -> AnnotatedValue[list[float]]:
    """Calculate annotated molarity values from a sequence of component moles.

    Parameters
    ----------
    component_moles : Sequence[float]
        A sequence of moles for each component.
    solution_volume : float
        The volume of the solution.
    name : str, optional
        The name stored in the annotated result. Defaults to ``"molarity"``.
    description : str, optional
        The description stored in the annotated result.
    unit : str, optional
        The unit stored in the annotated result.
    symbol : str, optional
        The symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[list[float]]
        The calculated molarity values with metadata.
    """
    return to_annotated_value(
        _calc_molarities(
            component_moles=component_moles,
            solution_volume=solution_volume,
            as_list=True
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molarities_from_sequence"
    )

# ::: annotated for mapping


@calculation_info(
    name="molarity",
    description="Calculate the molarity of each keyed component in a solution.",
    equation="molarity = component_moles / solution_volume",
    inputs={
        "component_moles": "Mapping of component identifiers to component moles in the solution.",
        "solution_volume": "Volume of the solution."
    },
    outputs={
        "molarity": "Mapping of component identifiers to molarity values."
    },
    aliases=(
        "molar concentration",
        "amount concentration",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solution_volume are already expressed on that molarity basis.",
    ),
    tags=(
        "numeric",
        "molarity",
        "mapping",
        "component_aware",
        "keyed",
    )
)
def calc_molarities_from_mapping(
        component_moles: Mapping[str, float | int],
        solution_volume: float,
        *,
        name: str = "molarity",
        description: str = "Calculate the molarity of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molarity values from a component-mole mapping.

    Parameters
    ----------
    component_moles : Mapping[str, float | int]
        A mapping of component names to their respective moles.
    solution_volume : float
        The volume of the solution.
    name : str, optional
        The name stored in the annotated result. Defaults to ``"molarity"``.
    description : str, optional
        The description stored in the annotated result.
    unit : str, optional
        The unit stored in the annotated result.
    symbol : str, optional
        The symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        The calculated molarity values with metadata.
    """
    return to_annotated_value(
        _calc_molarities_from_mapping(
            component_moles=component_moles,
            solution_volume=solution_volume
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molarities_from_mapping"
    )


# ::: annotated for mapping with custom properties

@calculation_info(
    name="component_molarity",
    description="Calculate unit-aware component molarities with optional component-key remapping and ordering.",
    equation="molarity = component_moles / solution_volume",
    inputs={
        "component_moles": "Mapping of component identifiers to unit-aware component mole amounts.",
        "solution_volume": "Unit-aware volume of the solution.",
        "output_unit": "Molarity unit used to normalize component moles and solution volume.",
        "components": "Optional component definitions used to resolve and order component identifiers.",
        "component_key": "Optional component key used for identifier matching."
    },
    outputs={
        "component_molarity": "Mapping of resolved component identifiers to molarity values in output_unit."
    },
    aliases=(
        "molar concentration",
        "amount concentration",
    ),
    notes=(
        "The output_unit must be a ratio such as mol/L with amount in the numerator and volume in the denominator.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "component_molarity",
        "molarity",
        "mapping",
        "component_wise",
        "unit_aware",
        "unit_conversion",
        "component_key",
        "component_ordering",
    )
)
def calc_molarities_from_props(
    component_moles: Mapping[str, CustomProp],
        solution_volume: CustomProp,
        output_unit: str = 'mol/L',
        components: Optional[Sequence[Component]] = None,
        component_key: Optional[ComponentKey] = None,
        case_sensitive: bool = True,
        sort_by_components_order: bool = True,
        *,
        name: str = "molarity",
        description: str = "Calculate the molarity of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molarity values for mapped components with units.

    Parameters
    ----------
    component_moles : Mapping[str, CustomProp]
        A mapping of component identifiers to their mole amounts with units.
    solution_volume : CustomProp
        The solution volume with units.
    output_unit : str, optional
        The output molarity unit. Defaults to ``"mol/L"``.
    components : Optional[Sequence[Component]], optional
        Components used to resolve and order the component identifiers.
    component_key : Optional[ComponentKey], optional
        The component key used for mapping identifiers.
    case_sensitive : bool, optional
        Whether component identifier matching is case sensitive.
    sort_by_components_order : bool, optional
        Whether to order results according to ``components``.
    name : str, optional
        The name stored in the annotated result. Defaults to ``"molarity"``.
    description : str, optional
        The description stored in the annotated result.
    unit : str, optional
        The unit stored in the annotated result. Must match ``output_unit``.
    symbol : str, optional
        The symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        The calculated molarity values with metadata.
    """
    # SECTION: set default unit for output if not provided
    if unit is None:
        unit = output_unit

    # check unit & output unit consistency
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    res = _calc_molarities_from_props(
        component_moles=component_moles,
        solution_volume=solution_volume,
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
        unit=output_unit,  # ! set output unit
        symbol=symbol,
        implementation="_calc_molarities_from_props"
    )


# all
__all__ = [
    "calc_molarities",
    "calc_molarities_from_sequence",
    "calc_molarities_from_mapping",
    "calc_molarities_from_props",
]
