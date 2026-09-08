# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional, cast
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp, AnnotatedValue
from pythermodb_settings.utils import (
    config_components_values,
    to_annotated_value,
)
from pythermodb_settings.decorators import calculation_info
# locals
from ..utils.conversions import _to_moles, _to_units, _to_volume

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Helper functions
# ======================================================================
def _validate_moles_and_volume(
    moles: NDArray[np.float64],
    volume: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Validate the shapes and values of moles and volume arrays.

    Parameters
    ----------
    moles : NDArray[np.number]
        Array of component moles.
    volume : NDArray[np.number]
        Array of solution volumes.

    Returns
    -------
    NDArray[np.number]
        Validated volume array. A 1-D per-state volume for 2-D moles is
        returned as ``(n_states, 1)`` so NumPy broadcasts across components.

    Raises
    ------
    ValueError
        If the shapes of moles and volume are incompatible or if any volume is zero.
    """
    if moles.ndim not in (1, 2):
        raise ValueError("component_moles must be a 1-D or 2-D array.")

    if volume.ndim > 2:
        raise ValueError(
            "solution_volume must be a scalar, 1-D array, or 2-D array."
        )

    if moles.ndim == 1:
        if volume.ndim not in (0, 1):
            raise ValueError(
                "For 1-D component_moles, solution_volume must be a scalar "
                "or have the same shape as component_moles."
            )
        if volume.ndim == 1 and volume.shape != moles.shape:
            raise ValueError(
                "For 1-D component_moles, solution_volume must be a scalar "
                "or have the same shape as component_moles."
            )

    elif moles.ndim == 2:
        if volume.ndim == 0:
            pass
        elif volume.ndim == 1:
            if volume.shape[0] != moles.shape[0]:
                raise ValueError(
                    "For 2-D component_moles, 1-D solution_volume must have "
                    "one entry per state."
                )
            volume = volume[:, None]

        elif volume.ndim == 2:
            if volume.shape not in (moles.shape, (moles.shape[0], 1)):
                raise ValueError(
                    "For 2-D component_moles, solution_volume must be a "
                    "scalar, have shape (n_states,), shape (n_states, 1), "
                    "or the same shape as component_moles."
                )

    # NOTE: volume must be finite and greater than zero
    # REVIEW
    if not np.all(np.isfinite(volume)):
        raise ValueError("solution_volume must contain finite values.")

    if np.any(volume <= 0):
        raise ValueError("solution_volume must be greater than zero.")

    return volume

# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

# ! ::: Molarity from array-like inputs


def _calc_molarities(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solution_volume: float | int | NDArray[np.number],
) -> NDArray[np.float64]:
    """
    Calculate molarities using NumPy vectorization.

    Parameters
    ----------
    component_moles : Sequence[float | int] | NDArray[np.number]
        Component mole amounts. May be a Python sequence or a
        1-D/2-D NumPy array.
    solution_volume : float | int | NDArray[np.number]
        Solution volume. Must be a scalar or have the same shape as
        ``component_moles``.

    Returns
    -------
    NDArray[np.floating]
        Molarities with the same shape as ``component_moles``.
    """
    # set
    moles: NDArray[np.float64] = np.asarray(component_moles, dtype=np.float64)
    volume: NDArray[np.float64] = np.asarray(solution_volume, dtype=np.float64)

    # validate
    volume = _validate_moles_and_volume(moles, volume)

    return cast(NDArray[np.float64], moles / volume)

# ! ::: Molarity from sequence


def _calc_molarities_from_sequence(
        component_moles: Sequence[float | int],
        solution_volume: float,
) -> list[float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume.

    Parameters
    ----------
    component_moles : Sequence[float | int]
        A sequence of moles for each component.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    list[float]
        A list of molarity values for each component.
    """
    # calc
    return _calc_molarities(
        component_moles,
        solution_volume
    ).tolist()


# ! ::: Molarity from mapping


def _calc_molarities_from_mapping(
    component_moles: Mapping[str, float | int],
    solution_volume: float,
) -> dict[str, float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume.

    Parameters
    ----------
    component_moles : Mapping[str, float | int]
        A mapping of component names to their respective moles.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molarity values.
    """
    # calc
    molarities_ = _calc_molarities(
        component_moles=list(component_moles.values()),
        solution_volume=solution_volume
    )

    # to dict
    component_molarity_dict = dict(
        zip(component_moles.keys(), molarities_.tolist())
    )

    return component_molarity_dict

# ! ::: Molarity with solution volume as CustomProp


def _calc_molarities_from_props(
    component_moles: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
) -> dict[str, float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume as a CustomProp. The default
    volume unit is litre (L).

    Parameters
    ----------
    component_moles : Mapping[str, CustomProp]
        A mapping of component names to their respective mole amounts.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output molarity values. Defaults to 'mol/L'.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molarity values.

    Notes
    -----
    - The solution volume is expected to be provided as a CustomProp object. If the output_unit is not specified, it defaults to mol/L.
    - Component mole values are expected to be CustomProp objects so their units can be converted to the mole unit from output_unit.
    """
    # SECTION: set default units for moles and volume
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    volume_unit = units_[1]

    # NOTE: component moles
    # ! convert component moles to the specified unit if necessary
    component_moles_dict: dict[str, float] = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
    )

    # NOTE: solution volume unit should match output unit denominator
    # ! convert solution volume to the specified unit
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit
    )

    # SECTION: calculate molarity for each component
    return _calc_molarities_from_mapping(
        component_moles=component_moles_dict,
        solution_volume=solution_volume_scalar,
    )

# ! ::: Molarity with component ID mapping and sorting


def _calc_component_molarities_from_props(
    component_moles: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume as a CustomProp. The default
    volume unit is litre (L). The component molarity list and dictionary will be ordered according to the components list if sort_by_components_order is True.

    Parameters
    ----------
    component_moles : Mapping[str, CustomProp]
        A mapping of component names to CustomProp objects representing the moles.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output molarity values. Defaults to 'mol/L'.
    components : Optional[Sequence[Component]], optional
        A sequence of Component objects to map the component moles to, by default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component moles to components, by default None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default True.
    sort_by_components_order : bool, optional
        Whether to sort the component molarities by the order of components, by default True.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molarity values, or None if the calculation could not be performed.
    """
    # SECTION: Unit validation
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    volume_unit = units_[1]

    # SECTION: convert component moles to moles if they are CustomProp objects
    # ! component moles
    component_moles_dict = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
    )

    # ! volume
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit
    )

    # SECTION: get component values
    if component_key is not None:
        if not components:
            logger.error(
                "Component key is provided but components list is empty."
            )
            raise ValueError(
                "Component key is provided but components list is empty."
            )

        component_values: tuple[
            dict[str, Any],
            list[Any]
        ] | None = config_components_values(
            values=component_moles_dict,
            components=list(components),
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order
        )
        # >> check
        if component_values is None:
            logger.error(
                "Failed to configure component values."
            )
            raise ValueError(
                "Failed to configure component values."
            )

        # unpack
        component_values_dict, _ = component_values
    else:
        component_values_dict = component_moles_dict

    # SECTION: calculate molarity for each component
    return _calc_molarities_from_mapping(
        component_moles=component_values_dict,
        solution_volume=solution_volume_scalar,
    )

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
def _molarity_annotated(
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
def _molarity_1_annotated(
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
        _calc_molarities_from_sequence(
            component_moles=component_moles,
            solution_volume=solution_volume
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
def _molarity_2_annotated(
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
    name="molarity",
    description="Calculate keyed molarity values from unit-aware component moles and solution volume.",
    equation="molarity = component_moles / solution_volume",
    inputs={
        "component_moles": "Mapping of component identifiers to unit-aware component mole amounts.",
        "solution_volume": "Unit-aware volume of the solution.",
        "output_unit": "Molarity unit used to normalize component moles and solution volume."
    },
    outputs={
        "molarity": "Mapping of component identifiers to molarity values in output_unit."
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
        "molarity",
        "mapping",
        "unit_aware",
        "component_aware",
        "unit_conversion",
    )
)
def _molarity_3_annotated(
    component_moles: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    *,
    name: str = "molarity",
    description: str = "Calculate the molarity of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molarity values from unit-aware component moles.

    Parameters
    ----------
    component_moles : Mapping[str, CustomProp]
        A mapping of component names to their mole amounts with units.
    solution_volume : CustomProp
        The solution volume with units.
    output_unit : str, optional
        The output molarity unit. Defaults to ``"mol/L"``.
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

    return to_annotated_value(
        _calc_molarities_from_props(
            component_moles=component_moles,
            solution_volume=solution_volume,
            output_unit=output_unit
        ),
        name=name,
        description=description,
        unit=output_unit,  # ! set output unit
        symbol=symbol,
        implementation="_calc_molarities_from_props"
    )


# ::: annotated for component mapping with custom properties


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
def _molarity_4_annotated(
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

    res = _calc_component_molarities_from_props(
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
        implementation="_calc_component_molarities_from_props"
    )


# ======================================================================
# *** Aliases
# ======================================================================
# >> molarities
calc_molarities = _molarity_annotated

# >> molarities from sequence
calc_molarities_from_sequence = _molarity_1_annotated

# >> molarities from mapping
calc_molarities_from_mapping = _molarity_2_annotated

# >> molarities with units
calc_molarities_from_props = _molarity_3_annotated

# >> component molarities
calc_component_molarities_from_props = _molarity_4_annotated

# all
__all__ = [
    "_calc_molarities",
    "calc_molarities",
    "_calc_molarities_from_sequence",
    "calc_molarities_from_sequence",
    "_calc_molarities_from_mapping",
    "calc_molarities_from_mapping",
    "_calc_molarities_from_props",
    "calc_molarities_from_props",
    "_calc_component_molarities_from_props",
    "calc_component_molarities_from_props",
]
