# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional
from pythermodb_settings.models import Component, ComponentKey, CustomProp, AnnotatedValue
from pythermodb_settings.utils import (
    config_components_values,
    to_annotated_value,
)
from pythermodb_settings.decorators import annotated_value
# locals
from ..utils.conversions import _to_moles, _to_units, _to_volume
# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

# ! ::: Molarity from sequence

def _calc_molarities_from_sequence(
        component_moles: Sequence[float],
        solution_volume: float,
) -> list[float]:
    """
    Calculate the molarity of each component in a solution given the component moles and the solution volume.

    Parameters
    ----------
    component_moles : Sequence[float]
        A sequence of moles for each component.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    list[float]
        A list of molarity values for each component.
    """
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")
    return [moles / solution_volume for moles in component_moles]

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
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")

    # ! dict
    component_molarity_dict = {
        key: value / solution_volume for key, value in component_moles.items()
    }

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

# ::: annotated for sequence


def _molarity_1_annotated(
        component_moles: Sequence[float],
        solution_volume: float,
        *,
        name: str = "molarity",
        description: str = "Calculate the molarity of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None
) -> AnnotatedValue[list[float]]:
    return to_annotated_value(
        _calc_molarities_from_sequence(
            component_moles=component_moles,
            solution_volume=solution_volume
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol
    )

# ::: annotated for mapping


def _molarity_2_annotated(
        component_moles: Mapping[str, float | int],
        solution_volume: float,
        *,
        name: str = "molarity",
        description: str = "Calculate the molarity of each component in a solution.",
        unit: str | None = None,
        symbol: str | None = None
) -> AnnotatedValue[dict[str, float]]:
    return to_annotated_value(
        _calc_molarities_from_mapping(
            component_moles=component_moles,
            solution_volume=solution_volume
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol
    )


# ::: annotated for mapping with custom properties

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
        symbol=symbol
    )


# ::: annotated for component mapping with custom properties


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
        symbol=symbol
    )


# ======================================================================
# *** Aliases
# ======================================================================
# >> molarities
calc_molarities_from_sequence = _molarity_1_annotated

# >> molarities from mapping
calc_molarities_from_mapping = _molarity_2_annotated

# >> molarities with units
calc_molarities_from_props = _molarity_3_annotated

# >> component molarities
calc_component_molarities_from_props = _molarity_4_annotated


# all
__all__ = [
    "_calc_molarities_from_sequence",
    "calc_molarities_from_sequence",
    "_calc_molarities_from_mapping",
    "calc_molarities_from_mapping",
    "_calc_molarities_from_props",
    "calc_molarities_from_props",
    "_calc_component_molarities_from_props",
    "calc_component_molarities_from_props",
]
