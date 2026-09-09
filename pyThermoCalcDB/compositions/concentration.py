# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional, Dict, List, Optional
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp, AnnotatedValue
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.utils import (
    to_annotated_value,
)

# old
from pythermodb_settings.models import Component, ComponentKey, CustomProp, UnitConversionFn
# locals
from .core.concentration import (
    _calc_concentrations,
    _calc_concentrations_from_sequence,
    _calc_concentrations_from_mapping,
    _calc_concentrations_from_props,
    _calc_component_concentrations_from_props,
)

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

# ! ::: Mass concentration
def calc_mass_concentrations(
        component_amounts: Sequence[float | int] | NDArray[np.number],
        solution_volume: float | int | NDArray[np.number],
) -> NDArray[np.float64]:
    """Calculate numeric mass concentrations from mass sequence and volume."""
    return _calc_concentrations(
        component_amounts=component_amounts,
        solution_volume=solution_volume,
    )


def calc_mass_concentrations_from_sequence(
        component_mass: Sequence[float | int],
        solution_volume: float | int,
) -> List[float]:
    """Calculate numeric mass concentrations from mass sequence and volume."""
    return _calc_concentrations_from_sequence(
        component_amounts=component_mass,
        solution_volume=solution_volume,
    )


def calc_mass_concentrations_from_mapping(
    component_mass: Dict[str, float | int],
    solution_volume: float,
) -> Dict[str, float]:
    """Calculate numeric mass concentrations from keyed masses and volume."""
    return _calc_concentrations_from_mapping(
        component_amounts=component_mass,
        solution_volume=solution_volume,
    )


def calc_mass_concentrations_from_props(
    component_mass: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
) -> Dict[str, float]:
    """Calculate unit-aware mass concentrations on the requested output basis."""
    return _calc_concentrations_from_props(
        component_amounts=component_mass,
        solution_volume=solution_volume,
        output_unit=output_unit,
    )


def calc_component_mass_concentration_from_props(
    component_mass: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Dict[str, float]:
    """Calculate component-aware unit-normalized mass concentrations."""
    return _calc_component_concentrations_from_props(
        component_amounts=component_mass,
        solution_volume=solution_volume,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
        unit_conversion_fn=unit_conversion_fn,
    )


# ! ::: Molar concentration
def calc_molar_concentrations(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solution_volume: float | int | NDArray[np.number],
) -> NDArray[np.float64]:
    """Calculate numeric molar concentrations from mole sequence and volume."""
    return _calc_concentrations(
        component_amounts=component_moles,
        solution_volume=solution_volume,
    )


def calc_molar_concentrations_from_sequence(
        component_moles: Sequence[float],
        solution_volume: float,
) -> List[float]:
    """Calculate numeric molar concentrations from mole sequence and volume."""
    return _calc_concentrations_from_sequence(
        component_amounts=component_moles,
        solution_volume=solution_volume,
    )


def calc_molar_concentrations_from_mapping(
    component_moles: Dict[str, float | int],
    solution_volume: float,
) -> Dict[str, float]:
    """Calculate numeric molar concentrations from keyed moles and volume."""
    return _calc_concentrations_from_mapping(
        component_amounts=component_moles,
        solution_volume=solution_volume,
    )


def calc_molar_concentrations_from_props(
    component_moles: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Dict[str, float]:
    """Calculate unit-aware molar concentrations on the requested output basis."""
    return _calc_concentrations_from_props(
        component_amounts=component_moles,
        solution_volume=solution_volume,
        output_unit=output_unit,
    )


def calc_component_molar_concentrations_from_props(
    component_moles: Mapping[str, CustomProp],
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Dict[str, float]:
    """Calculate component-aware unit-normalized molar concentrations."""
    return _calc_component_concentrations_from_props(
        component_amounts=component_moles,
        solution_volume=solution_volume,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
        unit_conversion_fn=unit_conversion_fn,
    )


# export functions
__all__ = [
    "calc_mass_concentrations",
    "calc_molar_concentrations_from_sequence",
    "calc_molar_concentrations_from_mapping",
    "calc_molar_concentrations_from_props",
    "calc_component_molar_concentrations_from_props",
    "calc_molar_concentrations",
    "calc_molar_concentrations_from_sequence",
    "calc_molar_concentrations_from_mapping",
    "calc_molar_concentrations_from_props",
    "calc_component_molar_concentrations_from_props",
]
