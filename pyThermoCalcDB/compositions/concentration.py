# import libs
import logging
from typing import Any, Dict, List, Optional, Tuple
from pythermodb_settings.models import Component, ComponentKey, CustomProp, UnitConversionFn
from pythermodb_settings.utils import config_components_values, to_amounts
# locals
from ..models import ComponentAmounts
from ..utils.conversions import _to_units, _to_volume, _resolve_unit_conversion_fn

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ! ::: Concentration [amount/volume]
def concentration_amount_volume_1(
        component_amounts: List[float],
        solution_volume: float,
) -> List[float]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume.

    Parameters
    ----------
    component_amounts : List[float]
        A list of amounts for each component.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    List[float]
        A list of concentration values for each component.
    """
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")

    return [amount / solution_volume for amount in component_amounts]


def concentration_amount_volume_2(
    component_amounts: Dict[str, float | int],
    solution_volume: float,
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume.

    Parameters
    ----------
    component_amounts : Dict[str, float | int]
        A dictionary mapping component names to their respective amounts.
    solution_volume : float
        The volume of the solution.

    Returns
    -------
    Tuple[Dict[str, float], List[float]]
        A tuple containing a dictionary of component concentrations and a list of concentration values.
    """
    # check
    if solution_volume == 0:
        logger.error("Volume of the solution cannot be zero.")
        raise ValueError("Volume of the solution cannot be zero.")

    # ! dict
    component_concentration_dict = {
        key: value / solution_volume for key, value in component_amounts.items()
    }
    # ! list
    component_concentration_list = list(component_concentration_dict.values())

    return component_concentration_dict, component_concentration_list


# ! ::: Concentration [amount/volume] with solution volume as CustomProp
def concentration_amount_volume_3(
    component_amounts: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Tuple[Dict[str, float], List[float]]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume as a CustomProp.

    Parameters
    ----------
    component_amounts : ComponentAmounts
        A dictionary mapping component names to their respective amounts. Numeric values are assumed to already be in the numerator unit from output_unit.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output concentration values. Defaults to 'kg/m^3'.
    unit_conversion_fn : UnitConversionFn, optional
        The function to use for unit conversion. Defaults to None. Then it will use the default conversion function `pycuc.convert_from_to`.

    Returns
    -------
    Tuple[Dict[str, float], List[float]]
        A tuple containing a dictionary of component concentrations and a list of concentration values.
    """
    # SECTION: set default units for amount and volume
    units_ = _to_units(output_unit)
    # >> set
    amount_unit = units_[0]
    volume_unit = units_[1]

    # NOTE: component amounts
    component_amounts_dict: Dict[str, float] = to_amounts(
        component_amounts=component_amounts,
        output_unit=amount_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )

    # NOTE: solution volume unit should match output unit denominator
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )

    # SECTION: calculate concentration for each component
    component_concentration = concentration_amount_volume_2(
        component_amounts=component_amounts_dict,
        solution_volume=solution_volume_scalar,
    )
    # unpack the result
    component_concentration_dict, component_concentration_list = component_concentration

    return component_concentration_dict, component_concentration_list


# ! ::: Concentration [amount/volume] with component ID mapping and sorting
def concentration_amount_volume_4(
    component_amounts: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    """
    Calculate the concentration of each component in a solution given component amounts and solution volume as a CustomProp. The component concentration list and dictionary will be ordered according to the components list if sort_by_components_order is True.

    Parameters
    ----------
    component_amounts : ComponentAmounts
        A dictionary mapping component names to their respective amounts or CustomProp objects representing the amounts.
    solution_volume : CustomProp
        The volume of the solution as a CustomProp object.
    output_unit : str, optional
        The unit for the output concentration values. Defaults to 'kg/m^3'.
    components : Optional[List[Component]], optional
        A list of Component objects to map the component amounts to, by default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component amounts to components, by default None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default True.
    sort_by_components_order : bool, optional
        Whether to sort the component concentrations by the order of components, by default True.
    unit_conversion_fn : UnitConversionFn, optional
        The function to use for unit conversion. Defaults to None. Then it will use the default conversion function `pycuc.convert_from_to`.

    Returns
    -------
    Optional[Tuple[Dict[str, float], List[float]]]
        A tuple containing a dictionary of component concentrations and a list of concentration values, or None if the calculation could not be performed.
    """
    # SECTION: Unit validation
    units_ = _to_units(output_unit)
    # >> set
    amount_unit = units_[0]
    volume_unit = units_[1]

    # SECTION: convert component amounts if they are CustomProp objects
    # ! component amounts
    component_amounts_dict = to_amounts(
        component_amounts=component_amounts,
        output_unit=amount_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )

    # ! volume
    solution_volume_scalar = _to_volume(
        solution_volume=solution_volume,
        output_unit=volume_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )

    # SECTION: get component values
    if component_key is not None:
        if not components:
            logger.warning(
                "Component key is provided but components list is empty.")
            return None

        component_values: Tuple[
            Dict[str, Any],
            List[Any]
        ] | None = config_components_values(
            values=component_amounts_dict,
            components=components,
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order
        )
        # >> check
        if component_values is None:
            logger.warning("Failed to configure component values")
            return None

        # unpack
        component_values_dict, _ = component_values
    else:
        component_values_dict = component_amounts_dict

    # SECTION: calculate concentration for each component
    component_concentration = concentration_amount_volume_2(
        component_amounts=component_values_dict,
        solution_volume=solution_volume_scalar,
    )
    # unpack
    component_concentration_dict, component_concentration_list = component_concentration

    return component_concentration_dict, component_concentration_list


# ! ::: Mass concentration wrappers
def mass_concentration1(
        component_mass: List[float],
        solution_volume: float,
) -> List[float]:
    return concentration_amount_volume_1(
        component_amounts=component_mass,
        solution_volume=solution_volume,
    )


def mass_concentration2(
    component_mass: Dict[str, float | int],
    solution_volume: float,
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_2(
        component_amounts=component_mass,
        solution_volume=solution_volume,
    )


def mass_concentration3(
    component_mass: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_3(
        component_amounts=component_mass,
        solution_volume=solution_volume,
        output_unit=output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )


def mass_concentration4(
    component_mass: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'kg/m^3',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    return concentration_amount_volume_4(
        component_amounts=component_mass,
        solution_volume=solution_volume,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )


# ! ::: Molar concentration wrappers
def molar_concentration1(
        component_moles: List[float],
        solution_volume: float,
) -> List[float]:
    return concentration_amount_volume_1(
        component_amounts=component_moles,
        solution_volume=solution_volume,
    )


def molar_concentration2(
    component_moles: Dict[str, float | int],
    solution_volume: float,
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_2(
        component_amounts=component_moles,
        solution_volume=solution_volume,
    )


def molar_concentration3(
    component_moles: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Tuple[Dict[str, float], List[float]]:
    return concentration_amount_volume_3(
        component_amounts=component_moles,
        solution_volume=solution_volume,
        output_unit=output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )


def molar_concentration4(
    component_moles: ComponentAmounts,
    solution_volume: CustomProp,
    output_unit: str = 'mol/L',
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Optional[Tuple[Dict[str, float], List[float]]]:
    return concentration_amount_volume_4(
        component_amounts=component_moles,
        solution_volume=solution_volume,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )


# SECTION: Aliases
# ! generic concentration
calculate_concentrations = concentration_amount_volume_1
calculate_keyed_concentrations = concentration_amount_volume_2
calculate_component_concentrations = concentration_amount_volume_4
