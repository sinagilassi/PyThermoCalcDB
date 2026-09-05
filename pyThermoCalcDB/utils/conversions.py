# import libs
import logging
from collections.abc import Mapping
from typing import List, Optional, Dict, Any
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.models import Temperature, CustomProp, ComponentMoles, UnitConversionFn, Component, ComponentKey, ScalarValue
from pythermodb_settings.utils.quantity import to_amounts, to_custom_props_mapping, to_custom_prop_scalar, pos, to_scalar
from pycuc import convert_from_to
# locals

# NOTE: logger setup
logger = logging.getLogger(__name__)


# SECTION: Unit conversion function resolver
def _resolve_unit_conversion_fn(
    unit_conversion_fn: UnitConversionFn | None,
) -> UnitConversionFn:
    """Return the provided converter or the module default converter."""
    return convert_from_to if unit_conversion_fn is None else unit_conversion_fn


# SECTION: Unit handling helpers

# ! ::: Helper function to split unit string into individual units


def _to_units(unit: str) -> List[str]:
    # NOTE: validation
    if not unit or not isinstance(unit, str):
        raise ValueError("Invalid unit string provided.")

    # NOTE: only support units separated by '/'
    if '/' not in unit:
        raise ValueError("Unit string must contain '/' to separate units.")

    # NOTE: more than one '/'
    if unit.count('/') != 1:
        raise ValueError(
            "Unit string must contain exactly one '/' to separate units."
        )

    # NOTE: split the unit string by '/' to get individual units
    return [unit_.strip() for unit_ in unit.strip().split('/')]

# SECTION: Internal helpers
# ! ::: Convert energy value to J/mol


def _to_J__mol(
    value: float,
    from_unit: str,
    **kwargs
) -> float:
    """
    Convert energy value to J/mol.

    Parameters
    ----------
    value : float
        The energy value to convert.
    from_unit : str
        The unit of the input energy value.

    Returns
    -------
    float
        The energy value in J/mol.
    """
    try:
        converted_value = convert_from_to(
            value=value,
            from_unit=from_unit,
            to_unit="J/mol",
        )
        return converted_value
    except Exception as e:
        logger.error(f"Error converting energy to J/mol: {e}")
        raise

# ! ::: Convert temperature value to K


def _to_kelvin(temperature: Temperature) -> float:
    """Return temperature value in K."""
    T_value = temperature.value
    T_unit = temperature.unit.strip()
    if T_unit != "K":
        T_value = convert_from_to(
            value=T_value,
            from_unit=T_unit,
            to_unit="K"
        )
    return float(T_value)

# ! ::: Convert energy value to g/mol


def to_g_mol(
    value: float,
    from_unit: str,
    **kwargs
) -> Optional[float]:
    """
    Convert energy value to g/mol.

    Parameters
    ----------
    value : float
        The energy value to convert.
    from_unit : str
        The unit of the input energy value.

    Returns
    -------
    float
        The energy value in g/mol.
    """
    try:
        converted_value = convert_from_to(
            value=value,
            from_unit=from_unit,
            to_unit="g/mol",
        )
        return converted_value
    except Exception as e:
        logger.error(f"Error converting energy to g/mol: {e}")
        return None


# ! ::: Convert component moles to the requested output unit
def _to_moles(
        component_moles: ComponentMoles,
        output_unit: Optional[str] = None,
        unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Dict[str, float]:
    """
    Convert a dictionary of component moles to float values.

    Parameters
    ----------
    component_moles : ComponentMoles
        A dictionary mapping component names to mole amounts. Numeric values are assumed to already be in output_unit.
    output_unit : str, optional
        The unit to which CustomProp component moles should be converted. Default is None.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping component names to their respective moles as floats.
    """
    # NOTE: resolver
    return to_custom_props_mapping(
        values=component_moles,
        to_unit=output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )


# ! ::: Convert to desired volume unit
def _to_volume(
        solution_volume: CustomProp,
        output_unit: Optional[str] = None,
        unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> float:
    """
    Convert a solution volume defined as a CustomProp object to a float value in the desired unit.

    Parameters
    ----------
    solution_volume : CustomProp
        The solution volume as a CustomProp object.
    output_unit : str, optional
        The unit to which the solution volume should be converted. Default is None.

    Returns
    -------
    float
        The solution volume in the desired unit as a float.
    """
    return to_custom_prop_scalar(
        prop=solution_volume,
        output_unit=output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )


# ! ::: Convert to desired mass unit
def _to_mass(
        solvent_mass: CustomProp,
        output_unit: Optional[str] = None,
        unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> float:
    """
    Convert a solvent mass defined as a CustomProp object to a float value in the desired unit.

    Parameters
    ----------
    solvent_mass : CustomProp
        The solvent mass as a CustomProp object.
    output_unit : str, optional
        The unit to which the solvent mass should be converted. Default is None.

    Returns
    -------
    float
        The solvent mass in the desired unit as a float.
    """
    return to_custom_prop_scalar(
        prop=solvent_mass,
        output_unit=output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )

# ! ::: COnvert to desired molecular weight unit


def _to_molecular_weight(
        molecular_weight: CustomProp,
        output_unit: Optional[str] = None,
        unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> float:
    """
    Convert a molecular weight defined as a CustomProp object to a float value in the desired unit.

    Parameters
    ----------
    molecular_weight : CustomProp
        The molecular weight as a CustomProp object.
    output_unit : str, optional
        The unit to which the molecular weight should be converted. Default is None.

    Returns
    -------
    float
        The molecular weight in the desired unit as a float.
    """
    return to_custom_prop_scalar(
        prop=molecular_weight,
        output_unit=output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )

# SECTION: Internal helpers


def _scalar(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a scalar input to float, optionally normalizing units."""
    return to_scalar(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


def _pos(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a scalar input to a positive float, optionally normalizing units."""
    return pos(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )

# ! ::: Configure component values with name


def _configure_component_values(
    values: Dict[str, Any] | Mapping[str, Any],
    components: Optional[List[Component]],
    component_key: Optional[ComponentKey],
    case_sensitive: bool,
    sort_by_components_order: bool,
    name: str,
) -> dict[str, float]:
    """Remap and order mapping values using component metadata when requested.

    Notes
    -----
    This helper does not change value units. Any unit normalization must happen
    before values are passed here.
    """
    # ! When no component key is requested, preserve caller mapping keys/order.
    if component_key is None:
        return dict(values)

    # NOTE: Component metadata is required only for key remapping.
    if not components:
        raise ValueError(
            f"component_key is provided but components is empty for {name}.")

    # SECTION: Remap values through pythermodb-settings utilities
    component_values = config_components_values(
        values=values,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    if component_values is None:
        raise ValueError(f"Failed to configure {name} component values.")
    component_values_dict, _ = component_values
    return component_values_dict

# ! ::: Validate matching component keys


def _validate_same_keys(left: Mapping[str, float], right: Mapping[str, float]) -> None:
    """Validate matching component keys.

    Notes
    -----
    This helper performs only key validation and returns ``None``. It has no
    calculated value or return unit.
    """
    # ? Mismatched keys usually indicate a missing concentration or charge.
    if set(left) != set(right):
        raise ValueError(
            "concentrations and charges must have the same component keys.")
