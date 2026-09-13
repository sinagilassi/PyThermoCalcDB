# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import List, Optional, Dict, Any, cast, TypeAlias, Literal, overload
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.utils import config_components_values, get_unit
from pythermodb_settings.models import Temperature, CustomProp, ComponentMoles, UnitConversionFn, Component, ComponentKey, ScalarValue
from pythermodb_settings.utils.quantity import to_amounts, to_custom_props_mapping, to_custom_prop_scalar, pos, to_scalar, to_dict, to_values
from pycuc import convert_from_to
# locals

# NOTE: logger setup
logger = logging.getLogger(__name__)

# SECTION: Numeric helper aliases
NumericArrayInput: TypeAlias = \
    float | int | Sequence[float | int] | NDArray[np.number]


# SECTION: Unit conversion function resolver
def _resolve_unit_conversion_fn(
    unit_conversion_fn: UnitConversionFn | None,
) -> UnitConversionFn:
    """Return the provided converter or the module default converter."""
    return convert_from_to if unit_conversion_fn is None else unit_conversion_fn


def _iter_values(
    values: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
):
    """Yield scalar values from mapping or sequence component input."""
    return values.values() if isinstance(values, Mapping) else values


def _contains_custom_prop(
    values: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
) -> bool:
    """Return True when component values carry explicit unit metadata."""
    return any(isinstance(value, CustomProp) for value in _iter_values(values))


def _all_custom_props(values: Mapping[str, float | int | CustomProp]) -> bool:
    """Return True when every mapping value is a ``CustomProp``."""
    # ? Public wrappers use this to choose props-only adapters.
    return all(isinstance(value, CustomProp) for value in values.values())


def _validate_custom_prop_mapping(values: Mapping[str, object], name: str) -> None:
    """Validate that all mapping values are ``CustomProp`` objects."""
    # ! Props adapters accept only fully unit-aware CustomProp mappings.
    if not all(isinstance(value, CustomProp) for value in values.values()):
        raise TypeError(f"{name} must be a mapping of CustomProp values.")


def _validate_custom_prop_scalar(value: object, name: str) -> None:
    """Validate that a scalar props-adapter input is a ``CustomProp`` object."""
    # ! Scalar props adapters accept only unit-aware CustomProp scalars.
    if not isinstance(value, CustomProp):
        raise TypeError(f"{name} must be a CustomProp value.")


def _resolve_result_unit(
    identifier: str,
    values: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_unit: str | None,
) -> str | None:
    """Resolve an annotated result unit from unit-carrying component values.

    Numeric values do not carry source units, so the function cannot prove or
    perform any conversion for them. Unit annotation is therefore managed only
    when at least one component value is a ``CustomProp``.
    """
    if not _contains_custom_prop(values):
        return None

    if output_unit is not None:
        return output_unit

    unit_info = get_unit(identifier=identifier, data=values)
    if not unit_info["consistent"]:
        raise ValueError(
            f"{identifier} CustomProp values must have consistent units when "
            "no output unit is provided."
        )

    return str(unit_info["unit"]) if unit_info["unit"] is not None else None


# SECTION: NumPy numeric helpers

def _as_float_array(
    values: NumericArrayInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert an array-like numeric input to a finite float64 NumPy array."""
    # NOTE: Core kernels operate on finite NumPy float arrays.
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 0:
        raise ValueError(f"{name} must contain component values.")
    if arr.ndim > 2:
        raise ValueError(f"{name} must be one- or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _return_scalar_if_zero_dim(
    value: float | np.float64 | NDArray[np.float64],
) -> float | NDArray[np.float64]:
    """Return a Python float for 0-D results and a float64 array otherwise."""
    # NOTE: Preserves scalar return behavior for component reductions.
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim == 0:
        return float(arr)
    return cast(NDArray[np.float64], arr)


def _validate_fraction_array(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative fractions that close along the last axis."""
    # ! Fractions are component-wise and must close along the last axis.
    if np.any(values < 0.0):
        raise ValueError(f"{name} must be non-negative.")
    totals = np.sum(values, axis=-1)
    if not np.allclose(totals, 1.0):
        raise ValueError(f"{name} must sum to 1.0 along the component axis.")


def _validate_non_negative_array(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative numeric array values."""
    # ! Negative component values are outside the supported physical domain.
    if np.any(values < 0.0):
        raise ValueError(f"{name} must be non-negative.")


def _validate_positive_array(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive numeric array values."""
    # ! Denominators and positive physical properties cannot be zero or negative.
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


def _validate_positive_scalar(value: float | int, name: str) -> float:
    """Validate and return a strictly positive finite scalar."""
    # ! Scalar physical inputs cannot be NaN, infinite, zero, or negative.
    scalar = float(value)
    if not np.isfinite(scalar) or scalar <= 0.0:
        raise ValueError(f"{name} must be greater than zero.")
    return scalar


def _validate_same_array_shape(
    left: NDArray[np.float64],
    right: NDArray[np.float64],
    left_name: str,
    right_name: str,
) -> None:
    """Validate identical array shapes for pairwise component calculations."""
    # ? Pairwise mixture rules need component arrays aligned by shape.
    if left.shape != right.shape:
        raise ValueError(
            f"{left_name} and {right_name} must have the same shape.")


def _validate_same_mapping_keys(
    left: Mapping[str, object],
    right: Mapping[str, object],
    left_name: str,
    right_name: str,
) -> None:
    """Validate matching component keys for pairwise mapping calculations."""
    # ? Mapping adapters must preserve component identity before array math.
    if set(left) != set(right):
        raise ValueError(
            f"{left_name} and {right_name} must have the same component keys.")


# SECTION: Unit handling helpers

# ! ::: get all CustomProp instances from a collection

@overload
def _get_all_custom_props(
    collection: object,
    return_type: Literal["list"],
) -> List[CustomProp]: ...


@overload
def _get_all_custom_props(
    collection: object,
    return_type: Literal["mapping"],
) -> Mapping[str, CustomProp]: ...


def _get_all_custom_props(
        collection: object,
        return_type: Literal["list", "mapping"],
) -> List[CustomProp] | Mapping[str, CustomProp]:
    """Retrieve all CustomProp instances from a collection."""
    custom_props = []
    if isinstance(collection, Mapping):
        for value in collection.values():
            if isinstance(value, CustomProp):
                custom_props.append(value)
    elif isinstance(collection, Sequence) and not isinstance(collection, str):
        for item in collection:
            if isinstance(item, CustomProp):
                custom_props.append(item)
    return custom_props

# ! ::: Helper function to split unit string into individual units


def _to_units(unit: str) -> List[str]:
    """
    Split a unit string into individual units.

    Parameters
    ----------
    unit : str
        The unit string to split, expected to contain exactly one '/'.

    Returns
    -------
    List[str]
        A list containing the individual units.
    """
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

def _to_amounts(
    component_amounts: Mapping[str, float | int | CustomProp],
    output_unit: str,
    unit_conversion_fn: Optional[UnitConversionFn] = None,
) -> Dict[str, float]:
    """
    Convert a dictionary of component amounts to float values in the desired output unit.

    Parameters
    ----------
    component_amounts : Mapping[str, float | int | CustomProp]
        A dictionary mapping component names to their respective amounts or CustomProp objects representing the amounts.
    output_unit : str
        The unit to which the component amounts should be converted.
    unit_conversion_fn : UnitConversionFn, optional
        The function to use for unit conversion. Defaults to None. Then it will use the default conversion function `pycuc.convert_from_to`.

    Returns
    -------
    Dict[str, float]
        A dictionary mapping component names to their respective amounts as floats in the desired output unit.
    """
    return to_amounts(
        component_amounts=component_amounts,
        output_unit=output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn)
    )


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
    value: float | int | CustomProp,
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
    extract_values: bool = True,
) -> dict[str, float]:
    """Remap and order mapping values using component metadata when requested.

    Parameters
    ----------
    values : Dict[str, Any] | Mapping[str, Any]
        The input values keyed by component names.
    components : Optional[List[Component]]
        The list of component metadata.
    component_key : Optional[ComponentKey]
        The key to use for mapping components.
    case_sensitive : bool
        Whether the component key matching should be case-sensitive.
    sort_by_components_order : bool
        Whether to sort the output by the order of components.
    name : str
        The name of the values being configured, used for error messages.
    extract_values : bool, optional
        Whether to extract the values from the configured component mapping. Default is True.

    Returns
    -------
    dict[str, float]
        The configured component values, optionally extracted from the component mapping.

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
        logger.warning(
            f"component_key is provided but components is empty for {name}.")
        components = []

    # SECTION: Remap values through pythermodb-settings utilities
    component_values = config_components_values(
        values=values,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
        extract_values=extract_values,
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

# ! ::: to values


def _to_values(
        data: Mapping[str, float | int | CustomProp] | Dict[str, float | int | CustomProp],
        name: str,
        output_unit: str | None = None,
        unit_conversion_fn: UnitConversionFn | None = None,
):
    try:
        return to_values(
            data=data,
            output_unit=output_unit,
            unit_conversion_fn=unit_conversion_fn,
        )
    except Exception as e:
        raise ValueError(f"Failed to convert {name} to values: {e}")

# ! ::: Scalar conversion helper


def _to_scalar(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert scalar input to float, optionally normalizing units."""
    return to_scalar(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )
