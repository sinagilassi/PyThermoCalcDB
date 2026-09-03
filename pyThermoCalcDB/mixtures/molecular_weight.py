"""Mixture molecular-weight rules."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import fractions, positive, same_shape
# locals
from ..utils.conversions import _resolve_unit_conversion_fn


# SECTION: Internal helpers

def _configure_component_values(
    values: Mapping[str, float],
    components: Optional[List[Component]],
    component_key: Optional[ComponentKey],
    case_sensitive: bool,
    sort_by_components_order: bool,
    name: str,
) -> dict[str, float]:
    """Remap and order mapping values using component metadata when requested."""
    # ! When no component key is requested, preserve caller mapping keys/order.
    if component_key is None:
        return dict(values)

    # NOTE: Component metadata is required only for key remapping.
    if not components:
        raise ValueError(
            f"component_key is provided but components is empty for {name}."
        )

    # SECTION: Remap values through pythermodb-settings utilities
    component_values = config_components_values(
        values=dict(values),
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    if component_values is None:
        raise ValueError(f"Failed to configure {name} component values.")

    # NOTE: config_components_values returns both dict and ordered list.
    component_values_dict, _ = component_values
    return component_values_dict


# SECTION: Mole-fraction basis

def calc_mixture_molecular_weight_from_mole_fractions(
    mole_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molecular_weights: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate mixture molecular weight from mole fractions.

    Parameters
    ----------
    mole_fractions : mapping or sequence of float | int | CustomProp
        Component mole fractions. Fractions must be non-negative and valid
        according to the shared ``fractions`` validator.
    molecular_weights : mapping or sequence of float | int | CustomProp
        Component molecular weights. Values must be positive. ``CustomProp``
        values are converted to ``output_molecular_weight_unit`` when provided.
    output_molecular_weight_unit : str, optional
        Unit used to normalize molecular weights before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.
    components : list[Component], optional
        Component metadata used for mapping-key remapping when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format used for mapping inputs.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether mapping values should follow ``components`` order.

    Returns
    -------
    float
        Mixture molecular weight in the normalized molecular-weight unit.

    Notes
    -----
    Equation: ``M_mix = sum_i(x_i*M_i)``. No silent fraction normalization is
    performed beyond the shared package validator behavior.
    """
    # SECTION: Validate inputs
    fractions(mole_fractions, "mole_fractions")
    positive(molecular_weights, "molecular_weights")
    same_shape(mole_fractions, molecular_weights)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(mole_fractions, Mapping) and isinstance(molecular_weights, Mapping):
        x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)
        mw = to_dict(
            molecular_weights,
            output_molecular_weight_unit,
            unit_conversion_fn=conversion_fn,
        )
        x = _configure_component_values(
            x, components, component_key, case_sensitive, sort_by_components_order, "mole_fractions"
        )
        mw = _configure_component_values(
            mw, components, component_key, case_sensitive, sort_by_components_order, "molecular_weights"
        )
        return sum(x[key] * mw[key] for key in x)

    # ! Mixed mapping/sequence input is ambiguous.
    if isinstance(mole_fractions, Mapping) or isinstance(molecular_weights, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    mw = to_list(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn=conversion_fn,
    )
    return sum(x_i * mw_i for x_i, mw_i in zip(x, mw))


# SECTION: Mass-fraction basis

def calc_mixture_molecular_weight_from_mass_fractions(
    mass_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molecular_weights: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate mixture molecular weight from mass fractions.

    Parameters
    ----------
    mass_fractions : mapping or sequence of float | int | CustomProp
        Component mass fractions. Fractions must be non-negative and valid
        according to the shared ``fractions`` validator.
    molecular_weights : mapping or sequence of float | int | CustomProp
        Component molecular weights. Values must be positive. ``CustomProp``
        values are converted to ``output_molecular_weight_unit`` when provided.
    output_molecular_weight_unit : str, optional
        Unit used to normalize molecular weights before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.
    components : list[Component], optional
        Component metadata used for mapping-key remapping when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format used for mapping inputs.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether mapping values should follow ``components`` order.

    Returns
    -------
    float
        Mixture molecular weight in the normalized molecular-weight unit.

    Notes
    -----
    Equation: ``M_mix = 1/sum_i(w_i/M_i)``. No silent fraction normalization is
    performed beyond the shared package validator behavior.

    Raises
    ------
    TypeError
        If one component input is a mapping and the other is a sequence.
    ValueError
        If fractions, molecular weights, keys, or lengths are invalid.
    """
    # SECTION: Validate inputs
    fractions(mass_fractions, "mass_fractions")
    positive(molecular_weights, "molecular_weights")
    same_shape(mass_fractions, molecular_weights)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(mass_fractions, Mapping) and isinstance(molecular_weights, Mapping):
        w = to_dict(mass_fractions, unit_conversion_fn=conversion_fn)
        mw = to_dict(
            molecular_weights,
            output_molecular_weight_unit,
            unit_conversion_fn=conversion_fn,
        )
        w = _configure_component_values(
            w, components, component_key, case_sensitive, sort_by_components_order, "mass_fractions"
        )
        mw = _configure_component_values(
            mw, components, component_key, case_sensitive, sort_by_components_order, "molecular_weights"
        )
        denominator = sum(w[key] / mw[key] for key in w)
    elif isinstance(mass_fractions, Mapping) or isinstance(molecular_weights, Mapping):
        # ! Mixed mapping/sequence input is ambiguous.
        raise TypeError("Both component inputs must be mappings or both sequences.")
    else:
        # SECTION: Sequence implementation
        w = to_list(mass_fractions, unit_conversion_fn=conversion_fn)
        mw = to_list(
            molecular_weights,
            output_molecular_weight_unit,
            unit_conversion_fn=conversion_fn,
        )
        denominator = sum(w_i / mw_i for w_i, mw_i in zip(w, mw))

    # NOTE: A zero denominator means all supplied fractions are zero.
    if denominator <= 0.0:
        raise ValueError("The reciprocal molecular-weight sum must be positive.")
    return 1.0 / denominator


# SECTION: Legacy compatibility wrappers

def calc_mixture_molecular_weight_1(
    component_mole_fractions,
    component_molecular_weights,
    input_unit: str | None = None,
    output_unit: str | None = None,
) -> float:
    """Compatibility wrapper for the original list-based mole-fraction API."""
    # NOTE: Original API converted the final mixture MW when both units were set.
    result = calc_mixture_molecular_weight_from_mole_fractions(
        component_mole_fractions,
        component_molecular_weights,
    )
    if input_unit and output_unit:
        return _resolve_unit_conversion_fn(None)(result, input_unit, output_unit)
    return result


def calc_mixture_molecular_weight_2(
    component_mole_fractions,
    component_molecular_weights,
    output_unit: str = "g/mol",
) -> float:
    """Compatibility wrapper for the original keyed mole-fraction API."""
    return calc_mixture_molecular_weight_from_mole_fractions(
        component_mole_fractions,
        component_molecular_weights,
        output_molecular_weight_unit=output_unit,
    )


def calc_mixture_molecular_weight_from_mass_fractions_1(
    component_mass_fractions,
    component_molecular_weights,
    input_unit: str | None = None,
    output_unit: str | None = None,
) -> float:
    """Compatibility wrapper for the original list-based mass-fraction API."""
    # NOTE: Original API converted the final mixture MW when both units were set.
    result = calc_mixture_molecular_weight_from_mass_fractions(
        component_mass_fractions,
        component_molecular_weights,
    )
    if input_unit and output_unit:
        return _resolve_unit_conversion_fn(None)(result, input_unit, output_unit)
    return result


def calc_mixture_molecular_weight_from_mass_fractions_2(
    component_mass_fractions,
    component_molecular_weights,
    output_unit: str = "g/mol",
) -> float:
    """Compatibility wrapper for the original keyed mass-fraction API."""
    return calc_mixture_molecular_weight_from_mass_fractions(
        component_mass_fractions,
        component_molecular_weights,
        output_molecular_weight_unit=output_unit,
    )


# SECTION: Public exports
__all__ = [
    "calc_mixture_molecular_weight_from_mole_fractions",
    "calc_mixture_molecular_weight_from_mass_fractions",
    "calc_mixture_molecular_weight_1",
    "calc_mixture_molecular_weight_2",
    "calc_mixture_molecular_weight_from_mass_fractions_1",
    "calc_mixture_molecular_weight_from_mass_fractions_2",
]

