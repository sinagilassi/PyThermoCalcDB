"""Electrolyte composition primitives."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, List

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import non_negative, same_shape
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
            f"component_key is provided but components is empty for {name}.")

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
    component_values_dict, _ = component_values
    return component_values_dict


def _validate_same_keys(left: Mapping[str, float], right: Mapping[str, float]) -> None:
    """Validate matching component keys."""
    # ? Mismatched keys usually indicate a missing concentration or charge.
    if set(left) != set(right):
        raise ValueError(
            "concentrations and charges must have the same component keys.")


# ! ::: Ionic strength

def _ionic_strength(
    concentrations: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_concentration_unit: str | None,
    unit_conversion_fn: UnitConversionFn | None,
    components: Optional[List[Component]],
    component_key: Optional[ComponentKey],
    case_sensitive: bool,
    sort_by_components_order: bool,
) -> float:
    """Calculate ionic strength from a supplied concentration basis."""
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(concentrations, Mapping) and isinstance(charges, Mapping):
        c = to_dict(
            concentrations,
            output_concentration_unit,
            unit_conversion_fn=conversion_fn,
        )
        z = to_dict(
            charges,
            unit_conversion_fn=conversion_fn
        )
        c = _configure_component_values(
            c,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "concentrations"
        )
        z = _configure_component_values(
            z,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "charges"
        )
        _validate_same_keys(c, z)
        return 0.5 * sum(c[key] * z[key] ** 2 for key in c)

    # ! Mixed mapping/sequence input is ambiguous.
    if isinstance(concentrations, Mapping) or isinstance(charges, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    c = to_list(
        concentrations,
        output_concentration_unit,
        unit_conversion_fn=conversion_fn,
    )
    z = to_list(charges, unit_conversion_fn=conversion_fn)
    if len(c) != len(z):
        raise ValueError(
            "concentrations and charges must have the same length.")
    return 0.5 * sum(c_i * z_i ** 2 for c_i, z_i in zip(c, z))


# ! ::: Calculate molality-based ionic strength

def calc_ionic_strength_molality(
    molalities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molality_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate molality-based ionic strength.

    Parameters
    ----------
    molalities : mapping or sequence of float | int | CustomProp
        Species molalities. Values must be non-negative.
    charges : mapping or sequence of float | int | CustomProp
        Species charges. Neutral species may be supplied with charge zero.
    output_molality_unit : str, optional
        Unit used to normalize molalities when supplied as ``CustomProp``.
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
        Molality-based ionic strength on the normalized molality basis.

    Notes
    -----
    Equation: ``I_m = 0.5*sum_i(m_i*z_i**2)``. The sum is over the species
    supplied by the caller; no dissociation/speciation assumptions are made.
    """
    # NOTE: No speciation or dissociation assumptions are made here.
    return _ionic_strength(
        molalities,
        charges,
        output_molality_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )

# ! ::: Molarity-based ionic strength using Component metadata


def calc_ionic_strength_molality_2(
    molalities: Mapping[str, float | int | CustomProp],
    components: List[Component],
    output_molality_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
):
    # SECTION: build a dictionary of charge for components
    # molarities_normalized = config_components_values(
    #     values=molalities,
    #     components=components,
    #     component_key=component_key,
    #     case_sensitive=case_sensitive,
    #     sort_by_components_order=sort_by_components_order,
    # )
    pass


# ! ::: Calculate molarity-based ionic strength


def calc_ionic_strength_molarity(
    molarities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molarity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate molarity-based ionic strength.

    Parameters
    ----------
    molarities : mapping or sequence of float | int | CustomProp
        Species molarities. Values must be non-negative.
    charges : mapping or sequence of float | int | CustomProp
        Species charges. Neutral species may be supplied with charge zero.
    output_molarity_unit : str, optional
        Unit used to normalize molarities when supplied as ``CustomProp``.
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
        Molarity-based ionic strength on the normalized molarity basis.

    Notes
    -----
    Equation: ``I_c = 0.5*sum_i(c_i*z_i**2)``. Neutral species with ``z = 0``
    automatically contribute zero.
    """
    # NOTE: Neutral species with z = 0 automatically contribute zero.
    return _ionic_strength(
        molarities,
        charges,
        output_molarity_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )


# SECTION: Charge balance

def calc_charge_balance(
    concentrations: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate charge-balance residual.

    Parameters
    ----------
    concentrations : mapping or sequence of float | int | CustomProp
        Species concentrations on a common basis. Values must be non-negative.
    charges : mapping or sequence of float | int | CustomProp
        Species charges.
    output_concentration_unit : str, optional
        Unit used to normalize concentrations when supplied as ``CustomProp``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.
    components : list[Component], optional
        Component metadata used for mapping-key remapping.
    component_key : ComponentKey, optional
        Component identifier format used for mapping inputs.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether mapping values should follow ``components`` order.

    Returns
    -------
    float
        Charge-balance residual ``sum_i(z_i*c_i)``. A neutral bulk solution has
        a residual of zero within numerical tolerance.

    Notes
    -----
    This helper does not adjust concentrations to impose electroneutrality.
    """
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(concentrations, Mapping) and isinstance(charges, Mapping):
        c = to_dict(concentrations, output_concentration_unit,
                    unit_conversion_fn=conversion_fn)
        z = to_dict(charges, unit_conversion_fn=conversion_fn)
        c = _configure_component_values(
            c, components, component_key, case_sensitive, sort_by_components_order, "concentrations"
        )
        z = _configure_component_values(
            z, components, component_key, case_sensitive, sort_by_components_order, "charges"
        )
        _validate_same_keys(c, z)
        return sum(z[key] * c[key] for key in c)

    # ! Mixed mapping/sequence input is ambiguous.
    if isinstance(concentrations, Mapping) or isinstance(charges, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    c = to_list(concentrations, output_concentration_unit,
                unit_conversion_fn=conversion_fn)
    z = to_list(charges, unit_conversion_fn=conversion_fn)
    if len(c) != len(z):
        raise ValueError(
            "concentrations and charges must have the same length.")
    return sum(z_i * c_i for c_i, z_i in zip(c, z))


# SECTION: Electroneutrality check

def check_electroneutrality(
    concentrations: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    tolerance: float = 1e-12,
    output_concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> bool:
    """Check whether a composition is electrically neutral.

    Parameters
    ----------
    concentrations : mapping or sequence of float | int | CustomProp
        Species concentrations on a common basis.
    charges : mapping or sequence of float | int | CustomProp
        Species charges.
    tolerance : float, optional
        Absolute tolerance for the charge-balance residual. Must be
        non-negative.
    output_concentration_unit : str, optional
        Unit used to normalize concentrations when supplied as ``CustomProp``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.
    components : list[Component], optional
        Component metadata used for mapping-key remapping.
    component_key : ComponentKey, optional
        Component identifier format used for mapping inputs.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether mapping values should follow ``components`` order.

    Returns
    -------
    bool
        ``True`` when ``abs(sum_i(z_i*c_i)) <= tolerance``.

    Notes
    -----
    This is only a diagnostic helper. It does not modify concentrations or
    solve a speciation/electroneutrality problem.
    """
    # SECTION: Validate tolerance
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate and compare charge-balance residual
    return abs(
        calc_charge_balance(
            concentrations,
            charges,
            output_concentration_unit,
            unit_conversion_fn,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
        )
    ) <= tolerance


# SECTION: Public exports
__all__ = [
    "calc_ionic_strength_molality",
    "calc_ionic_strength_molarity",
    "calc_charge_balance",
    "check_electroneutrality",
]
