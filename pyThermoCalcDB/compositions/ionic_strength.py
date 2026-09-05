"""Electrolyte composition primitives."""

# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Optional, List

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, Component, ComponentKey, AnnotatedValue
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.components import extract_components_values
from pythermodb_settings.utils.validators import non_negative, same_shape
# locals
from ..utils.conversions import _resolve_unit_conversion_fn, _configure_component_values, _validate_same_keys
from ..utils.tools import to_annotated_value, get_unit

# NOTE: setup logger
logger = logging.getLogger(__name__)

# SECTION: Internal helpers


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
    """Calculate ionic strength from a supplied concentration basis.

    Notes
    -----
    The ion charge ``z_i`` is dimensionless, so the calculated ionic strength
    has the same unit as the normalized concentration values. For molality input
    this is the normalized molality unit; for molarity input this is the
    normalized molarity unit. If no output concentration unit is supplied, the
    returned float follows the numeric basis of the supplied values.
    """
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if (
        isinstance(concentrations, Mapping) and
        isinstance(charges, Mapping)
    ):
        # ! Convert concentrations and charges to dictionaries with normalized units
        c = to_dict(
            concentrations,
            output_concentration_unit,
            unit_conversion_fn=conversion_fn,
        )
        z = to_dict(
            charges,
            unit_conversion_fn=conversion_fn
        )

        # ! Configure component values using metadata to dictionary
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

        # validate dicts
        _validate_same_keys(c, z)

        # calculate ionic strength
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

def _calc_ionic_strength_molality(
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
        Molality-based ionic strength on the normalized molality basis. Because
        charge is dimensionless, the return unit is the normalized molality
        unit, for example ``mol/kg`` or the unit selected by
        ``output_molality_unit``. If no unit is supplied, the returned float
        follows the numeric basis of ``molalities``.

    Notes
    -----
    Equation: ``I_m = 0.5*sum_i(m_i*z_i**2)``. The sum is over the species
    supplied by the caller; no dissociation/speciation assumptions are made.
    The squared charge term is dimensionless, so the result keeps the same unit
    as ``m_i`` after optional normalization.
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

# ! ::: return AnnotatedValue


def calc_ionic_strength_molality(
    molalities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molality_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength on the normalized molality basis.",
    symbol: str = "I_m",
) -> AnnotatedValue:
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
        name : str, optional
            Name of the property.
        description : str, optional
            Description of the property.
        symbol : str, optional
            Symbol representing the property.

        Returns
        -------
        AnnotatedValue
            Molality-based ionic strength on the normalized molality basis. Because
            charge is dimensionless, the return unit is the normalized molality
            unit, for example ``mol/kg`` or the unit selected by
            ``output_molality_unit``. If no unit is supplied, the return will have the same unit as the input molalities.

        Notes
        -----
        Equation:

            I_m = 0.5*sum_i(m_i*z_i**2)

        The sum is over the species
        supplied by the caller; no dissociation/speciation assumptions are made.
        The squared charge term is dimensionless, so the result keeps the same unit
        as ``m_i`` after optional normalization.
    """
    # NOTE: Calculate the ionic strength based on molality using the helper function.
    res = _calc_ionic_strength_molality(
        molalities=molalities,
        charges=charges,
        output_molality_unit=output_molality_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )

    # NOTE: set return unit
    get_unit_res = get_unit(
        identifier="molalities",
        data=molalities
    )
    # >> unit
    unit = str(
        get_unit_res["unit"]
    ) if get_unit_res["unit"] is not None else output_molality_unit

    # res contains the calculated ionic strength value
    return to_annotated_value(
        name=name,
        description=description,
        value=res,
        unit=unit,
        symbol=symbol,
    )


# ! ::: Molarity-based ionic strength using Component metadata


def calc_ionic_strength_molality_2(
    molalities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    components: List[Component],
    output_molality_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    name: str = "ionic_strength",
    description: str = "Molality-based ionic strength on the normalized molality basis.",
    symbol: str = "I_m",
) -> AnnotatedValue:
    """
    Calculate molality-based ionic strength using component metadata.

    Parameters
    ----------
    molalities : Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp]
        Molalities of the species.
    components : list[Component]
        Component metadata used for mapping-key remapping.
    output_molality_unit : str, optional
        Unit used to normalize molalities when supplied as ``CustomProp``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.
    component_key : ComponentKey, optional
        Component identifier format used for mapping inputs.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    name : str, optional
        Name of the returned annotated value. Defaults to ``"ionic_strength"``.
    description : str, optional
        Description of the returned annotated value. Defaults to
        ``"Molality-based ionic strength on the normalized molality basis."``.
    symbol : str, optional
        Symbol representing the ionic strength. Defaults to ``"I_m"``.

    Returns
    -------
    AnnotatedValue
        Molality-based ionic strength on the normalized molality basis. The
        ``unit`` field is ``output_molality_unit`` when supplied; otherwise it
        is the string ``"None"`` because no explicit unit normalization was
        requested.

    Notes
    -----
    Component charges are read from ``Component.net_charge`` and treated as
    dimensionless. The calculated value therefore has the same unit as the
    normalized molalities used in the calculation. No dissociation/speciation
    assumptions are made beyond the species represented by ``components``.
    """
    # SECTION: normalize molalities
    if isinstance(molalities, Sequence):
        # ? to list
        molalities_list = to_list(
            values=list(molalities),
            output_unit=output_molality_unit,
            unit_conversion_fn=unit_conversion_fn,
        )
    elif isinstance(molalities, Mapping):
        # ? mapping to dict
        molalities_dict = dict(molalities)

        # ! >>> Normalize molalities according to component metadata
        normalized_molalities = config_components_values(
            values=molalities_dict,
            components=components,
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=True,
        )
        # >> check
        if normalized_molalities is None:
            raise ValueError(
                "Failed to normalize molalities component values.")
        # >> unpack
        _, molalities_list = normalized_molalities

    # NOTE: build charge dictionary for components
    charges = extract_components_values(
        attribute_name='net_charge',
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
    )
    # >> check
    if charges is None:
        raise ValueError("Failed to extract component charges.")

    # >>> unpack
    _, charges_list = charges

    # SECTION: calculate ionic strength
    res = _ionic_strength(
        concentrations=molalities_list,
        charges=charges_list,
        output_concentration_unit=output_molality_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=True,
    )

    # res
    return to_annotated_value(
        name=name,
        description=description,
        value=res,
        unit=output_molality_unit if output_molality_unit is not None else 'None',
        symbol=symbol,
    )


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
    name: str = "ionic_strength",
    description: str = "Molarity-based ionic strength on the normalized molarity basis.",
    symbol: str = "I",
) -> AnnotatedValue:
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
    name : str, optional
        Name of the returned annotated value. Defaults to ``"ionic_strength"``.
    description : str, optional
        Description of the returned annotated value. Defaults to
        ``"Molarity-based ionic strength on the normalized molarity basis."``.
    symbol : str, optional
        Symbol representing the ionic strength. Defaults to ``"I"``.

    Returns
    -------
    AnnotatedValue
        Molarity-based ionic strength on the normalized molarity basis. Because
        charge is dimensionless, the return unit is the normalized molarity
        unit, for example ``mol/L`` or the unit selected by
        ``output_molarity_unit``. If no unit is supplied, the returned AnnotatedValue
        follows the numeric basis of ``molarities``.

    Notes
    -----
    Equation: ``I_c = 0.5*sum_i(c_i*z_i**2)``. Neutral species with ``z = 0``
    automatically contribute zero. The squared charge term is dimensionless, so
    the result keeps the same unit as ``c_i`` after optional normalization.
    """
    # NOTE: Neutral species with z = 0 automatically contribute zero.
    res = _ionic_strength(
        molarities,
        charges,
        output_molarity_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )

    # ! to AnnotatedValue
    return to_annotated_value(
        value=res,
        unit=output_molarity_unit,
        name=name,
        description=description,
        symbol=symbol,
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
    name: str = "charge_balance",
    description: str = "Charge-balance residual on the normalized concentration basis.",
    symbol: str = "ChBa_RES"
) -> AnnotatedValue:
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
    AnnotatedValue
        Charge-balance residual ``sum_i(z_i*c_i)``. A neutral bulk solution has
        a residual of zero within numerical tolerance. Because charge is
        dimensionless, the residual has the same unit as the normalized
        concentration values selected by ``output_concentration_unit``. If no
        unit is supplied, the returned AnnotatedValue follows the numeric basis of
        ``concentrations``.

    Notes
    -----
    This helper does not adjust concentrations to impose electroneutrality. It
    only sums signed charge concentration, so its unit is inherited from
    ``c_i`` after optional normalization.
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
        res_ = sum(z[key] * c[key] for key in c)
        # Ensure the result is a float before creating AnnotatedValue
        res_ = float(res_)

        # res
        return to_annotated_value(
            value=res_,
            unit=output_concentration_unit,
            name=name,
            description=description,
            symbol=symbol,
        )

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
    res_ = sum(z_i * c_i for c_i, z_i in zip(c, z))
    res_ = float(res_)

    # res
    return to_annotated_value(
        value=res_,
        unit=output_concentration_unit,
        name=name,
        description=description,
        symbol=symbol,
    )


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
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str = "ElNe",
) -> AnnotatedValue:
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
    name : str, optional
        Name of the returned AnnotatedValue.
    description : str, optional
        Description of the returned AnnotatedValue.
    symbol : str, optional
        Symbol representing the returned AnnotatedValue.

    Returns
    -------
    bool
        ``True`` when ``abs(sum_i(z_i*c_i)) <= tolerance``. The boolean return
        value has no unit.

    Notes
    -----
    This is only a diagnostic helper. It does not modify concentrations or
    solve a speciation/electroneutrality problem. The charge-balance residual
    compared against ``tolerance`` has the same unit as the normalized
    concentration values, so ``tolerance`` should be provided on that same
    concentration basis.
    """
    # SECTION: Validate tolerance
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate and compare charge-balance residual
    res_ = calc_charge_balance(
        concentrations,
        charges,
        output_concentration_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )
    res_ = abs(
        res_.value
    ) <= tolerance

    # return
    return to_annotated_value(
        value=res_,
        unit=None,
        name=name,
        description=description,
        symbol=symbol,
    )


# SECTION: Public exports
__all__ = [
    "calc_ionic_strength_molality",
    "calc_ionic_strength_molality_2",
    "calc_ionic_strength_molarity",
    "calc_charge_balance",
    "check_electroneutrality",
]
