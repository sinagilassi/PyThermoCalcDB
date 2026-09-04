"""Composition basis conversion functions.

These helpers perform deterministic arithmetic conversions between common
composition bases. Inputs may be plain numeric values or ``CustomProp`` values.
Mapping inputs are normalized through ``pythermodb_settings.utils.to_amounts``.

Unit convention
---------------
Arguments named ``output_*_unit`` are optional and define the unit used to
normalize the matching input before calculation. When an ``output_*_unit`` is
``None``, the input value is used as-is. When it is provided, numeric values are
assumed to already use that unit and ``CustomProp`` values are converted to that
unit with ``pycuc.convert_from_to``. The functions return plain floats,
dictionaries of floats, or lists of floats.
"""

from collections.abc import Mapping, Sequence
from typing import Optional, List
# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, ScalarValue, Component, ComponentKey
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.validators import (
    non_empty,
    fractions,
    positive,
    non_negative,
    same_shape
)
from pythermodb_settings.utils.quantity import (
    to_dict,
    to_list,
    to_scalar,
    pos
)
# locals
from ..utils.conversions import _resolve_unit_conversion_fn

# SECTION: Unit conversion helpers


def _dict(
    values: Mapping[str, float | int | CustomProp],
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> dict[str, float]:
    return to_dict(
        values,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


def _list(
    values: Sequence[float | int | CustomProp],
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> list[float]:
    return to_list(
        values,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


def _scalar(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
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
    return pos(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


_non_empty = non_empty
_fractions = fractions
_positive = positive
_non_negative = non_negative
_same_shape = same_shape


def _configure_component_values(
    values: Mapping[str, float],
    components: Optional[List[Component]],
    component_key: Optional[ComponentKey],
    case_sensitive: bool,
    sort_by_components_order: bool,
    name: str,
) -> dict[str, float]:
    """Remap and order mapping values using component metadata when requested."""
    # ! If no component key is provided, return the original values as a dictionary.
    if component_key is None:
        return dict(values)

    if not components:
        raise ValueError(
            f"component_key is provided but components is empty for {name}."
        )

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


# SECTION: Mole fraction and mass fraction conversions

# ! mole fraction to mass fraction conversion


def mole_fraction_to_mass_fraction(
    mole_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molecular_weights: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float] | list[float]:
    """Convert component mole fractions to mass fractions.

    Parameters
    ----------
    mole_fractions : Mapping[str, float | int | CustomProp] or Sequence[float | int | CustomProp]
        Component mole fractions. Values must be between zero and one. Mapping
        keys must match ``molecular_weights`` keys; sequence order must match
        ``molecular_weights`` order.
    molecular_weights : Mapping[str, float | int | CustomProp] or Sequence[float | int | CustomProp]
        Component molecular weights. Numeric values are assumed to already be in
        ``output_molecular_weight_unit``. ``CustomProp`` values are converted to
        ``output_molecular_weight_unit`` when it is provided.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weights`` before calculation. Leave as ``None`` to use input values as-is. Leave as
        ``None`` when all molecular weights are already numerically consistent.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.
    components : list[Component], optional
        Components used to remap and order mapping-input results when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format for remapping mapping-input keys. When
        ``None``, original mapping keys and insertion order are preserved.
    case_sensitive : bool, optional
        Whether component ID matching is case-sensitive. Defaults to ``True``.
    sort_by_components_order : bool, optional
        Whether mapping results should follow the order of ``components``.
        Defaults to ``True``.

    Returns
    -------
    dict[str, float] or list[float]
        Component mass fractions. Mapping input returns a dictionary; sequence
        input returns a list.

    Notes
    -----
    Equation defines as:
    w_i = x_i*M_i / sum_j(x_j*M_j)
    """
    # SECTION: Validate inputs
    fractions(mole_fractions, "mole_fractions")
    positive(molecular_weights, "molecular_weights")
    same_shape(mole_fractions, molecular_weights)

    # SECTION: Mapping implementation
    if isinstance(mole_fractions, Mapping) and isinstance(molecular_weights, Mapping):
        x = _dict(mole_fractions)
        mw = _dict(
            molecular_weights,
            output_molecular_weight_unit,
            unit_conversion_fn
        )
        x = _configure_component_values(
            x,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "mole_fractions",
        )
        mw = _configure_component_values(
            mw,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "molecular_weights",
        )
        denom = sum(x[key] * mw[key] for key in x)
        if denom <= 0.0:
            raise ValueError(
                "The weighted molecular-weight sum must be positive.")
        return {key: x[key] * mw[key] / denom for key in x}

    if isinstance(mole_fractions, Mapping) or isinstance(molecular_weights, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    x = _list(mole_fractions)
    mw = _list(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn
    )
    denom = sum(x_i * mw_i for x_i, mw_i in zip(x, mw))
    if denom <= 0.0:
        raise ValueError("The weighted molecular-weight sum must be positive.")
    return [x_i * mw_i / denom for x_i, mw_i in zip(x, mw)]

# ! ::: Convert Mass Fraction to Mole Fraction


def mass_fraction_to_mole_fraction(
    mass_fractions: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molecular_weights: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float] | list[float]:
    """Convert component mass fractions to mole fractions.

    Parameters
    ----------
    mass_fractions : Mapping[str, float | int | CustomProp] or Sequence[float | int | CustomProp]
        Component mass fractions. Values must be between zero and one. Mapping
        keys must match ``molecular_weights`` keys; sequence order must match
        ``molecular_weights`` order.
    molecular_weights : Mapping[str, float | int | CustomProp] or Sequence[float | int | CustomProp]
        Component molecular weights. Numeric values are assumed to already be in
        ``output_molecular_weight_unit``. ``CustomProp`` values are converted to
        ``output_molecular_weight_unit`` when it is provided.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weights`` before calculation. Leave as ``None`` to use input values as-is. Leave as
        ``None`` when all molecular weights are already numerically consistent.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.
    components : list[Component], optional
        Components used to remap and order mapping-input results when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format for remapping mapping-input keys. When
        ``None``, original mapping keys and insertion order are preserved.
    case_sensitive : bool, optional
        Whether component ID matching is case-sensitive. Defaults to ``True``.
    sort_by_components_order : bool, optional
        Whether mapping results should follow the order of ``components``.
        Defaults to ``True``.

    Returns
    -------
    dict[str, float] or list[float]
        Component mole fractions. Mapping input returns a dictionary; sequence
        input returns a list.

    Notes
    -----
    Equation defines as:
    x_i = (w_i/M_i) / sum_j(w_j/M_j)
    """
    # SECTION: Validate inputs
    _fractions(mass_fractions, "mass_fractions")
    _positive(molecular_weights, "molecular_weights")
    _same_shape(mass_fractions, molecular_weights)

    # SECTION: Mapping implementation
    if isinstance(mass_fractions, Mapping) and isinstance(molecular_weights, Mapping):
        w = _dict(mass_fractions)
        mw = _dict(
            molecular_weights,
            output_molecular_weight_unit,
            unit_conversion_fn
        )
        w = _configure_component_values(
            w,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "mass_fractions",
        )
        mw = _configure_component_values(
            mw,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "molecular_weights",
        )
        denom = sum(w[key] / mw[key] for key in w)
        if denom <= 0.0:
            raise ValueError(
                "The reciprocal molecular-weight sum must be positive.")
        return {key: (w[key] / mw[key]) / denom for key in w}

    if isinstance(mass_fractions, Mapping) or isinstance(molecular_weights, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    w = _list(mass_fractions)
    mw = _list(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn
    )
    denom = sum(w_i / mw_i for w_i, mw_i in zip(w, mw))
    if denom <= 0.0:
        raise ValueError(
            "The reciprocal molecular-weight sum must be positive.")
    return [(w_i / mw_i) / denom for w_i, mw_i in zip(w, mw)]


# SECTION: Molarity and molality conversions
# ! ::: Convert Molarity to Molality
def molarity_to_molality(
    molarity: ScalarValue,
    molecular_weight: ScalarValue,
    solution_density: ScalarValue,
    output_molarity_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    output_solution_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert single-solute molarity to molality.

    Parameters
    ----------
    molarity : float | int | CustomProp
        Solute molarity. Numeric values are assumed to already be in
        ``output_molarity_unit`` when it is provided.
    molecular_weight : float | int | CustomProp
        Solute molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    solution_density : float | int | CustomProp
        Solution density. Numeric values are assumed to already be in
        ``output_solution_density_unit`` when it is provided.
    output_molarity_unit : str, optional
        Unit used to normalize ``molarity`` before calculation. Leave as ``None`` to use the input value as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    output_solution_density_unit : str, optional
        Unit used to normalize ``solution_density`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Solute molality, consistent with the normalized input units.

    Notes
    -----
    Equation defines as:
    b = C / (rho - C*M)
    """
    # SECTION: Validate inputs
    c = _pos(
        molarity,
        "molarity",
        output_molarity_unit,
        unit_conversion_fn
    )
    mw = _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn
    )
    rho = _pos(
        solution_density,
        "solution_density",
        output_solution_density_unit,
        unit_conversion_fn
    )

    # ! The solvent mass on a 1-volume basis must remain positive.
    solvent_mass = rho - c * mw
    if solvent_mass <= 0.0:
        raise ValueError(
            "solution_density - molarity*molecular_weight must be positive.")

    # SECTION: Calculate molality
    return c / solvent_mass

# ! ::: Convert Molarity to Molality


def molality_to_molarity(
    molality: ScalarValue,
    molecular_weight: ScalarValue,
    solution_density: ScalarValue,
    output_molality_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    output_solution_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert single-solute molality to molarity.

    Parameters
    ----------
    molality : float | int | CustomProp
        Solute molality. Numeric values are assumed to already be in
        ``output_molality_unit`` when it is provided.
    molecular_weight : float | int | CustomProp
        Solute molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    solution_density : float | int | CustomProp
        Solution density. Numeric values are assumed to already be in
        ``output_solution_density_unit`` when it is provided.
    output_molality_unit : str, optional
        Unit used to normalize ``molality`` before calculation. Leave as ``None`` to use the input value as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    output_solution_density_unit : str, optional
        Unit used to normalize ``solution_density`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Solute molarity, consistent with the normalized input units.

    Notes
    -----
    Equation defines as:
    C = b*rho / (1 + b*M)
    """
    # SECTION: Validate inputs
    b = _pos(
        molality,
        "molality",
        output_molality_unit,
        unit_conversion_fn
    )
    mw = _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn
    )
    rho = _pos(
        solution_density,
        "solution_density",
        output_solution_density_unit,
        unit_conversion_fn
    )

    # SECTION: Calculate molarity
    return b * rho / (1.0 + b * mw)

# ! ::: Convert Molality to Molarity


def molarities_to_molalities(
    molarities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molecular_weights: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    solution_density: ScalarValue,
    output_molarity_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    output_solution_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float] | list[float]:
    """Convert multisolute molarities to molalities.

    Parameters
    ----------
    molarities : Mapping[str, float | int | CustomProp] or Sequence[float | int | CustomProp]
        Solute molarities. Numeric values are assumed to already be in
        ``output_molarity_unit`` when it is provided. Mapping keys must match ``molecular_weights``
        keys; sequence order must match ``molecular_weights`` order.
    molecular_weights : Mapping[str, float | int | CustomProp] or Sequence[float | int | CustomProp]
        Solute molecular weights. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    solution_density : float | int | CustomProp
        Solution density. Numeric values are assumed to already be in
        ``output_solution_density_unit`` when it is provided.
    output_molarity_unit : str, optional
        Unit used to normalize ``molarities`` before calculation. Leave as ``None`` to use input values as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weights`` before calculation. Leave as ``None`` to use input values as-is.
    output_solution_density_unit : str, optional
        Unit used to normalize ``solution_density`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.
    components : list[Component], optional
        Components used to remap and order mapping-input results when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format for remapping mapping-input keys. When
        ``None``, original mapping keys and insertion order are preserved.
    case_sensitive : bool, optional
        Whether component ID matching is case-sensitive. Defaults to ``True``.
    sort_by_components_order : bool, optional
        Whether mapping results should follow the order of ``components``.
        Defaults to ``True``.

    Returns
    -------
    dict[str, float] or list[float]
        Solute molalities. Mapping input returns a dictionary; sequence input
        returns a list.

    Notes
    -----
    Equation defines as:
    b_i = C_i / (rho - sum_j(C_j*M_j))
    """
    # SECTION: Validate inputs
    _non_empty(molarities, "molarities")
    _non_negative(molarities, "molarities")
    _positive(molecular_weights, "molecular_weights")
    _same_shape(molarities, molecular_weights)
    rho = _pos(
        solution_density,
        "solution_density",
        output_solution_density_unit,
        unit_conversion_fn
    )

    # SECTION: Mapping implementation
    if isinstance(molarities, Mapping) and isinstance(molecular_weights, Mapping):
        c = _dict(molarities, output_molarity_unit, unit_conversion_fn)
        mw = _dict(
            molecular_weights,
            output_molecular_weight_unit,
            unit_conversion_fn
        )
        c = _configure_component_values(
            c,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "molarities",
        )
        mw = _configure_component_values(
            mw,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "molecular_weights",
        )
        solvent_mass = rho - sum(c[key] * mw[key] for key in c)
        if solvent_mass <= 0.0:
            raise ValueError(
                "solution_density - sum(molarity*molecular_weight) must be positive.")
        return {key: c[key] / solvent_mass for key in c}

    if isinstance(molarities, Mapping) or isinstance(molecular_weights, Mapping):
        raise TypeError(
            "Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    c = _list(molarities, output_molarity_unit, unit_conversion_fn)
    mw = _list(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn
    )
    solvent_mass = rho - sum(c_i * mw_i for c_i, mw_i in zip(c, mw))
    if solvent_mass <= 0.0:
        raise ValueError(
            "solution_density - sum(molarity*molecular_weight) must be positive.")
    return [c_i / solvent_mass for c_i in c]


# SECTION: Molality and mole fraction conversions

# ! ::: Convert Molality to Mole Fraction
def molality_to_mole_fraction(
    molalities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    solvent_molecular_weight: ScalarValue,
    solvent_key: str = "solvent",
    output_molality_unit: str | None = None,
    output_solvent_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float] | list[float]:
    """Convert solute molalities to mole fractions, including solvent.

    Parameters
    ----------
    molalities : Mapping[str, float | int | CustomProp] or Sequence[float | int | CustomProp]
        Solute molalities. Numeric values are assumed to already be in
        ``output_molality_unit`` when it is provided.
    solvent_molecular_weight : float | int | CustomProp
        Solvent molecular weight. Numeric values are assumed to already be in
        ``output_solvent_molecular_weight_unit`` when it is provided.
    solvent_key : str, default="solvent"
        Key used for the solvent entry when ``molalities`` is a mapping.
    output_molality_unit : str, optional
        Unit used to normalize ``molalities`` before calculation. Leave as ``None`` to use input values as-is.
    output_solvent_molecular_weight_unit : str, optional
        Unit used to normalize ``solvent_molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.
    components : list[Component], optional
        Components used to remap and order mapping-input results when
        ``component_key`` is provided.
    component_key : ComponentKey, optional
        Component identifier format for remapping mapping-input keys. When
        ``None``, original mapping keys and insertion order are preserved.
    case_sensitive : bool, optional
        Whether component ID matching is case-sensitive. Defaults to ``True``.
    sort_by_components_order : bool, optional
        Whether mapping results should follow the order of ``components``.
        Defaults to ``True``.

    Returns
    -------
    dict[str, float] or list[float]
        Mole fractions for solutes plus solvent. Sequence output appends the
        solvent mole fraction as the final item.

    Notes
    -----
    Uses a 1 kg solvent basis, so n_solvent = 1/M_solvent.
    """
    # SECTION: Validate inputs
    _non_empty(molalities, "molalities")
    _non_negative(molalities, "molalities")
    solvent_moles = 1.0 / _pos(
        solvent_molecular_weight,
        "solvent_molecular_weight",
        output_solvent_molecular_weight_unit,
        unit_conversion_fn,
    )

    # SECTION: Mapping implementation
    if isinstance(molalities, Mapping):
        b = _dict(
            molalities,
            output_molality_unit,
            unit_conversion_fn
        )
        b = _configure_component_values(
            b,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "molalities",
        )
        total = solvent_moles + sum(b.values())
        res = {key: value / total for key, value in b.items()}
        res[solvent_key] = solvent_moles / total
        return res

    # SECTION: Sequence implementation
    b = _list(molalities, output_molality_unit, unit_conversion_fn)
    total = solvent_moles + sum(b)
    return [value / total for value in b] + [solvent_moles / total]

# ! ::: Convert Mole Fraction to Molality


def mole_fraction_to_molality(
    solute_mole_fraction: ScalarValue,
    solvent_mole_fraction: ScalarValue,
    solvent_molecular_weight: ScalarValue,
    output_solvent_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a solute mole fraction to molality.

    Parameters
    ----------
    solute_mole_fraction : float | int | CustomProp
        Mole fraction of the solute. Must be between zero and one.
    solvent_mole_fraction : float | int | CustomProp
        Mole fraction of the solvent. Must be greater than zero and no greater
        than one.
    solvent_molecular_weight : float | int | CustomProp
        Solvent molecular weight. Numeric values are assumed to already be in
        ``output_solvent_molecular_weight_unit`` when it is provided.
    output_solvent_molecular_weight_unit : str, optional
        Unit used to normalize ``solvent_molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Solute molality, consistent with the normalized solvent molecular weight.

    Notes
    -----
    Equation defines as:
    b_i = x_i / (x_s*M_s)
    """
    # SECTION: Validate inputs
    x_i = _scalar(solute_mole_fraction, "solute_mole_fraction")
    x_s = _scalar(solvent_mole_fraction, "solvent_mole_fraction")
    mw_s = _pos(
        solvent_molecular_weight,
        "solvent_molecular_weight",
        output_solvent_molecular_weight_unit,
        unit_conversion_fn,
    )

    # ! Solvent fraction must be positive because it appears in the denominator.
    if x_i < 0.0 or x_i > 1.0:
        raise ValueError("solute_mole_fraction must be between zero and one.")
    if x_s <= 0.0 or x_s > 1.0:
        raise ValueError(
            "solvent_mole_fraction must be greater than zero and no greater than one.")

    # SECTION: Calculate molality
    return x_i / (x_s * mw_s)


# SECTION: Molarity and mass fraction conversions
def molarity_to_mass_fraction(
    molarity: ScalarValue,
    molecular_weight: ScalarValue,
    solution_density: ScalarValue,
    output_molarity_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    output_solution_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert single-solute molarity to mass fraction.

    Parameters
    ----------
    molarity : float | int | CustomProp
        Solute molarity. Numeric values are assumed to already be in
        ``output_molarity_unit`` when it is provided.
    molecular_weight : float | int | CustomProp
        Solute molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    solution_density : float | int | CustomProp
        Solution density. Numeric values are assumed to already be in
        ``output_solution_density_unit`` when it is provided.
    output_molarity_unit : str, optional
        Unit used to normalize ``molarity`` before calculation. Leave as ``None`` to use the input value as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    output_solution_density_unit : str, optional
        Unit used to normalize ``solution_density`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Solute mass fraction.

    Notes
    -----
    Equation defines as:
    w_i = C_i*M_i / rho
    """
    # SECTION: Calculate and validate result
    result = _pos(molarity, "molarity", output_molarity_unit, unit_conversion_fn) * _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn,
    )
    result = result / _pos(
        solution_density,
        "solution_density",
        output_solution_density_unit,
        unit_conversion_fn
    )
    if result > 1.0:
        raise ValueError("Calculated mass fraction is greater than one.")
    return result

# ! ::: Convert Mass Fraction to Molarity


def mass_fraction_to_molarity(
    mass_fraction: ScalarValue,
    solution_density: ScalarValue,
    molecular_weight: ScalarValue,
    output_solution_density_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert single-solute mass fraction to molarity.

    Parameters
    ----------
    mass_fraction : float | int | CustomProp
        Solute mass fraction. Must be between zero and one.
    solution_density : float | int | CustomProp
        Solution density. Numeric values are assumed to already be in
        ``output_solution_density_unit`` when it is provided.
    molecular_weight : float | int | CustomProp
        Solute molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    output_solution_density_unit : str, optional
        Unit used to normalize ``solution_density`` before calculation. Leave as ``None`` to use the input value as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Solute molarity, consistent with the normalized input units.

    Notes
    -----
    Equation defines as:
    C_i = w_i*rho / M_i
    """
    # SECTION: Validate inputs
    w = _scalar(mass_fraction, "mass_fraction")
    if w < 0.0 or w > 1.0:
        raise ValueError("mass_fraction must be between zero and one.")

    # SECTION: Calculate molarity
    return w * _pos(
        solution_density,
        "solution_density",
        output_solution_density_unit,
        unit_conversion_fn
    ) / _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn,
    )


# SECTION: Molality and mass fraction conversions

# ! ::: Convert Molarity to Mass Fraction
def molality_to_mass_fraction(
    molality: ScalarValue,
    molecular_weight: ScalarValue,
    output_molality_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert single-solute molality to mass fraction.

    Parameters
    ----------
    molality : float | int | CustomProp
        Solute molality. Numeric values are assumed to already be in
        ``output_molality_unit`` when it is provided.
    molecular_weight : float | int | CustomProp
        Solute molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    output_molality_unit : str, optional
        Unit used to normalize ``molality`` before calculation. Leave as ``None`` to use the input value as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Solute mass fraction.

    Notes
    -----
    Equation defines as:
    w_i = b_i*M_i / (1 + b_i*M_i)
    """
    # SECTION: Calculate from a 1 kg solvent basis
    solute_mass = _pos(molality, "molality", output_molality_unit, unit_conversion_fn) * _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn,
    )
    return solute_mass / (1.0 + solute_mass)

# ! ::: Convert Mass Fraction to Molality


def mass_fraction_to_molality(
    mass_fraction: ScalarValue,
    molecular_weight: ScalarValue,
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert single-solute mass fraction to molality.

    Parameters
    ----------
    mass_fraction : float | int | CustomProp
        Solute mass fraction. Must be at least zero and less than one.
    molecular_weight : float | int | CustomProp
        Solute molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Solute molality, consistent with the normalized molecular weight.

    Notes
    -----
    Equation defines as:
    b_i = w_i / (M_i*(1 - w_i))
    """
    # SECTION: Validate inputs
    w = _scalar(mass_fraction, "mass_fraction")
    if w < 0.0 or w >= 1.0:
        raise ValueError(
            "mass_fraction must be at least zero and less than one.")

    # SECTION: Calculate molality
    return w / (
        _pos(molecular_weight, "molecular_weight",
             output_molecular_weight_unit, unit_conversion_fn) * (1.0 - w)
    )


# SECTION: Molarity and mass concentration conversions

# ! ::: Convert Molarity to Mass Concentration
def molarity_to_mass_concentration(
    molarity: ScalarValue,
    molecular_weight: ScalarValue,
    output_molarity_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert molarity to mass concentration.

    Parameters
    ----------
    molarity : float | int | CustomProp
        Component molarity. Numeric values are assumed to already be in
        ``output_molarity_unit`` when it is provided.
    molecular_weight : float | int | CustomProp
        Component molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    output_molarity_unit : str, optional
        Unit used to normalize ``molarity`` before calculation. Leave as ``None`` to use the input value as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Component mass concentration, consistent with the normalized input units.

    Notes
    -----
    Equation defines as:
    c_m,i = C_i*M_i
    """
    return _pos(molarity, "molarity", output_molarity_unit, unit_conversion_fn) * _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn,
    )


# ! ::: Convert Mass Concentration to Molarity
def mass_concentration_to_molarity(
    mass_concentration: ScalarValue,
    molecular_weight: ScalarValue,
    output_mass_concentration_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert mass concentration to molarity.

    Parameters
    ----------
    mass_concentration : float | int | CustomProp
        Component mass concentration. Numeric values are assumed to already be in
        ``output_mass_concentration_unit`` when it is provided.
    molecular_weight : float | int | CustomProp
        Component molecular weight. Numeric values are assumed to already be in
        ``output_molecular_weight_unit`` when it is provided.
    output_mass_concentration_unit : str, optional
        Unit used to normalize ``mass_concentration`` before calculation. Leave as ``None`` to use the input value as-is.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation. Leave as ``None`` to use the input value as-is.
    unit_conversion_fn : UnitConversionFn, optional
        Function used to convert ``CustomProp`` values when a matching
        ``output_*_unit`` is provided. Defaults to ``pycuc.convert_from_to``
        when ``None``.

    Returns
    -------
    float
        Component molarity, consistent with the normalized input units.

    Notes
    -----
    Equation defines as:
    C_i = c_m,i / M_i
    """
    return _pos(mass_concentration, "mass_concentration", output_mass_concentration_unit, unit_conversion_fn) / _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn,
    )


# SECTION: Percent conversions
# ! ::: Convert Mass Fraction to Weight Percent
def mass_fraction_to_weight_percent(mass_fraction: ScalarValue) -> float:
    """Convert mass fraction to weight percent.

    Parameters
    ----------
    mass_fraction : float | int | CustomProp
        Mass fraction on [0, 1]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for fraction inputs.

    Returns
    -------
    float
        Weight percent on [0, 100].
    """
    # SECTION: Validate and convert
    w = _scalar(mass_fraction, "mass_fraction")
    if w < 0.0 or w > 1.0:
        raise ValueError("mass_fraction must be between zero and one.")
    return 100.0 * w

# ! ::: Convert Weight Percent to Mass Fraction


def weight_percent_to_mass_fraction(weight_percent: ScalarValue) -> float:
    """Convert weight percent to mass fraction.

    Parameters
    ----------
    weight_percent : float | int | CustomProp
        Weight percent on [0, 100]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for percent inputs.

    Returns
    -------
    float
        Mass fraction on [0, 1].
    """
    # SECTION: Validate and convert
    value = _scalar(weight_percent, "weight_percent")
    if value < 0.0 or value > 100.0:
        raise ValueError("weight_percent must be between zero and 100.")
    return value / 100.0

# ! ::: Convert Mole Fraction to Mole Percent


def mole_fraction_to_mole_percent(mole_fraction: ScalarValue) -> float:
    """Convert mole fraction to mole percent.

    Parameters
    ----------
    mole_fraction : float | int | CustomProp
        Mole fraction on [0, 1]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for fraction inputs.

    Returns
    -------
    float
        Mole percent on [0, 100].
    """
    # SECTION: Validate and convert
    x = _scalar(mole_fraction, "mole_fraction")
    if x < 0.0 or x > 1.0:
        raise ValueError("mole_fraction must be between zero and one.")
    return 100.0 * x

# ! ::: Convert Mole Percent to Mole Fraction


def mole_percent_to_mole_fraction(mole_percent: ScalarValue) -> float:
    """Convert mole percent to mole fraction.

    Parameters
    ----------
    mole_percent : float | int | CustomProp
        Mole percent on [0, 100]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for percent inputs.

    Returns
    -------
    float
        Mole fraction on [0, 1].
    """
    # SECTION: Validate and convert
    value = _scalar(mole_percent, "mole_percent")
    if value < 0.0 or value > 100.0:
        raise ValueError("mole_percent must be between zero and 100.")
    return value / 100.0


# SECTION: ppm conversions
# ! ::: Convert Mass Fraction to Mass-based ppm
def mass_fraction_to_ppm(mass_fraction: ScalarValue) -> float:
    """Convert mass fraction to mass-based ppm.

    Parameters
    ----------
    mass_fraction : float | int | CustomProp
        Mass fraction on [0, 1]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for fraction inputs.

    Returns
    -------
    float
        Mass-based ppm.
    """
    # NOTE: ppm_mass = mass_fraction * 1e6
    return mass_fraction_to_weight_percent(mass_fraction) * 10000.0

# ! ::: Convert Mass-based ppm to Mass Fraction


def ppm_mass_to_mass_fraction(ppm: ScalarValue) -> float:
    """Convert mass-based ppm to mass fraction.

    Parameters
    ----------
    ppm : float | int | CustomProp
        Mass-based parts per million. Must be non-negative. ``CustomProp`` values
        are used as-is; no unit conversion is performed.

    Returns
    -------
    float
        Mass fraction.
    """
    # SECTION: Validate and convert
    value = _scalar(ppm, "ppm")
    if value < 0.0:
        raise ValueError("ppm must be non-negative.")
    return value * 1e-6


# ! ::: Convert Mole Fraction to Mole-based ppm

def mole_fraction_to_ppm(mole_fraction: ScalarValue) -> float:
    """Convert mole fraction to mole-based ppm.

    Parameters
    ----------
    mole_fraction : float | int | CustomProp
        Mole fraction on [0, 1]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for fraction inputs.

    Returns
    -------
    float
        Mole-based ppm.
    """
    # NOTE: ppm_mole = mole_fraction * 1e6
    return mole_fraction_to_mole_percent(mole_fraction) * 10000.0

# ! ::: Convert Mole-based ppm to Mole Fraction


def ppm_mole_to_mole_fraction(ppm: ScalarValue) -> float:
    """Convert mole-based ppm to mole fraction.

    Parameters
    ----------
    ppm : float | int | CustomProp
        Mole-based parts per million. Must be non-negative. ``CustomProp`` values
        are used as-is; no unit conversion is performed.

    Returns
    -------
    float
        Mole fraction.
    """
    # SECTION: Validate and convert
    value = _scalar(ppm, "ppm")
    if value < 0.0:
        raise ValueError("ppm must be non-negative.")
    return value * 1e-6


# SECTION: ppb conversions

# ! ::: Convert Mass Fraction to Mass-based ppb
def mass_fraction_to_ppb(mass_fraction: ScalarValue) -> float:
    """Convert mass fraction to mass-based ppb.

    Parameters
    ----------
    mass_fraction : float | int | CustomProp
        Mass fraction on [0, 1]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for fraction inputs.

    Returns
    -------
    float
        Mass-based ppb.
    """
    # NOTE: ppb_mass = mass_fraction * 1e9
    return mass_fraction_to_weight_percent(mass_fraction) * 10000000.0

# ! ::: Convert Mass-based ppb to Mass Fraction


def ppb_mass_to_mass_fraction(ppb: ScalarValue) -> float:
    """Convert mass-based ppb to mass fraction.

    Parameters
    ----------
    ppb : float | int | CustomProp
        Mass-based parts per billion. Must be non-negative. ``CustomProp`` values
        are used as-is; no unit conversion is performed.

    Returns
    -------
    float
        Mass fraction.
    """
    # SECTION: Validate and convert
    value = _scalar(ppb, "ppb")
    if value < 0.0:
        raise ValueError("ppb must be non-negative.")
    return value * 1e-9

# ! ::: Convert Mole Fraction to Mole-based ppb


def mole_fraction_to_ppb(mole_fraction: ScalarValue) -> float:
    """Convert mole fraction to mole-based ppb.

    Parameters
    ----------
    mole_fraction : float | int | CustomProp
        Mole fraction on [0, 1]. ``CustomProp`` values are used as-is; no unit
        conversion is performed for fraction inputs.

    Returns
    -------
    float
        Mole-based ppb.
    """
    # NOTE: ppb_mole = mole_fraction * 1e9
    return mole_fraction_to_mole_percent(mole_fraction) * 10000000.0

# ! ::: Convert Mole-based ppb to Mole Fraction


def ppb_mole_to_mole_fraction(ppb: ScalarValue) -> float:
    """Convert mole-based ppb to mole fraction.

    Parameters
    ----------
    ppb : float | int | CustomProp
        Mole-based parts per billion. Must be non-negative. ``CustomProp`` values
        are used as-is; no unit conversion is performed.

    Returns
    -------
    float
        Mole fraction.
    """
    # SECTION: Validate and convert
    value = _scalar(ppb, "ppb")
    if value < 0.0:
        raise ValueError("ppb must be non-negative.")
    return value * 1e-9


# SECTION: Public exports
__all__ = [
    "mole_fraction_to_mass_fraction",
    "mass_fraction_to_mole_fraction",
    "molarity_to_molality",
    "molality_to_molarity",
    "molarities_to_molalities",
    "molality_to_mole_fraction",
    "mole_fraction_to_molality",
    "molarity_to_mass_fraction",
    "mass_fraction_to_molarity",
    "molality_to_mass_fraction",
    "mass_fraction_to_molality",
    "molarity_to_mass_concentration",
    "mass_concentration_to_molarity",
    "mass_fraction_to_weight_percent",
    "weight_percent_to_mass_fraction",
    "mole_fraction_to_mole_percent",
    "mole_percent_to_mole_fraction",
    "mass_fraction_to_ppm",
    "ppm_mass_to_mass_fraction",
    "mole_fraction_to_ppm",
    "ppm_mole_to_mole_fraction",
    "mass_fraction_to_ppb",
    "ppb_mass_to_mass_fraction",
    "mole_fraction_to_ppb",
    "ppb_mole_to_mole_fraction",
]
