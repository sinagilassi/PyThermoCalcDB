"""Charge-balance calculations."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional

# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import non_negative, same_shape

# locals
from ..utils.conversions import (
    _configure_component_values,
    _resolve_result_unit,
    _resolve_unit_conversion_fn,
    _validate_same_keys,
)
from ..utils.tools import to_annotated_value


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

def _calc_charge_balance_v1(
    concentrations: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate charge-balance residual ``sum_i(z_i * c_i)``."""
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)

    # ? Normalize unit-aware concentration values only when needed.
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Mapping implementation
    if isinstance(concentrations, Mapping) and isinstance(charges, Mapping):
        concentration_values = to_dict(
            concentrations,
            output_concentration_unit,
            unit_conversion_fn=conversion_fn,
        )
        charge_values = to_dict(
            charges,
            unit_conversion_fn=conversion_fn,
        )

        # ? Component metadata can remap keys and enforce component order.
        concentration_values = _configure_component_values(
            concentration_values,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "concentrations",
        )
        charge_values = _configure_component_values(
            charge_values,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
            "charges",
        )

        # ! Mapping inputs must pair the same species.
        _validate_same_keys(concentration_values, charge_values)

        return float(
            sum(
                charge_values[key] * concentration_values[key]
                for key in concentration_values
            )
        )

    # ! Mixed mapping/sequence inputs cannot be paired unambiguously.
    if isinstance(concentrations, Mapping) or isinstance(charges, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")

    # SECTION: Sequence implementation
    concentration_values = to_list(
        concentrations,
        output_concentration_unit,
        unit_conversion_fn=conversion_fn,
    )
    charge_values = to_list(
        charges,
        unit_conversion_fn=conversion_fn,
    )

    if len(concentration_values) != len(charge_values):
        raise ValueError("concentrations and charges must have the same length.")

    return float(
        sum(
            z_i * c_i
            for c_i, z_i in zip(concentration_values, charge_values)
        )
    )


# =====================================================================
# *** Public annotated API
# =====================================================================

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from mapping inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Keyed species concentrations.",
        "charges": "Keyed species charges.",
    },
    outputs={"charge_balance": "Charge-balance residual."},
    tags=("charge_balance", "mapping"),
)
def calc_charge_balance_from_mapping(
    concentrations: Mapping[str, float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp],
    output_concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    symbol: str = "CB",
) -> AnnotatedValue[float]:
    """Return annotated charge-balance residual from mapping inputs."""
    # SECTION: Calculate value
    value = _calc_charge_balance_v1(
        concentrations,
        charges,
        output_concentration_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        unit=_resolve_result_unit(
            "concentrations",
            concentrations,
            output_concentration_unit,
        ),
        name=name,
        description=description,
        symbol=symbol,
    )


@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from sequence inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
    },
    outputs={"charge_balance": "Charge-balance residual."},
    tags=("charge_balance", "sequence"),
)
def calc_charge_balance_from_sequence(
    concentrations: Sequence[float | int | CustomProp],
    charges: Sequence[float | int | CustomProp],
    output_concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    symbol: str = "CB",
) -> AnnotatedValue[float]:
    """Return annotated charge-balance residual from sequence inputs."""
    # SECTION: Calculate value
    value = _calc_charge_balance_v1(
        concentrations,
        charges,
        output_concentration_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        unit=_resolve_result_unit(
            "concentrations",
            concentrations,
            output_concentration_unit,
        ),
        name=name,
        description=description,
        symbol=symbol,
    )


@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
    },
    outputs={"charge_balance": "Charge-balance residual."},
    tags=("charge_balance",),
)
def calc_charge_balance(
    concentrations: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    symbol: str = "CB",
) -> AnnotatedValue[float]:
    """Return annotated charge-balance residual."""
    # SECTION: Calculate value
    value = _calc_charge_balance_v1(
        concentrations,
        charges,
        output_concentration_unit,
        unit_conversion_fn,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        unit=_resolve_result_unit(
            "concentrations",
            concentrations,
            output_concentration_unit,
        ),
        name=name,
        description=description,
        symbol=symbol,
    )


@calculation_info(
    name="electroneutrality",
    description="Check whether a composition is electrically neutral.",
    equation="abs(sum_i(z_i * c_i)) <= tolerance",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
        "tolerance": "Absolute residual tolerance.",
    },
    outputs={"electroneutrality": "Whether the composition is electrically neutral."},
    tags=("charge_balance", "electroneutrality"),
)
def check_electroneutrality(
    concentrations: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    charges: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    tolerance: float = 1e-12,
    output_concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str = "ElNe",
) -> AnnotatedValue[bool]:
    """Check whether ``abs(charge_balance) <= tolerance``."""
    # ! Tolerance must be a valid non-negative residual magnitude.
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate charge-balance residual
    result = calc_charge_balance(
        concentrations,
        charges,
        output_concentration_unit,
        unit_conversion_fn,
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=abs(result.value) <= tolerance,
        unit=None,
        name=name,
        description=description,
        symbol=symbol,
    )


# =====================================================================
# *** Aliases
# =====================================================================

calc_mapping_charge_balance = calc_charge_balance_from_mapping
calc_sequence_charge_balance = calc_charge_balance_from_sequence


# SECTION: Public exports
__all__ = [
    "_calc_charge_balance_v1",
    "calc_charge_balance_from_mapping",
    "calc_charge_balance_from_sequence",
    "calc_mapping_charge_balance",
    "calc_sequence_charge_balance",
    "calc_charge_balance",
    "check_electroneutrality",
]
