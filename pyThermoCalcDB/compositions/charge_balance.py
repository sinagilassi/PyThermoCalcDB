"""Charge-balance calculations."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional, cast

import numpy as np
from numpy.typing import NDArray

# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict
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
# ! ::: Numeric array-like core

def _calc_charge_balance(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate charge-balance residual: ``sum_i(z_i * c_i)``."""
    # SECTION: Convert inputs to numeric arrays
    concentration_values: NDArray[np.float64] = np.asarray(
        concentrations,
        dtype=np.float64,
    )
    charge_values: NDArray[np.float64] = np.asarray(
        charges,
        dtype=np.float64,
    )

    # SECTION: Validate inputs
    if concentration_values.ndim > 2:
        raise ValueError(
            "concentrations must be a scalar, 1-D array, or 2-D array."
        )

    if charge_values.ndim > 2:
        raise ValueError("charges must be a scalar, 1-D array, or 2-D array.")

    if not np.all(np.isfinite(concentration_values)):
        raise ValueError("concentrations must contain finite values.")

    if not np.all(np.isfinite(charge_values)):
        raise ValueError("charges must contain finite values.")

    # ! Concentrations cannot be negative for charge-balance composition input.
    if np.any(concentration_values < 0):
        raise ValueError("concentrations must be non-negative.")

    # ? NumPy handles scalar, 1-D, and 2-D broadcasting consistently here.
    try:
        concentration_values, charge_values = np.broadcast_arrays(
            concentration_values,
            charge_values,
        )
    except ValueError as exc:
        raise ValueError(
            "concentrations and charges must have broadcast-compatible shapes."
        ) from exc

    # SECTION: Calculate residual contributions
    values = concentration_values * charge_values

    # NOTE: Scalar and 1-D inputs return a scalar; 2-D inputs return row sums.
    if values.ndim == 0:
        return float(values)

    if values.ndim == 1:
        return float(np.sum(values))

    return cast(NDArray[np.float64], np.sum(values, axis=1))


# ! ::: Sequence adapter

def _calc_charge_balance_from_sequence(
    concentrations: Sequence[float | int],
    charges: Sequence[float | int],
) -> float:
    """Calculate charge-balance residual from numeric sequence inputs."""
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)

    # SECTION: Calculate charge-balance residual
    return cast(
        float,
        _calc_charge_balance(
            concentrations,
            charges,
        ),
    )


# ! ::: Mapping adapter

def _calc_charge_balance_from_mapping(
    concentrations: Mapping[str, float | int],
    charges: Mapping[str, float | int],
) -> float:
    """Calculate charge-balance residual from numeric mapping inputs."""
    # SECTION: Validate inputs
    non_negative(concentrations, "concentrations")
    same_shape(concentrations, charges)
    _validate_same_keys(concentrations, charges)

    # SECTION: Calculate charge-balance residual
    return cast(
        float,
        _calc_charge_balance(
            list(concentrations.values()),
            [charges[key] for key in concentrations],
        ),
    )


# ! ::: Unit-aware mapping adapter

def _calc_charge_balance_from_props(
    concentrations: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_concentration_unit: str,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate charge-balance residual from unit-aware concentration inputs."""
    # SECTION: Normalize unit-aware concentrations
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    normalized_concentrations = to_dict(
        concentrations,
        output_concentration_unit,
        unit_conversion_fn=conversion_fn,
    )

    # ? Component metadata can remap keys and enforce component order.
    normalized_concentrations = _configure_component_values(
        normalized_concentrations,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "concentrations",
    )
    normalized_charges = _configure_component_values(
        dict(charges),
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "charges",
    )

    # ! Mapping inputs must pair the same species.
    _validate_same_keys(normalized_concentrations, normalized_charges)

    # SECTION: Calculate charge-balance residual
    return cast(
        float,
        _calc_charge_balance(
            list(normalized_concentrations.values()),
            [normalized_charges[key] for key in normalized_concentrations],
        ),
    )


# =====================================================================
# *** Public annotated API
# =====================================================================
# ::: annotated for numpy array

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "array_like",
        "numpy",
        "numeric",
    ),
)
def calc_charge_balance(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float | NDArray[np.float64]]:
    """Return annotated charge-balance residual from numeric inputs."""
    # SECTION: Calculate value
    value = _calc_charge_balance(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_charge_balance",
    )


# ::: annotated for sequence

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from sequence inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "sequence",
        "numeric",
    ),
)
def calc_charge_balance_from_sequence(
    concentrations: Sequence[float | int],
    charges: Sequence[float | int],
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated charge-balance residual from numeric sequences."""
    # SECTION: Calculate value
    value = _calc_charge_balance_from_sequence(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_charge_balance_from_sequence",
    )


# ::: annotated for mapping

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from mapping inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Keyed species concentrations.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "mapping",
        "numeric",
    ),
)
def calc_charge_balance_from_mapping(
    concentrations: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated charge-balance residual from numeric mappings."""
    # SECTION: Calculate value
    value = _calc_charge_balance_from_mapping(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_charge_balance_from_mapping",
    )


# ::: annotated for CustomProp mapping

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from unit-aware mapping inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Keyed unit-aware species concentrations.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "mapping",
        "unit_aware",
    ),
)
def calc_charge_balance_from_props(
    concentrations: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_concentration_unit: str,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated charge-balance residual from unit-aware mappings."""
    # SECTION: Calculate value
    value = _calc_charge_balance_from_props(
        concentrations=concentrations,
        charges=charges,
        output_concentration_unit=output_concentration_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=_resolve_result_unit(
            "concentrations",
            concentrations,
            output_concentration_unit,
        ),
        symbol=symbol,
        implementation="_calc_charge_balance_from_props",
    )


# =====================================================================
# *** Electroneutrality checks
# =====================================================================
# ! ::: Numeric array-like check

def check_electroneutrality(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
    tolerance: float = 1e-12,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool | NDArray[np.bool_]]:
    """Check whether numeric charge-balance residuals are within tolerance."""
    # ! Tolerance must be a valid non-negative residual magnitude.
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate charge-balance residual
    result = _calc_charge_balance(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=np.abs(result) <= tolerance,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# ! ::: Sequence check

def check_electroneutrality_from_sequence(
    concentrations: Sequence[float | int],
    charges: Sequence[float | int],
    tolerance: float = 1e-12,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool]:
    """Check whether a numeric sequence composition is electrically neutral."""
    # ! Tolerance must be a valid non-negative residual magnitude.
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate charge-balance residual
    result = _calc_charge_balance_from_sequence(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=abs(result) <= tolerance,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# ! ::: Mapping check

def check_electroneutrality_from_mapping(
    concentrations: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    tolerance: float = 1e-12,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool]:
    """Check whether a numeric mapping composition is electrically neutral."""
    # ! Tolerance must be a valid non-negative residual magnitude.
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate charge-balance residual
    result = _calc_charge_balance_from_mapping(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=abs(result) <= tolerance,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# ! ::: Unit-aware mapping check

def check_electroneutrality_from_props(
    concentrations: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_concentration_unit: str,
    tolerance: float = 1e-12,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool]:
    """Check whether a unit-aware mapping composition is electrically neutral."""
    # ! Tolerance must be a valid non-negative residual magnitude.
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative.")

    # SECTION: Calculate charge-balance residual
    result = _calc_charge_balance_from_props(
        concentrations=concentrations,
        charges=charges,
        output_concentration_unit=output_concentration_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=abs(result) <= tolerance,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# =====================================================================
# *** Aliases
# =====================================================================

calc_mapping_charge_balance = calc_charge_balance_from_mapping
calc_sequence_charge_balance = calc_charge_balance_from_sequence


# SECTION: Public exports
__all__ = [
    # internal
    "_calc_charge_balance",
    "_calc_charge_balance_from_sequence",
    "_calc_charge_balance_from_mapping",
    "_calc_charge_balance_from_props",

    # public
    "calc_charge_balance",
    "calc_charge_balance_from_sequence",
    "calc_charge_balance_from_mapping",
    "calc_charge_balance_from_props",
    "check_electroneutrality",
    "check_electroneutrality_from_sequence",
    "check_electroneutrality_from_mapping",
    "check_electroneutrality_from_props",

    # compatibility
    "calc_mapping_charge_balance",
    "calc_sequence_charge_balance",
]
