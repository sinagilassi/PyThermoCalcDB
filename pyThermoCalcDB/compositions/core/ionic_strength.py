"""Shared numeric ionic-strength core."""

# import libs
from collections.abc import Mapping, Sequence
from typing import cast, Optional
from pythermodb_settings.models import Component, ComponentKey, CustomProp, UnitConversionFn
from pythermodb_settings.utils import config_components_values
from pythermodb_settings.utils.components import extract_components_values
from pythermodb_settings.utils.quantity import to_dict
from pythermodb_settings.utils.validators import non_negative, same_shape

import numpy as np
from numpy.typing import NDArray
from ...utils.conversions import (
    _resolve_result_unit,
    _resolve_unit_conversion_fn,
    _validate_same_keys,
    _configure_component_values
)


# =====================================================================
# *** Helper functions
# =====================================================================

def _validate_concentrations_and_charges(
    concentrations: NDArray[np.float64],
    charges: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Validate and broadcast concentration and charge arrays."""
    # SECTION: Validate array dimensions
    if concentrations.ndim > 2:
        raise ValueError(
            "concentrations must be a scalar, 1-D array, or 2-D array."
        )

    if charges.ndim > 2:
        raise ValueError("charges must be a scalar, 1-D array, or 2-D array.")

    # SECTION: Validate finite numeric values
    if not np.all(np.isfinite(concentrations)):
        raise ValueError("concentrations must contain finite values.")

    if not np.all(np.isfinite(charges)):
        raise ValueError("charges must contain finite values.")

    # ! Concentrations cannot be negative for ionic strength.
    if np.any(concentrations < 0):
        raise ValueError("concentrations must be non-negative.")

    # ? NumPy handles scalar, 1-D, and 2-D broadcasting consistently here.
    try:
        concentrations, charges = np.broadcast_arrays(concentrations, charges)
    except ValueError as exc:
        raise ValueError(
            "concentrations and charges must have broadcast-compatible shapes."
        ) from exc

    return (
        cast(NDArray[np.float64], concentrations),
        cast(NDArray[np.float64], charges),
    )


# ======================================================================
# *** Core calculation
# ======================================================================

def _calc_ionic_strength(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate ionic strength: ``I = 0.5 * sum_i(c_i * z_i**2)``."""
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
    concentration_values, charge_values = _validate_concentrations_and_charges(
        concentration_values,
        charge_values,
    )

    # SECTION: Calculate component ionic-strength contributions
    values = 0.5 * concentration_values * charge_values**2

    # NOTE: Scalar and 1-D inputs return a scalar; 2-D inputs return row sums.
    if values.ndim == 0:
        return float(values)

    if values.ndim == 1:
        return float(np.sum(values))

    return cast(NDArray[np.float64], np.sum(values, axis=1))


# ======================================================================
# *** Internal deterministic calculations (molality)
# ======================================================================
# ! ::: Numeric array-like core

def _calc_ionic_strength_molality(
    molalities: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate molality-based ionic strength from numeric array-like inputs."""
    # SECTION: Calculate ionic strength
    return _calc_ionic_strength(
        concentrations=molalities,
        charges=charges,
    )


# ! ::: Sequence adapter

def _calc_ionic_strength_molality_from_sequence(
    molalities: Sequence[float | int],
    charges: Sequence[float | int],
) -> float:
    """Calculate molality-based ionic strength from numeric sequence inputs."""
    # SECTION: Validate inputs
    non_negative(molalities, "molalities")
    same_shape(molalities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            molalities,
            charges,
        ),
    )


# ! ::: Mapping adapter

def _calc_ionic_strength_molality_from_mapping(
    molalities: Mapping[str, float | int],
    charges: Mapping[str, float | int],
) -> float:
    """Calculate molality-based ionic strength from numeric mapping inputs."""
    # SECTION: Validate inputs
    non_negative(molalities, "molalities")
    same_shape(molalities, charges)
    _validate_same_keys(molalities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            list(molalities.values()),
            [charges[key] for key in molalities],
        ),
    )


# ! ::: Unit-aware mapping adapter

def _calc_ionic_strength_molality_from_props(
    molalities: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_molality_unit: str = "mol/kg",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate molality-based ionic strength from unit-aware molality inputs."""
    # SECTION: Normalize unit-aware molalities
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    normalized_molalities = to_dict(
        molalities,
        output_molality_unit,
        unit_conversion_fn=conversion_fn,
    )

    # ! Mapping inputs must pair the same species.
    _validate_same_keys(normalized_molalities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            list(normalized_molalities.values()),
            [charges[key] for key in normalized_molalities],
        ),
    )


# ! ::: Mapping input with component metadata

def _calc_ionic_strength_molality_with_components_from_mapping(
    molalities: Mapping[str, float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
) -> float:
    """Calculate molality ionic strength from mapping input and components."""
    # SECTION: Normalize molality values by component metadata
    configured = config_components_values(
        values=dict(molalities),
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=True,
    )

    # ! Component configuration must succeed before charges are extracted.
    if configured is None:
        raise ValueError("Failed to normalize molalities component values.")

    _, molality_values = configured

    # SECTION: Extract component charges
    extracted = extract_components_values(
        attribute_name="net_charge",
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    if extracted is None:
        raise ValueError("Failed to extract component charges.")

    _, charge_values = extracted

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            molality_values,
            charge_values,
        ),
    )


# ! ::: Sequence input with component metadata

def _calc_ionic_strength_molality_with_components_from_sequence(
    molalities: Sequence[float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
) -> float:
    """Calculate molality ionic strength from sequence input and components."""
    # SECTION: Extract component charges
    extracted = extract_components_values(
        attribute_name="net_charge",
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    if extracted is None:
        raise ValueError("Failed to extract component charges.")

    _, charge_values = extracted

    # SECTION: Calculate ionic strength
    return _calc_ionic_strength_molality_from_sequence(
        molalities=molalities,
        charges=charge_values,
    )

# ======================================================================
# *** Internal deterministic calculations
# ======================================================================
# ! ::: Numeric array-like core


def _calc_ionic_strength_molarity(
    molarities: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
) -> float | NDArray[np.float64]:
    """Calculate molarity-based ionic strength from numeric array-like inputs."""
    # SECTION: Calculate ionic strength
    return _calc_ionic_strength(
        concentrations=molarities,
        charges=charges,
    )


# ! ::: Sequence adapter

def _calc_ionic_strength_molarity_from_sequence(
    molarities: Sequence[float | int],
    charges: Sequence[float | int],
) -> float:
    """Calculate molarity-based ionic strength from numeric sequence inputs."""
    # SECTION: Validate inputs
    non_negative(molarities, "molarities")
    same_shape(molarities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            molarities,
            charges,
        ),
    )


# ! ::: Mapping adapter

def _calc_ionic_strength_molarity_from_mapping(
    molarities: Mapping[str, float | int],
    charges: Mapping[str, float | int],
) -> float:
    """Calculate molarity-based ionic strength from numeric mapping inputs."""
    # SECTION: Validate inputs
    non_negative(molarities, "molarities")
    same_shape(molarities, charges)
    _validate_same_keys(molarities, charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            list(molarities.values()),
            [charges[key] for key in molarities],
        ),
    )


# ! ::: Unit-aware mapping adapter

def _calc_ionic_strength_molarity_from_props(
    molarities: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_molarity_unit: str = "mol/L",
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate molarity-based ionic strength from unit-aware molarity inputs."""
    # SECTION: Normalize unit-aware molarities
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    normalized_molarities = to_dict(
        molarities,
        output_molarity_unit,
        unit_conversion_fn=conversion_fn,
    )

    # ? Component metadata can remap keys and enforce component order.
    normalized_molarities = _configure_component_values(
        normalized_molarities,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molarities",
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
    _validate_same_keys(normalized_molarities, normalized_charges)

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            list(normalized_molarities.values()),
            [normalized_charges[key] for key in normalized_molarities],
        ),
    )

# ! ::: Mapping input with component metadata


def _calc_ionic_strength_molarity_with_components_from_mapping(
    molarities: Mapping[str, float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
) -> float:
    """Calculate molarity ionic strength from mapping input and components."""
    # SECTION: Normalize molarity values by component metadata
    configured = config_components_values(
        values=dict(molarities),
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=True,
    )

    # ! Component configuration must succeed before charges are extracted.
    if configured is None:
        raise ValueError("Failed to normalize molarities component values.")

    _, molarity_values = configured

    # SECTION: Extract component charges
    extracted = extract_components_values(
        attribute_name="net_charge",
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    if extracted is None:
        raise ValueError("Failed to extract component charges.")

    _, charge_values = extracted

    # SECTION: Calculate ionic strength
    return cast(
        float,
        _calc_ionic_strength(
            molarity_values,
            charge_values,
        ),
    )


# ! ::: Sequence input with component metadata

def _calc_ionic_strength_molarity_with_components_from_sequence(
    molarities: Sequence[float | int],
    components: Sequence[Component],
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
) -> float:
    """Calculate molarity ionic strength from sequence input and components."""
    # SECTION: Extract component charges
    extracted = extract_components_values(
        attribute_name="net_charge",
        components=list(components),
        component_key=component_key,
        case_sensitive=case_sensitive,
    )

    if extracted is None:
        raise ValueError("Failed to extract component charges.")

    _, charge_values = extracted

    # SECTION: Calculate ionic strength
    return _calc_ionic_strength_molarity_from_sequence(
        molarities=molarities,
        charges=charge_values,
    )


# SECTION: Public exports
__all__ = [
    "_calc_ionic_strength",
    # molality
    "_calc_ionic_strength_molality",
    "_calc_ionic_strength_molality_from_sequence",
    "_calc_ionic_strength_molality_from_mapping",
    "_calc_ionic_strength_molality_from_props",
    "_calc_ionic_strength_molality_with_components_from_sequence",
    "_calc_ionic_strength_molality_with_components_from_mapping",
    # molarity
    "_calc_ionic_strength_molarity",
    "_calc_ionic_strength_molarity_from_sequence",
    "_calc_ionic_strength_molarity_from_mapping",
    "_calc_ionic_strength_molarity_from_props",
    "_calc_ionic_strength_molarity_with_components_from_sequence",
    "_calc_ionic_strength_molarity_with_components_from_mapping",
]
