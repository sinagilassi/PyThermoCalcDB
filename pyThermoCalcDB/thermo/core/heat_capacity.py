"""Core heat-capacity relationship calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, ScalarValue, UnitConversionFn

# locals
from pythermocalcdb.utils.conversions import (
    NumericArrayInput,
    _pos,
    _return_scalar_if_zero_dim,
    _validate_custom_prop_scalar,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators

def _as_positive_heat_capacity_array(
    values: NumericInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert heat-capacity inputs to finite positive float64 arrays."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, 1-D, or 2-D values.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return cast(NDArray[np.float64], arr)


# SECTION: Core numeric calculations

def _calc_ideal_gas_cv_from_cp(
    cp: NumericInput,
    gas_constant: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ideal-gas ``Cv = Cp - R``.

    Inputs may be scalars, 1-D arrays, or 2-D arrays when their shapes are
    broadcast-compatible. For 2-D inputs, axis 0 is states and axis 1 is
    component/property values.
    """
    # SECTION: Normalize and validate
    cp_arr = _as_positive_heat_capacity_array(cp, "cp")
    r_arr = _as_positive_heat_capacity_array(gas_constant, "gas_constant")

    try:
        cv = cp_arr - r_arr
    except ValueError as exc:
        raise ValueError(
            "cp and gas_constant must be broadcast-compatible."
        ) from exc

    # ! Ideal-gas Cv must remain physically positive after subtracting R.
    if np.any(cv <= 0.0):
        raise ValueError(
            "calculated cv must be greater than zero; cp must be greater than R."
        )
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], cv))


def _calc_ideal_gas_cp_from_cv(
    cv: NumericInput,
    gas_constant: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ideal-gas ``Cp = Cv + R``."""
    # SECTION: Normalize and validate
    cv_arr = _as_positive_heat_capacity_array(cv, "cv")
    r_arr = _as_positive_heat_capacity_array(gas_constant, "gas_constant")

    try:
        cp = cv_arr + r_arr
    except ValueError as exc:
        raise ValueError(
            "cv and gas_constant must be broadcast-compatible."
        ) from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], cp))


def _calc_heat_capacity_ratio(
    cp: NumericInput,
    cv: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate heat-capacity ratio ``gamma = Cp/Cv``."""
    # SECTION: Normalize and validate
    cp_arr = _as_positive_heat_capacity_array(cp, "cp")
    cv_arr = _as_positive_heat_capacity_array(cv, "cv")

    try:
        ratio = cp_arr / cv_arr
    except ValueError as exc:
        raise ValueError("cp and cv must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], ratio))


# SECTION: Props adapters

def _calc_ideal_gas_cv_from_cp_from_props(
    cp: CustomProp,
    gas_constant: CustomProp,
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal-gas Cv from unit-aware scalar heat capacities."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(cp, "cp")
    _validate_custom_prop_scalar(gas_constant, "gas_constant")

    # SECTION: Normalize and calculate
    cp_value = _pos(cp, "cp", output_heat_capacity_unit, unit_conversion_fn)
    r_value = _pos(
        gas_constant,
        "gas_constant",
        output_heat_capacity_unit,
        unit_conversion_fn,
    )
    return float(_calc_ideal_gas_cv_from_cp(cp_value, r_value))


def _calc_ideal_gas_cp_from_cv_from_props(
    cv: CustomProp,
    gas_constant: CustomProp,
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal-gas Cp from unit-aware scalar heat capacities."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(cv, "cv")
    _validate_custom_prop_scalar(gas_constant, "gas_constant")

    # SECTION: Normalize and calculate
    cv_value = _pos(cv, "cv", output_heat_capacity_unit, unit_conversion_fn)
    r_value = _pos(
        gas_constant,
        "gas_constant",
        output_heat_capacity_unit,
        unit_conversion_fn,
    )
    return float(_calc_ideal_gas_cp_from_cv(cv_value, r_value))


def _calc_heat_capacity_ratio_from_props(
    cp: CustomProp,
    cv: CustomProp,
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate heat-capacity ratio from unit-aware scalar heat capacities."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(cp, "cp")
    _validate_custom_prop_scalar(cv, "cv")

    # SECTION: Normalize and calculate
    cp_value = _pos(cp, "cp", output_heat_capacity_unit, unit_conversion_fn)
    cv_value = _pos(cv, "cv", output_heat_capacity_unit, unit_conversion_fn)
    return float(_calc_heat_capacity_ratio(cp_value, cv_value))


# SECTION: Scalar adapters

def _calc_ideal_gas_cv_from_cp_from_scalars(
    cp: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate ideal-gas Cv."""
    cp_value = _pos(cp, "cp", output_heat_capacity_unit, unit_conversion_fn)
    r_value = _pos(gas_constant, "gas_constant")
    return float(_calc_ideal_gas_cv_from_cp(cp_value, r_value))


def _calc_ideal_gas_cp_from_cv_from_scalars(
    cv: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate ideal-gas Cp."""
    cv_value = _pos(cv, "cv", output_heat_capacity_unit, unit_conversion_fn)
    r_value = _pos(gas_constant, "gas_constant")
    return float(_calc_ideal_gas_cp_from_cv(cv_value, r_value))


def _calc_heat_capacity_ratio_from_scalars(
    cp: ScalarValue,
    cv: ScalarValue,
    output_heat_capacity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate heat-capacity ratio."""
    cp_value = _pos(cp, "cp", output_heat_capacity_unit, unit_conversion_fn)
    cv_value = _pos(cv, "cv", output_heat_capacity_unit, unit_conversion_fn)
    return float(_calc_heat_capacity_ratio(cp_value, cv_value))


# SECTION: Core exports
__all__ = [
    "_calc_ideal_gas_cv_from_cp",
    "_calc_ideal_gas_cp_from_cv",
    "_calc_heat_capacity_ratio",
    "_calc_ideal_gas_cv_from_cp_from_props",
    "_calc_ideal_gas_cp_from_cv_from_props",
    "_calc_heat_capacity_ratio_from_props",
    "_calc_ideal_gas_cv_from_cp_from_scalars",
    "_calc_ideal_gas_cp_from_cv_from_scalars",
    "_calc_heat_capacity_ratio_from_scalars",
]
