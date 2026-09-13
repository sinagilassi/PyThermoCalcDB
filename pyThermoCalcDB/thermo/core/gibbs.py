"""Core Gibbs-energy identity calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp, Temperature, UnitConversionFn

# locals
from pythermocalcdb.utils.conversions import (
    NumericArrayInput,
    _generic_temperature,
    _resolve_unit_conversion_fn,
    _return_scalar_if_zero_dim,
    _scalar,
    _validate_custom_prop_scalar,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Core numeric calculations

def _calc_gibbs_energy(
    enthalpy: NumericInput,
    temperature: NumericInput,
    entropy: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``G = H - T*S`` using finite numeric inputs.

    Inputs may be scalars, 1-D arrays, or 2-D arrays when their shapes are
    broadcast-compatible. For 2-D inputs, axis 0 is states and axis 1 is
    component/property values.
    """
    # SECTION: Normalize and validate
    h = np.asarray(enthalpy, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    s = np.asarray(entropy, dtype=np.float64)

    if any(arr.ndim > 2 for arr in (h, t, s)):
        raise ValueError(
            "enthalpy, temperature, and entropy must be scalar, 1-D, or 2-D values."
        )
    if not all(np.all(np.isfinite(arr)) for arr in (h, t, s)):
        raise ValueError(
            "enthalpy, temperature, and entropy values must be finite."
        )

    try:
        gibbs_energy = h - t * s
    except ValueError as exc:
        raise ValueError(
            "enthalpy, temperature, and entropy must be broadcast-compatible."
        ) from exc

    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], gibbs_energy))


# SECTION: Props adapters

def _calc_gibbs_energy_from_props(
    enthalpy: CustomProp,
    temperature: Temperature,
    entropy: CustomProp,
    output_enthalpy_unit: str | None = None,
    output_entropy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate Gibbs energy from unit-aware scalar inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(enthalpy, "enthalpy")
    _validate_custom_prop_scalar(entropy, "entropy")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    h = _scalar(enthalpy, "enthalpy", output_enthalpy_unit, conversion_fn)
    t = _generic_temperature(
        temperature,
        output_temperature_unit,
        conversion_fn,
    )
    s = _scalar(entropy, "entropy", output_entropy_unit, conversion_fn)

    # SECTION: Calculate Gibbs energy
    return float(_calc_gibbs_energy(h, t, s))


def _calc_gibbs_energy_from_scalars(
    enthalpy: float | int | CustomProp,
    temperature: Temperature,
    entropy: float | int | CustomProp,
    output_enthalpy_unit: str | None = None,
    output_entropy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate Gibbs energy."""
    # SECTION: Normalize scalar inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    h = _scalar(enthalpy, "enthalpy", output_enthalpy_unit, conversion_fn)
    t = _generic_temperature(
        temperature,
        output_temperature_unit,
        conversion_fn,
    )
    s = _scalar(entropy, "entropy", output_entropy_unit, conversion_fn)

    # SECTION: Calculate Gibbs energy
    return float(_calc_gibbs_energy(h, t, s))


def _calc_gibbs_energy_change(
    enthalpy_change: NumericInput,
    entropy_change: NumericInput,
    temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``dG = dH - T*dS`` using finite numeric inputs."""
    # SECTION: Delegate to the generic identity
    return _calc_gibbs_energy(
        enthalpy=enthalpy_change,
        temperature=temperature,
        entropy=entropy_change,
    )


def _calc_gibbs_energy_change_from_props(
    enthalpy_change: CustomProp,
    entropy_change: CustomProp,
    temperature: Temperature,
    output_enthalpy_change_unit: str | None = None,
    output_entropy_change_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate Gibbs energy change from unit-aware scalar inputs."""
    return _calc_gibbs_energy_from_props(
        enthalpy=enthalpy_change,
        temperature=temperature,
        entropy=entropy_change,
        output_enthalpy_unit=output_enthalpy_change_unit,
        output_entropy_unit=output_entropy_change_unit,
        output_temperature_unit=output_temperature_unit,
        unit_conversion_fn=unit_conversion_fn,
    )


def _calc_gibbs_energy_change_from_scalars(
    enthalpy_change: float | int | CustomProp,
    entropy_change: float | int | CustomProp,
    temperature: Temperature,
    output_enthalpy_change_unit: str | None = None,
    output_entropy_change_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate Gibbs energy change."""
    return _calc_gibbs_energy_from_scalars(
        enthalpy=enthalpy_change,
        temperature=temperature,
        entropy=entropy_change,
        output_enthalpy_unit=output_enthalpy_change_unit,
        output_entropy_unit=output_entropy_change_unit,
        output_temperature_unit=output_temperature_unit,
        unit_conversion_fn=unit_conversion_fn,
    )


# SECTION: Core exports
__all__ = [
    "_calc_gibbs_energy",
    "_calc_gibbs_energy_from_props",
    "_calc_gibbs_energy_from_scalars",
    "_calc_gibbs_energy_change",
    "_calc_gibbs_energy_change_from_props",
    "_calc_gibbs_energy_change_from_scalars",
]
