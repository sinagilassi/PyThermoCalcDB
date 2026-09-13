"""Core Helmholtz-energy identity calculations."""

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

def _calc_helmholtz_energy(
    internal_energy: NumericInput,
    temperature: NumericInput,
    entropy: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate ``A = U - T*S`` using finite numeric inputs.

    Inputs may be scalars, 1-D arrays, or 2-D arrays when their shapes are
    broadcast-compatible. For 2-D inputs, axis 0 is states and axis 1 is
    component/property values.
    """
    # SECTION: Normalize and validate
    u = np.asarray(internal_energy, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    s = np.asarray(entropy, dtype=np.float64)

    if any(arr.ndim > 2 for arr in (u, t, s)):
        raise ValueError(
            "internal_energy, temperature, and entropy must be scalar, "
            "1-D, or 2-D values."
        )
    if not all(np.all(np.isfinite(arr)) for arr in (u, t, s)):
        raise ValueError(
            "internal_energy, temperature, and entropy values must be finite."
        )

    try:
        helmholtz_energy = u - t * s
    except ValueError as exc:
        raise ValueError(
            "internal_energy, temperature, and entropy must be "
            "broadcast-compatible."
        ) from exc

    return _return_scalar_if_zero_dim(
        cast(NDArray[np.float64], helmholtz_energy)
    )


# SECTION: Props adapters

def _calc_helmholtz_energy_from_props(
    internal_energy: CustomProp,
    temperature: Temperature,
    entropy: CustomProp,
    output_internal_energy_unit: str | None = None,
    output_entropy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate Helmholtz energy from unit-aware scalar inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(internal_energy, "internal_energy")
    _validate_custom_prop_scalar(entropy, "entropy")

    # SECTION: Normalize units
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    u = _scalar(
        internal_energy,
        "internal_energy",
        output_internal_energy_unit,
        conversion_fn,
    )
    t = _generic_temperature(
        temperature,
        output_temperature_unit,
        conversion_fn,
    )
    s = _scalar(entropy, "entropy", output_entropy_unit, conversion_fn)

    # SECTION: Calculate Helmholtz energy
    return float(_calc_helmholtz_energy(u, t, s))


def _calc_helmholtz_energy_from_scalars(
    internal_energy: float | int | CustomProp,
    temperature: Temperature,
    entropy: float | int | CustomProp,
    output_internal_energy_unit: str | None = None,
    output_entropy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Normalize scalar public inputs and calculate Helmholtz energy."""
    # SECTION: Normalize scalar inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    u = _scalar(
        internal_energy,
        "internal_energy",
        output_internal_energy_unit,
        conversion_fn,
    )
    t = _generic_temperature(
        temperature,
        output_temperature_unit,
        conversion_fn,
    )
    s = _scalar(entropy, "entropy", output_entropy_unit, conversion_fn)

    # SECTION: Calculate Helmholtz energy
    return float(_calc_helmholtz_energy(u, t, s))


# SECTION: Core exports
__all__ = [
    "_calc_helmholtz_energy",
    "_calc_helmholtz_energy_from_props",
    "_calc_helmholtz_energy_from_scalars",
]
