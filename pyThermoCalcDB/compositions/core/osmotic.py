"""Core osmotic composition primitives."""

# import libs
from collections.abc import Mapping
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import (
    NumericArrayInput,
    _return_scalar_if_zero_dim,
    _validate_non_negative_array,
    _validate_positive_array,
)
from ...configs.constants import R_J_molK

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_component_array(
    values: NumericInput,
    name: str,
) -> NDArray[np.float64]:
    """Convert component values to a finite 1-D/2-D float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim not in (1, 2):
        raise ValueError(f"{name} must be a one- or two-dimensional array.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    _validate_non_negative_array(arr, name)
    return cast(NDArray[np.float64], arr)


# SECTION: Core numeric calculations
def _calc_osmolarity(
    species_molarities: NumericInput,
) -> float | NDArray[np.float64]:
    """Sum supplied dissolved particle molarities along the component axis."""
    c = _as_component_array(species_molarities, "species_molarities")
    return _return_scalar_if_zero_dim(np.sum(c, axis=-1))


def _calc_osmolality(
    species_molalities: NumericInput,
) -> float | NDArray[np.float64]:
    """Sum supplied dissolved particle molalities along the component axis."""
    b = _as_component_array(species_molalities, "species_molalities")
    return _return_scalar_if_zero_dim(np.sum(b, axis=-1))


def _calc_osmolarity_from_mapping(
    species_molarities: Mapping[str, float | int],
) -> float:
    """Calculate osmolarity from a species-keyed molarity mapping."""
    return float(_calc_osmolarity(list(species_molarities.values())))


def _calc_osmolality_from_mapping(
    species_molalities: Mapping[str, float | int],
) -> float:
    """Calculate osmolality from a species-keyed molality mapping."""
    return float(_calc_osmolality(list(species_molalities.values())))


def _calc_ideal_osmotic_pressure(
    molar_concentration: NumericInput,
    temperature: NumericInput,
    vant_hoff_factor: NumericInput = 1.0,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate ideal osmotic pressure ``Pi = i*c*R*T``."""
    c = np.asarray(molar_concentration, dtype=np.float64)
    t = np.asarray(temperature, dtype=np.float64)
    i = np.asarray(vant_hoff_factor, dtype=np.float64)
    r = np.asarray(gas_constant, dtype=np.float64)
    for name, arr in (
        ("molar_concentration", c),
        ("temperature", t),
        ("vant_hoff_factor", i),
        ("gas_constant", r),
    ):
        if arr.ndim > 2:
            raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} values must be finite.")
    _validate_non_negative_array(c, "molar_concentration")
    _validate_positive_array(t, "temperature")
    _validate_positive_array(i, "vant_hoff_factor")
    _validate_positive_array(r, "gas_constant")
    return _return_scalar_if_zero_dim(c * r * t * i)


__all__ = [
    "_calc_osmolarity",
    "_calc_osmolality",
    "_calc_osmolarity_from_mapping",
    "_calc_osmolality_from_mapping",
    "_calc_ideal_osmotic_pressure",
]
