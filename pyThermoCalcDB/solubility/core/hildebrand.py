"""Core Hildebrand solubility-parameter calculations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_finite_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric scalar, 1-D, or 2-D input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _as_positive_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to finite positive float64 array."""
    arr = _as_finite_array(values, name)
    # ! Denominators and absolute thermodynamic quantities must be positive.
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return arr


# SECTION: Core numeric calculations
def _calc_internal_energy_of_vaporization(
    heat_of_vaporization: NumericInput,
    temperature: NumericInput,
    gas_constant: NumericInput = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate internal energy of vaporization.

    Equation: ``delta_U_vap = delta_H_vap - R*T``. Heat of vaporization is in
    J/mol, temperature is K, and gas constant is J/(mol.K), so the output is
    J/mol. This exact ideal-vapor correction supports scalar and array inputs.
    """
    hvap = _as_positive_array(heat_of_vaporization, "heat_of_vaporization")
    t = _as_positive_array(temperature, "temperature")
    r = _as_positive_array(gas_constant, "gas_constant")
    try:
        du = hvap - r * t
    except ValueError as exc:
        raise ValueError("heat_of_vaporization, temperature, and gas_constant must be broadcast-compatible.") from exc
    # ! Hildebrand CED requires positive cohesive energy.
    if np.any(du <= 0.0):
        raise ValueError("calculated internal energy of vaporization must be greater than zero.")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], du))


def _calc_cohesive_energy_density(
    internal_energy_of_vaporization: NumericInput,
    molar_volume: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate cohesive energy density.

    Equation: ``CED = delta_U_vap / Vm``. Internal energy is J/mol and molar
    volume is m3/mol, so output is J/m3. This exact definition preserves
    scalar/array behavior.
    """
    du = _as_positive_array(internal_energy_of_vaporization, "internal_energy_of_vaporization")
    vm = _as_positive_array(molar_volume, "molar_volume")
    try:
        ced = du / vm
    except ValueError as exc:
        raise ValueError("internal_energy_of_vaporization and molar_volume must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], ced))


def _calc_solubility_parameter(
    cohesive_energy_density: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate Hildebrand solubility parameter from cohesive energy density.

    Equation: ``delta = sqrt(CED)``. CED is J/m3, so the output is
    ``sqrt(J/m3)`` or Pa**0.5. This exact transform preserves scalar/array
    behavior.
    """
    ced = _as_positive_array(cohesive_energy_density, "cohesive_energy_density")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], np.sqrt(ced)))


def _calc_solubility_parameter_from_hvap_volume(
    heat_of_vaporization: NumericInput,
    temperature: NumericInput,
    molar_volume: NumericInput,
    gas_constant: NumericInput = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate Hildebrand parameter from heat of vaporization and molar volume.

    Equation: ``delta = sqrt((delta_H_vap - R*T) / Vm)``. Inputs use J/mol, K,
    m3/mol, and J/(mol.K). The output is Pa**0.5. This empirical solubility
    parameter transform assumes the simple Hildebrand cohesive-energy model.
    """
    du = _calc_internal_energy_of_vaporization(heat_of_vaporization, temperature, gas_constant)
    ced = _calc_cohesive_energy_density(du, molar_volume)
    return _calc_solubility_parameter(ced)


def _calc_solubility_parameter_from_hvap_density(
    heat_of_vaporization: NumericInput,
    temperature: NumericInput,
    molar_density: NumericInput,
    gas_constant: NumericInput = 8.31446261815324,
) -> float | NDArray[np.float64]:
    """Calculate Hildebrand parameter from heat of vaporization and molar density.

    Equation: ``delta = sqrt((delta_H_vap - R*T) * rho_m)``. Inputs use J/mol,
    K, mol/m3, and J/(mol.K). The output is Pa**0.5 and preserves scalar/array
    behavior under the Hildebrand cohesive-energy model.
    """
    du = np.asarray(
        _calc_internal_energy_of_vaporization(heat_of_vaporization, temperature, gas_constant),
        dtype=np.float64,
    )
    rho_m = _as_positive_array(molar_density, "molar_density")
    try:
        ced = du * rho_m
    except ValueError as exc:
        raise ValueError("internal energy and molar_density must be broadcast-compatible.") from exc
    return _calc_solubility_parameter(cast(NDArray[np.float64], ced))


__all__ = [
    "_calc_internal_energy_of_vaporization",
    "_calc_cohesive_energy_density",
    "_calc_solubility_parameter",
    "_calc_solubility_parameter_from_hvap_volume",
    "_calc_solubility_parameter_from_hvap_density",
]
