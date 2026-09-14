"""Core density definitions and low-level mixing relations."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Validators
def _as_finite_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert scalar, sequence, or array input to finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _as_positive_array(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert numeric input to finite positive float64 array."""
    arr = _as_finite_array(values, name)
    # ! Densities, molar masses, volumes, pressure scales, and temperatures are positive here.
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")
    return arr


def _validate_fractions(values: NDArray[np.float64], name: str) -> None:
    """Validate non-negative fractions that close along the last axis."""
    # ! Phase volume fractions must form a physical partition of volume.
    if np.any(values < 0.0):
        raise ValueError(f"{name} values must be non-negative.")
    if not np.allclose(values.sum(axis=-1), 1.0):
        raise ValueError(f"{name} values must sum to 1.0 along the last axis.")


# SECTION: Core numeric calculations
def _calc_density_from_molar_volume(
    molar_mass: NumericInput,
    molar_volume: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate mass density from molar mass and molar volume.

    Equation: ``rho = M / Vm``. Molar mass is kg/mol and molar volume is
    m3/mol, so output density is kg/m3. This is a fundamental definition and
    supports scalar, 1-D, and 2-D broadcast-compatible inputs.
    """
    m = _as_positive_array(molar_mass, "molar_mass")
    vm = _as_positive_array(molar_volume, "molar_volume")
    try:
        rho = m / vm
    except ValueError as exc:
        raise ValueError("molar_mass and molar_volume must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], rho))


def _calc_gas_density_from_compressibility(
    molar_mass: NumericInput,
    pressure: NumericInput,
    temperature: NumericInput,
    compressibility_factor: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate gas density from a supplied compressibility factor.

    Equation: ``rho = M*P/(Z*R*T)``. Inputs use kg/mol, Pa, K, dimensionless
    ``Z``, and J/(mol.K), producing kg/m3. This is EOS-independent; the caller
    supplies the compressibility factor from another model or measurement.
    """
    m = _as_positive_array(molar_mass, "molar_mass")
    p = _as_positive_array(pressure, "pressure")
    t = _as_positive_array(temperature, "temperature")
    z = _as_positive_array(compressibility_factor, "compressibility_factor")
    r = _as_positive_array(gas_constant, "gas_constant")
    try:
        rho = m * p / (z * r * t)
    except ValueError as exc:
        raise ValueError("gas-density inputs must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], rho))


def _calc_multiphase_density_from_volume_fractions(
    phase_densities: NumericInput,
    phase_volume_fractions: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate mixture density from phase volume fractions.

    Equation: ``rho_m = sum_k f_k*rho_k``. Phase densities are kg/m3 and phase
    volume fractions are dimensionless fractions that sum to one along the last
    axis. This is a low-level n-phase mixture definition.
    """
    rho = _as_positive_array(phase_densities, "phase_densities")
    f = _as_finite_array(phase_volume_fractions, "phase_volume_fractions")
    _validate_fractions(f, "phase_volume_fractions")
    if rho.shape != f.shape:
        raise ValueError("phase_densities and phase_volume_fractions must have the same shape.")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], (f * rho).sum(axis=-1)))


__all__ = [
    "_calc_density_from_molar_volume",
    "_calc_gas_density_from_compressibility",
    "_calc_multiphase_density_from_volume_fractions",
]

