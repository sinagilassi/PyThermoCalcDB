"""Core reduced thermodynamic property calculations."""

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
    """Convert numeric scalar, 1-D, or 2-D input to a finite float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _validate_positive(values: NDArray[np.float64], name: str) -> None:
    """Validate strictly positive denominator-like property values."""
    # ! Critical properties are denominators in reduced-property definitions.
    if np.any(values <= 0.0):
        raise ValueError(f"{name} values must be greater than zero.")


# SECTION: Core numeric calculations
def _calc_reduced_temperature(
    temperature: NumericInput,
    critical_temperature: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reduced temperature.

    The reduced temperature is defined as ``Tr = T / Tc``.

    Parameters
    ----------
    temperature : NumericInput
        Absolute temperature, K.
    critical_temperature : NumericInput
        Critical temperature, K.

    Returns
    -------
    float | NDArray[np.float64]
        Reduced temperature, dimensionless.

    Raises
    ------
    ValueError
        If inputs are non-finite, if ``critical_temperature`` is not positive,
        or if the inputs are not broadcast-compatible.

    Notes
    -----
    This is an exact thermodynamic definition, not an empirical correlation.
    Scalar inputs return ``float``; array inputs return ``float64`` arrays.
    For 2-D inputs, axis 0 is states and axis 1 is components/properties.
    """
    t = _as_finite_array(temperature, "temperature")
    tc = _as_finite_array(critical_temperature, "critical_temperature")
    _validate_positive(tc, "critical_temperature")
    try:
        result = t / tc
    except ValueError as exc:
        raise ValueError("temperature and critical_temperature must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], result))


def _calc_reduced_pressure(
    pressure: NumericInput,
    critical_pressure: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reduced pressure.

    The reduced pressure is defined as ``Pr = P / Pc``.

    Parameters
    ----------
    pressure : NumericInput
        Pressure, Pa.
    critical_pressure : NumericInput
        Critical pressure, Pa.

    Returns
    -------
    float | NDArray[np.float64]
        Reduced pressure, dimensionless.

    Notes
    -----
    This exact definition assumes a shared pressure unit basis. Scalar inputs
    return ``float``; array inputs return ``float64`` arrays.
    """
    p = _as_finite_array(pressure, "pressure")
    pc = _as_finite_array(critical_pressure, "critical_pressure")
    _validate_positive(pc, "critical_pressure")
    try:
        result = p / pc
    except ValueError as exc:
        raise ValueError("pressure and critical_pressure must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], result))


def _calc_reduced_volume(
    molar_volume: NumericInput,
    critical_molar_volume: NumericInput,
) -> float | NDArray[np.float64]:
    """Calculate reduced molar volume.

    The reduced volume is defined as ``Vr = Vm / Vc``.

    Parameters
    ----------
    molar_volume : NumericInput
        Molar volume, m3/mol.
    critical_molar_volume : NumericInput
        Critical molar volume, m3/mol.

    Returns
    -------
    float | NDArray[np.float64]
        Reduced volume, dimensionless.

    Notes
    -----
    This exact definition assumes a shared molar-volume unit basis. Scalar
    inputs return ``float``; array inputs return ``float64`` arrays.
    """
    v = _as_finite_array(molar_volume, "molar_volume")
    vc = _as_finite_array(critical_molar_volume, "critical_molar_volume")
    _validate_positive(vc, "critical_molar_volume")
    try:
        result = v / vc
    except ValueError as exc:
        raise ValueError("molar_volume and critical_molar_volume must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], result))


__all__ = [
    "_calc_reduced_temperature",
    "_calc_reduced_pressure",
    "_calc_reduced_volume",
]
