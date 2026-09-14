"""Core heat-capacity polynomial calculations."""

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
    """Convert numeric input to finite scalar, 1-D, or 2-D float64 array."""
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim > 2:
        raise ValueError(f"{name} must be scalar, one-dimensional, or two-dimensional.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} values must be finite.")
    return cast(NDArray[np.float64], arr)


def _as_positive_temperature(values: NumericInput, name: str) -> NDArray[np.float64]:
    """Convert absolute temperature input to finite positive float64 array."""
    arr = _as_finite_array(values, name)
    # ! Absolute temperatures in heat-capacity correlations must be positive.
    if np.any(arr <= 0.0):
        raise ValueError(f"{name} values must be greater than zero K.")
    return arr


# SECTION: Core numeric calculations
def _calc_heat_capacity_polynomial(
    temperature: NumericInput,
    A: NumericInput,
    B: NumericInput = 0.0,
    C: NumericInput = 0.0,
    D: NumericInput = 0.0,
    E: NumericInput = 0.0,
    F: NumericInput = 0.0,
) -> float | NDArray[np.float64]:
    """Calculate heat capacity from a sixth-term temperature polynomial.

    Equation: ``Cp = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5``.
    Temperature is in K. Coefficients must be on a consistent heat-capacity
    basis, commonly J/(mol.K). The output uses that same heat-capacity basis.
    This empirical correlation accepts scalar, 1-D, and 2-D broadcastable inputs.
    """
    t = _as_positive_temperature(temperature, "temperature")
    a = _as_finite_array(A, "A")
    b = _as_finite_array(B, "B")
    c = _as_finite_array(C, "C")
    d = _as_finite_array(D, "D")
    e = _as_finite_array(E, "E")
    f = _as_finite_array(F, "F")
    try:
        cp = a + b * t + c * t**2 + d * t**3 + e * t**4 + f * t**5
    except ValueError as exc:
        raise ValueError("temperature and polynomial coefficients must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], cp))


def _calc_enthalpy_change_from_cp_polynomial(
    initial_temperature: NumericInput,
    final_temperature: NumericInput,
    A: NumericInput,
    B: NumericInput = 0.0,
    C: NumericInput = 0.0,
    D: NumericInput = 0.0,
    E: NumericInput = 0.0,
    F: NumericInput = 0.0,
) -> float | NDArray[np.float64]:
    """Calculate analytical enthalpy change from a Cp polynomial.

    Equation: ``dH = integral(Cp dT)`` from ``T1`` to ``T2``. Temperatures are
    in K; coefficients are on a consistent heat-capacity basis. The output is
    energy on the matching molar or mass basis, for example J/mol. This is the
    exact integral of the empirical polynomial and preserves scalar/array shape.
    """
    t1 = _as_positive_temperature(initial_temperature, "initial_temperature")
    t2 = _as_positive_temperature(final_temperature, "final_temperature")
    a = _as_finite_array(A, "A")
    b = _as_finite_array(B, "B")
    c = _as_finite_array(C, "C")
    d = _as_finite_array(D, "D")
    e = _as_finite_array(E, "E")
    f = _as_finite_array(F, "F")
    try:
        dh = (
            a * (t2 - t1)
            + b / 2.0 * (t2**2 - t1**2)
            + c / 3.0 * (t2**3 - t1**3)
            + d / 4.0 * (t2**4 - t1**4)
            + e / 5.0 * (t2**5 - t1**5)
            + f / 6.0 * (t2**6 - t1**6)
        )
    except ValueError as exc:
        raise ValueError("temperatures and polynomial coefficients must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], dh))


def _calc_entropy_change_from_cp_polynomial(
    initial_temperature: NumericInput,
    final_temperature: NumericInput,
    A: NumericInput,
    B: NumericInput = 0.0,
    C: NumericInput = 0.0,
    D: NumericInput = 0.0,
    E: NumericInput = 0.0,
    F: NumericInput = 0.0,
) -> float | NDArray[np.float64]:
    """Calculate analytical entropy change from a Cp polynomial.

    Equation: ``dS = integral(Cp/T dT)`` from ``T1`` to ``T2``. Temperatures are
    positive absolute temperatures in K. The output uses the same entropy basis
    as the coefficients, commonly J/(mol.K). This exact integral preserves
    scalar/array behavior.
    """
    t1 = _as_positive_temperature(initial_temperature, "initial_temperature")
    t2 = _as_positive_temperature(final_temperature, "final_temperature")
    a = _as_finite_array(A, "A")
    b = _as_finite_array(B, "B")
    c = _as_finite_array(C, "C")
    d = _as_finite_array(D, "D")
    e = _as_finite_array(E, "E")
    f = _as_finite_array(F, "F")
    try:
        ds = (
            a * np.log(t2 / t1)
            + b * (t2 - t1)
            + c / 2.0 * (t2**2 - t1**2)
            + d / 3.0 * (t2**3 - t1**3)
            + e / 4.0 * (t2**4 - t1**4)
            + f / 5.0 * (t2**5 - t1**5)
        )
    except ValueError as exc:
        raise ValueError("temperatures and polynomial coefficients must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], ds))


__all__ = [
    "_calc_heat_capacity_polynomial",
    "_calc_enthalpy_change_from_cp_polynomial",
    "_calc_entropy_change_from_cp_polynomial",
]
