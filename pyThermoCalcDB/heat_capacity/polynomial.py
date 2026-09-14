"""Heat-capacity polynomial public wrappers."""

# import libs
from pythermodb_settings.models import ScalarValue

# locals
from ..utils.conversions import _pos, _scalar
from .core.polynomial import (
    _calc_enthalpy_change_from_cp_polynomial,
    _calc_entropy_change_from_cp_polynomial,
    _calc_heat_capacity_polynomial,
)


# SECTION: Public wrappers
def calc_heat_capacity_polynomial(
    temperature: ScalarValue,
    A: ScalarValue,
    B: ScalarValue = 0.0,
    C: ScalarValue = 0.0,
    D: ScalarValue = 0.0,
    E: ScalarValue = 0.0,
    F: ScalarValue = 0.0,
) -> float:
    """Calculate ``Cp = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5``."""
    return float(_calc_heat_capacity_polynomial(
        _pos(temperature, "temperature"),
        _scalar(A, "A"),
        _scalar(B, "B"),
        _scalar(C, "C"),
        _scalar(D, "D"),
        _scalar(E, "E"),
        _scalar(F, "F"),
    ))


def calc_enthalpy_change_from_cp_polynomial(
    initial_temperature: ScalarValue,
    final_temperature: ScalarValue,
    A: ScalarValue,
    B: ScalarValue = 0.0,
    C: ScalarValue = 0.0,
    D: ScalarValue = 0.0,
    E: ScalarValue = 0.0,
    F: ScalarValue = 0.0,
) -> float:
    """Calculate analytical ``integral(Cp dT)`` for a heat-capacity polynomial."""
    return float(_calc_enthalpy_change_from_cp_polynomial(
        _pos(initial_temperature, "initial_temperature"),
        _pos(final_temperature, "final_temperature"),
        _scalar(A, "A"),
        _scalar(B, "B"),
        _scalar(C, "C"),
        _scalar(D, "D"),
        _scalar(E, "E"),
        _scalar(F, "F"),
    ))


def calc_entropy_change_from_cp_polynomial(
    initial_temperature: ScalarValue,
    final_temperature: ScalarValue,
    A: ScalarValue,
    B: ScalarValue = 0.0,
    C: ScalarValue = 0.0,
    D: ScalarValue = 0.0,
    E: ScalarValue = 0.0,
    F: ScalarValue = 0.0,
) -> float:
    """Calculate analytical ``integral(Cp/T dT)`` for a heat-capacity polynomial."""
    return float(_calc_entropy_change_from_cp_polynomial(
        _pos(initial_temperature, "initial_temperature"),
        _pos(final_temperature, "final_temperature"),
        _scalar(A, "A"),
        _scalar(B, "B"),
        _scalar(C, "C"),
        _scalar(D, "D"),
        _scalar(E, "E"),
        _scalar(F, "F"),
    ))


__all__ = [
    "calc_heat_capacity_polynomial",
    "calc_enthalpy_change_from_cp_polynomial",
    "calc_entropy_change_from_cp_polynomial",
]
