"""Heat-capacity polynomial public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, ScalarValue
from pythermodb_settings.utils import to_annotated_value

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
    *,
    name: str = "heat_capacity",
    description: str = "Calculate heat capacity from a temperature polynomial.",
    unit: str | None = "J/(mol.K)",
    symbol: str | None = "Cp",
) -> AnnotatedValue[float]:
    """Calculate annotated ``Cp = A + B*T + C*T**2 + D*T**3 + E*T**4 + F*T**5``."""
    return to_annotated_value(
        float(_calc_heat_capacity_polynomial(
            _pos(temperature, "temperature"),
            _scalar(A, "A"),
            _scalar(B, "B"),
            _scalar(C, "C"),
            _scalar(D, "D"),
            _scalar(E, "E"),
            _scalar(F, "F"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_heat_capacity_polynomial",
    )


def calc_enthalpy_change_from_cp_polynomial(
    initial_temperature: ScalarValue,
    final_temperature: ScalarValue,
    A: ScalarValue,
    B: ScalarValue = 0.0,
    C: ScalarValue = 0.0,
    D: ScalarValue = 0.0,
    E: ScalarValue = 0.0,
    F: ScalarValue = 0.0,
    *,
    name: str = "enthalpy_change",
    description: str = "Calculate enthalpy change from a heat-capacity polynomial.",
    unit: str | None = "J/mol",
    symbol: str | None = "delta_H",
) -> AnnotatedValue[float]:
    """Calculate annotated analytical ``integral(Cp dT)`` for a Cp polynomial."""
    return to_annotated_value(
        float(_calc_enthalpy_change_from_cp_polynomial(
            _pos(initial_temperature, "initial_temperature"),
            _pos(final_temperature, "final_temperature"),
            _scalar(A, "A"),
            _scalar(B, "B"),
            _scalar(C, "C"),
            _scalar(D, "D"),
            _scalar(E, "E"),
            _scalar(F, "F"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_enthalpy_change_from_cp_polynomial",
    )


def calc_entropy_change_from_cp_polynomial(
    initial_temperature: ScalarValue,
    final_temperature: ScalarValue,
    A: ScalarValue,
    B: ScalarValue = 0.0,
    C: ScalarValue = 0.0,
    D: ScalarValue = 0.0,
    E: ScalarValue = 0.0,
    F: ScalarValue = 0.0,
    *,
    name: str = "entropy_change",
    description: str = "Calculate entropy change from a heat-capacity polynomial.",
    unit: str | None = "J/(mol.K)",
    symbol: str | None = "delta_S",
) -> AnnotatedValue[float]:
    """Calculate annotated analytical ``integral(Cp/T dT)`` for a Cp polynomial."""
    return to_annotated_value(
        float(_calc_entropy_change_from_cp_polynomial(
            _pos(initial_temperature, "initial_temperature"),
            _pos(final_temperature, "final_temperature"),
            _scalar(A, "A"),
            _scalar(B, "B"),
            _scalar(C, "C"),
            _scalar(D, "D"),
            _scalar(E, "E"),
            _scalar(F, "F"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_entropy_change_from_cp_polynomial",
    )


__all__ = [
    "calc_heat_capacity_polynomial",
    "calc_enthalpy_change_from_cp_polynomial",
    "calc_entropy_change_from_cp_polynomial",
]
