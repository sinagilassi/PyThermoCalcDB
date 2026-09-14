"""Rackett saturated-liquid density public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, ScalarValue
from pythermodb_settings.utils import to_annotated_value

# locals
from ..utils.conversions import _pos, _scalar
from .core.rackett import (
    _calc_rackett_constant_from_acentric_factor,
    _calc_saturated_liquid_molar_volume_rackett,
)


# SECTION: Public wrappers
def calc_saturated_liquid_molar_volume_rackett(
    temperature: ScalarValue,
    critical_temperature: ScalarValue,
    critical_pressure: ScalarValue,
    rackett_constant: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "saturated_liquid_molar_volume",
    description: str = "Calculate saturated-liquid molar volume with the Rackett equation.",
    unit: str | None = "m3/mol",
    symbol: str | None = "Vs",
) -> AnnotatedValue[float]:
    """Calculate annotated saturated-liquid molar volume with Rackett."""
    return to_annotated_value(
        float(_calc_saturated_liquid_molar_volume_rackett(
            _pos(temperature, "temperature"),
            _pos(critical_temperature, "critical_temperature"),
            _pos(critical_pressure, "critical_pressure"),
            _pos(rackett_constant, "rackett_constant"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_saturated_liquid_molar_volume_rackett",
    )


def calc_rackett_constant_from_acentric_factor(
    acentric_factor: ScalarValue,
    *,
    name: str = "rackett_constant",
    description: str = "Estimate Rackett constant from acentric factor.",
    symbol: str | None = "Z_RA",
) -> AnnotatedValue[float]:
    """Calculate annotated empirical Rackett constant estimate from acentric factor."""
    return to_annotated_value(
        float(_calc_rackett_constant_from_acentric_factor(_scalar(acentric_factor, "acentric_factor"))),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_rackett_constant_from_acentric_factor",
    )


__all__ = [
    "calc_saturated_liquid_molar_volume_rackett",
    "calc_rackett_constant_from_acentric_factor",
]
