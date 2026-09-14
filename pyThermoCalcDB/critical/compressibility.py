"""Critical compressibility public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, ScalarValue
from pythermodb_settings.utils import to_annotated_value

# locals
from ..utils.conversions import _pos
from .core.compressibility import (
    _calc_critical_compressibility_factor,
    _calc_critical_pressure_from_zc,
    _calc_critical_temperature_from_zc,
    _calc_critical_volume_from_zc,
)


# SECTION: Public wrappers
def calc_critical_compressibility_factor(
    critical_pressure: ScalarValue,
    critical_molar_volume: ScalarValue,
    critical_temperature: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "critical_compressibility_factor",
    description: str = "Calculate critical compressibility factor.",
    symbol: str | None = "Zc",
) -> AnnotatedValue[float]:
    """Calculate annotated ``Zc = Pc*Vc/(R*Tc)`` from critical properties."""
    return to_annotated_value(
        float(_calc_critical_compressibility_factor(
            _pos(critical_pressure, "critical_pressure"),
            _pos(critical_molar_volume, "critical_molar_volume"),
            _pos(critical_temperature, "critical_temperature"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_critical_compressibility_factor",
    )


def calc_critical_volume_from_zc(
    critical_compressibility_factor: ScalarValue,
    critical_temperature: ScalarValue,
    critical_pressure: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "critical_molar_volume",
    description: str = "Calculate critical molar volume from critical compressibility factor.",
    unit: str | None = "m3/mol",
    symbol: str | None = "Vc",
) -> AnnotatedValue[float]:
    """Calculate annotated critical molar volume from ``Zc`` in SI-consistent units."""
    return to_annotated_value(
        float(_calc_critical_volume_from_zc(
            _pos(critical_compressibility_factor, "critical_compressibility_factor"),
            _pos(critical_temperature, "critical_temperature"),
            _pos(critical_pressure, "critical_pressure"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_critical_volume_from_zc",
    )


def calc_critical_pressure_from_zc(
    critical_compressibility_factor: ScalarValue,
    critical_temperature: ScalarValue,
    critical_molar_volume: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "critical_pressure",
    description: str = "Calculate critical pressure from critical compressibility factor.",
    unit: str | None = "Pa",
    symbol: str | None = "Pc",
) -> AnnotatedValue[float]:
    """Calculate annotated critical pressure from ``Zc`` in SI-consistent units."""
    return to_annotated_value(
        float(_calc_critical_pressure_from_zc(
            _pos(critical_compressibility_factor, "critical_compressibility_factor"),
            _pos(critical_temperature, "critical_temperature"),
            _pos(critical_molar_volume, "critical_molar_volume"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_critical_pressure_from_zc",
    )


def calc_critical_temperature_from_zc(
    critical_pressure: ScalarValue,
    critical_molar_volume: ScalarValue,
    critical_compressibility_factor: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "critical_temperature",
    description: str = "Calculate critical temperature from critical compressibility factor.",
    unit: str | None = "K",
    symbol: str | None = "Tc",
) -> AnnotatedValue[float]:
    """Calculate annotated critical temperature from ``Zc`` in SI-consistent units."""
    return to_annotated_value(
        float(_calc_critical_temperature_from_zc(
            _pos(critical_pressure, "critical_pressure"),
            _pos(critical_molar_volume, "critical_molar_volume"),
            _pos(critical_compressibility_factor, "critical_compressibility_factor"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_critical_temperature_from_zc",
    )


__all__ = [
    "calc_critical_compressibility_factor",
    "calc_critical_volume_from_zc",
    "calc_critical_pressure_from_zc",
    "calc_critical_temperature_from_zc",
]
