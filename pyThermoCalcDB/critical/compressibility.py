"""Critical compressibility public wrappers."""

# import libs
from pythermodb_settings.models import ScalarValue

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
) -> float:
    """Calculate ``Zc = Pc*Vc/(R*Tc)`` from critical properties.

    Inputs are positive finite scalars in consistent SI units. The output is a
    dimensionless exact critical-property definition.
    """
    return float(_calc_critical_compressibility_factor(
        _pos(critical_pressure, "critical_pressure"),
        _pos(critical_molar_volume, "critical_molar_volume"),
        _pos(critical_temperature, "critical_temperature"),
        _pos(gas_constant, "gas_constant"),
    ))


def calc_critical_volume_from_zc(
    critical_compressibility_factor: ScalarValue,
    critical_temperature: ScalarValue,
    critical_pressure: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> float:
    """Calculate critical molar volume from ``Zc`` in SI-consistent units."""
    return float(_calc_critical_volume_from_zc(
        _pos(critical_compressibility_factor, "critical_compressibility_factor"),
        _pos(critical_temperature, "critical_temperature"),
        _pos(critical_pressure, "critical_pressure"),
        _pos(gas_constant, "gas_constant"),
    ))


def calc_critical_pressure_from_zc(
    critical_compressibility_factor: ScalarValue,
    critical_temperature: ScalarValue,
    critical_molar_volume: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> float:
    """Calculate critical pressure from ``Zc`` in SI-consistent units."""
    return float(_calc_critical_pressure_from_zc(
        _pos(critical_compressibility_factor, "critical_compressibility_factor"),
        _pos(critical_temperature, "critical_temperature"),
        _pos(critical_molar_volume, "critical_molar_volume"),
        _pos(gas_constant, "gas_constant"),
    ))


def calc_critical_temperature_from_zc(
    critical_pressure: ScalarValue,
    critical_molar_volume: ScalarValue,
    critical_compressibility_factor: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> float:
    """Calculate critical temperature from ``Zc`` in SI-consistent units."""
    return float(_calc_critical_temperature_from_zc(
        _pos(critical_pressure, "critical_pressure"),
        _pos(critical_molar_volume, "critical_molar_volume"),
        _pos(critical_compressibility_factor, "critical_compressibility_factor"),
        _pos(gas_constant, "gas_constant"),
    ))


__all__ = [
    "calc_critical_compressibility_factor",
    "calc_critical_volume_from_zc",
    "calc_critical_pressure_from_zc",
    "calc_critical_temperature_from_zc",
]
