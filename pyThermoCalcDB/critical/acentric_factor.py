"""Acentric-factor public wrappers."""

# import libs
from pythermodb_settings.models import ScalarValue

# locals
from ..utils.conversions import _pos, _scalar
from .core.acentric_factor import (
    _calc_acentric_factor_from_reduced_vapor_pressure,
    _calc_acentric_factor_from_vapor_pressure,
    _calc_reduced_vapor_pressure_from_acentric_factor,
)


# SECTION: Public wrappers
def calc_acentric_factor_from_reduced_vapor_pressure(
    reduced_vapor_pressure: ScalarValue,
) -> float:
    """Calculate Pitzer acentric factor from ``Pr_sat`` at ``Tr = 0.7``."""
    return float(_calc_acentric_factor_from_reduced_vapor_pressure(
        _pos(reduced_vapor_pressure, "reduced_vapor_pressure")
    ))


def calc_acentric_factor_from_vapor_pressure(
    saturation_pressure: ScalarValue,
    critical_pressure: ScalarValue,
) -> float:
    """Calculate Pitzer acentric factor from ``Psat/Pc`` at ``Tr = 0.7``."""
    return float(_calc_acentric_factor_from_vapor_pressure(
        _pos(saturation_pressure, "saturation_pressure"),
        _pos(critical_pressure, "critical_pressure"),
    ))


def calc_reduced_vapor_pressure_from_acentric_factor(
    acentric_factor: ScalarValue,
) -> float:
    """Calculate ``Pr_sat = 10**(-omega - 1)`` from acentric factor."""
    return float(_calc_reduced_vapor_pressure_from_acentric_factor(
        _scalar(acentric_factor, "acentric_factor")
    ))


__all__ = [
    "calc_acentric_factor_from_reduced_vapor_pressure",
    "calc_acentric_factor_from_vapor_pressure",
    "calc_reduced_vapor_pressure_from_acentric_factor",
]
