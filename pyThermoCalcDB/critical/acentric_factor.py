"""Acentric-factor public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, ScalarValue
from pythermodb_settings.utils import to_annotated_value

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
    *,
    name: str = "acentric_factor",
    description: str = "Calculate Pitzer acentric factor from reduced vapor pressure.",
    symbol: str | None = "omega",
) -> AnnotatedValue[float]:
    """Calculate annotated Pitzer acentric factor from ``Pr_sat`` at ``Tr = 0.7``."""
    return to_annotated_value(
        float(_calc_acentric_factor_from_reduced_vapor_pressure(
            _pos(reduced_vapor_pressure, "reduced_vapor_pressure")
        )),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_acentric_factor_from_reduced_vapor_pressure",
    )


def calc_acentric_factor_from_vapor_pressure(
    saturation_pressure: ScalarValue,
    critical_pressure: ScalarValue,
    *,
    name: str = "acentric_factor",
    description: str = "Calculate Pitzer acentric factor from saturation and critical pressure.",
    symbol: str | None = "omega",
) -> AnnotatedValue[float]:
    """Calculate annotated Pitzer acentric factor from ``Psat/Pc`` at ``Tr = 0.7``."""
    return to_annotated_value(
        float(_calc_acentric_factor_from_vapor_pressure(
            _pos(saturation_pressure, "saturation_pressure"),
            _pos(critical_pressure, "critical_pressure"),
        )),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_acentric_factor_from_vapor_pressure",
    )


def calc_reduced_vapor_pressure_from_acentric_factor(
    acentric_factor: ScalarValue,
    *,
    name: str = "reduced_vapor_pressure",
    description: str = "Calculate reduced vapor pressure from Pitzer acentric factor.",
    symbol: str | None = "Pr_sat",
) -> AnnotatedValue[float]:
    """Calculate annotated ``Pr_sat = 10**(-omega - 1)`` from acentric factor."""
    return to_annotated_value(
        float(_calc_reduced_vapor_pressure_from_acentric_factor(
            _scalar(acentric_factor, "acentric_factor")
        )),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_reduced_vapor_pressure_from_acentric_factor",
    )


__all__ = [
    "calc_acentric_factor_from_reduced_vapor_pressure",
    "calc_acentric_factor_from_vapor_pressure",
    "calc_reduced_vapor_pressure_from_acentric_factor",
]
