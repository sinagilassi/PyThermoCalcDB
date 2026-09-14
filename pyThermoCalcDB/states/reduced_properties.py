"""Reduced thermodynamic property public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, CustomProp, ScalarValue, Temperature
from pythermodb_settings.utils import to_annotated_value

# locals
from ..utils.conversions import _pos, _scalar, _to_kelvin
from .core.reduced_properties import (
    _calc_reduced_pressure,
    _calc_reduced_temperature,
    _calc_reduced_volume,
)


# SECTION: Public wrappers
def calc_reduced_temperature(
    temperature: Temperature | ScalarValue,
    critical_temperature: Temperature | ScalarValue,
    unit_conversion_fn=None,
    *,
    name: str = "reduced_temperature",
    description: str = "Calculate reduced temperature.",
    symbol: str | None = "Tr",
) -> AnnotatedValue[float]:
    """Calculate annotated reduced temperature ``Tr = T / Tc``.

    Temperatures are normalized to K for ``Temperature`` or unit-aware
    ``CustomProp`` inputs. The output is dimensionless and exact.
    """
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _scalar(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        unit_conversion_fn,
    )
    tc = _to_kelvin(critical_temperature) if isinstance(critical_temperature, Temperature) else _pos(
        critical_temperature,
        "critical_temperature",
        "K" if isinstance(critical_temperature, CustomProp) else None,
        unit_conversion_fn,
    )
    return to_annotated_value(
        float(_calc_reduced_temperature(t, tc)),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_reduced_temperature",
    )


def calc_reduced_pressure(
    pressure: ScalarValue,
    critical_pressure: ScalarValue,
    output_pressure_unit: str | None = None,
    unit_conversion_fn=None,
    *,
    name: str = "reduced_pressure",
    description: str = "Calculate reduced pressure.",
    symbol: str | None = "Pr",
) -> AnnotatedValue[float]:
    """Calculate annotated reduced pressure ``Pr = P / Pc``.

    Unit-aware inputs are normalized to ``output_pressure_unit`` when provided.
    Numeric inputs are assumed to already share a pressure unit basis.
    """
    p = _scalar(pressure, "pressure", output_pressure_unit, unit_conversion_fn)
    pc = _pos(critical_pressure, "critical_pressure", output_pressure_unit, unit_conversion_fn)
    return to_annotated_value(
        float(_calc_reduced_pressure(p, pc)),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_reduced_pressure",
    )


def calc_reduced_volume(
    molar_volume: ScalarValue,
    critical_molar_volume: ScalarValue,
    output_volume_unit: str | None = None,
    unit_conversion_fn=None,
    *,
    name: str = "reduced_volume",
    description: str = "Calculate reduced molar volume.",
    symbol: str | None = "Vr",
) -> AnnotatedValue[float]:
    """Calculate annotated reduced molar volume ``Vr = Vm / Vc``.

    Unit-aware inputs are normalized to ``output_volume_unit`` when provided.
    Numeric inputs are assumed to already share a molar-volume unit basis.
    """
    v = _scalar(molar_volume, "molar_volume", output_volume_unit, unit_conversion_fn)
    vc = _pos(critical_molar_volume, "critical_molar_volume", output_volume_unit, unit_conversion_fn)
    return to_annotated_value(
        float(_calc_reduced_volume(v, vc)),
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
        implementation="_calc_reduced_volume",
    )


__all__ = [
    "calc_reduced_temperature",
    "calc_reduced_pressure",
    "calc_reduced_volume",
]
