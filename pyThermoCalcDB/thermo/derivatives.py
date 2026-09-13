"""Thermodynamic derivative-property identities."""

# import libs
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pycuc import convert_from_to

# locals
from ..utils.conversions import _pos, _scalar, _to_kelvin
from .core.derivatives import (
    _calc_isothermal_compressibility,
    _calc_joule_thomson_coefficient,
    _calc_thermal_expansion_coefficient,
)


# SECTION: Public wrappers
def calc_thermal_expansion_coefficient(
    volume: ScalarValue,
    dvolume_dtemperature_at_pressure: ScalarValue,
    output_unit: str = "1/K",
    unit_conversion_fn=None,
) -> float:
    """Calculate the thermal expansion coefficient from supplied derivative data.

    Equation
    --------
    alpha = (1/V) * (dV/dT)_P

    This function consumes an already-known derivative; it does not derive the
    derivative from an equation of state.
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    v = _pos(volume, "volume", "m3/mol" if isinstance(volume, CustomProp) else None, conversion_fn)
    dvd_t = _scalar(
        dvolume_dtemperature_at_pressure,
        "dvolume_dtemperature_at_pressure",
        "m3/(mol.K)" if isinstance(dvolume_dtemperature_at_pressure, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_thermal_expansion_coefficient(v, dvd_t))
    if output_unit != "1/K":
        result = float(conversion_fn(result, "1/K", output_unit))
    return result


def calc_isothermal_compressibility(
    volume: ScalarValue,
    dvolume_dpressure_at_temperature: ScalarValue,
    output_unit: str = "1/Pa",
    unit_conversion_fn=None,
) -> float:
    """Calculate isothermal compressibility from supplied derivative data.

    Equation
    --------
    kappa_T = -(1/V) * (dV/dP)_T

    This function does not calculate the pressure derivative from an EOS.
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    v = _pos(volume, "volume", "m3/mol" if isinstance(volume, CustomProp) else None, conversion_fn)
    dvd_p = _scalar(
        dvolume_dpressure_at_temperature,
        "dvolume_dpressure_at_temperature",
        "m3/(mol.Pa)" if isinstance(dvolume_dpressure_at_temperature, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_isothermal_compressibility(v, dvd_p))
    if output_unit != "1/Pa":
        result = float(conversion_fn(result, "1/Pa", output_unit))
    return result


def calc_joule_thomson_coefficient(
    temperature: Temperature | ScalarValue,
    molar_volume: ScalarValue,
    heat_capacity_cp: ScalarValue,
    dvolume_dtemperature_at_pressure: ScalarValue,
    output_unit: str = "K/Pa",
    unit_conversion_fn=None,
) -> float:
    """Calculate the Joule-Thomson coefficient from supplied property data.

    Equation
    --------
    mu_JT = [T * (dV/dT)_P - V] / Cp

    The required Cp and volume derivative must be provided by the caller or a
    higher-level model; this function does not fit or evaluate an EOS.
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    if isinstance(temperature, Temperature):
        t = _to_kelvin(temperature)
    else:
        t = _pos(
            temperature,
            "temperature",
            "K" if isinstance(temperature, CustomProp) else None,
            conversion_fn,
        )
    v = _pos(
        molar_volume,
        "molar_volume",
        "m3/mol" if isinstance(molar_volume, CustomProp) else None,
        conversion_fn,
    )
    cp = _scalar(
        heat_capacity_cp,
        "heat_capacity_cp",
        "J/(mol.K)" if isinstance(heat_capacity_cp, CustomProp) else None,
        conversion_fn,
    )
    dvd_t = _scalar(
        dvolume_dtemperature_at_pressure,
        "dvolume_dtemperature_at_pressure",
        "m3/(mol.K)" if isinstance(dvolume_dtemperature_at_pressure, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_joule_thomson_coefficient(t, v, cp, dvd_t))
    if output_unit != "K/Pa":
        result = float(conversion_fn(result, "K/Pa", output_unit))
    return result


__all__ = [
    "calc_thermal_expansion_coefficient",
    "calc_isothermal_compressibility",
    "calc_joule_thomson_coefficient",
]
