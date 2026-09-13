"""Thermodynamic derivative-property identities."""

# import libs
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pycuc import convert_from_to

# locals
from ..utils.conversions import _pos, _scalar, _to_kelvin
from .core.derivatives import (
    _calc_isothermal_compressibility,
    _calc_isothermal_compressibility_from_density,
    _calc_isentropic_compressibility,
    _calc_joule_thomson_coefficient,
    _calc_joule_thomson_coefficient_from_alpha,
    _calc_speed_of_sound,
    _calc_speed_of_sound_from_isentropic_compressibility,
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


def calc_isothermal_compressibility_from_density(
    density: ScalarValue,
    ddensity_dpressure_at_temperature: ScalarValue,
    output_unit: str = "1/Pa",
    unit_conversion_fn=None,
) -> float:
    """Calculate isothermal compressibility from density derivative data.

    Equation
    --------
    kappa_T = (1/rho) * (drho/dP)_T
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    rho = _pos(density, "density", "kg/m3" if isinstance(density, CustomProp) else None, conversion_fn)
    drho_dp = _scalar(
        ddensity_dpressure_at_temperature,
        "ddensity_dpressure_at_temperature",
        "kg/(m3.Pa)" if isinstance(ddensity_dpressure_at_temperature, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_isothermal_compressibility_from_density(rho, drho_dp))
    if output_unit != "1/Pa":
        result = float(conversion_fn(result, "1/Pa", output_unit))
    return result


def calc_isentropic_compressibility(
    isothermal_compressibility: ScalarValue,
    heat_capacity_cv: ScalarValue,
    heat_capacity_cp: ScalarValue,
    output_unit: str = "1/Pa",
    unit_conversion_fn=None,
) -> float:
    """Calculate isentropic compressibility from ``kappa_T``, ``Cv``, and ``Cp``.

    Equation
    --------
    kappa_S = kappa_T * Cv/Cp
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    kappa_t = _pos(
        isothermal_compressibility,
        "isothermal_compressibility",
        "1/Pa" if isinstance(isothermal_compressibility, CustomProp) else None,
        conversion_fn,
    )
    cv = _pos(
        heat_capacity_cv,
        "heat_capacity_cv",
        "J/(mol.K)" if isinstance(heat_capacity_cv, CustomProp) else None,
        conversion_fn,
    )
    cp = _pos(
        heat_capacity_cp,
        "heat_capacity_cp",
        "J/(mol.K)" if isinstance(heat_capacity_cp, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_isentropic_compressibility(kappa_t, cv, cp))
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


def calc_joule_thomson_coefficient_from_alpha(
    temperature: Temperature | ScalarValue,
    molar_volume: ScalarValue,
    heat_capacity_cp: ScalarValue,
    thermal_expansion_coefficient: ScalarValue,
    output_unit: str = "K/Pa",
    unit_conversion_fn=None,
) -> float:
    """Calculate Joule-Thomson coefficient from thermal expansion coefficient.

    Equation
    --------
    mu_JT = V*(alpha*T - 1)/Cp
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
    cp = _pos(
        heat_capacity_cp,
        "heat_capacity_cp",
        "J/(mol.K)" if isinstance(heat_capacity_cp, CustomProp) else None,
        conversion_fn,
    )
    alpha = _scalar(
        thermal_expansion_coefficient,
        "thermal_expansion_coefficient",
        "1/K" if isinstance(thermal_expansion_coefficient, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_joule_thomson_coefficient_from_alpha(t, v, cp, alpha))
    if output_unit != "K/Pa":
        result = float(conversion_fn(result, "K/Pa", output_unit))
    return result


def calc_speed_of_sound_from_isentropic_compressibility(
    density: ScalarValue,
    isentropic_compressibility: ScalarValue,
    output_unit: str = "m/s",
    unit_conversion_fn=None,
) -> float:
    """Calculate speed of sound from density and isentropic compressibility."""
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    rho = _pos(density, "density", "kg/m3" if isinstance(density, CustomProp) else None, conversion_fn)
    kappa_s = _pos(
        isentropic_compressibility,
        "isentropic_compressibility",
        "1/Pa" if isinstance(isentropic_compressibility, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_speed_of_sound_from_isentropic_compressibility(rho, kappa_s))
    if output_unit != "m/s":
        result = float(conversion_fn(result, "m/s", output_unit))
    return result


def calc_speed_of_sound(
    density: ScalarValue,
    isothermal_compressibility: ScalarValue,
    heat_capacity_cp: ScalarValue,
    heat_capacity_cv: ScalarValue,
    output_unit: str = "m/s",
    unit_conversion_fn=None,
) -> float:
    """Calculate speed of sound from ``rho``, ``kappa_T``, ``Cp``, and ``Cv``."""
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    rho = _pos(density, "density", "kg/m3" if isinstance(density, CustomProp) else None, conversion_fn)
    kappa_t = _pos(
        isothermal_compressibility,
        "isothermal_compressibility",
        "1/Pa" if isinstance(isothermal_compressibility, CustomProp) else None,
        conversion_fn,
    )
    cp = _pos(
        heat_capacity_cp,
        "heat_capacity_cp",
        "J/(mol.K)" if isinstance(heat_capacity_cp, CustomProp) else None,
        conversion_fn,
    )
    cv = _pos(
        heat_capacity_cv,
        "heat_capacity_cv",
        "J/(mol.K)" if isinstance(heat_capacity_cv, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_speed_of_sound(rho, kappa_t, cp, cv))
    if output_unit != "m/s":
        result = float(conversion_fn(result, "m/s", output_unit))
    return result


__all__ = [
    "calc_thermal_expansion_coefficient",
    "calc_isothermal_compressibility",
    "calc_isothermal_compressibility_from_density",
    "calc_isentropic_compressibility",
    "calc_joule_thomson_coefficient",
    "calc_joule_thomson_coefficient_from_alpha",
    "calc_speed_of_sound_from_isentropic_compressibility",
    "calc_speed_of_sound",
]
