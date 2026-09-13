"""Phase-change thermodynamic identities."""

# import libs
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pycuc import convert_from_to

# locals
from ..utils.conversions import _pos, _scalar, _to_kelvin
from .core.phase_change import (
    _calc_clapeyron_slope,
    _calc_enthalpy_of_sublimation,
    _calc_enthalpy_vaporization_watson,
    _calc_phase_transition_entropy,
    _calc_transition_enthalpy_from_constant_delta_cp,
    _calc_transition_enthalpy_from_cp_integral,
)


# SECTION: Public wrappers
def calc_phase_transition_entropy(
    transition_enthalpy: ScalarValue,
    transition_temperature: Temperature | ScalarValue,
    output_unit: str = "J/(mol.K)",
    unit_conversion_fn=None,
) -> float:
    """Calculate reversible phase-transition entropy.

    Equation
    --------
    delta_S_tr = delta_H_tr / T_tr

    The enthalpy is normalized to J/mol and temperature to K when unit-aware
    inputs are supplied. The caller is responsible for using a transition
    enthalpy and transition temperature for the same equilibrium condition.
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    h_tr = _scalar(
        transition_enthalpy,
        "transition_enthalpy",
        "J/mol",
        conversion_fn,
    )
    if isinstance(transition_temperature, Temperature):
        t_tr = _to_kelvin(transition_temperature)
    else:
        t_tr = _pos(
            transition_temperature,
            "transition_temperature",
            "K" if isinstance(transition_temperature, CustomProp) else None,
            conversion_fn,
        )
    result = float(_calc_phase_transition_entropy(h_tr, t_tr))
    if output_unit != "J/(mol.K)":
        result = float(conversion_fn(result, "J/(mol.K)", output_unit))
    return result


def calc_enthalpy_of_sublimation(
    enthalpy_fusion: ScalarValue,
    enthalpy_vaporization: ScalarValue,
    output_unit: str = "J/mol",
    unit_conversion_fn=None,
) -> float:
    """Calculate sublimation enthalpy from fusion and vaporization enthalpies.

    Equation
    --------
    delta_H_sub = delta_H_fus + delta_H_vap

    The two enthalpy values must refer to a thermodynamically consistent
    reference path; this helper does not reconcile incompatible temperatures or
    reference states.
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    h_fus = _scalar(enthalpy_fusion, "enthalpy_fusion", "J/mol", conversion_fn)
    h_vap = _scalar(
        enthalpy_vaporization,
        "enthalpy_vaporization",
        "J/mol",
        conversion_fn,
    )
    result = float(_calc_enthalpy_of_sublimation(h_fus, h_vap))
    if output_unit != "J/mol":
        result = float(conversion_fn(result, "J/mol", output_unit))
    return result


def calc_clapeyron_slope(
    transition_enthalpy: ScalarValue,
    temperature: Temperature | ScalarValue,
    delta_molar_volume: ScalarValue,
    output_unit: str = "Pa/K",
    unit_conversion_fn=None,
) -> float:
    """Calculate the equilibrium phase-boundary slope.

    Equation
    --------
    dP/dT = delta_H_tr / (T * delta_V_tr)

    The input enthalpy and molar-volume difference must refer to the same
    transition state. This function does not calculate phase molar volumes.
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    h_tr = _scalar(
        transition_enthalpy,
        "transition_enthalpy",
        "J/mol",
        conversion_fn,
    )
    if isinstance(temperature, Temperature):
        t = _to_kelvin(temperature)
    else:
        t = _pos(
            temperature,
            "temperature",
            "K" if isinstance(temperature, CustomProp) else None,
            conversion_fn,
        )
    dv = _scalar(
        delta_molar_volume,
        "delta_molar_volume",
        "m3/mol",
        conversion_fn,
    )
    result = float(_calc_clapeyron_slope(h_tr, t, dv))
    if output_unit != "Pa/K":
        result = float(conversion_fn(result, "Pa/K", output_unit))
    return result


def calc_enthalpy_vaporization_watson(
    enthalpy_vaporization_reference: ScalarValue,
    temperature_reference: Temperature | ScalarValue,
    temperature: Temperature | ScalarValue,
    critical_temperature: Temperature | ScalarValue,
    exponent: float = 0.38,
    output_unit: str = "J/mol",
    unit_conversion_fn=None,
) -> float:
    """Correct vaporization enthalpy from one temperature to another.

    Equation
    --------
    delta_H_vap(T2) = delta_H_vap(T1) * [(1 - Tr2)/(1 - Tr1)]**n

    The default exponent is the classical Watson value ``0.38``.
    """
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    h_ref = _pos(
        enthalpy_vaporization_reference,
        "enthalpy_vaporization_reference",
        "J/mol",
        conversion_fn,
    )
    t_ref = _to_kelvin(temperature_reference) if isinstance(temperature_reference, Temperature) else _pos(
        temperature_reference,
        "temperature_reference",
        "K" if isinstance(temperature_reference, CustomProp) else None,
        conversion_fn,
    )
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    tc = _to_kelvin(critical_temperature) if isinstance(critical_temperature, Temperature) else _pos(
        critical_temperature,
        "critical_temperature",
        "K" if isinstance(critical_temperature, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_enthalpy_vaporization_watson(h_ref, t_ref, t, tc, exponent))
    if output_unit != "J/mol":
        result = float(conversion_fn(result, "J/mol", output_unit))
    return result


def calc_transition_enthalpy_from_constant_delta_cp(
    transition_enthalpy_reference: ScalarValue,
    temperature_reference: Temperature | ScalarValue,
    temperature: Temperature | ScalarValue,
    delta_heat_capacity: ScalarValue,
    output_unit: str = "J/mol",
    unit_conversion_fn=None,
) -> float:
    """Apply Kirchhoff transition-enthalpy correction with constant ``delta_Cp``."""
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    h_ref = _scalar(
        transition_enthalpy_reference,
        "transition_enthalpy_reference",
        "J/mol",
        conversion_fn,
    )
    t_ref = _to_kelvin(temperature_reference) if isinstance(temperature_reference, Temperature) else _pos(
        temperature_reference,
        "temperature_reference",
        "K" if isinstance(temperature_reference, CustomProp) else None,
        conversion_fn,
    )
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    delta_cp = _scalar(
        delta_heat_capacity,
        "delta_heat_capacity",
        "J/(mol.K)" if isinstance(delta_heat_capacity, CustomProp) else None,
        conversion_fn,
    )
    result = float(_calc_transition_enthalpy_from_constant_delta_cp(h_ref, t_ref, t, delta_cp))
    if output_unit != "J/mol":
        result = float(conversion_fn(result, "J/mol", output_unit))
    return result


def calc_transition_enthalpy_from_cp_integral(
    transition_enthalpy_reference: ScalarValue,
    delta_cp_integral: ScalarValue,
    output_unit: str = "J/mol",
    unit_conversion_fn=None,
) -> float:
    """Apply Kirchhoff correction from supplied ``integral(delta_Cp dT)``."""
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn
    h_ref = _scalar(
        transition_enthalpy_reference,
        "transition_enthalpy_reference",
        "J/mol",
        conversion_fn,
    )
    integral = _scalar(
        delta_cp_integral,
        "delta_cp_integral",
        "J/mol",
        conversion_fn,
    )
    result = float(_calc_transition_enthalpy_from_cp_integral(h_ref, integral))
    if output_unit != "J/mol":
        result = float(conversion_fn(result, "J/mol", output_unit))
    return result


__all__ = [
    "calc_phase_transition_entropy",
    "calc_enthalpy_of_sublimation",
    "calc_clapeyron_slope",
    "calc_enthalpy_vaporization_watson",
    "calc_transition_enthalpy_from_constant_delta_cp",
    "calc_transition_enthalpy_from_cp_integral",
]
