"""Excess-property and Gibbs-Duhem mixture identities."""

# import libs
from collections.abc import Sequence

from pythermodb_settings.models import CustomProp, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_list
from pythermodb_settings.utils.validators import fractions, positive, same_shape
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import _pos, _resolve_unit_conversion_fn, _scalar, _to_kelvin
from .core.excess import (
    _calc_excess_gibbs_energy_from_activity_coefficients,
    _calc_excess_property,
    _calc_excess_entropy_from_gibbs_enthalpy,
    _calc_gibbs_duhem_residual,
    _check_gibbs_duhem_consistency,
)


# SECTION: Public excess-property identities

def calc_excess_property(
    real_property,
    ideal_property,
) -> float:
    """Calculate excess property ``M^E = M - M_ideal``."""
    return float(_calc_excess_property(real_property, ideal_property))


def calc_excess_gibbs_energy_from_activity_coefficients(
    mole_fractions: Sequence[float | int | CustomProp],
    activity_coefficients: Sequence[float | int | CustomProp],
    temperature,
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate molar excess Gibbs energy from activity coefficients."""
    fractions(mole_fractions, "mole_fractions")
    positive(activity_coefficients, "activity_coefficients")
    same_shape(mole_fractions, activity_coefficients)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    gamma = to_list(activity_coefficients, unit_conversion_fn=conversion_fn)
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    r = _pos(gas_constant, "gas_constant")
    return float(_calc_excess_gibbs_energy_from_activity_coefficients(x, gamma, t, r))


def calc_excess_entropy_from_gibbs_enthalpy(
    excess_gibbs_energy,
    excess_enthalpy,
    temperature,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate excess entropy from excess Gibbs energy and excess enthalpy."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    g_e = _scalar(
        excess_gibbs_energy,
        "excess_gibbs_energy",
        "J/mol" if isinstance(excess_gibbs_energy, CustomProp) else None,
        conversion_fn,
    )
    h_e = _scalar(
        excess_enthalpy,
        "excess_enthalpy",
        "J/mol" if isinstance(excess_enthalpy, CustomProp) else None,
        conversion_fn,
    )
    t = _to_kelvin(temperature) if isinstance(temperature, Temperature) else _pos(
        temperature,
        "temperature",
        "K" if isinstance(temperature, CustomProp) else None,
        conversion_fn,
    )
    return float(_calc_excess_entropy_from_gibbs_enthalpy(g_e, h_e, t))


# SECTION: Public Gibbs-Duhem helpers

def calc_gibbs_duhem_residual(
    mole_fractions: Sequence[float | int | CustomProp],
    dlog_activity_coefficients: Sequence[float | int | CustomProp],
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate Gibbs-Duhem residual ``sum_i x_i*dln(gamma_i)``."""
    fractions(mole_fractions, "mole_fractions")
    same_shape(mole_fractions, dlog_activity_coefficients)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    dln_gamma = to_list(dlog_activity_coefficients, unit_conversion_fn=conversion_fn)
    return float(_calc_gibbs_duhem_residual(x, dln_gamma))


def check_gibbs_duhem_consistency(
    mole_fractions: Sequence[float | int | CustomProp],
    dlog_activity_coefficients: Sequence[float | int | CustomProp],
    tolerance: float = 1.0e-8,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> bool:
    """Return whether the Gibbs-Duhem residual is within tolerance."""
    residual = calc_gibbs_duhem_residual(
        mole_fractions,
        dlog_activity_coefficients,
        unit_conversion_fn,
    )
    return bool(_check_gibbs_duhem_consistency([1.0], [residual], tolerance))


# SECTION: Public exports
__all__ = [
    "calc_excess_property",
    "calc_excess_gibbs_energy_from_activity_coefficients",
    "calc_excess_entropy_from_gibbs_enthalpy",
    "calc_gibbs_duhem_residual",
    "check_gibbs_duhem_consistency",
]
