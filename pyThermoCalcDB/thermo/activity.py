"""Activity transformations independent of activity-coefficient models."""

# import libs
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
# locals
from ..utils.conversions import _pos, _resolve_unit_conversion_fn, _scalar
from .core.activity import (
    _calc_activity_from_concentration,
    _calc_activity_from_mole_fraction,
    _calc_activity_coefficient_from_fugacity,
    _calc_effective_concentration,
)


# SECTION: Public helpers

def calc_activity_from_mole_fraction(
    mole_fraction,
    activity_coefficient,
) -> float | NDArray[np.float64]:
    """Calculate dimensionless activity from mole fraction.

    Parameters
    ----------
    mole_fraction : float | sequence | ndarray
        Component mole fraction(s), dimensionless and non-negative.
    activity_coefficient : float | sequence | ndarray
        Supplied activity coefficient(s), dimensionless and positive.

    Returns
    -------
    float | NDArray[np.float64]
        Dimensionless activity.

    Notes
    -----
    Equation: ``a_i = gamma_i*x_i``. This function does not calculate
    ``gamma_i`` from an activity-coefficient model.
    """
    return _calc_activity_from_mole_fraction(mole_fraction, activity_coefficient)


def calc_activity_from_concentration(
    concentration,
    activity_coefficient,
    reference_concentration,
    concentration_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | NDArray[np.float64]:
    """Calculate dimensionless activity from concentration and a reference concentration.

    Numeric concentrations are assumed to already share the same basis. When
    ``CustomProp`` values are supplied, both concentrations are normalized to
    ``concentration_unit`` before calculation.
    """
    # SECTION: Normalize scalar unit-bearing inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    c = _scalar(concentration, "concentration", concentration_unit, conversion_fn)
    gamma = _pos(activity_coefficient, "activity_coefficient")
    c_ref = _pos(reference_concentration, "reference_concentration", concentration_unit, conversion_fn)
    return _calc_activity_from_concentration(c, gamma, c_ref)


def calc_effective_concentration(
    activity,
    reference_concentration,
    concentration_unit: str | None = None,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | CustomProp | NDArray[np.float64]:
    """Calculate effective concentration from activity.

    Parameters
    ----------
    activity : float | int
        Dimensionless activity.
    reference_concentration : float | int | CustomProp
        Reference concentration. Unit-bearing values can be normalized before
        calculation.
    concentration_unit : str, optional
        Unit used to normalize a ``CustomProp`` reference concentration.
    output_unit : str, optional
        Unit annotation for a ``CustomProp`` result. Defaults to
        ``concentration_unit`` when supplied.

    Returns
    -------
    float | CustomProp | NDArray[np.float64]
        Effective concentration ``c_eff = a_i*c_ref``.
    """
    # SECTION: Normalize scalar unit-bearing inputs
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    a = _scalar(activity, "activity")
    c_ref = _pos(reference_concentration, "reference_concentration", concentration_unit, conversion_fn)
    value = _calc_effective_concentration(a, c_ref)
    # NOTE: Preserve unit metadata when the reference concentration carried it.
    if isinstance(reference_concentration, CustomProp) or output_unit is not None:
        return CustomProp(value=float(value), unit=output_unit or concentration_unit or "")
    return value



def calc_activity_coefficient_from_fugacity(
    liquid_fugacity,
    mole_fraction,
    standard_state_fugacity,
) -> float | NDArray[np.float64]:
    """Calculate activity coefficient from fugacity definition.

    Equation: ``gamma_i = f_i^L/(x_i*f_i^0)``. Numeric fugacities are assumed
    to already share a pressure unit basis.
    """
    f_l = _pos(liquid_fugacity, "liquid_fugacity")
    x = _pos(mole_fraction, "mole_fraction")
    f0 = _pos(standard_state_fugacity, "standard_state_fugacity")
    return _calc_activity_coefficient_from_fugacity(f_l, x, f0)
# SECTION: Public exports
__all__ = [
    "calc_activity_from_mole_fraction",
    "calc_activity_from_concentration",
    "calc_effective_concentration",
    "calc_activity_coefficient_from_fugacity",
]


