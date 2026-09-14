"""Core Rackett saturated-liquid density helpers."""

# import libs
from typing import cast

import numpy as np
from numpy.typing import NDArray

# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import NumericArrayInput, _return_scalar_if_zero_dim
from .definitions import _as_finite_array, _as_positive_array

# SECTION: Type aliases
NumericInput = NumericArrayInput


# SECTION: Core numeric calculations
def _calc_saturated_liquid_molar_volume_rackett(
    temperature: NumericInput,
    critical_temperature: NumericInput,
    critical_pressure: NumericInput,
    rackett_constant: NumericInput,
    gas_constant: NumericInput = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate saturated-liquid molar volume with the modified Rackett equation.

    Equation: ``Vs = (R*Tc/Pc) * Z_RA**(1 + (1 - Tr)**(2/7))`` where
    ``Tr = T/Tc``. Inputs use K, Pa, dimensionless Rackett constant, and
    J/(mol.K), producing m3/mol. This empirical general correlation is used
    below the critical temperature.
    """
    t = _as_positive_array(temperature, "temperature")
    tc = _as_positive_array(critical_temperature, "critical_temperature")
    pc = _as_positive_array(critical_pressure, "critical_pressure")
    zra = _as_positive_array(rackett_constant, "rackett_constant")
    r = _as_positive_array(gas_constant, "gas_constant")
    # ! The saturated-liquid Rackett form is intended for subcritical temperatures.
    if np.any(t >= tc):
        raise ValueError("temperature must be less than critical_temperature.")
    try:
        tr = t / tc
        exponent = 1.0 + np.power(1.0 - tr, 2.0 / 7.0)
        volume = (r * tc / pc) * np.power(zra, exponent)
    except ValueError as exc:
        raise ValueError("Rackett inputs must be broadcast-compatible.") from exc
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], volume))


def _calc_rackett_constant_from_acentric_factor(
    acentric_factor: NumericInput,
) -> float | NDArray[np.float64]:
    """Estimate Rackett constant from acentric factor.

    Equation: ``Z_RA = 0.2956 - 0.08775*omega``. The acentric factor and output
    are dimensionless. This is an empirical estimate, not an identity.
    """
    omega = _as_finite_array(acentric_factor, "acentric_factor")
    zra = 0.2956 - 0.08775 * omega
    if np.any(zra <= 0.0):
        raise ValueError("estimated rackett_constant must be greater than zero.")
    return _return_scalar_if_zero_dim(cast(NDArray[np.float64], zra))


__all__ = [
    "_calc_saturated_liquid_molar_volume_rackett",
    "_calc_rackett_constant_from_acentric_factor",
]
