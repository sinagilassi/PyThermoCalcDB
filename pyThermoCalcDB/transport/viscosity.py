"""Viscosity correlations and mixing-rule public wrappers."""

# import libs
import numpy as np
from numpy.typing import NDArray

# locals
from .core.viscosity import (
    _calc_liquid_mixture_viscosity_log_rule,
    _calc_viscosity_exponential_correlation,
)


# SECTION: Public wrappers
def calc_liquid_mixture_viscosity_log_rule(
    mole_fractions,
    component_viscosities,
) -> float | NDArray[np.float64]:
    """Calculate liquid-mixture viscosity with ``eta = exp(sum(x_i ln eta_i))``."""
    return _calc_liquid_mixture_viscosity_log_rule(mole_fractions, component_viscosities)


def calc_viscosity_exponential_correlation(
    temperature,
    A,
    B,
    C,
    D,
    E,
) -> float | NDArray[np.float64]:
    """Evaluate ``eta = exp(A + B/T + C*ln(T) + D*T**E)`` with supplied coefficients."""
    return _calc_viscosity_exponential_correlation(temperature, A, B, C, D, E)


__all__ = [
    "calc_liquid_mixture_viscosity_log_rule",
    "calc_viscosity_exponential_correlation",
]
