"""Heat-capacity correlations and analytical integrals."""

# locals
from .polynomial import (
    calc_enthalpy_change_from_cp_polynomial,
    calc_entropy_change_from_cp_polynomial,
    calc_heat_capacity_polynomial,
)

__all__ = [
    "calc_heat_capacity_polynomial",
    "calc_enthalpy_change_from_cp_polynomial",
    "calc_entropy_change_from_cp_polynomial",
]
