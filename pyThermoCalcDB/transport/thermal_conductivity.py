"""Public thermal-conductivity correlations and mixing rules."""

# import libs
from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_list
from pythermodb_settings.utils.validators import fractions, positive, same_shape
# locals
from ..configs.constants import R_J_molK
from ..utils.conversions import _pos, _resolve_unit_conversion_fn
from .core.thermal_conductivity import (
    _calc_ideal_liquid_thermal_conductivity,
    _calc_stiel_thodos_gas_thermal_conductivity,
)


# SECTION: Public thermal conductivity calculations

def calc_stiel_thodos_gas_thermal_conductivity(
    viscosity,
    molecular_weight,
    molar_heat_capacity,
    gas_constant: float = R_J_molK,
    output_unit: str = "W/m/K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float | CustomProp | NDArray[np.float64]:
    """Calculate Stiel-Thodos low-pressure gas thermal conductivity."""
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    mu = _pos(viscosity, "viscosity", "Pa*s", conversion_fn)
    mw = _pos(molecular_weight, "molecular_weight", "kg/mol", conversion_fn)
    # FIXME: add to pycuc
    cp = _pos(molar_heat_capacity, "molar_heat_capacity",
              "J/mol/K", conversion_fn)
    r = _pos(gas_constant, "gas_constant")
    value = _calc_stiel_thodos_gas_thermal_conductivity(mu, mw, cp, r)
    if output_unit != "W/m/K":
        value = conversion_fn(value=float(
            value), from_unit="W/m/K", to_unit=output_unit)
    if any(isinstance(item, CustomProp) for item in (viscosity, molecular_weight, molar_heat_capacity)) or output_unit != "W/m/K":
        return CustomProp(value=float(value), unit=output_unit)
    return value


def calc_ideal_liquid_thermal_conductivity(
    mole_fractions: Sequence[float | int | CustomProp],
    thermal_conductivities: Sequence[float | int | CustomProp],
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal liquid thermal conductivity by mole-fraction averaging."""
    # SECTION: Validate and normalize sequence inputs
    fractions(mole_fractions, "mole_fractions")
    positive(thermal_conductivities, "thermal_conductivities")
    same_shape(mole_fractions, thermal_conductivities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_list(mole_fractions, unit_conversion_fn=conversion_fn)
    k = to_list(thermal_conductivities, output_unit,
                unit_conversion_fn=conversion_fn)
    return float(_calc_ideal_liquid_thermal_conductivity(x, k))


# SECTION: Public exports
__all__ = [
    "calc_stiel_thodos_gas_thermal_conductivity",
    "calc_ideal_liquid_thermal_conductivity",
]
