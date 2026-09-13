"""Public vapor-liquid surface-tension mixing rules."""

# import libs
from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_list
from pythermodb_settings.utils.validators import fractions, positive, same_shape
# locals
from ..utils.conversions import _resolve_unit_conversion_fn
from .core.surface_tension import (
    _calc_ideal_vapor_liquid_surface_tension,
    _calc_winterfeld_vapor_liquid_surface_tension,
)


# SECTION: Public surface-tension mixing rules

def calc_ideal_vapor_liquid_surface_tension(
    liquid_mole_fractions: Sequence[float | int | CustomProp],
    pure_surface_tensions: Sequence[float | int | CustomProp],
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal vapor-liquid surface tension by liquid mole-fraction averaging."""
    fractions(liquid_mole_fractions, "liquid_mole_fractions")
    positive(pure_surface_tensions, "pure_surface_tensions")
    same_shape(liquid_mole_fractions, pure_surface_tensions)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_list(liquid_mole_fractions, unit_conversion_fn=conversion_fn)
    sigma = to_list(pure_surface_tensions, output_unit, unit_conversion_fn=conversion_fn)
    return float(_calc_ideal_vapor_liquid_surface_tension(x, sigma))


def calc_winterfeld_vapor_liquid_surface_tension(
    liquid_mole_fractions: Sequence[float | int | CustomProp],
    pure_surface_tensions: Sequence[float | int | CustomProp],
    liquid_molar_densities: Sequence[float | int | CustomProp],
    output_surface_tension_unit: str | None = None,
    output_molar_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate Winterfeld-Scriven-Davis vapor-liquid surface tension.

    The density input is pure-component liquid molar density, not mass density.
    """
    fractions(liquid_mole_fractions, "liquid_mole_fractions")
    positive(pure_surface_tensions, "pure_surface_tensions")
    positive(liquid_molar_densities, "liquid_molar_densities")
    same_shape(liquid_mole_fractions, pure_surface_tensions)
    same_shape(liquid_mole_fractions, liquid_molar_densities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    x = to_list(liquid_mole_fractions, unit_conversion_fn=conversion_fn)
    sigma = to_list(pure_surface_tensions, output_surface_tension_unit, unit_conversion_fn=conversion_fn)
    rho = to_list(liquid_molar_densities, output_molar_density_unit, unit_conversion_fn=conversion_fn)
    return float(_calc_winterfeld_vapor_liquid_surface_tension(x, sigma, rho))


# SECTION: Public exports
__all__ = [
    "calc_ideal_vapor_liquid_surface_tension",
    "calc_winterfeld_vapor_liquid_surface_tension",
]
