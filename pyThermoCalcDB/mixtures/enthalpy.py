"""Ideal mixture enthalpy helpers."""

# import libs
from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

# locals
from .core.enthalpy import _calc_ideal_enthalpy_of_mixing


def calc_ideal_enthalpy_of_mixing(
    mole_fractions: Sequence[float | int] | NDArray[np.number] | None = None,
) -> float | NDArray[np.float64]:
    """Return the ideal-solution enthalpy of mixing.

    For an ideal solution, ``delta_H_mix^ideal = 0``. The optional composition
    argument exists for API symmetry and validation only; no dissociation or
    non-ideal interaction contribution is inferred.
    """
    return _calc_ideal_enthalpy_of_mixing(mole_fractions)


calc_ideal_molar_enthalpy_of_mixing = calc_ideal_enthalpy_of_mixing


__all__ = [
    "calc_ideal_enthalpy_of_mixing",
    "calc_ideal_molar_enthalpy_of_mixing",
]
