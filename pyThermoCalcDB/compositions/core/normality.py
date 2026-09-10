"""Normality and equivalent-concentration helpers."""

from collections.abc import Sequence
from typing import Literal, overload, cast

import numpy as np
from numpy.typing import NDArray
# >> pythermodb-settings
from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_custom_prop_scalar

# locals
from ...utils.conversions import _resolve_unit_conversion_fn, _to_units

# ======================================================================
# *** Helper functions
# ======================================================================


def _validate_positive_values(
    values: NDArray[np.float64],
    name: str,
) -> None:
    """Validate finite positive numeric scalar or array values."""
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain finite values.")

    if np.any(values <= 0):
        raise ValueError(f"{name} must be greater than zero.")


def _molarity_unit_from_normality_unit(
    output_unit: str,
) -> str:
    """Return the molarity unit matching a normality unit denominator."""
    units_ = _to_units(output_unit)
    return f"mol/{units_[1]}"


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================
@overload
def _calc_normality(
    molarity: float | int,
    equivalence_factor: float | int,
    *,
    as_list: Literal[False] = False,
) -> float:
    ...


@overload
def _calc_normality(
    molarity: Sequence[float | int] | NDArray[np.number],
    equivalence_factor: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    as_list: Literal[False] = False,
) -> NDArray[np.float64]:
    ...


@overload
def _calc_normality(
    molarity: Sequence[float | int],
    equivalence_factor: float | int | Sequence[float | int],
    *,
    as_list: Literal[True],
) -> list[float]:
    ...


def _calc_normality(
    molarity: float | int | Sequence[float | int] | NDArray[np.number],
    equivalence_factor: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    as_list: bool = False,
) -> float | NDArray[np.float64] | list[float]:
    """Calculate normality from numeric molarity and equivalence factor.

    Parameters
    ----------
    molarity : float | int | Sequence[float | int] | NDArray[np.number]
        Molar concentration of the solute on the desired basis.
    equivalence_factor : float | int | Sequence[float | int] | NDArray[np.number]
        Reaction-context equivalence factor. Must be supplied by the caller and
        must be greater than zero.
    as_list : bool, optional
        Return sequence calculations as a Python list instead of a NumPy array.

    Returns
    -------
    float | NDArray[np.float64] | list[float]
        Normality on the same volume basis as ``molarity``. Scalar inputs return
        a float, array-like inputs return a NumPy array by default, and
        ``as_list=True`` returns a Python list.

    Notes
    -----
    Equation: ``normality = molarity * equivalence_factor``. Numeric inputs do
    not carry unit metadata, so no unit conversion or unit inference is done.
    """
    # SECTION: Validate inputs
    molarity_values: NDArray[np.float64] = np.asarray(
        molarity,
        dtype=np.float64,
    )
    equivalence_factor_values: NDArray[np.float64] = np.asarray(
        equivalence_factor,
        dtype=np.float64,
    )
    _validate_positive_values(molarity_values, "molarity")
    _validate_positive_values(equivalence_factor_values, "equivalence_factor")

    # SECTION: Calculate normality
    normality = cast(
        NDArray[np.float64],
        molarity_values * equivalence_factor_values,
    )
    if as_list:
        if normality.ndim == 0:
            return [float(normality)]
        return cast(list[float], normality.tolist())
    if normality.ndim == 0:
        return float(normality)
    return normality


def _calc_normality_from_props(
    molarity: CustomProp,
    equivalence_factor: float | int,
    output_unit: str = "eq/L",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate normality from unit-aware molarity and equivalence factor.

    Parameters
    ----------
    molarity : CustomProp
        Unit-aware molarity value.
    equivalence_factor : float | int
        Reaction-context equivalence factor. Must be greater than zero.
    output_unit : str, optional
        Normality unit used for the result. Defaults to ``"eq/L"``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Normality in ``output_unit``.

    Notes
    -----
    The denominator of ``output_unit`` defines the molarity basis used before
    applying the equivalence factor. For example, ``eq/L`` normalizes molarity
    to ``mol/L``.
    """
    # SECTION: Normalize molarity to the denominator basis of output_unit
    molarity_value = to_custom_prop_scalar(
        prop=molarity,
        output_unit=_molarity_unit_from_normality_unit(output_unit),
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )

    # SECTION: Calculate normality with the numeric core
    return _calc_normality(
        molarity=molarity_value,
        equivalence_factor=equivalence_factor,
    )


# SECTION: Public exports
__all__ = [
    "_calc_normality",
    "_calc_normality_from_props",
]
