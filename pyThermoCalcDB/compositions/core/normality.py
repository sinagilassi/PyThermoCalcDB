"""Normality and equivalent-concentration helpers."""

# import libs
import math

# >> pythermodb-settings
from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_custom_prop_scalar

# locals
from ...utils.conversions import _resolve_unit_conversion_fn, _to_units

# ======================================================================
# *** Helper functions
# ======================================================================


def _validate_positive_scalar(
    value: float | int,
    name: str,
) -> float:
    """Validate and normalize a positive numeric scalar."""
    value_ = float(value)

    if not math.isfinite(value_):
        raise ValueError(f"{name} must be finite.")

    if value_ <= 0:
        raise ValueError(f"{name} must be greater than zero.")

    return value_


def _molarity_unit_from_normality_unit(
    output_unit: str,
) -> str:
    """Return the molarity unit matching a normality unit denominator."""
    units_ = _to_units(output_unit)
    return f"mol/{units_[1]}"


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================
def _calc_normality(
    molarity: float | int,
    equivalence_factor: float | int,
) -> float:
    """Calculate normality from numeric molarity and equivalence factor.

    Parameters
    ----------
    molarity : float | int
        Molar concentration of the solute on the desired basis.
    equivalence_factor : float | int
        Reaction-context equivalence factor. Must be supplied by the caller and
        must be greater than zero.

    Returns
    -------
    float
        Normality on the same volume basis as ``molarity``.

    Notes
    -----
    Equation: ``normality = molarity * equivalence_factor``. Numeric inputs do
    not carry unit metadata, so no unit conversion or unit inference is done.
    """
    # SECTION: Validate inputs
    molarity_value = _validate_positive_scalar(molarity, "molarity")
    equivalence_factor_value = _validate_positive_scalar(
        equivalence_factor,
        "equivalence_factor",
    )

    # SECTION: Calculate normality
    return molarity_value * equivalence_factor_value


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
