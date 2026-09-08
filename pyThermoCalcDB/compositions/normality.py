"""Normality and equivalent-concentration helpers."""

# import libs
import math

# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils import to_annotated_value

# locals
from ..utils.conversions import _resolve_unit_conversion_fn, _to_units
from pythermodb_settings.utils.quantity import to_custom_prop_scalar


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


# ======================================================================
# *** Public annotated API
# ======================================================================
@calculation_info(
    name="normality",
    description="Calculate normality from numeric molarity and an equivalence factor.",
    equation="normality = molarity * equivalence_factor",
    inputs={
        "molarity": "Molar concentration of the solute on the desired basis.",
        "equivalence_factor": "Reaction-context equivalence factor.",
    },
    outputs={
        "normality": "Normality on the same volume basis as molarity."
    },
    aliases=(
        "equivalent concentration",
        "normal concentration",
    ),
    notes=(
        "The equivalence factor must be supplied by the caller; it is not inferred from formula, valence, charge, or stoichiometry.",
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when molarity and equivalence_factor are already expressed on that normality basis.",
    ),
    tags=(
        "normality",
        "equivalent_concentration",
        "numeric",
        "scalar",
    )
)
def _normality_annotated(
    molarity: float | int,
    equivalence_factor: float | int,
    *,
    name: str = "normality",
    description: str = "Calculate normality from molarity and an equivalence factor.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate annotated normality from numeric inputs."""
    return to_annotated_value(
        _calc_normality(
            molarity=molarity,
            equivalence_factor=equivalence_factor,
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_normality",
    )


@calculation_info(
    name="normality",
    description="Calculate normality from unit-aware molarity and an equivalence factor.",
    equation="normality = molarity * equivalence_factor",
    inputs={
        "molarity": "Unit-aware molarity of the solute.",
        "equivalence_factor": "Reaction-context equivalence factor.",
        "output_unit": "Normality unit used to normalize molarity and annotate the result.",
    },
    outputs={
        "normality": "Normality in output_unit."
    },
    aliases=(
        "equivalent concentration",
        "normal concentration",
    ),
    notes=(
        "The output_unit must be a ratio such as eq/L.",
        "The denominator of output_unit defines the molarity conversion basis, for example eq/L uses mol/L.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "normality",
        "equivalent_concentration",
        "unit_aware",
        "unit_conversion",
        "scalar",
    )
)
def _normality_from_props_annotated(
    molarity: CustomProp,
    equivalence_factor: float | int,
    output_unit: str = "eq/L",
    unit_conversion_fn: UnitConversionFn | None = None,
    *,
    name: str = "normality",
    description: str = "Calculate normality from molarity and an equivalence factor.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Calculate annotated normality from unit-aware molarity."""
    # SECTION: set default unit for output if not provided
    if unit is None:
        unit = output_unit

    # check unit & output unit consistency
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    return to_annotated_value(
        _calc_normality_from_props(
            molarity=molarity,
            equivalence_factor=equivalence_factor,
            output_unit=output_unit,
            unit_conversion_fn=unit_conversion_fn,
        ),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_normality_from_props",
    )


# ======================================================================
# *** Aliases
# ======================================================================
calc_normality = _normality_annotated

calc_normality_from_props = _normality_from_props_annotated


# SECTION: Public exports
__all__ = [
    "_calc_normality",
    "calc_normality",
    "_calc_normality_from_props",
    "calc_normality_from_props",
]
