"""Normality and equivalent-concentration helpers."""

# >> pythermodb-settings
from pythermodb_settings.models import ScalarValue
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos
# locals
from ..utils.conversions import _resolve_unit_conversion_fn


# SECTION: Normality

def calc_normality(
    molarity: ScalarValue,
    equivalence_factor: ScalarValue,
    output_molarity_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate normality from molarity and an equivalence factor.

    Parameters
    ----------
    molarity : float | int | CustomProp
        Molar concentration of the solute. Numeric values are assumed to already
        use ``output_molarity_unit`` when provided.
    equivalence_factor : float | int | CustomProp
        Reaction-context equivalence factor. Must be supplied by the caller and
        must be greater than zero.
    output_molarity_unit : str, optional
        Unit used to normalize ``molarity`` when supplied as ``CustomProp``.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Normality in equivalents per volume on the same volume basis as the
        normalized molarity.

    Notes
    -----
    Equation: ``N = M*f_eq``. The equivalence factor may depend on acid-base,
    redox, precipitation, or other reaction context. This function does not
    infer it from formula, valence, charge, or stoichiometry.

    Raises
    ------
    ValueError
        If ``molarity`` or ``equivalence_factor`` is not positive.
    """
    # SECTION: Normalize inputs
    molarity_value = pos(
        molarity,
        "molarity",
        output_molarity_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )
    # ! The equivalence factor must be supplied by caller/model context.
    equivalence_factor_value = pos(equivalence_factor, "equivalence_factor")

    # SECTION: Calculate normality
    return molarity_value * equivalence_factor_value


# SECTION: Public exports
__all__ = ["calc_normality"]
