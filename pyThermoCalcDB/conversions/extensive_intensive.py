"""Extensive and intensive property conversion helpers.

# SECTION: Scope
These helpers convert between total properties and molar or mass-specific
intensive properties using explicitly supplied amount or mass.
"""

# import libs
from pythermodb_settings.models import ScalarValue
from pythermodb_settings.models.units import UnitConversionFn
# locals
from ..utils.conversions import _scalar, _pos


# SECTION: Molar extensive/intensive conversions

# ! ::: Conversions between molar and total properties
def molar_property_to_total(
    moles: ScalarValue,
    molar_property: ScalarValue,
    output_moles_unit: str | None = None,
    output_molar_property_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a molar property to a total extensive property.

    Parameters
    ----------
    moles : float | int | CustomProp
        Amount of substance.
    molar_property : float | int | CustomProp
        Property on a molar basis.
    output_moles_unit : str, optional
        Unit used to normalize ``moles`` before calculation.
    output_molar_property_unit : str, optional
        Unit used to normalize ``molar_property`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Total extensive property.

    Notes
    -----
    Equation
        `Y_total = n * Y_molar`
    """
    # SECTION: Validate and calculate
    n = _pos(moles, "moles", output_moles_unit, unit_conversion_fn)
    y_molar = _scalar(
        molar_property,
        "molar_property",
        output_molar_property_unit,
        unit_conversion_fn,
    )
    return n * y_molar

# ! ::: Conversions between total and molar properties


def total_to_molar_property(
    total_property: ScalarValue,
    moles: ScalarValue,
    output_total_property_unit: str | None = None,
    output_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a total extensive property to a molar property.

    Parameters
    ----------
    total_property : float | int | CustomProp
        Total extensive property.
    moles : float | int | CustomProp
        Amount of substance.
    output_total_property_unit : str, optional
        Unit used to normalize ``total_property`` before calculation.
    output_moles_unit : str, optional
        Unit used to normalize ``moles`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Property on a molar basis.

    Notes
    -----
    Equation
        `Y_molar = Y_total / n`
    """
    # SECTION: Validate and calculate
    y_total = _scalar(
        total_property,
        "total_property",
        output_total_property_unit,
        unit_conversion_fn,
    )
    n = _pos(moles, "moles", output_moles_unit, unit_conversion_fn)
    return y_total / n


# SECTION: Mass extensive/intensive conversions

# ! ::: Conversions between mass-specific and total properties
def specific_property_to_total(
    mass: ScalarValue,
    specific_property: ScalarValue,
    output_mass_unit: str | None = None,
    output_specific_property_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a mass-specific property to a total extensive property.

    Parameters
    ----------
    mass : float | int | CustomProp
        Mass of material.
    specific_property : float | int | CustomProp
        Property on a mass basis.
    output_mass_unit : str, optional
        Unit used to normalize ``mass`` before calculation.
    output_specific_property_unit : str, optional
        Unit used to normalize ``specific_property`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Total extensive property.

    Notes
    -----
    Equation
        `Y_total = m * Y_specific`
    """
    # SECTION: Validate and calculate
    mass_value = _pos(mass, "mass", output_mass_unit, unit_conversion_fn)
    y_specific = _scalar(
        specific_property,
        "specific_property",
        output_specific_property_unit,
        unit_conversion_fn,
    )
    return mass_value * y_specific

# ! ::: Conversions between total and mass-specific properties


def total_to_specific_property(
    total_property: ScalarValue,
    mass: ScalarValue,
    output_total_property_unit: str | None = None,
    output_mass_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a total extensive property to a mass-specific property.

    Parameters
    ----------
    total_property : float | int | CustomProp
        Total extensive property.
    mass : float | int | CustomProp
        Mass of material.
    output_total_property_unit : str, optional
        Unit used to normalize ``total_property`` before calculation.
    output_mass_unit : str, optional
        Unit used to normalize ``mass`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Property on a mass-specific basis.

    Notes
    -----
    Equation
        `Y_specific = Y_total / m`
    """
    # SECTION: Validate and calculate
    y_total = _scalar(
        total_property,
        "total_property",
        output_total_property_unit,
        unit_conversion_fn,
    )
    mass_value = _pos(mass, "mass", output_mass_unit, unit_conversion_fn)
    return y_total / mass_value


# SECTION: Public exports
__all__ = [
    "molar_property_to_total",
    "specific_property_to_total",
    "total_to_molar_property",
    "total_to_specific_property",
]
