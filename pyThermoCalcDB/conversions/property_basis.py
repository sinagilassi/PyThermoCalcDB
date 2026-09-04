"""Property basis conversion helpers.

# SECTION: Scope
These helpers convert thermodynamic properties between molar and mass-specific
bases using explicitly supplied molecular weight.

# NOTE: Unit convention
Numeric values are assumed to already use the requested ``output_*_unit``. When
inputs are ``CustomProp`` objects and a matching ``output_*_unit`` is supplied,
they are converted with ``pycuc.convert_from_to`` before calculation.
"""

# import libs
from pythermodb_settings.models import ScalarValue
from pythermodb_settings.models.units import UnitConversionFn
# locals
from ..utils.conversions import _scalar, _pos

# SECTION: Molar and mass-specific property conversions


def molar_to_mass_specific(
    molar_property: ScalarValue,
    molecular_weight: ScalarValue,
    output_molar_property_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a molar property to a mass-specific property.

    Parameters
    ----------
    molar_property : float | int | CustomProp
        Property on a molar basis, for example J/mol or J/mol.K.
    molecular_weight : float | int | CustomProp
        Molecular weight, for example kg/mol.
    output_molar_property_unit : str, optional
        Unit used to normalize ``molar_property`` before calculation.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Property on a mass-specific basis.

    Notes
    -----
    Equation
        `Y_mass = Y_molar / M`
    """
    # SECTION: Validate and normalize inputs
    y_molar = _scalar(
        molar_property,
        "molar_property",
        output_molar_property_unit,
        unit_conversion_fn,
    )
    mw = _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn,
    )

    # SECTION: Calculate mass-specific property
    return y_molar / mw


def mass_specific_to_molar(
    mass_specific_property: ScalarValue,
    molecular_weight: ScalarValue,
    output_mass_specific_property_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a mass-specific property to a molar property.

    Parameters
    ----------
    mass_specific_property : float | int | CustomProp
        Property on a mass basis, for example J/kg or J/kg.K.
    molecular_weight : float | int | CustomProp
        Molecular weight, for example kg/mol.
    output_mass_specific_property_unit : str, optional
        Unit used to normalize ``mass_specific_property`` before calculation.
    output_molecular_weight_unit : str, optional
        Unit used to normalize ``molecular_weight`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Property on a molar basis.

    Notes
    -----
    Equation
        `Y_molar = Y_mass * M`
    """
    # SECTION: Validate and normalize inputs
    y_mass = _scalar(
        mass_specific_property,
        "mass_specific_property",
        output_mass_specific_property_unit,
        unit_conversion_fn,
    )
    mw = _pos(
        molecular_weight,
        "molecular_weight",
        output_molecular_weight_unit,
        unit_conversion_fn,
    )

    # SECTION: Calculate molar property
    return y_mass * mw


# SECTION: Heat capacity aliases
# ! Cp follows the same basis conversion equations as any other property.
molar_cp_to_mass_cp = molar_to_mass_specific
mass_cp_to_molar_cp = mass_specific_to_molar


# SECTION: Public exports
__all__ = [
    "molar_to_mass_specific",
    "mass_specific_to_molar",
    "molar_cp_to_mass_cp",
    "mass_cp_to_molar_cp",
]
