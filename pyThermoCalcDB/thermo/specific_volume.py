"""Specific volume and density conversion helpers."""

# import libs
from pycuc import convert_from_to
from pythermodb_settings.models import ScalarValue
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos


# SECTION: Internal helpers
def _resolve_unit_conversion_fn(
    unit_conversion_fn: UnitConversionFn | None,
) -> UnitConversionFn:
    """Return the provided converter or the module default converter."""
    return convert_from_to if unit_conversion_fn is None else unit_conversion_fn


def _pos(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a scalar input to a positive float, optionally normalizing units."""
    return pos(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


# SECTION: Density/specific-volume conversions

def density_to_specific_volume(
    density: ScalarValue,
    output_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert density to specific volume.

    Parameters
    ----------
    density : float | int | CustomProp
        Density value, for example kg/m^3.
    output_density_unit : str, optional
        Unit used to normalize ``density`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Specific volume reciprocal to the normalized density unit.

    Notes
    -----
    Equation
        `v = 1/rho`
    """
    # SECTION: Validate and calculate
    rho = _pos(density, "density", output_density_unit, unit_conversion_fn)
    return 1.0 / rho


def specific_volume_to_density(
    specific_volume: ScalarValue,
    output_specific_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert specific volume to density.

    Parameters
    ----------
    specific_volume : float | int | CustomProp
        Specific volume value, for example m^3/kg.
    output_specific_volume_unit : str, optional
        Unit used to normalize ``specific_volume`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Density reciprocal to the normalized specific-volume unit.

    Notes
    -----
    Equation
        `rho = 1/v`
    """
    # SECTION: Validate and calculate
    v = _pos(
        specific_volume,
        "specific_volume",
        output_specific_volume_unit,
        unit_conversion_fn,
    )
    return 1.0 / v


# SECTION: Public exports
__all__ = ["density_to_specific_volume", "specific_volume_to_density"]
