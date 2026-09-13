"""Specific volume and density conversion helpers."""

# import libs
from pythermodb_settings.models import CustomProp, ScalarValue
from pythermodb_settings.models.units import UnitConversionFn
# locals
from .core.specific_volume import (
    _calc_density_to_specific_volume_from_props,
    _calc_specific_volume_to_density_from_props,
    _calc_density_to_specific_volume_from_scalars,
    _calc_specific_volume_to_density_from_scalars,
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
    # SECTION: Delegate unit-aware inputs to the props adapter
    if isinstance(density, CustomProp):
        return _calc_density_to_specific_volume_from_props(
            density,
            output_density_unit,
            unit_conversion_fn,
        )

    # SECTION: Normalize numeric scalar inputs
    return _calc_density_to_specific_volume_from_scalars(
        density,
        output_density_unit,
        unit_conversion_fn,
    )


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
    # SECTION: Delegate unit-aware inputs to the props adapter
    if isinstance(specific_volume, CustomProp):
        return _calc_specific_volume_to_density_from_props(
            specific_volume,
            output_specific_volume_unit,
            unit_conversion_fn,
        )

    # SECTION: Normalize numeric scalar inputs
    return _calc_specific_volume_to_density_from_scalars(
        specific_volume,
        output_specific_volume_unit,
        unit_conversion_fn,
    )


# SECTION: Public exports
__all__ = [
    "density_to_specific_volume",
    "specific_volume_to_density"
]
