"""Internal energy identity helpers."""

# import libs
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
# locals
from .core.internal_energy import (
    _calc_internal_energy_from_props,
    _calc_ideal_gas_internal_energy_from_props,
    _calc_internal_energy_from_scalars,
    _calc_ideal_gas_internal_energy_from_scalars,
)


# SECTION: Internal energy calculations

def calc_internal_energy(
    enthalpy: ScalarValue,
    pressure: ScalarValue,
    volume: ScalarValue,
    output_enthalpy_unit: str | None = None,
    output_pressure_unit: str | None = None,
    output_volume_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate internal energy from enthalpy, pressure, and volume.

    Parameters
    ----------
    enthalpy : float | int | CustomProp
        Enthalpy on the same amount basis as the desired internal energy.
    pressure : float | int | CustomProp
        Pressure.
    volume : float | int | CustomProp
        Volume on the same amount basis implied by ``enthalpy``.
    output_enthalpy_unit : str, optional
        Unit used to normalize ``enthalpy`` before calculation.
    output_pressure_unit : str, optional
        Unit used to normalize ``pressure`` before calculation.
    output_volume_unit : str, optional
        Unit used to normalize ``volume`` before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Internal energy.

    Notes
    -----
    # NOTE: Equation
    U = H - P*V
    """
    # SECTION: Delegate unit-aware inputs to the props adapter
    if (
        isinstance(enthalpy, CustomProp)
        and isinstance(pressure, CustomProp)
        and isinstance(volume, CustomProp)
    ):
        return _calc_internal_energy_from_props(
            enthalpy=enthalpy,
            pressure=pressure,
            volume=volume,
            output_enthalpy_unit=output_enthalpy_unit,
            output_pressure_unit=output_pressure_unit,
            output_volume_unit=output_volume_unit,
            unit_conversion_fn=unit_conversion_fn,
        )

    # SECTION: Normalize mixed/numeric scalar inputs
    return _calc_internal_energy_from_scalars(
        enthalpy=enthalpy,
        pressure=pressure,
        volume=volume,
        output_enthalpy_unit=output_enthalpy_unit,
        output_pressure_unit=output_pressure_unit,
        output_volume_unit=output_volume_unit,
        unit_conversion_fn=unit_conversion_fn,
    )


def calc_ideal_gas_internal_energy(
    molar_enthalpy: ScalarValue,
    temperature: Temperature,
    output_molar_enthalpy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    universal_gas_constant: float = 8.31446261815324,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal-gas molar internal energy from molar enthalpy.

    Parameters
    ----------
    molar_enthalpy : float | int | CustomProp
        Ideal-gas molar enthalpy.
    temperature : Temperature
        Temperature value used in the ``R*T`` product. When
        ``output_temperature_unit`` is ``None``, ``temperature.value`` is used
        as-is. When ``output_temperature_unit`` is provided, ``temperature`` is
        converted to that unit before calculation.
    output_molar_enthalpy_unit : str, optional
        Unit used to normalize ``molar_enthalpy`` before calculation.
    output_temperature_unit : str, optional
        Unit used to normalize ``temperature`` before calculation. Leave as
        ``None`` to use ``temperature.value`` as supplied.
    universal_gas_constant : float, optional
        Gas constant in units consistent with ``molar_enthalpy`` per the
        temperature unit used in the ``R*T`` product.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Ideal-gas molar internal energy.

    Notes
    -----
    # NOTE: Equation
    U_molar = H_molar - R*T
    """
    # SECTION: Delegate unit-aware inputs to the props adapter
    if isinstance(molar_enthalpy, CustomProp):
        return _calc_ideal_gas_internal_energy_from_props(
            molar_enthalpy=molar_enthalpy,
            temperature=temperature,
            output_molar_enthalpy_unit=output_molar_enthalpy_unit,
            output_temperature_unit=output_temperature_unit,
            universal_gas_constant=universal_gas_constant,
            unit_conversion_fn=unit_conversion_fn,
        )

    # SECTION: Normalize mixed/numeric scalar inputs
    return _calc_ideal_gas_internal_energy_from_scalars(
        molar_enthalpy=molar_enthalpy,
        temperature=temperature,
        output_molar_enthalpy_unit=output_molar_enthalpy_unit,
        output_temperature_unit=output_temperature_unit,
        universal_gas_constant=universal_gas_constant,
        unit_conversion_fn=unit_conversion_fn,
    )


# SECTION: Public exports
__all__ = ["calc_internal_energy", "calc_ideal_gas_internal_energy"]
