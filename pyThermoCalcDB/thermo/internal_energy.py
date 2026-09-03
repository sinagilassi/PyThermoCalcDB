"""Internal energy identity helpers."""

# import libs
from pycuc import convert_from_to
from pythermodb_settings.models import ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos, to_scalar


# SECTION: Internal helpers
def _resolve_unit_conversion_fn(
    unit_conversion_fn: UnitConversionFn | None,
) -> UnitConversionFn:
    """Return the provided converter or the module default converter."""
    return convert_from_to if unit_conversion_fn is None else unit_conversion_fn


def _scalar(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert a scalar input to float, optionally normalizing units."""
    return to_scalar(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


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


def _temperature(
    temperature: Temperature,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Return temperature value, optionally converted to the requested unit."""
    # SECTION: Resolve conversion function
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    # SECTION: Normalize temperature only when requested
    temperature_value = float(temperature.value)
    temperature_unit = temperature.unit.strip()
    if output_temperature_unit is not None and temperature_unit != output_temperature_unit:
        temperature_value = float(
            conversion_fn(
                temperature_value,
                temperature_unit,
                output_temperature_unit,
            )
        )

    # ! Thermodynamic identities require a physically valid absolute temperature.
    validation_unit = output_temperature_unit or temperature_unit
    try:
        temperature_k = temperature_value
        if validation_unit != "K":
            temperature_k = float(conversion_fn(temperature_value, validation_unit, "K"))
        if temperature_k <= 0.0:
            raise ValueError("temperature must be greater than zero K after conversion.")
    except Exception:
        # NOTE: If a custom unit cannot be converted to K, validate the numeric
        # value used in the T term directly.
        if temperature_value <= 0.0:
            raise ValueError(f"temperature must be greater than zero {validation_unit}.")

    return temperature_value


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
    # SECTION: Validate and normalize inputs
    h = _scalar(enthalpy, "enthalpy", output_enthalpy_unit, unit_conversion_fn)
    p = _pos(pressure, "pressure", output_pressure_unit, unit_conversion_fn)
    v = _pos(volume, "volume", output_volume_unit, unit_conversion_fn)

    # SECTION: Calculate internal energy
    return h - p * v


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
    # SECTION: Validate and normalize inputs
    h = _scalar(
        molar_enthalpy,
        "molar_enthalpy",
        output_molar_enthalpy_unit,
        unit_conversion_fn,
    )
    t = _temperature(temperature, output_temperature_unit, unit_conversion_fn)
    r = _pos(universal_gas_constant, "universal_gas_constant")

    # SECTION: Calculate ideal-gas internal energy
    return h - r * t


# SECTION: Public exports
__all__ = ["calc_internal_energy", "calc_ideal_gas_internal_energy"]
