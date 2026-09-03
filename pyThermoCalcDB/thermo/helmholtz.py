"""Helmholtz energy identity helpers."""

# import libs
from pycuc import convert_from_to
from pythermodb_settings.models import ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_scalar


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
                value=temperature_value,
                from_unit=temperature_unit,
                to_unit=output_temperature_unit,
            )
        )

    # ! Thermodynamic identities require a physically valid absolute temperature.
    validation_unit = output_temperature_unit or temperature_unit
    try:
        temperature_k = temperature_value
        if validation_unit != "K":
            temperature_k = float(conversion_fn(
                value=temperature_value,
                from_unit=validation_unit,
                to_unit="K",
            ))
        if temperature_k <= 0.0:
            raise ValueError(
                "temperature must be greater than zero K after conversion.")
    except Exception:
        # NOTE: If a custom unit cannot be converted to K, validate the numeric
        # value used in the T*S product directly.
        if temperature_value <= 0.0:
            raise ValueError(
                f"temperature must be greater than zero {validation_unit}."
            )

    return temperature_value


# SECTION: Helmholtz energy calculations

def calc_helmholtz_energy(
    internal_energy: ScalarValue,
    temperature: Temperature,
    entropy: ScalarValue,
    output_internal_energy_unit: str | None = None,
    output_entropy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate Helmholtz energy from internal energy and entropy.

    Parameters
    ----------
    internal_energy : float | int | CustomProp
        Internal energy on the desired amount basis.
    temperature : Temperature
        Temperature value used in the ``T*S`` product. When
        ``output_temperature_unit`` is ``None``, ``temperature.value`` is used
        as-is. When ``output_temperature_unit`` is provided, ``temperature`` is
        converted to that unit before calculation. Supported unit labels are
        determined by pycuc, commonly ``C``, ``K``, ``R``, and ``F``.
    entropy : float | int | CustomProp
        Entropy on the same amount basis as ``internal_energy`` and per the
        temperature unit used in the ``T*S`` product.
    output_internal_energy_unit : str, optional
        Unit used to normalize ``internal_energy`` before calculation.
    output_entropy_unit : str, optional
        Unit used to normalize ``entropy`` before calculation.
    output_temperature_unit : str, optional
        Unit used to normalize ``temperature`` before calculation. Leave as
        ``None`` to use ``temperature.value`` as supplied.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Helmholtz energy.

    Notes
    -----
    Equation
        `A = U - T*S`
    """
    # SECTION: Validate and normalize inputs
    u = _scalar(
        internal_energy,
        "internal_energy",
        output_internal_energy_unit,
        unit_conversion_fn,
    )
    t = _temperature(
        temperature,
        output_temperature_unit,
        unit_conversion_fn,
    )
    s = _scalar(entropy, "entropy", output_entropy_unit, unit_conversion_fn)

    # SECTION: Calculate Helmholtz energy
    return u - t * s


# SECTION: Public exports
__all__ = ["calc_helmholtz_energy"]
