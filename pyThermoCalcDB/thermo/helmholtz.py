"""Helmholtz energy identity helpers."""

# import libs
from pythermodb_settings.models import CustomProp, ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
# locals
from .core.helmholtz import (
    _calc_helmholtz_energy_from_props,
    _calc_helmholtz_energy_from_scalars,
)


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
    # SECTION: Delegate unit-aware inputs to the props adapter
    if (
        isinstance(internal_energy, CustomProp)
        and isinstance(entropy, CustomProp)
    ):
        return _calc_helmholtz_energy_from_props(
            internal_energy=internal_energy,
            temperature=temperature,
            entropy=entropy,
            output_internal_energy_unit=output_internal_energy_unit,
            output_entropy_unit=output_entropy_unit,
            output_temperature_unit=output_temperature_unit,
            unit_conversion_fn=unit_conversion_fn,
        )

    # SECTION: Normalize mixed/numeric scalar inputs
    return _calc_helmholtz_energy_from_scalars(
        internal_energy=internal_energy,
        temperature=temperature,
        entropy=entropy,
        output_internal_energy_unit=output_internal_energy_unit,
        output_entropy_unit=output_entropy_unit,
        output_temperature_unit=output_temperature_unit,
        unit_conversion_fn=unit_conversion_fn,
    )


# SECTION: Public exports
__all__ = ["calc_helmholtz_energy"]
