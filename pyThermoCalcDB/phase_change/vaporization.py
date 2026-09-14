"""Vaporization-property public wrappers."""

# locals
from ..thermo.phase_change import calc_enthalpy_vaporization_watson
from .core.vaporization import _calc_heat_of_vaporization_watson


# SECTION: Public wrappers
def calc_heat_of_vaporization_watson(
    heat_of_vaporization_ref,
    temperature,
    reference_temperature,
    critical_temperature,
    exponent: float = 0.38,
    output_unit: str = "J/mol",
    unit_conversion_fn=None,
) -> float:
    """Correct heat of vaporization from a reference temperature with Watson.

    Equation: ``Hvap2 = Hvap1 * ((1 - Tr2) / (1 - Tr1))**exponent``. Inputs
    use the same units as ``calc_enthalpy_vaporization_watson`` and the result
    is returned in ``output_unit``. The correlation is empirical and intended
    for temperatures below ``critical_temperature``.
    """
    return calc_enthalpy_vaporization_watson(
        heat_of_vaporization_ref,
        reference_temperature,
        temperature,
        critical_temperature,
        exponent=exponent,
        output_unit=output_unit,
        unit_conversion_fn=unit_conversion_fn,
    )


__all__ = [
    "calc_heat_of_vaporization_watson",
    "_calc_heat_of_vaporization_watson",
]
