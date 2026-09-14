"""Core vaporization-property corrections."""

# locals
from ...thermo.core.phase_change import _calc_enthalpy_vaporization_watson


# SECTION: Core numeric calculations
def _calc_heat_of_vaporization_watson(
    heat_of_vaporization_ref,
    temperature,
    reference_temperature,
    critical_temperature,
    exponent=0.38,
):
    """Calculate Watson heat-of-vaporization temperature correction.

    Equation: ``Hvap2 = Hvap1 * ((1 - Tr2) / (1 - Tr1))**exponent``.
    Temperatures are absolute temperatures in K and ``critical_temperature`` is
    the critical temperature in K. Heat of vaporization is commonly J/mol, and
    the output keeps that same basis. This empirical Watson correlation is valid
    below the critical temperature and preserves scalar/array behavior.
    """
    # NOTE: Reuse the existing thermo core kernel to keep one numeric source.
    return _calc_enthalpy_vaporization_watson(
        heat_of_vaporization_ref,
        reference_temperature,
        temperature,
        critical_temperature,
        exponent,
    )


__all__ = ["_calc_heat_of_vaporization_watson"]
