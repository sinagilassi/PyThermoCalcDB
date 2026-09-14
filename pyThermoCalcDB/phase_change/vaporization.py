"""Vaporization-property public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue
from pythermodb_settings.utils import to_annotated_value

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
    *,
    name: str = "heat_of_vaporization",
    description: str = "Correct heat of vaporization with the Watson correlation.",
    symbol: str | None = "delta_H_vap",
) -> AnnotatedValue[float]:
    """Correct annotated heat of vaporization from a reference temperature with Watson."""
    # NOTE: Existing thermo wrapper returns a scalar; this package adds annotation metadata.
    value = calc_enthalpy_vaporization_watson(
        heat_of_vaporization_ref,
        reference_temperature,
        temperature,
        critical_temperature,
        exponent=exponent,
        output_unit=output_unit,
        unit_conversion_fn=unit_conversion_fn,
    )
    return to_annotated_value(
        value,
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_heat_of_vaporization_watson",
    )


__all__ = [
    "calc_heat_of_vaporization_watson",
    "_calc_heat_of_vaporization_watson",
]
