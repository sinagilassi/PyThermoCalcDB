"""Hildebrand solubility-parameter public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, ScalarValue
from pythermodb_settings.utils import to_annotated_value

# locals
from ..utils.conversions import _pos
from .core.hildebrand import (
    _calc_cohesive_energy_density,
    _calc_internal_energy_of_vaporization,
    _calc_solubility_parameter,
    _calc_solubility_parameter_from_hvap_density,
    _calc_solubility_parameter_from_hvap_volume,
)


# SECTION: Public wrappers
def calc_internal_energy_of_vaporization(
    heat_of_vaporization: ScalarValue,
    temperature: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "internal_energy_of_vaporization",
    description: str = "Calculate internal energy of vaporization.",
    unit: str | None = "J/mol",
    symbol: str | None = "delta_U_vap",
) -> AnnotatedValue[float]:
    """Calculate annotated ``delta_U_vap = delta_H_vap - R*T``."""
    return to_annotated_value(
        float(_calc_internal_energy_of_vaporization(
            _pos(heat_of_vaporization, "heat_of_vaporization"),
            _pos(temperature, "temperature"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_internal_energy_of_vaporization",
    )


def calc_cohesive_energy_density(
    internal_energy_of_vaporization: ScalarValue,
    molar_volume: ScalarValue,
    *,
    name: str = "cohesive_energy_density",
    description: str = "Calculate cohesive energy density.",
    unit: str | None = "J/m3",
    symbol: str | None = "CED",
) -> AnnotatedValue[float]:
    """Calculate annotated cohesive energy density ``CED = delta_U_vap / Vm``."""
    return to_annotated_value(
        float(_calc_cohesive_energy_density(
            _pos(internal_energy_of_vaporization, "internal_energy_of_vaporization"),
            _pos(molar_volume, "molar_volume"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_cohesive_energy_density",
    )


def calc_solubility_parameter(
    cohesive_energy_density: ScalarValue,
    *,
    name: str = "solubility_parameter",
    description: str = "Calculate Hildebrand solubility parameter.",
    unit: str | None = "Pa^0.5",
    symbol: str | None = "delta",
) -> AnnotatedValue[float]:
    """Calculate annotated Hildebrand solubility parameter ``delta = sqrt(CED)``."""
    return to_annotated_value(
        float(_calc_solubility_parameter(_pos(cohesive_energy_density, "cohesive_energy_density"))),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_solubility_parameter",
    )


def calc_solubility_parameter_from_hvap_volume(
    heat_of_vaporization: ScalarValue,
    temperature: ScalarValue,
    molar_volume: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "solubility_parameter",
    description: str = "Calculate Hildebrand solubility parameter from heat of vaporization and molar volume.",
    unit: str | None = "Pa^0.5",
    symbol: str | None = "delta",
) -> AnnotatedValue[float]:
    """Calculate annotated Hildebrand parameter from ``delta_H_vap``, ``T``, and ``Vm``."""
    return to_annotated_value(
        float(_calc_solubility_parameter_from_hvap_volume(
            _pos(heat_of_vaporization, "heat_of_vaporization"),
            _pos(temperature, "temperature"),
            _pos(molar_volume, "molar_volume"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_solubility_parameter_from_hvap_volume",
    )


def calc_solubility_parameter_from_hvap_density(
    heat_of_vaporization: ScalarValue,
    temperature: ScalarValue,
    molar_density: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "solubility_parameter",
    description: str = "Calculate Hildebrand solubility parameter from heat of vaporization and molar density.",
    unit: str | None = "Pa^0.5",
    symbol: str | None = "delta",
) -> AnnotatedValue[float]:
    """Calculate annotated Hildebrand parameter from ``delta_H_vap``, ``T``, and molar density."""
    return to_annotated_value(
        float(_calc_solubility_parameter_from_hvap_density(
            _pos(heat_of_vaporization, "heat_of_vaporization"),
            _pos(temperature, "temperature"),
            _pos(molar_density, "molar_density"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_solubility_parameter_from_hvap_density",
    )


__all__ = [
    "calc_internal_energy_of_vaporization",
    "calc_cohesive_energy_density",
    "calc_solubility_parameter",
    "calc_solubility_parameter_from_hvap_volume",
    "calc_solubility_parameter_from_hvap_density",
]
