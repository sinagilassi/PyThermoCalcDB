"""Hildebrand solubility-parameter public wrappers."""

# import libs
from pythermodb_settings.models import ScalarValue

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
) -> float:
    """Calculate ``delta_U_vap = delta_H_vap - R*T`` in J/mol."""
    return float(_calc_internal_energy_of_vaporization(
        _pos(heat_of_vaporization, "heat_of_vaporization"),
        _pos(temperature, "temperature"),
        _pos(gas_constant, "gas_constant"),
    ))


def calc_cohesive_energy_density(
    internal_energy_of_vaporization: ScalarValue,
    molar_volume: ScalarValue,
) -> float:
    """Calculate cohesive energy density ``CED = delta_U_vap / Vm`` in J/m3."""
    return float(_calc_cohesive_energy_density(
        _pos(internal_energy_of_vaporization, "internal_energy_of_vaporization"),
        _pos(molar_volume, "molar_volume"),
    ))


def calc_solubility_parameter(cohesive_energy_density: ScalarValue) -> float:
    """Calculate Hildebrand solubility parameter ``delta = sqrt(CED)``."""
    return float(_calc_solubility_parameter(_pos(cohesive_energy_density, "cohesive_energy_density")))


def calc_solubility_parameter_from_hvap_volume(
    heat_of_vaporization: ScalarValue,
    temperature: ScalarValue,
    molar_volume: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> float:
    """Calculate Hildebrand parameter from ``delta_H_vap``, ``T``, and ``Vm``."""
    return float(_calc_solubility_parameter_from_hvap_volume(
        _pos(heat_of_vaporization, "heat_of_vaporization"),
        _pos(temperature, "temperature"),
        _pos(molar_volume, "molar_volume"),
        _pos(gas_constant, "gas_constant"),
    ))


def calc_solubility_parameter_from_hvap_density(
    heat_of_vaporization: ScalarValue,
    temperature: ScalarValue,
    molar_density: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> float:
    """Calculate Hildebrand parameter from ``delta_H_vap``, ``T``, and molar density."""
    return float(_calc_solubility_parameter_from_hvap_density(
        _pos(heat_of_vaporization, "heat_of_vaporization"),
        _pos(temperature, "temperature"),
        _pos(molar_density, "molar_density"),
        _pos(gas_constant, "gas_constant"),
    ))


__all__ = [
    "calc_internal_energy_of_vaporization",
    "calc_cohesive_energy_density",
    "calc_solubility_parameter",
    "calc_solubility_parameter_from_hvap_volume",
    "calc_solubility_parameter_from_hvap_density",
]
