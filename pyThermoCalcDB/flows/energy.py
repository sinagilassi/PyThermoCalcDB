"""Stream energy flow calculations."""

# import libs
from collections.abc import Mapping, Sequence

from pythermodb_settings.models import CustomProp
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict, to_list
from pythermodb_settings.utils.validators import positive, same_shape

# locals
from ..utils.conversions import _resolve_unit_conversion_fn
from .core.energy import (
    _calc_enthalpy_flow_rate,
    _calc_enthalpy_flow_rate_from_mapping,
    _calc_flowing_heat_capacity,
    _calc_flowing_heat_capacity_from_mapping,
)


def calc_flowing_heat_capacity(
    molar_flow_rates: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molar_heat_capacities: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molar_flow_unit: str | None = "mol/s",
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate stream flowing heat capacity: ``Cpdot = sum_i(F_i*Cp_i)``."""
    positive(molar_heat_capacities, "molar_heat_capacities")
    same_shape(molar_flow_rates, molar_heat_capacities)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    if isinstance(molar_flow_rates, Mapping) and isinstance(molar_heat_capacities, Mapping):
        f = to_dict(molar_flow_rates, output_molar_flow_unit, unit_conversion_fn=conversion_fn)
        cp = to_dict(molar_heat_capacities, output_heat_capacity_unit, unit_conversion_fn=conversion_fn)
        return _calc_flowing_heat_capacity_from_mapping(f, cp)
    if isinstance(molar_flow_rates, Mapping) or isinstance(molar_heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")

    f = to_list(molar_flow_rates, output_molar_flow_unit, unit_conversion_fn=conversion_fn)
    cp = to_list(molar_heat_capacities, output_heat_capacity_unit, unit_conversion_fn=conversion_fn)
    return float(_calc_flowing_heat_capacity(f, cp))


def calc_flowing_heat_capacity_from_mapping(
    molar_flow_rates: Mapping[str, float | int],
    molar_heat_capacities: Mapping[str, float | int],
    output_molar_flow_unit: str | None = "mol/s",
    output_heat_capacity_unit: str | None = "J/mol.K",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate flowing heat capacity from mapping inputs."""
    if not isinstance(molar_flow_rates, Mapping) or not isinstance(molar_heat_capacities, Mapping):
        raise TypeError("Both component inputs must be mappings.")
    return calc_flowing_heat_capacity(
        molar_flow_rates,
        molar_heat_capacities,
        output_molar_flow_unit,
        output_heat_capacity_unit,
        unit_conversion_fn,
    )


calc_flowing_heat_capacity_from_props = calc_flowing_heat_capacity_from_mapping


def calc_enthalpy_flow_rate(
    molar_flow_rates: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    molar_enthalpies: Mapping[str, float | int | CustomProp] | Sequence[float | int | CustomProp],
    output_molar_flow_unit: str | None = "mol/s",
    output_enthalpy_unit: str | None = "J/mol",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate stream enthalpy flow rate: ``Hdot = sum_i(F_i*h_i)``."""
    same_shape(molar_flow_rates, molar_enthalpies)
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)

    if isinstance(molar_flow_rates, Mapping) and isinstance(molar_enthalpies, Mapping):
        f = to_dict(molar_flow_rates, output_molar_flow_unit, unit_conversion_fn=conversion_fn)
        h = to_dict(molar_enthalpies, output_enthalpy_unit, unit_conversion_fn=conversion_fn)
        return _calc_enthalpy_flow_rate_from_mapping(f, h)
    if isinstance(molar_flow_rates, Mapping) or isinstance(molar_enthalpies, Mapping):
        raise TypeError("Both component inputs must be mappings or both sequences.")

    f = to_list(molar_flow_rates, output_molar_flow_unit, unit_conversion_fn=conversion_fn)
    h = to_list(molar_enthalpies, output_enthalpy_unit, unit_conversion_fn=conversion_fn)
    return float(_calc_enthalpy_flow_rate(f, h))


def calc_enthalpy_flow_rate_from_mapping(
    molar_flow_rates: Mapping[str, float | int],
    molar_enthalpies: Mapping[str, float | int],
    output_molar_flow_unit: str | None = "mol/s",
    output_enthalpy_unit: str | None = "J/mol",
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate enthalpy flow rate from mapping inputs."""
    if not isinstance(molar_flow_rates, Mapping) or not isinstance(molar_enthalpies, Mapping):
        raise TypeError("Both component inputs must be mappings.")
    return calc_enthalpy_flow_rate(
        molar_flow_rates,
        molar_enthalpies,
        output_molar_flow_unit,
        output_enthalpy_unit,
        unit_conversion_fn,
    )


calc_enthalpy_flow_rate_from_props = calc_enthalpy_flow_rate_from_mapping


__all__ = [
    "calc_flowing_heat_capacity",
    "calc_flowing_heat_capacity_from_mapping",
    "calc_flowing_heat_capacity_from_props",
    "calc_enthalpy_flow_rate",
    "calc_enthalpy_flow_rate_from_mapping",
    "calc_enthalpy_flow_rate_from_props",
]
