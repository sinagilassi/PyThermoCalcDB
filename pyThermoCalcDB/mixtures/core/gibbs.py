"""Core ideal Gibbs-energy-of-mixing calculations."""

# import libs
from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, CustomProp, Temperature, UnitConversionFn
from pythermodb_settings.utils.quantity import pos, to_dict

# locals
from ...configs.constants import R_J_molK
from ...utils.conversions import (
    NumericArrayInput,
    _as_float_array,
    _configure_component_values,
    _resolve_unit_conversion_fn,
    _return_scalar_if_zero_dim,
    _validate_custom_prop_mapping,
    _validate_custom_prop_scalar,
    _validate_fraction_array,
    _validate_positive_scalar,
)

# SECTION: Type aliases
NumericInput = NumericArrayInput


def _temperature_k(
    temperature: Temperature,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Return absolute temperature in K."""
    # SECTION: Normalize temperature
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    value = float(temperature.value)
    unit = temperature.unit.strip()
    if unit != "K":
        value = float(conversion_fn(value=value, from_unit=unit, to_unit="K"))

    # ! Gibbs energy identities require T > 0 K.
    return _validate_positive_scalar(value, "temperature")


# SECTION: Ideal-mixing shared term

def _x_log_x_sum(mole_fractions: NDArray[np.float64]) -> np.float64 | NDArray[np.float64]:
    # NOTE: Apply the x*ln(x) -> 0 limiting value at zero mole fraction.
    terms = np.zeros_like(mole_fractions, dtype=np.float64)
    mask = mole_fractions > 0.0
    terms[mask] = mole_fractions[mask] * np.log(mole_fractions[mask])
    return np.sum(terms, axis=-1)


# SECTION: Ideal molar Gibbs energy of mixing

def _calc_ideal_molar_gibbs_energy_of_mixing(
    mole_fractions: NumericInput,
    temperature: float | int,
    gas_constant: float = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate Delta G_mix = R * T * sum_i(x_i * ln(x_i)).

    For 2-D inputs, axis 0 is states and axis 1 is components.
    """
    # SECTION: Normalize and validate
    x = _as_float_array(mole_fractions, "mole_fractions")
    t = _validate_positive_scalar(temperature, "temperature")
    r = _validate_positive_scalar(gas_constant, "gas_constant")
    _validate_fraction_array(x, "mole_fractions")

    # SECTION: Calculate ideal molar Gibbs energy of mixing
    return _return_scalar_if_zero_dim(r * t * _x_log_x_sum(x))


# SECTION: Total ideal Gibbs energy of mixing

def _calc_ideal_gibbs_energy_of_mixing(
    total_moles: float | int,
    mole_fractions: NumericInput,
    temperature: float | int,
    gas_constant: float = R_J_molK,
) -> float | NDArray[np.float64]:
    """Calculate total ideal Gibbs energy of mixing."""
    # SECTION: Validate extensive amount and scale molar result
    n_total = _validate_positive_scalar(total_moles, "total_moles")
    return _return_scalar_if_zero_dim(
        n_total * np.asarray(
            _calc_ideal_molar_gibbs_energy_of_mixing(
                mole_fractions, temperature, gas_constant),
            dtype=np.float64,
        )
    )


# SECTION: Mapping adapter

def _calc_ideal_molar_gibbs_energy_of_mixing_from_mapping(
    mole_fractions: Mapping[str, float | int],
    temperature: float | int,
    gas_constant: float = R_J_molK,
) -> float:
    """Calculate ideal molar Gibbs energy of mixing from keyed mole fractions."""
    # NOTE: Component identity is not needed after public remapping/order is complete.
    return float(
        _calc_ideal_molar_gibbs_energy_of_mixing(
            list(mole_fractions.values()),
            temperature,
            gas_constant,
        )
    )


# SECTION: Props adapters

def _calc_ideal_molar_gibbs_energy_of_mixing_from_props(
    mole_fractions: Mapping[str, CustomProp],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate ideal molar Gibbs energy of mixing from unit-aware keyed inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_mapping(mole_fractions, "mole_fractions")

    # SECTION: Normalize composition and temperature
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    temperature_value = _temperature_k(temperature, conversion_fn)
    x = to_dict(mole_fractions, unit_conversion_fn=conversion_fn)

    # SECTION: Remap component keys
    x = _configure_component_values(
        x,
        list(components) if components is not None else None,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "mole_fractions",
    )

    # SECTION: Calculate from mapping
    return _calc_ideal_molar_gibbs_energy_of_mixing_from_mapping(
        x,
        temperature_value,
        gas_constant,
    )


def _calc_ideal_gibbs_energy_of_mixing_from_props(
    total_moles: CustomProp,
    mole_fractions: Mapping[str, CustomProp],
    temperature: Temperature,
    gas_constant: float = R_J_molK,
    output_total_moles_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Sequence[Component] | None = None,
    component_key: ComponentKey | None = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> float:
    """Calculate total ideal Gibbs energy of mixing from unit-aware keyed inputs."""
    # SECTION: Validate props input contract
    _validate_custom_prop_scalar(total_moles, "total_moles")
    _validate_custom_prop_mapping(mole_fractions, "mole_fractions")

    # SECTION: Normalize total amount
    conversion_fn = _resolve_unit_conversion_fn(unit_conversion_fn)
    n_total = pos(
        total_moles,
        "total_moles",
        output_total_moles_unit,
        unit_conversion_fn=conversion_fn,
    )

    # SECTION: Calculate molar Gibbs energy from props and scale
    return float(
        n_total * _calc_ideal_molar_gibbs_energy_of_mixing_from_props(
            mole_fractions,
            temperature,
            gas_constant,
            conversion_fn,
            components,
            component_key,
            case_sensitive,
            sort_by_components_order,
        )
    )


# SECTION: Core exports
__all__ = [
    "_calc_ideal_molar_gibbs_energy_of_mixing",
    "_calc_ideal_gibbs_energy_of_mixing",
    "_calc_ideal_molar_gibbs_energy_of_mixing_from_mapping",
    "_calc_ideal_molar_gibbs_energy_of_mixing_from_props",
    "_calc_ideal_gibbs_energy_of_mixing_from_props",
]
