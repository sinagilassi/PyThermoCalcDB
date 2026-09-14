"""Thermodynamic departure-property public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, ScalarValue
from pythermodb_settings.utils import to_annotated_value

# locals
from ..utils.conversions import _pos, _scalar
from .core.departure import (
    _calc_cp_departure_from_eos_derivatives,
    _calc_dimensionless_enthalpy_departure,
    _calc_dimensionless_entropy_departure,
    _calc_enthalpy_departure,
    _calc_enthalpy_from_ideal_and_departure,
    _calc_entropy_departure,
    _calc_entropy_from_ideal_and_departure,
)


# SECTION: Public wrappers
def calc_enthalpy_departure(enthalpy: ScalarValue, ideal_gas_enthalpy: ScalarValue) -> AnnotatedValue[float]:
    """Calculate annotated enthalpy departure ``H - H_id``."""
    return to_annotated_value(
        float(_calc_enthalpy_departure(_scalar(enthalpy, "enthalpy"), _scalar(ideal_gas_enthalpy, "ideal_gas_enthalpy"))),
        name="enthalpy_departure",
        description="Calculate enthalpy departure.",
        unit="J/mol",
        symbol="H^R",
        implementation="_calc_enthalpy_departure",
    )


def calc_entropy_departure(entropy: ScalarValue, ideal_gas_entropy: ScalarValue) -> AnnotatedValue[float]:
    """Calculate annotated entropy departure ``S - S_id``."""
    return to_annotated_value(
        float(_calc_entropy_departure(_scalar(entropy, "entropy"), _scalar(ideal_gas_entropy, "ideal_gas_entropy"))),
        name="entropy_departure",
        description="Calculate entropy departure.",
        unit="J/(mol.K)",
        symbol="S^R",
        implementation="_calc_entropy_departure",
    )


def calc_dimensionless_enthalpy_departure(
    enthalpy_departure: ScalarValue,
    temperature: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> AnnotatedValue[float]:
    """Calculate annotated dimensionless enthalpy departure ``H^R/(R*T)``."""
    return to_annotated_value(
        float(_calc_dimensionless_enthalpy_departure(
            _scalar(enthalpy_departure, "enthalpy_departure"),
            _pos(temperature, "temperature"),
            _pos(gas_constant, "gas_constant"),
        )),
        name="dimensionless_enthalpy_departure",
        description="Calculate dimensionless enthalpy departure.",
        unit=None,
        symbol="H^R/(RT)",
        implementation="_calc_dimensionless_enthalpy_departure",
    )


def calc_dimensionless_entropy_departure(
    entropy_departure: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> AnnotatedValue[float]:
    """Calculate annotated dimensionless entropy departure ``S^R/R``."""
    return to_annotated_value(
        float(_calc_dimensionless_entropy_departure(
            _scalar(entropy_departure, "entropy_departure"),
            _pos(gas_constant, "gas_constant"),
        )),
        name="dimensionless_entropy_departure",
        description="Calculate dimensionless entropy departure.",
        unit=None,
        symbol="S^R/R",
        implementation="_calc_dimensionless_entropy_departure",
    )


def calc_enthalpy_from_ideal_and_departure(
    ideal_gas_enthalpy: ScalarValue,
    enthalpy_departure: ScalarValue,
) -> AnnotatedValue[float]:
    """Calculate annotated real-fluid enthalpy from ideal and departure terms."""
    return to_annotated_value(
        float(_calc_enthalpy_from_ideal_and_departure(
            _scalar(ideal_gas_enthalpy, "ideal_gas_enthalpy"),
            _scalar(enthalpy_departure, "enthalpy_departure"),
        )),
        name="enthalpy",
        description="Calculate enthalpy from ideal-gas and departure contributions.",
        unit="J/mol",
        symbol="H",
        implementation="_calc_enthalpy_from_ideal_and_departure",
    )


def calc_entropy_from_ideal_and_departure(
    ideal_gas_entropy: ScalarValue,
    entropy_departure: ScalarValue,
) -> AnnotatedValue[float]:
    """Calculate annotated real-fluid entropy from ideal and departure terms."""
    return to_annotated_value(
        float(_calc_entropy_from_ideal_and_departure(
            _scalar(ideal_gas_entropy, "ideal_gas_entropy"),
            _scalar(entropy_departure, "entropy_departure"),
        )),
        name="entropy",
        description="Calculate entropy from ideal-gas and departure contributions.",
        unit="J/(mol.K)",
        symbol="S",
        implementation="_calc_entropy_from_ideal_and_departure",
    )


def calc_cp_departure_from_eos_derivatives(
    temperature: ScalarValue,
    integral_d2p_dt2_dv: ScalarValue,
    dpressure_dtemperature_at_volume: ScalarValue,
    dpressure_dvolume_at_temperature: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
) -> AnnotatedValue[float]:
    """Calculate annotated generic Cp departure from EOS derivative data."""
    return to_annotated_value(
        float(_calc_cp_departure_from_eos_derivatives(
            _pos(temperature, "temperature"),
            _scalar(integral_d2p_dt2_dv, "integral_d2p_dt2_dv"),
            _scalar(dpressure_dtemperature_at_volume, "dpressure_dtemperature_at_volume"),
            _scalar(dpressure_dvolume_at_temperature, "dpressure_dvolume_at_temperature"),
            _pos(gas_constant, "gas_constant"),
        )),
        name="heat_capacity_departure",
        description="Calculate heat-capacity departure from generic EOS derivative data.",
        unit="J/(mol.K)",
        symbol="Cp^R",
        implementation="_calc_cp_departure_from_eos_derivatives",
    )


__all__ = [
    "calc_enthalpy_departure",
    "calc_entropy_departure",
    "calc_dimensionless_enthalpy_departure",
    "calc_dimensionless_entropy_departure",
    "calc_enthalpy_from_ideal_and_departure",
    "calc_entropy_from_ideal_and_departure",
    "calc_cp_departure_from_eos_derivatives",
]
