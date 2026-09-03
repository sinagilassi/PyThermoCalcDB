"""Chemical-potential identity helpers."""

# import libs
import math

# >> pythermodb-settings
from pythermodb_settings.models import ScalarValue, Temperature
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import pos, to_scalar
# locals
from ..configs.constants import R_J_molK, P_ref_Pa
from ..utils.conversions import _resolve_unit_conversion_fn


# SECTION: Internal scalar helpers

def _scalar(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert scalar input to float, optionally normalizing units."""
    return to_scalar(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


def _pos(
    value: ScalarValue,
    name: str,
    output_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Convert scalar input to a positive float, optionally normalizing units."""
    return pos(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


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
        value = float(conversion_fn(value, unit, "K"))

    # ! Chemical-potential identities require T > 0 K.
    if value <= 0.0:
        raise ValueError("temperature must be greater than zero K.")
    return value


# SECTION: Activity-based chemical potential

def calc_chemical_potential_from_activity(
    chemical_potential_std: ScalarValue,
    activity: ScalarValue,
    temperature: Temperature,
    output_chemical_potential_unit: str | None = "J/mol",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate chemical potential from a dimensionless activity.

    Parameters
    ----------
    chemical_potential_std : float | int | CustomProp
        Standard chemical potential. Converted to
        ``output_chemical_potential_unit`` when supplied as ``CustomProp``.
    activity : float | int | CustomProp
        Dimensionless thermodynamic activity. Must be greater than zero.
    temperature : Temperature
        System temperature. Converted to K before calculation.
    output_chemical_potential_unit : str, optional
        Unit used to normalize ``chemical_potential_std``. Defaults to
        ``J/mol``.
    gas_constant : float, optional
        Gas constant in units consistent with chemical potential per mol per K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Chemical potential in the normalized energy unit.

    Notes
    -----
    Equation: ``mu_i = mu_i_std + R*T*ln(a_i)``. The activity must already be
    dimensionless relative to the selected standard state.

    Raises
    ------
    ValueError
        If activity, temperature, or gas constant is not positive.
    """
    # SECTION: Normalize inputs
    mu_std = _scalar(
        chemical_potential_std,
        "chemical_potential_std",
        output_chemical_potential_unit,
        unit_conversion_fn,
    )
    # ! Activity must already be dimensionless and strictly positive.
    a = _pos(activity, "activity")
    temperature_k = _temperature_k(temperature, unit_conversion_fn)
    r = _pos(gas_constant, "gas_constant")

    # SECTION: Calculate chemical potential
    return mu_std + r * temperature_k * math.log(a)


# SECTION: Ideal-gas chemical potential

def calc_ideal_gas_chemical_potential(
    chemical_potential_std: ScalarValue,
    partial_pressure: ScalarValue,
    temperature: Temperature,
    standard_pressure: ScalarValue = P_ref_Pa,
    output_chemical_potential_unit: str | None = "J/mol",
    output_pressure_unit: str | None = "Pa",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate ideal-gas chemical potential from partial pressure.

    Parameters
    ----------
    chemical_potential_std : float | int | CustomProp
        Standard chemical potential at the same temperature/reference state.
    partial_pressure : float | int | CustomProp
        Species partial pressure. Must be positive.
    temperature : Temperature
        System temperature. Converted to K before calculation.
    standard_pressure : float | int | CustomProp, optional
        Standard pressure used to form the dimensionless pressure ratio.
        Defaults to ``P_ref_Pa``.
    output_chemical_potential_unit : str, optional
        Unit used to normalize ``chemical_potential_std``. Defaults to
        ``J/mol``.
    output_pressure_unit : str, optional
        Unit used to normalize pressure values. Defaults to ``Pa``.
    gas_constant : float, optional
        Gas constant in units consistent with chemical potential per mol per K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Ideal-gas chemical potential in the normalized energy unit.

    Notes
    -----
    Equation: ``mu_i = mu_i_std + R*T*ln(P_i/P_std)``. The logarithm is applied
    only to the dimensionless pressure ratio.
    """
    # SECTION: Build dimensionless pressure activity
    p_i = _pos(partial_pressure, "partial_pressure", output_pressure_unit, unit_conversion_fn)
    p_std = _pos(standard_pressure, "standard_pressure", output_pressure_unit, unit_conversion_fn)

    # NOTE: The logarithm is applied to P_i/P_std, not a dimensional pressure.
    return calc_chemical_potential_from_activity(
        chemical_potential_std,
        p_i / p_std,
        temperature,
        output_chemical_potential_unit,
        gas_constant,
        unit_conversion_fn,
    )


# SECTION: Fugacity-based chemical potential

def calc_chemical_potential_from_fugacity(
    chemical_potential_std: ScalarValue,
    fugacity: ScalarValue,
    temperature: Temperature,
    standard_fugacity: ScalarValue,
    output_chemical_potential_unit: str | None = "J/mol",
    output_fugacity_unit: str | None = None,
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate chemical potential from a fugacity ratio.

    Parameters
    ----------
    chemical_potential_std : float | int | CustomProp
        Standard chemical potential at the same temperature/reference state.
    fugacity : float | int | CustomProp
        Species fugacity. Must be positive.
    temperature : Temperature
        System temperature. Converted to K before calculation.
    standard_fugacity : float | int | CustomProp
        Standard-state fugacity used to form a dimensionless ratio. Must be
        positive.
    output_chemical_potential_unit : str, optional
        Unit used to normalize ``chemical_potential_std``. Defaults to
        ``J/mol``.
    output_fugacity_unit : str, optional
        Unit used to normalize fugacity values when supplied as ``CustomProp``.
    gas_constant : float, optional
        Gas constant in units consistent with chemical potential per mol per K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Real-gas chemical potential in the normalized energy unit.

    Notes
    -----
    Equation: ``mu_i = mu_i_std + R*T*ln(f_i/f_std)``. This function does not
    calculate fugacity or fugacity coefficients.
    """
    # SECTION: Build dimensionless fugacity activity
    f_i = _pos(fugacity, "fugacity", output_fugacity_unit, unit_conversion_fn)
    f_std = _pos(standard_fugacity, "standard_fugacity", output_fugacity_unit, unit_conversion_fn)

    # NOTE: Fugacity coefficients are model outputs and are not calculated here.
    return calc_chemical_potential_from_activity(
        chemical_potential_std,
        f_i / f_std,
        temperature,
        output_chemical_potential_unit,
        gas_constant,
        unit_conversion_fn,
    )


# SECTION: Solution chemical potential

def calc_solution_chemical_potential(
    chemical_potential_std: ScalarValue,
    mole_fraction: ScalarValue,
    activity_coefficient: ScalarValue,
    temperature: Temperature,
    output_chemical_potential_unit: str | None = "J/mol",
    gas_constant: float = R_J_molK,
    unit_conversion_fn: UnitConversionFn | None = None,
) -> float:
    """Calculate solution chemical potential from mole fraction and activity coefficient.

    Parameters
    ----------
    chemical_potential_std : float | int | CustomProp
        Standard chemical potential at the same temperature/reference state.
    mole_fraction : float | int | CustomProp
        Species mole fraction. Must be greater than zero and no greater than
        one.
    activity_coefficient : float | int | CustomProp
        Supplied activity coefficient. Must be greater than zero.
    temperature : Temperature
        System temperature. Converted to K before calculation.
    output_chemical_potential_unit : str, optional
        Unit used to normalize ``chemical_potential_std``. Defaults to
        ``J/mol``.
    gas_constant : float, optional
        Gas constant in units consistent with chemical potential per mol per K.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function.

    Returns
    -------
    float
        Solution chemical potential in the normalized energy unit.

    Notes
    -----
    Uses ``a_i = x_i*gamma_i`` followed by
    ``mu_i = mu_i_std + R*T*ln(a_i)``. Activity coefficients are supplied by the
    caller/model layer and are not calculated here.
    """
    # SECTION: Build dimensionless solution activity
    x_i = _pos(mole_fraction, "mole_fraction")
    if x_i > 1.0:
        raise ValueError("mole_fraction must be no greater than one.")
    gamma_i = _pos(activity_coefficient, "activity_coefficient")

    # NOTE: gamma_i is supplied by a model layer; this function only uses it.
    return calc_chemical_potential_from_activity(
        chemical_potential_std,
        x_i * gamma_i,
        temperature,
        output_chemical_potential_unit,
        gas_constant,
        unit_conversion_fn,
    )


# SECTION: Public exports
__all__ = [
    "calc_chemical_potential_from_activity",
    "calc_ideal_gas_chemical_potential",
    "calc_chemical_potential_from_fugacity",
    "calc_solution_chemical_potential",
]

