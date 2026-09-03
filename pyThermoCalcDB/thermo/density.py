# import libs
import logging
from typing import Dict, Optional, Any
from pythermodb_settings.models import (
    Temperature,
    Pressure,
    CustomProp,
)
import pycuc
# locals


# NOTE: set up logger
logger = logging.getLogger(__name__)


def rackett(
        temperature: Temperature,
        critical_temperature: Temperature,
        critical_pressure: Pressure,
        molecular_weight: CustomProp,
        critical_compressibility: CustomProp,
        message: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Calculate the density of a substance using the Rackett equation.

    Parameters
    ----------
    temperature : Temperature
        The temperature at which to calculate the density.
    critical_temperature : Temperature
        The critical temperature of the substance.
    critical_pressure : Pressure
        The critical pressure of the substance.
    molecular_weight : CustomProp
        The molecular weight of the substance, in kg/mol.
    critical_compressibility : CustomProp
        The critical compressibility factor of the substance, dimensionless.
    message : Optional[Dict], optional
        A dictionary to store messages or warnings.

    Returns
    -------
    Optional[Dict[str, Any]]
        A dictionary containing the calculation results, or None if an error occurs.
    """
    try:
        # SECTION: Input validation
        # NOTE: check for valid temperature
        if temperature.value <= 0 or critical_temperature.value <= 0:
            logger.warning("Temperature values must be greater than zero.")
            return None

        # NOTE: check for valid pressure
        if critical_pressure.value <= 0:
            logger.warning("Critical pressure must be greater than zero.")
            return None

        # SECTION: Unit conversions
        T_value = temperature.value
        T_unit = temperature.unit

        # >> to K
        if T_unit != "K":
            T_value = pycuc.convert_from_to(
                value=T_value,
                from_unit=T_unit,
                to_unit="K",
            )

        # NOTE: critical temperature to [K]
        Tc_value = critical_temperature.value
        Tc_unit = critical_temperature.unit

        if Tc_unit != "K":
            Tc_value = pycuc.convert_from_to(
                value=Tc_value,
                from_unit=Tc_unit,
                to_unit="K",
            )

        # NOTE: critical pressure to [bar]
        Pc_value = critical_pressure.value
        Pc_unit = critical_pressure.unit
        if Pc_unit != "bar":
            Pc_value = pycuc.convert_from_to(
                value=Pc_value,
                from_unit=Pc_unit,
                to_unit="bar",
            )

        # NOTE: molecular weight to [kg/mol]
        MW_value = float(molecular_weight.value)  # kg/mol
        MW_unit = molecular_weight.unit

        # NOTE: compressibility factor (dimensionless)
        Z_value = float(critical_compressibility.value)
        Z_unit = critical_compressibility.unit

        # ! >> to kg/mol
        if MW_unit != "kg/mol":
            if MW_unit == "g/mol":
                MW_value = MW_value / 1000  # kg/mol
            elif MW_unit == "kg/kmol":
                MW_value = MW_value / 1000  # kg/mol
            else:
                logger.warning(
                    f"Unsupported molecular weight unit: {MW_unit}. Expected 'kg/mol'."
                )
                return None

        # SECTION: Rackett equation calculation
        # ! Gas constant R = 8.314 J/mol·K = 8.314e-5 bar·m³/mol·K
        R = 8.31446261815324e-5  # bar·m³/mol·K

        Tr = T_value / Tc_value  # Reduced temperature
        exponent = 1 + (1 - Tr)**(2/7)

        # ! >> Saturated molar volume in m³/mol
        V_sat = (R * Tc_value / Pc_value) * Z_value**exponent

        # NOTE: Density: ρ = M / V_sat
        # ! >> kg/m³
        rho = MW_value / V_sat

        res = {
            "result": {
                "value": rho,
                "unit": "kg/m3",
                "symbol": 'rho_LIQ'
            },
            "message": message,
        }

        # res
        return res
    except Exception as e:
        logger.error(f"Error in rackett calculation: {e}")
        return None


# SECTION: Ideal-gas and reciprocal-density helpers

def calc_ideal_gas_density(
        pressure: Pressure,
        molecular_weight: CustomProp,
        temperature: Temperature,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "kg/m3",
) -> Optional[CustomProp]:
    """Calculate ideal-gas density from pressure, molecular weight, and temperature.

    Parameters
    ----------
    pressure : Pressure
        Gas pressure. Converted to Pa before calculation.
    molecular_weight : CustomProp
        Molecular weight. Converted to kg/mol before calculation.
    temperature : Temperature
        Gas temperature. Converted to K before calculation.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to 8.31446261815324.
    output_unit : str, optional
        Desired density unit. Defaults to kg/m3.

    Returns
    -------
    Optional[CustomProp]
        Ideal-gas density or None if conversion/calculation fails.

    Notes
    -----
    Equation
        `rho = P*M / (R*T)`
    """
    try:
        # SECTION: Validate and normalize inputs
        if pressure.value <= 0 or temperature.value <= 0 or molecular_weight.value <= 0:
            logger.warning(
                "Pressure, temperature, and molecular weight must be positive.")
            return None

        p_value = pressure.value
        if pressure.unit != "Pa":
            p_value = pycuc.convert_from_to(p_value, pressure.unit, "Pa")

        t_value = temperature.value
        if temperature.unit != "K":
            t_value = pycuc.convert_from_to(t_value, temperature.unit, "K")

        mw_value = float(molecular_weight.value)
        if molecular_weight.unit != "kg/mol":
            mw_value = pycuc.convert_from_to(
                mw_value, molecular_weight.unit, "kg/mol")

        # SECTION: Calculate ideal-gas density
        density_value = p_value * mw_value / (universal_gas_constant * t_value)
        density_unit = "kg/m3"
        if output_unit != density_unit:
            density_value = pycuc.convert_from_to(
                density_value, density_unit, output_unit)
            density_unit = output_unit

        return CustomProp(value=density_value, unit=density_unit)
    except Exception as e:
        logger.error(f"Error in ideal gas density calculation: {e}")
        return None


def calc_ideal_gas_molar_volume(
        temperature: Temperature,
        pressure: Pressure,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "m3/mol",
) -> Optional[CustomProp]:
    """Calculate ideal-gas molar volume from temperature and pressure.

    Parameters
    ----------
    temperature : Temperature
        Gas temperature. Converted to K before calculation.
    pressure : Pressure
        Gas pressure. Converted to Pa before calculation.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to 8.31446261815324.
    output_unit : str, optional
        Desired molar-volume unit. Defaults to m3/mol.

    Returns
    -------
    Optional[CustomProp]
        Ideal-gas molar volume or None if conversion/calculation fails.

    Notes
    -----
    Equation
        `V_m = R*T / P`
    """
    try:
        # SECTION: Validate and normalize inputs
        if pressure.value <= 0 or temperature.value <= 0:
            logger.warning("Pressure and temperature must be positive.")
            return None

        p_value = pressure.value
        if pressure.unit != "Pa":
            p_value = pycuc.convert_from_to(p_value, pressure.unit, "Pa")

        t_value = temperature.value
        if temperature.unit != "K":
            t_value = pycuc.convert_from_to(t_value, temperature.unit, "K")

        # SECTION: Calculate ideal-gas molar volume
        volume_value = universal_gas_constant * t_value / p_value
        volume_unit = "m3/mol"
        if output_unit != volume_unit:
            volume_value = pycuc.convert_from_to(
                volume_value, volume_unit, output_unit)
            volume_unit = output_unit

        return CustomProp(value=volume_value, unit=volume_unit)
    except Exception as e:
        logger.error(f"Error in ideal gas molar volume calculation: {e}")
        return None

# SECTION: Compressibility-factor gas relations

def _scalar_value(
        value,
        name: str,
        output_unit: str | None = None,
        unit_conversion_fn=None,
) -> float:
    """Convert scalar input to float, optionally normalizing units."""
    # NOTE: Import locally to keep legacy module imports stable.
    from pythermodb_settings.utils.quantity import to_scalar
    from ..utils.conversions import _resolve_unit_conversion_fn

    return to_scalar(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


def _positive_scalar_value(
        value,
        name: str,
        output_unit: str | None = None,
        unit_conversion_fn=None,
) -> float:
    """Convert scalar input to a positive float, optionally normalizing units."""
    # NOTE: Import locally to keep legacy module imports stable.
    from pythermodb_settings.utils.quantity import pos
    from ..utils.conversions import _resolve_unit_conversion_fn

    return pos(
        value,
        name,
        output_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )


def calc_gas_molar_volume_from_z(
        temperature: Temperature,
        pressure: Pressure,
        compressibility_factor,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "m3/mol",
        unit_conversion_fn=None,
) -> CustomProp:
    """Calculate real-gas molar volume from a supplied compressibility factor.

    Parameters
    ----------
    temperature : Temperature
        Gas temperature. Converted to K before calculation.
    pressure : Pressure
        Gas pressure. Converted to Pa before calculation.
    compressibility_factor : float | int | CustomProp
        Supplied compressibility factor ``Z``. Must be greater than zero.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to ``8.31446261815324``.
    output_unit : str, optional
        Desired molar-volume unit. Defaults to ``m3/mol``.
    unit_conversion_fn : callable, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    CustomProp
        Real-gas molar volume in ``output_unit``.

    Notes
    -----
    Equation: ``V_m = Z*R*T/P``. This function does not calculate ``Z``; it
    only uses a caller/model supplied compressibility factor.

    Raises
    ------
    ValueError
        If pressure, temperature, gas constant, or ``Z`` is not positive.
    """
    # SECTION: Validate and normalize inputs
    conversion_fn = pycuc.convert_from_to if unit_conversion_fn is None else unit_conversion_fn

    p_value = float(pressure.value)
    if pressure.unit != "Pa":
        p_value = float(conversion_fn(p_value, pressure.unit, "Pa"))
    if p_value <= 0.0:
        raise ValueError("pressure must be greater than zero.")

    t_value = float(temperature.value)
    if temperature.unit != "K":
        t_value = float(conversion_fn(t_value, temperature.unit, "K"))
    if t_value <= 0.0:
        raise ValueError("temperature must be greater than zero K.")

    # ! Z is supplied by a model/source; this function does not calculate it.
    z_value = _positive_scalar_value(compressibility_factor, "compressibility_factor")
    r_value = _positive_scalar_value(universal_gas_constant, "universal_gas_constant")

    # SECTION: Calculate molar volume
    volume_value = z_value * r_value * t_value / p_value
    volume_unit = "m3/mol"
    if output_unit != volume_unit:
        volume_value = conversion_fn(volume_value, volume_unit, output_unit)
        volume_unit = output_unit

    return CustomProp(value=volume_value, unit=volume_unit)


def calc_gas_density_from_z(
        pressure: Pressure,
        molecular_weight: CustomProp,
        temperature: Temperature,
        compressibility_factor,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "kg/m3",
        unit_conversion_fn=None,
) -> CustomProp:
    """Calculate real-gas density from a supplied compressibility factor.

    Parameters
    ----------
    pressure : Pressure
        Gas pressure. Converted to Pa before calculation.
    molecular_weight : CustomProp
        Molecular weight. Converted to ``kg/mol`` before calculation.
    temperature : Temperature
        Gas temperature. Converted to K before calculation.
    compressibility_factor : float | int | CustomProp
        Supplied compressibility factor ``Z``. Must be greater than zero.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to ``8.31446261815324``.
    output_unit : str, optional
        Desired density unit. Defaults to ``kg/m3``.
    unit_conversion_fn : callable, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    CustomProp
        Real-gas density in ``output_unit``.

    Notes
    -----
    Equation: ``rho = P*M/(Z*R*T)``. For ``Z = 1`` this reduces to the ideal-gas
    density equation. This function does not calculate ``Z`` from an EOS.

    Raises
    ------
    ValueError
        If pressure, temperature, molecular weight, gas constant, or ``Z`` is
        not positive.
    """
    # SECTION: Validate and normalize inputs
    conversion_fn = pycuc.convert_from_to if unit_conversion_fn is None else unit_conversion_fn

    p_value = float(pressure.value)
    if pressure.unit != "Pa":
        p_value = float(conversion_fn(p_value, pressure.unit, "Pa"))
    if p_value <= 0.0:
        raise ValueError("pressure must be greater than zero.")

    t_value = float(temperature.value)
    if temperature.unit != "K":
        t_value = float(conversion_fn(t_value, temperature.unit, "K"))
    if t_value <= 0.0:
        raise ValueError("temperature must be greater than zero K.")

    mw_value = _positive_scalar_value(
        molecular_weight,
        "molecular_weight",
        "kg/mol",
        unit_conversion_fn,
    )
    z_value = _positive_scalar_value(compressibility_factor, "compressibility_factor")
    r_value = _positive_scalar_value(universal_gas_constant, "universal_gas_constant")

    # SECTION: Calculate density
    density_value = p_value * mw_value / (z_value * r_value * t_value)
    density_unit = "kg/m3"
    if output_unit != density_unit:
        density_value = conversion_fn(density_value, density_unit, output_unit)
        density_unit = output_unit

    return CustomProp(value=density_value, unit=density_unit)

