# import libs
import logging
from pycuc import convert_from_to
from ...utils.conversions import _scalar
# locals


# NOTE: set up logger
logger = logging.getLogger(__name__)


# ! ::: Ideal-gas and reciprocal-density helpers

def _calc_ideal_gas_density(
        pressure: float,
        molecular_weight: float,
        temperature: float,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "kg/m3",
) -> float:
    """
    Calculate ideal-gas density from pressure, molecular weight, and temperature.

    Parameters
    ----------
    pressure : float
        Gas pressure in Pa.
    molecular_weight : float
        Molecular weight in kg/mol.
    temperature : float
        Gas temperature in K.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to 8.31446261815324.
    output_unit : str, optional
        Desired density unit. Defaults to kg/m3.

    Returns
    -------
    float
        Ideal-gas density in kg/m3 or None if conversion/calculation fails.

    Notes
    -----
    Equation
        `rho = P*M / (R*T)`
    """
    try:
        # SECTION: Validate and normalize inputs
        if pressure <= 0 or temperature <= 0 or molecular_weight <= 0:
            logger.error(
                "Pressure, temperature, and molecular weight must be positive.")
            raise ValueError(
                "Pressure, temperature, and molecular weight must be positive.")

        # NOTE: pressure in Pa
        p_value = pressure

        # NOTE: temperature in K
        t_value = temperature

        # NOTE: molecular weight in kg/mol
        mw_value = float(molecular_weight)

        # SECTION: Calculate ideal-gas density
        density_value = p_value * mw_value / (universal_gas_constant * t_value)
        density_unit = "kg/m3"

        # NOTE: output unit conversion if necessary
        if output_unit != density_unit:
            density_value = convert_from_to(
                density_value, density_unit, output_unit)
            density_unit = output_unit

        return density_value
    except Exception as e:
        logger.error(f"Error in ideal gas density calculation: {e}")
        raise e


# ! ::: Ideal-gas molar volume calculation
def _calc_ideal_gas_molar_volume(
        temperature: float,
        pressure: float,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "m3/mol",
) -> float:
    """
    Calculate ideal-gas molar volume from temperature and pressure.

    Parameters
    ----------
    temperature : float
        Gas temperature in K.
    pressure : float
        Gas pressure in Pa.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to 8.31446261815324.
    output_unit : str, optional
        Desired molar-volume unit. Defaults to m3/mol.

    Returns
    -------
    float
        Ideal-gas molar volume in m3/mol or None if conversion/calculation fails.

    Notes
    -----
    Equation
        `V_m = R*T / P`
    """
    try:
        # SECTION: Validate and normalize inputs
        if pressure <= 0 or temperature <= 0:
            logger.error("Pressure and temperature must be positive.")
            raise ValueError("Pressure and temperature must be positive.")

        # NOTE: pressure in Pa
        p_value = pressure

        # NOTE: temperature in K
        t_value = temperature

        # SECTION: Calculate ideal-gas molar volume
        volume_value = universal_gas_constant * t_value / p_value
        volume_unit = "m3/mol"

        # NOTE: output unit conversion if necessary
        if output_unit != volume_unit:
            volume_value = convert_from_to(
                volume_value,
                volume_unit,
                output_unit
            )

        return volume_value
    except Exception as e:
        logger.error(f"Error in ideal gas molar volume calculation: {e}")
        raise e

# ! ::: Compressibility-factor gas relations


def _calc_gas_molar_volume_from_z(
        temperature: float,
        pressure: float,
        compressibility_factor: float,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "m3/mol",
        unit_conversion_fn=None,
) -> float:
    """
    Calculate real-gas molar volume from a supplied compressibility factor.

    Parameters
    ----------
    temperature : float
        Gas temperature in K.
    pressure : float
        Gas pressure in Pa.
    compressibility_factor : float
        Supplied compressibility factor ``Z``. Must be greater than zero.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to ``8.31446261815324``.
    output_unit : str, optional
        Desired molar-volume unit. Defaults to ``m3/mol``.
    unit_conversion_fn : callable, optional
        Unit conversion function. Defaults to ``convert_from_to``.

    Returns
    -------
    float
        Real-gas molar volume in m3/mol.

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
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn

    # NOTE: pressure in Pa
    p_value = float(pressure)
    # >> check
    if p_value <= 0.0:
        raise ValueError("pressure must be greater than zero.")

    # NOTE: temperature in K
    t_value = float(temperature)
    # >> check
    if t_value <= 0.0:
        raise ValueError("temperature must be greater than zero K.")

    # ! Z is supplied by a model/source; this function does not calculate it.
    z_value = _scalar(
        compressibility_factor,
        "compressibility_factor"
    )
    r_value = _scalar(
        universal_gas_constant,
        "universal_gas_constant"
    )

    # SECTION: Calculate molar volume
    volume_value = z_value * r_value * t_value / p_value
    volume_unit = "m3/mol"
    if output_unit != volume_unit:
        volume_value = conversion_fn(volume_value, volume_unit, output_unit)
        volume_unit = output_unit

    return volume_value

# ! ::: Calculate real-gas density from a supplied compressibility factor


def _calc_gas_density_from_z(
        pressure: float,
        molecular_weight: float,
        temperature: float,
        compressibility_factor: float | int,
        universal_gas_constant: float = 8.31446261815324,
        output_unit: str = "kg/m3",
        unit_conversion_fn=None,
) -> float:
    """
    Calculate real-gas density from a supplied compressibility factor.

    Parameters
    ----------
    pressure : float
        Gas pressure in Pa.
    molecular_weight : float
        Molecular weight in kg/mol.
    temperature : float
        Gas temperature in K.
    compressibility_factor : float | int
        Supplied compressibility factor ``Z``. Must be greater than zero.
    universal_gas_constant : float, optional
        Gas constant in J/mol/K. Defaults to ``8.31446261815324``.
    output_unit : str, optional
        Desired density unit. Defaults to ``kg/m3``.
    unit_conversion_fn : callable, optional
        Unit conversion function. Defaults to ``convert_from_to``.

    Returns
    -------
    float
        Real-gas density in kg/m3.

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
    conversion_fn = convert_from_to if unit_conversion_fn is None else unit_conversion_fn

    # NOTE: pressure in Pa
    p_value = float(pressure)
    # >> check
    if p_value <= 0.0:
        raise ValueError("pressure must be greater than zero.")

    # NOTE: temperature in K
    t_value = float(temperature)
    # >> check
    if t_value <= 0.0:
        raise ValueError("temperature must be greater than zero K.")

    # NOTE: molecular weight in kg/mol
    mw_value = float(molecular_weight)
    # >> check
    if mw_value <= 0.0:
        raise ValueError("molecular_weight must be greater than zero.")

    z_value = _scalar(
        compressibility_factor,
        "compressibility_factor"
    )
    r_value = _scalar(
        universal_gas_constant,
        "universal_gas_constant"
    )

    # SECTION: Calculate density
    density_value = p_value * mw_value / (z_value * r_value * t_value)
    density_unit = "kg/m3"
    if output_unit != density_unit:
        density_value = conversion_fn(density_value, density_unit, output_unit)

    return density_value


# all
__all__ = [
    "_calc_ideal_gas_density",
    "_calc_ideal_gas_molar_volume",
    "_calc_gas_molar_volume_from_z",
    "_calc_gas_density_from_z",
]
