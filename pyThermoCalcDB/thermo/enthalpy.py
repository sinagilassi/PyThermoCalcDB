# import libs
import logging
import numpy as np
from typing import Optional, Dict, Any
from pythermodb_settings.models import Temperature
import pycuc
# locals

# NOTE: logger setup
logger = logging.getLogger(__name__)


# SECTION: enthalpy calculations using NASA 9-coefficient polynomial
# NOTE: calculate ideal gas enthalpy using NASA 9-coefficient polynomial
def En_IG_NASA9_polynomial(
    a1: float,
    a2: float,
    a3: float,
    a4: float,
    a5: float,
    a6: float,
    a7: float,
    b1: float,
    temperature: Temperature,
    temperature_range: Optional[tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,  # J/mol.K
    message: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Calculate the ideal gas enthalpy (En_IG) using NASA 9-coefficient polynomial coefficients. The NASA polynomial equation is defined as:

    En_IG(T) = R * ( -a1/T + a2 ln T + a3 T + a4 T^2/2 + a5 T^3/3 + a6 T^4/4 + a7 T^5/5 + b1 )

    where En is the enthalpy at temperature T.

    Parameters
    ----------
    a1 : float
        NASA polynomial coefficient a1
    a2 : float
        NASA polynomial coefficient a2
    a3 : float
        NASA polynomial coefficient a3
    a4 : float
        NASA polynomial coefficient a4
    a5 : float
        NASA polynomial coefficient a5
    a6 : float
        NASA polynomial coefficient a6
    a7 : float
        NASA polynomial coefficient a7
    b1 : float
        NASA polynomial coefficient b1
    temperature : Temperature
        Temperature at which to calculate enthalpy defined in pythermodb_settings.models.Temperature, should be in Kelvin
    output_unit : Optional[str], optional
        Desired output unit for enthalpy (default is None, which returns J/mol)
    universal_gas_constant : float, optional
        Universal gas constant in J/mol.K (default is 8.31446261815324 J/mol.K)
    message : Optional[str], optional
        Optional message regarding the calculation

    Returns
    -------
    Optional[NASA9PolynomialIdealGasEnthalpyResult]
        Result model containing the calculated ideal gas enthalpy and related information, or None if an error occurs.
    """
    try:
        # SECTION: Input Validation
        # NOTE: extract temperature value in [K]
        T_value = temperature.value  # assuming temperature is in Kelvin
        T_unit = temperature.unit

        # >> check temperature unit
        if T_unit != "K":
            T_value = pycuc.convert_from_to(
                value=T_value, from_unit=T_unit, to_unit="K"
            )

        # NOTE: check temperature range validity
        if temperature_range is not None:
            T_low = temperature_range[0].value
            T_high = temperature_range[1].value
            T_low_unit = temperature_range[0].unit.strip()
            T_high_unit = temperature_range[1].unit.strip()

            # >> convert to K if necessary
            if T_low_unit != "K":
                T_low = pycuc.convert_from_to(
                    value=T_low, from_unit=T_low_unit, to_unit="K"
                )
            if T_high_unit != "K":
                T_high = pycuc.convert_from_to(
                    value=T_high, from_unit=T_high_unit, to_unit="K"
                )

            # >> check validity
            if not (T_low <= T_value <= T_high):
                logger.warning(
                    f"Temperature {T_value} K is out of the specified range [{T_low} K, {T_high} K]."
                )
                return None

        # NOTE: check coefficients
        coeffs = [a1, a2, a3, a4, a5, a6, a7, b1]
        for i, coeff in enumerate(coeffs):
            if not isinstance(coeff, (int, float)):
                logger.error(
                    f"Coefficient a{i+1} is not a valid number: {coeff}")
                return None

        # SECTION: calculate En using NASA 9-coefficient polynomial equation, temperature in K, and R in J/mol.K
        # En in J/mol
        En_value = universal_gas_constant * (
            -a1 / T_value
            + a2 * np.log(T_value)
            + a3 * T_value
            + (a4/2.0) * T_value**2
            + (a5/3.0) * T_value**3
            + (a6/4.0) * T_value**4
            + (a7/5.0) * T_value**5
            + b1
        )
        En_value = float(En_value)

        # set unit
        if output_unit is None:
            En_unit = "J/mol"
        else:
            # new unit
            En_unit = output_unit
            # convert
            En_value = pycuc.convert_from_to(
                value=En_value,
                from_unit="J/mol",
                to_unit=En_unit
            )

        # return result model
        res = {
            "result": {
                "value": En_value,
                "unit": En_unit,
                "symbol": "En_IG"
            },
            "message": message if message is not None else "Ideal gas enthalpy calculation using NASA-9 successful"
        }

        return res
    except Exception as e:
        logger.error(f"Error in ideal gas enthalpy calculation: {e}")
        return None

# NOTE: calculate enthalpy using NASA 9-coefficient integral polynomial at a given temperature range


def En_IG_NASA9_polynomial_range(
    a1: float,
    a2: float,
    a3: float,
    a4: float,
    a5: float,
    a6: float,
    a7: float,
    b1: float,
    T_low: Temperature,
    T_high: Temperature,
    T_points: int = 10,
    temperature_range: Optional[tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,  #
    message: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Calculate the integral ideal gas enthalpy (En_IG) over a temperature range using NASA 9-coefficient polynomial coefficients.

    Parameters
    ----------
    a1 : float
        NASA polynomial coefficient a1
    a2 : float
        NASA polynomial coefficient a2
    a3 : float
        NASA polynomial coefficient a3
    a4 : float
        NASA polynomial coefficient a4
    a5 : float
        NASA polynomial coefficient a5
    a6 : float
        NASA polynomial coefficient a6
    a7 : float
        NASA polynomial coefficient a7
    b1 : float
        NASA polynomial coefficient b1
    T_low : Temperature
        Lower temperature bound defined in pythermodb_settings.models.Temperature, should be in Kelvin
    T_high : Temperature
        Upper temperature bound defined in pythermodb_settings.models.Temperature, should be in Kelvin
    T_points : int, optional
        Number of temperature points to calculate within the range (default is 10)
    output_unit : Optional[str], optional
        Desired output unit for enthalpy (default is None, which returns J/mol)
    universal_gas_constant : float, optional
        Universal gas constant in J/mol.K (default is 8.31446261815324 J/mol.K)
    message : Optional[str], optional
        Optional message regarding the calculation

    Returns
    -------
    Optional[Dict[str, Any]]
        A dictionary containing the calculated integral ideal gas enthalpy values over the temperature range, or None if an error occurs.

    Notes
    -----
    - The NASA polynomial is commonly used to represent thermodynamic properties of substances over a range of temperatures.
    - The calculated enthalpy values are in J/mol by default, but can be converted to other units if specified.
    - In case of errors during calculation, enthalpy values are set to zero for the respective temperature points.
    """
    try:
        # SECTION: generate temperature points between T_low and T_high
        T_low_value = T_low.value
        T_high_value = T_high.value
        T_low_unit = T_low.unit.strip()
        T_high_unit = T_high.unit.strip()

        # >> convert to K if necessary
        if T_low_unit != "K":
            T_low_value = pycuc.convert_from_to(
                value=T_low_value, from_unit=T_low_unit, to_unit="K"
            )
        if T_high_unit != "K":
            T_high_value = pycuc.convert_from_to(
                value=T_high_value, from_unit=T_high_unit, to_unit="K"
            )

        temperatures = np.linspace(T_low_value, T_high_value, T_points)

        # SECTION: calculate En at each temperature point
        enthalpy_results = []

        # loop over temperatures
        for T in temperatures:
            temp_obj = Temperature(value=T, unit="K")
            res = En_IG_NASA9_polynomial(
                a1, a2, a3, a4, a5, a6, a7, b1,
                temperature=temp_obj,
                temperature_range=temperature_range,
                output_unit=output_unit,
                universal_gas_constant=universal_gas_constant,
                message=None
            )

            # check result
            if res is None:
                enthalpy_results.append(0.0)  # append zero if error
                continue

            # set
            enthalpy_results.append(res["result"]['value'])

        # return result model
        res = {
            "result": {
                "values": {
                    "x": temperatures.tolist(),
                    "y": enthalpy_results,
                },
                "unit": output_unit if output_unit is not None else "J/mol",
                "symbol": "En_IG"
            },
            "message": message if message is not None else "Integral ideal gas enthalpy calculation over range using NASA-9 successful"
        }
        return res
    except Exception as e:
        logger.error(
            f"Error in integral ideal gas enthalpy range calculation: {e}")
        return None


# NOTE: delta En using NASA 9-coefficient polynomial between two temperatures
def dEn_IG_NASA9_polynomial(
    a1: float,
    a2: float,
    a3: float,
    a4: float,
    a5: float,
    a6: float,
    a7: float,
    b1: float,
    T_initial: Temperature,
    T_final: Temperature,
    temperature_range: Optional[tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,  #
    message: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Calculate enthalpy change (dEn_IG) between two temperatures using NASA 9-coefficient polynomial coefficients.

    Parameters
    ----------
    a1 : float
        NASA polynomial coefficient a1
    a2 : float
        NASA polynomial coefficient a2
    a3 : float
        NASA polynomial coefficient a3
    a4 : float
        NASA polynomial coefficient a4
    a5 : float
        NASA polynomial coefficient a5
    a6 : float
        NASA polynomial coefficient a6
    a7 : float
        NASA polynomial coefficient a7
    b1 : float
        NASA polynomial coefficient b1
    T_initial : Temperature
        Initial temperature defined in pythermodb_settings.models.Temperature, should be in Kelvin
    T_final : Temperature
        Final temperature defined in pythermodb_settings.models.Temperature, should be in Kelvin
    output_unit : Optional[str], optional
        Desired output unit for enthalpy change (default is None, which returns J/mol)
    universal_gas_constant : float, optional
        Universal gas constant in J/mol.K (default is 8.31446261815324 J/mol.K)
    message : Optional[str], optional
        Optional message regarding the calculation

    Returns
    -------
    Optional[Dict[str, Any]]
        A dictionary containing the calculated sensible heat effect and related information, or None if an error occurs.
    """
    try:
        # SECTION: calculate En at initial temperature
        res_initial = En_IG_NASA9_polynomial(
            a1, a2, a3, a4, a5, a6, a7, b1,
            temperature=T_initial,
            temperature_range=temperature_range,
            output_unit=output_unit,
            universal_gas_constant=universal_gas_constant,
            message=None
        )
        if res_initial is None:
            return None

        En_initial = res_initial["result"]["value"]

        # SECTION: calculate En at final temperature
        res_final = En_IG_NASA9_polynomial(
            a1, a2, a3, a4, a5, a6, a7, b1,
            temperature=T_final,
            temperature_range=temperature_range,
            output_unit=output_unit,
            universal_gas_constant=universal_gas_constant,
            message=None
        )
        if res_final is None:
            return None

        En_final = res_final["result"]["value"]

        # SECTION: calculate sensible heat effect
        delta_En = En_final - En_initial
        delta_En = float(delta_En)

        # set unit
        if output_unit is None:
            delta_En_unit = "J/mol"
        else:
            delta_En_unit = output_unit

        # return result model
        res = {
            "result": {
                "value": delta_En,
                "unit": delta_En_unit,
                "symbol": "dEn_IG"
            },
            "message": message if message is not None else "Sensible heat effect calculation using NASA-9 successful"
        }
        return res

    except Exception as e:
        logger.error(f"Error in sensible heat effect calculation: {e}")
        return None

# SECTION: enthalpy calculations using Shomate equation


def En_IG_shomate(
    A: float,
    B: float,
    C: float,
    D: float,
    E: float,
    F: float,
    G: float,
    temperature: Temperature,
    temperature_range: Optional[tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,  # J/mol.K
    message: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Calculate the ideal gas enthalpy (En_IG) using Shomate equation coefficients. The Shomate equation is defined as:

    En_IG(T) = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F, where t = T/1000

    where En is the enthalpy [kJ/mol] at temperature T [K].

    Parameters
    ----------
    A : float
        Shomate equation coefficient A
    B : float
        Shomate equation coefficient B
    C : float
        Shomate equation coefficient C
    D : float
        Shomate equation coefficient D
    E : float
        Shomate equation coefficient E
    F : float
        Shomate equation coefficient F
    G : float
        Shomate equation coefficient G
    temperature : Temperature
        Temperature at which to calculate enthalpy defined in pythermodb_settings.models.Temperature, should be in Kelvin
    output_unit : Optional[str], optional
        Desired output unit for enthalpy (default is None, which returns J/mol)
    universal_gas_constant : float, optional
        Universal gas constant in J/mol.K (default is 8.31446261815324 J/mol.K)
    message : Optional[str], optional
        Optional message regarding the calculation

    Returns
    -------
    Optional[Dict[str, Any]]
        A dictionary containing the calculated ideal gas enthalpy and related information, or None if an error occurs.
    """
    try:
        # SECTION: Input Validation
        # NOTE: extract temperature value in [K]
        T_value = temperature.value  # assuming temperature is in Kelvin
        T_unit = temperature.unit

        # >> check temperature unit
        if T_unit != "K":
            T_value = pycuc.convert_from_to(
                value=T_value, from_unit=T_unit, to_unit="K"
            )
            T_unit = "K"

        # NOTE: check temperature range validity
        if temperature_range is not None:
            T_low = temperature_range[0].value
            T_high = temperature_range[1].value
            T_low_unit = temperature_range[0].unit.strip()
            T_high_unit = temperature_range[1].unit.strip()

            # >> convert to K if necessary
            if T_low_unit != "K":
                T_low = pycuc.convert_from_to(
                    value=T_low, from_unit=T_low_unit, to_unit="K"
                )
            if T_high_unit != "K":
                T_high = pycuc.convert_from_to(
                    value=T_high, from_unit=T_high_unit, to_unit="K"
                )

            # >> check validity
            if not (T_low <= T_value <= T_high):
                logger.warning(
                    f"Temperature {T_value} K is out of the specified range [{T_low} K, {T_high} K]."
                )
                return None

        # NOTE: check coefficients
        coeffs = [A, B, C, D, E, F, G]
        for i, coeff in enumerate(coeffs):
            if not isinstance(coeff, (int, float)):
                logger.error(
                    f"Coefficient {['A','B','C','D','E','F','G'][i]} is not a valid number: {coeff}")
                return None

        # SECTION: calculate En using Shomate equation, temperature in K
        # En in kJ/mol
        t = T_value / 1000.0

        En_value = (
            A * t +
            B * t**2 / 2 +
            C * t**3 / 3 +
            D * t**4 / 4 -
            E / t +
            F
        )
        En_value = float(En_value)

        # set unit
        if output_unit is None:
            En_unit = "kJ/mol"
        else:
            # set new unit
            En_unit = output_unit
            # convert
            En_value = pycuc.convert_from_to(
                value=En_value,
                from_unit="kJ/mol",
                to_unit=En_unit
            )

        # return result model
        res = {
            "result": {
                "value": En_value,
                "unit": En_unit,
                "symbol": "En_IG"
            },
            "message": message if message is not None else "Ideal gas enthalpy calculation using Shomate equation successful"
        }

        return res
    except Exception as e:
        logger.error(f"Error in ideal gas enthalpy calculation: {e}")
        return None


# NOTE: calculate enthalpy using Shomate equation at a given temperature range
def En_IG_shomate_range(
    A: float,
    B: float,
    C: float,
    D: float,
    E: float,
    F: float,
    G: float,
    T_low: Temperature,
    T_high: Temperature,
    T_points: int = 10,
    temperature_range: Optional[tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,  #
    message: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """
    Calculate the integral ideal gas enthalpy (En_IG) over a temperature range using Shomate equation coefficients.

    Parameters
    ----------
    A : float
        Shomate equation coefficient A
    B : float
        Shomate equation coefficient B
    C : float
        Shomate equation coefficient C
    D : float
        Shomate equation coefficient D
    E : float
        Shomate equation coefficient E
    F : float
        Shomate equation coefficient F
    G : float
        Shomate equation coefficient G
    T_low : Temperature
        Lower temperature bound defined in pythermodb_settings.models.Temperature, should be in Kelvin
    T_high : Temperature
        Upper temperature bound defined in pythermodb_settings.models.Temperature, should be in Kelvin
    T_points : int, optional
        Number of temperature points to calculate within the range (default is 10)
    output_unit : Optional[str], optional
        Desired output unit for enthalpy (default is None, which returns J/mol)
    universal_gas_constant : float, optional
        Universal gas constant in J/mol.K (default is 8.31446261815324 J/mol.K)
    message : Optional[str], optional
        Optional message regarding the calculation

    Returns
    -------
    Optional[Dict[str, Any]]
        A dictionary containing the calculated integral ideal gas enthalpy values over the temperature range, or None if an error occurs.

    Notes
    -----
    - The Shomate equation is commonly used to represent thermodynamic properties of substances over a range of temperatures.
    - The calculated enthalpy values are in kJ/mol by default, but can be converted to other units if specified.
    - In case of errors during calculation, enthalpy values are set to zero for the respective temperature points.
    """
    try:
        # SECTION: generate temperature points between T_low and T_high
        T_low_value = T_low.value
        T_high_value = T_high.value
        T_low_unit = T_low.unit.strip()
        T_high_unit = T_high.unit.strip()

        # >> convert to K if necessary
        if T_low_unit != "K":
            T_low_value = pycuc.convert_from_to(
                value=T_low_value, from_unit=T_low_unit, to_unit="K"
            )

        if T_high_unit != "K":
            T_high_value = pycuc.convert_from_to(
                value=T_high_value, from_unit=T_high_unit, to_unit="K"
            )

        temperatures = np.linspace(T_low_value, T_high_value, T_points)

        # SECTION: calculate En at each temperature point
        enthalpy_results = []

        # loop over temperatures
        for T in temperatures:
            temp_obj = Temperature(value=T, unit="K")
            res = En_IG_shomate(
                A, B, C, D, E, F, G,
                temperature=temp_obj,
                temperature_range=temperature_range,
                output_unit=output_unit,
                universal_gas_constant=universal_gas_constant,
                message=None
            )

            # check result
            if res is None:
                enthalpy_results.append(0.0)  # append zero if error
                continue

            # set
            enthalpy_results.append(res["result"]['value'])

        # return result model
        res = {
            "result": {
                "values": {
                    "x": temperatures.tolist(),
                    "y": enthalpy_results,
                },
                "unit": output_unit if output_unit is not None else "kJ/mol",
                "symbol": "En_IG"
            },
            "message": message if message is not None else "Integral ideal gas enthalpy calculation over range using Shomate equation successful"
        }

        return res
    except Exception as e:
        logger.error(
            f"Error in integral ideal gas enthalpy range calculation: {e}")
        return None
