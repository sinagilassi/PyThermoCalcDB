# import libs
import logging
from typing import Optional, Dict, Any, Tuple, Literal, List
from pythermodb_settings.models import Temperature, CustomProp
import pycuc
# local
from ..utils.conversions import _to_kelvin
from .enthalpy import calc_En_IG
from .entropy import calc_Ent_IG

# NOTE: logger setup
logger = logging.getLogger(__name__)


# SECTION: Gibbs free energy calculations (ideal gas)
# NOTE:
#   G_IG(T) = H_IG(T) - T * S_IG(T)
# where:
#   - H_IG is ideal-gas enthalpy [J/mol]
#   - S_IG is ideal-gas entropy  [J/mol.K]
#   - T is temperature [K]
#
# This module focuses on a *single component* Gibbs free energy.


G_IG_Method = Literal["NASA7", "NASA9", "Shomate"]


# NOTE: calculate Gibbs free energy at a temperature

def GiFrEn_IG(
    method: G_IG_Method,
    *,
    temperature: Temperature,
    temperature_range: Optional[Tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,  # J/mol.K
    message: Optional[str] = None,
    **coeffs: Any,
) -> Optional[CustomProp]:
    """
    Calculate ideal-gas Gibbs free energy of a component at a temperature.

    Parameters
    ----------
    method:
        One of: "NASA7", "NASA9", "Shomate".
        This must match the coefficient set provided in **coeffs.

    temperature:
        Temperature object (will be converted to K if needed).

    temperature_range:
        Optional validity range. If provided and T is outside, returns None.

    output_unit:
        If None, returns "J/mol".
        Otherwise uses pycuc.convert_from_to for conversion.

    universal_gas_constant:
        Passed through to En_IG / S_IG where applicable.

    coeffs:
        Same coefficient names expected by En_IG and S_IG:
        - NASA9: a1..a7, b1 (enthalpy), b2 (entropy)
        - NASA7: a1..a7
        - Shomate: A..G

    Returns
    -------
    CustomProp
        Object containing the calculated Gibbs free energy value and unit, or None if an error occurs.

    Notes
    -----
    - Internally forces H to J/mol and S to J/mol.K before computing G.
    - Gibbs here is the *ideal-gas* reference value (useful for ΔG° and Keq workflows).
    """
    try:
        # T in K
        T_k = _to_kelvin(temperature)

        # Optional: range check (convert bounds to K)
        if temperature_range is not None:
            T_low_k = _to_kelvin(temperature_range[0])
            T_high_k = _to_kelvin(temperature_range[1])
            if not (T_low_k <= T_k <= T_high_k):
                logger.warning(
                    f"Temperature {T_k} K is out of the specified range [{T_low_k} K, {T_high_k} K]."
                )
                return None

        # H in J/mol (force for unit consistency)
        H_res = calc_En_IG(
            method=method,
            temperature=Temperature(value=T_k, unit="K"),
            temperature_range=temperature_range,
            output_unit="J/mol",
            universal_gas_constant=universal_gas_constant,
            message=None,
            **coeffs,
        )
        if H_res is None:
            return None
        H = float(H_res.value)  # J/mol

        # S in J/mol.K (force for unit consistency)
        S_res = calc_Ent_IG(
            method=method,
            temperature=Temperature(value=T_k, unit="K"),
            temperature_range=temperature_range,
            output_unit="J/mol.K",
            universal_gas_constant=universal_gas_constant,
            message=None,
            **coeffs,
        )
        if S_res is None:
            return None
        S = float(S_res.value)  # J/mol.K

        # Gibbs in J/mol
        G_value = H - T_k * S
        G_value = float(G_value)

        # Output unit conversion
        if output_unit is None:
            G_unit = "J/mol"
        else:
            G_unit = output_unit
            G_value = pycuc.convert_from_to(
                value=G_value, from_unit="J/mol", to_unit=G_unit)

        # NOTE: prepare result
        result = CustomProp(
            value=G_value,
            unit=G_unit,
        )

        # NOTE: message
        if message is not None:
            message = "Ideal gas Gibbs free energy calculation successful" + message
            print(message)

        return result

    except Exception as e:
        logger.error(f"Error in ideal gas Gibbs free energy calculation: {e}")
        return None

# NOTE: calculate Gibbs free energy at different temperatures


def GiFrEn_IG_ranges(
    method: G_IG_Method,
    *,
    temperatures: List[Temperature],
    temperature_range: Optional[Tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,  # J/mol.K
    message: Optional[str] = None,
    **coeffs: Any,
) -> Optional[Dict[str, Any]]:
    """

    Calculate ideal-gas Gibbs free energy of a component at multiple temperatures.

    Parameters
    ----------
    method:
        One of: "NASA7", "NASA9", "Shomate".
        This must match the coefficient set provided in **coeffs.
    temperatures:
        List of Temperature objects (will be converted to K if needed).
    temperature_range:
        Optional validity range. If provided and T is outside, returns None for that T.
    output_unit:
        If None, returns "J/mol".
        Otherwise uses pycuc.convert_from_to for conversion.
    universal_gas_constant:
        Passed through to En_IG / S_IG where applicable.
    coeffs:
        Same coefficient names expected by En_IG and S_IG:
        - NASA9: a1..a7, b1 (enthalpy), b2 (entropy)
        - NASA7: a1..a7
        - Shomate: A..G

    Returns
    -------
    Dict with keys:
        result: List of {value, unit, symbol} for each temperature
        message

    Notes
    -----
    - Internally forces H to J/mol and S to J/mol.K before computing G.
    - Gibbs here is the *ideal-gas* reference value (useful for ΔG° and Keq workflows).
    """
    try:
        results = []

        # NOTE: looping over temperatures
        for temperature in temperatures:
            G_res = GiFrEn_IG(
                method=method,
                temperature=temperature,
                temperature_range=temperature_range,
                output_unit=output_unit,
                universal_gas_constant=universal_gas_constant,
                message=None,
                **coeffs,
            )
            if G_res is None:
                results.append(0.0)  # or None, depending on preference
            else:
                results.append(G_res.value)

        # NOTE: collect temperature values for reference
        temperature_values = [temp.value for temp in temperatures]

        # NOTE: message
        if message is not None:
            message = "Ideal gas Gibbs free energy range calculation successful" + message
            print(message)

        return {
            "values": {
                "x": temperature_values,
                "y": results,
            },
            "unit": output_unit if output_unit is not None else "J/mol",
        }

    except Exception as e:
        logger.error(
            f"Error in ideal gas Gibbs free energy range calculation: {e}")
        return None

# NOTE: calculate Gibbs free energy change between two temperatures


def dGiFrEn_IG(
    method: G_IG_Method,
    *,
    T_initial: Temperature,
    T_final: Temperature,
    temperature_range: Optional[Tuple[Temperature, Temperature]] = None,
    output_unit: Optional[str] = None,
    universal_gas_constant: float = 8.31446261815324,
    message: Optional[str] = None,
    **coeffs: Any,
) -> Optional[CustomProp]:
    """
    Calculate ideal-gas Gibbs free energy change between two temperatures.

    dG_IG = G_IG(T_final) - G_IG(T_initial)

    Default unit: J/mol.
    """
    try:
        G_i = GiFrEn_IG(
            method,
            temperature=T_initial,
            temperature_range=temperature_range,
            output_unit=output_unit,
            universal_gas_constant=universal_gas_constant,
            message=None,
            **coeffs,
        )
        if G_i is None:
            return None

        G_f = GiFrEn_IG(
            method,
            temperature=T_final,
            temperature_range=temperature_range,
            output_unit=output_unit,
            universal_gas_constant=universal_gas_constant,
            message=None,
            **coeffs,
        )
        if G_f is None:
            return None

        dG_value = float(G_f.value - G_i.value)
        dG_unit = output_unit if output_unit is not None else "J/mol"

        # NOTE: prepare result
        result = CustomProp(
            value=dG_value,
            unit=dG_unit,
        )

        # NOTE: message
        if message is not None:
            message = "Ideal gas Gibbs free energy change calculation successful" + message
            print(message)

        return result
    except Exception as e:
        logger.error(
            f"Error in ideal gas Gibbs free energy change calculation: {e}")
        return None


# SECTION: Generic Gibbs energy identities

def _generic_temperature(
    temperature: Temperature,
    output_temperature_unit: str | None = None,
    unit_conversion_fn=None,
) -> float:
    """Return temperature value, optionally converted to the requested unit."""
    # SECTION: Resolve conversion function
    conversion_fn = pycuc.convert_from_to if unit_conversion_fn is None else unit_conversion_fn

    # SECTION: Normalize temperature only when requested
    temperature_value = float(temperature.value)
    temperature_unit = temperature.unit.strip()
    if output_temperature_unit is not None and temperature_unit != output_temperature_unit:
        temperature_value = float(
            conversion_fn(
                temperature_value,
                temperature_unit,
                output_temperature_unit,
            )
        )

    # ! Thermodynamic identities require a physically valid absolute temperature.
    validation_unit = output_temperature_unit or temperature_unit
    try:
        temperature_k = temperature_value
        if validation_unit != "K":
            temperature_k = float(conversion_fn(
                temperature_value, validation_unit, "K"))
        if temperature_k <= 0.0:
            raise ValueError(
                "temperature must be greater than zero K after conversion.")
    except Exception:
        # NOTE: If a custom unit cannot be converted to K, validate the numeric
        # value used in the T*S product directly.
        if temperature_value <= 0.0:
            raise ValueError(
                f"temperature must be greater than zero {validation_unit}.")

    return temperature_value


def calc_gibbs_energy(
    enthalpy: CustomProp | float | int,
    temperature: Temperature,
    entropy: CustomProp | float | int,
    output_enthalpy_unit: str | None = None,
    output_entropy_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn=None,
) -> float:
    """Calculate Gibbs energy from enthalpy, temperature, and entropy.

    Parameters
    ----------
    enthalpy : float | int | CustomProp
        Enthalpy on the desired amount basis.
    temperature : Temperature
        Temperature value used in the ``T*S`` product. When
        ``output_temperature_unit`` is ``None``, ``temperature.value`` is used
        as-is. When ``output_temperature_unit`` is provided, ``temperature`` is
        converted to that unit before calculation. Supported unit labels are
        determined by pycuc, commonly ``C``, ``K``, ``R``, and ``F``.
    entropy : float | int | CustomProp
        Entropy on the same amount basis as ``enthalpy`` and per the
        temperature unit used in the ``T*S`` product.
    output_enthalpy_unit : str, optional
        Unit used to normalize ``enthalpy`` before calculation.
    output_entropy_unit : str, optional
        Unit used to normalize ``entropy`` before calculation.
    output_temperature_unit : str, optional
        Unit used to normalize ``temperature`` before calculation. Leave as
        ``None`` to use ``temperature.value`` as supplied.
    unit_conversion_fn : callable, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Gibbs energy.

    Notes
    -----
    Equation
        `G = H - T*S`
    """
    # SECTION: Resolve conversion function
    conversion_fn = pycuc.convert_from_to if unit_conversion_fn is None else unit_conversion_fn

    # SECTION: Normalize enthalpy
    h = enthalpy.value if isinstance(enthalpy, CustomProp) else enthalpy
    if isinstance(enthalpy, CustomProp) and output_enthalpy_unit and enthalpy.unit != output_enthalpy_unit:
        h = conversion_fn(h, enthalpy.unit, output_enthalpy_unit)

    # SECTION: Normalize temperature
    t = _generic_temperature(
        temperature,
        output_temperature_unit,
        unit_conversion_fn,
    )

    # SECTION: Normalize entropy
    s = entropy.value if isinstance(entropy, CustomProp) else entropy
    if isinstance(entropy, CustomProp) and output_entropy_unit and entropy.unit != output_entropy_unit:
        s = conversion_fn(s, entropy.unit, output_entropy_unit)

    # SECTION: Calculate Gibbs energy
    return float(h) - t * float(s)


def calc_gibbs_energy_change(
    enthalpy_change: CustomProp | float | int,
    entropy_change: CustomProp | float | int,
    temperature: Temperature,
    output_enthalpy_change_unit: str | None = None,
    output_entropy_change_unit: str | None = None,
    output_temperature_unit: str | None = None,
    unit_conversion_fn=None,
) -> float:
    """Calculate Gibbs energy change at a common temperature.

    Parameters
    ----------
    enthalpy_change : float | int | CustomProp
        Enthalpy change on the desired amount basis.
    entropy_change : float | int | CustomProp
        Entropy change on the same amount basis and per the temperature unit
        used in the ``T*dS`` product.
    temperature : Temperature
        Temperature value used in the ``T*dS`` product. When
        ``output_temperature_unit`` is ``None``, ``temperature.value`` is used
        as-is. When provided, ``temperature`` is converted to that unit before
        calculation.
    output_enthalpy_change_unit : str, optional
        Unit used to normalize ``enthalpy_change`` before calculation.
    output_entropy_change_unit : str, optional
        Unit used to normalize ``entropy_change`` before calculation.
    output_temperature_unit : str, optional
        Unit used to normalize ``temperature`` before calculation. Leave as
        ``None`` to use ``temperature.value`` as supplied.
    unit_conversion_fn : callable, optional
        Unit conversion function. Defaults to ``pycuc.convert_from_to``.

    Returns
    -------
    float
        Gibbs energy change.

    Notes
    -----
    Equation
        `dG = dH - T*dS`
    """
    # SECTION: Delegate to the generic identity
    return calc_gibbs_energy(
        enthalpy=enthalpy_change,
        temperature=temperature,
        entropy=entropy_change,
        output_enthalpy_unit=output_enthalpy_change_unit,
        output_entropy_unit=output_entropy_change_unit,
        output_temperature_unit=output_temperature_unit,
        unit_conversion_fn=unit_conversion_fn,
    )
