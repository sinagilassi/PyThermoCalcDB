# import libs
from typing import Any
import pycuc
from pythermodb_settings.models import CustomProp, Temperature
# locals
from pythermocalcdb.utils.conversions import _generic_temperature


# ! ::: Calculate Gibbs energy


def _calc_gibbs_energy(
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


def _calc_gibbs_energy_change(
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
    return _calc_gibbs_energy(
        enthalpy=enthalpy_change,
        temperature=temperature,
        entropy=entropy_change,
        output_enthalpy_unit=output_enthalpy_change_unit,
        output_entropy_unit=output_entropy_change_unit,
        output_temperature_unit=output_temperature_unit,
        unit_conversion_fn=unit_conversion_fn,
    )
