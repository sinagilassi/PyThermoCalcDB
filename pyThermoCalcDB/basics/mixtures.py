# import libs
import logging
from typing import Dict, List, Optional
from pythermodb_settings.models import ComponentAmounts, Component, ComponentKey
from pythermodb_settings.utils import to_amounts
from pycuc import convert_from_to
# ! locals
from ..utils.conversions import _to_molecular_weight


# NOTE: logger setup
logger = logging.getLogger(__name__)


# ! Mixture Molecular Weight
def calc_mixture_molecular_weight_1(
        component_mole_fractions: List[float],
        component_molecular_weights: List[float],
        input_unit: Optional[str] = None,
        output_unit: Optional[str] = None,
) -> float:
    """Calculate a mixture molecular weight from mole fractions and weights.

    Parameters
    ----------
    component_mole_fractions : list of float
        Mole fractions for the components.
    component_molecular_weights : list of float
        Molecular weights for the components.
    input_unit : str, optional
        Unit of the molecular weights, if conversion is needed.
    output_unit : str, optional
        Desired unit for the returned molecular weight.

    Returns
    -------
    float
        The mixture molecular weight, optionally converted to ``output_unit``.
    """
    mw_mix = sum(
        x * y for x, y in zip(component_mole_fractions, component_molecular_weights)
    )

    # >> check if conversion is needed
    if input_unit and output_unit:
        mw_mix = convert_from_to(
            value=mw_mix,
            from_unit=input_unit,
            to_unit=output_unit
        )

    return mw_mix

# ! Mixture Molecular Weight with desired output unit


def calc_mixture_molecular_weight_2(
        component_mole_fractions: ComponentAmounts,
        component_molecular_weights: ComponentAmounts,
        output_unit: str = "g/mol"
):
    """Calculate a mixture molecular weight from keyed component amounts.

    Parameters
    ----------
    component_mole_fractions : ComponentAmounts
        Keyed mole fractions for the components.
    component_molecular_weights : ComponentAmounts
        Keyed molecular weights for the components.
    output_unit : str, optional
        Unit to which the molecular weights are converted.

    Returns
    -------
    float
        The mixture molecular weight in ``output_unit``.
    """
    # SECTION: Convert component amounts and molecular weights to float values
    # >>> mole fraction
    component_mole_fractions_float = to_amounts(component_mole_fractions)
    # >>> molecular weight
    component_molecular_weights_float = to_amounts(
        component_amounts=component_molecular_weights,
        output_unit=output_unit,
        unit_conversion_fn=convert_from_to
    )

    return sum(
        component_mole_fractions_float[key] *
        component_molecular_weights_float[key]
        for key in component_mole_fractions_float
    )
