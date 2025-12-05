# import libs
import logging
from typing import Literal, Optional, Tuple
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component
from pythermodb_settings.utils import set_component_id
# local
from ..thermo import Source
from ..utils.tools import measure_time
from ..models import ComponentEquation
from ..core.x_prop import XProp

# NOTE: Logger
logger = logging.getLogger(__name__)

# SECTION: Equation Maker


@measure_time
def mkeq(
    prop_name: str,
    component: Component,
    model_source: ModelSource,
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = 'Name-State',
) -> Optional[ComponentEquation]:
    """
    Make an equation object for property calculations.

    Parameters
    ----------
    prop_name : str
        The ID of the equation to be used for calculations, e.g., 'VaPr', 'Cp_IG'.
    component : Component
        The chemical component for which properties are to be calculated.
    model_source : ModelSource
        The source containing data for calculations.
    component_key : Literal
        The key to identify the component in the source data. Defaults to 'Name-State'.

    Returns
    -------
    XProp
        An instance of the XProp class initialized with the provided parameters.
    """
    try:
        # SECTION: Validate inputs
        if not prop_name:
            logger.error("prop_name, equation ID is required.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        # NOTE: set component id
        component_id = set_component_id(
            component=component,
            component_key=component_key
        )

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # NOTE: search for equation id in source
        if Source_.is_prop_eq_available(
            component_id=component_id,
            prop_name=prop_name,
        ) is False:
            logger.error(
                f"Equation ID '{prop_name}' not found for component '{component_id}'.")
            return None

        # SECTION: Create XProp object
        XProp_ = XProp(
            prop_name=prop_name,
            component=component,
            source=Source_,
            component_key=component_key,
        )

        # SECTION: extract info
        equation = XProp_.equation
        body = XProp_.body
        args = XProp_.args
        arg_symbols = XProp_.arg_symbols
        returns = XProp_.returns
        return_symbol = XProp_.return_symbol

        # SECTION: create ComponentEquation object
        component_equation = ComponentEquation(
            equation=equation,
            body=body,
            args=args,
            arg_symbols=arg_symbols,
            returns=returns,
            return_symbol=return_symbol
        )

        return component_equation
    except Exception as e:
        logger.error(f"Error creating equation: {e}")
        return None
