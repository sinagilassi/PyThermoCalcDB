# import libs
import logging
from typing import Dict, Any, Literal, Optional
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, Temperature
# local
from ..core.hsg_properties import HSGProperties
from ..thermo import Source

# NOTE: Logger
logger = logging.getLogger(__name__)


def calc_enthalpy_of_formation_at_temperature(
    component: Component,
    model_source: ModelSource,
    temperature: Temperature,
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = 'Name-State',
) -> Optional[Dict[str, Any]]:
    """
    Calculate the enthalpy of formation at a given temperature for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the enthalpy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperature : Temperature
        The temperature at which to calculate the enthalpy of formation.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.

    Returns
    -------
    Optional[Dict[str, Any]]
        A dictionary containing the enthalpy of formation value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the enthalpy of formation.
    - The enthalpy of formation symbol must be consistent with the defined unit in the configuration which is EnFo_IG.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(temperature, Temperature):
            logger.error("Invalid temperature provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        EnFo_result = hsg_props.calc_enthalpy_of_formation(
            temperature=temperature
        )

        return EnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy of formation for component '{component.name}': {e}")
        return None


def calc_gibbs_energy_of_formation_at_temperature(
        component: Component,
        model_source: ModelSource,
        temperature: Temperature,
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
) -> Optional[Dict[str, Any]]:
    """
    Calculate the Gibbs free energy of formation at a given temperature for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the Gibbs free energy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperature : Temperature
        The temperature at which to calculate the Gibbs free energy of formation.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.

    Returns
    -------
    Optional[Dict[str, Any]]
        A dictionary containing the Gibbs free energy of formation value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the Gibbs free energy of formation.
    - The Gibbs free energy of formation symbol must be consistent with the defined unit in the configuration which is GiEnFo_IG.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(temperature, Temperature):
            logger.error("Invalid temperature provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        GiEnFo_result = hsg_props.calc_gibbs_free_energy_of_formation(
            temperature=temperature
        )

        return GiEnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating Gibbs free energy of formation for component '{component.name}': {e}")
        return None
