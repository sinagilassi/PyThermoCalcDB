# import libs
import logging
from typing import Literal, Optional, List
from pyreactlab_core.models.reaction import Reaction
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, Temperature, Pressure, CustomProp
from pyThermoLinkDB.thermo import Source
import pycuc
# locals
from ..core.hsg_reaction import HSGReaction

# NOTE: logger
logger = logging.getLogger(__name__)


def dH_rxn_STD(
    reaction: Reaction,
    temperature: Temperature,
    model_source: ModelSource,
) -> Optional[CustomProp]:
    """
    Calculate the standard enthalpy of reaction using HSG properties.

    Parameters
    ----------
    reaction : Reaction
        The Reaction object representing the chemical reaction.
    temperature : Temperature
        The Temperature object representing the temperature at which to calculate the enthalpy of reaction.
    model_source : ModelSource
        The ModelSource object representing the data source for the components.

    Returns
    -------
    Optional[CustomProp]
        The standard enthalpy of reaction in J/mol.K, or None if calculation fails.
    """
    try:
        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # NOTE: initialize HSGReaction
        hsg_reaction = HSGReaction(
            reaction=reaction,
            source=Source_,
        )

        # >> calculate standard enthalpy of reaction
        dH_rxn_std = hsg_reaction.calc_standard_rxn_enthalpy(
            temperature=temperature,
        )

        return dH_rxn_std
    except Exception as e:
        logger.error(f"Error calculating standard enthalpy of reaction: {e}")
        return None


def dG_rxn_STD(
    reaction: Reaction,
    temperature: Temperature,
    model_source: ModelSource,
) -> Optional[CustomProp]:
    """
    Calculate the standard Gibbs free energy of reaction using HSG properties.

    Parameters
    ----------
    reaction : Reaction
        The Reaction object representing the chemical reaction.
    temperature : Temperature
        The Temperature object representing the temperature at which to calculate the Gibbs free energy of reaction.
    model_source : ModelSource
        The ModelSource object representing the data source for the components.

    Returns
    -------
    Optional[CustomProp]
        The standard Gibbs free energy of reaction in J/mol.K, or None if calculation fails.
    """
    try:
        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # NOTE: initialize HSGReaction
        hsg_reaction = HSGReaction(
            reaction=reaction,
            source=Source_,
        )

        # >> calculate standard Gibbs free energy of reaction
        dG_rxn_std = hsg_reaction.calc_standard_rxn_gibbs_energy(
            temperature=temperature,
        )

        return dG_rxn_std
    except Exception as e:
        logger.error(
            f"Error calculating standard Gibbs free energy of reaction: {e}")
        return None
