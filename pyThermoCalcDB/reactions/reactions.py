# import libs
import logging
from typing import Literal, Optional, List, Dict, Any
from pyreactlab_core.models.reaction import Reaction
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Temperature, CustomProp, ComponentKey
from pythermodb_settings.utils import measure_time
from pyThermoLinkDB.thermo import Source
import pycuc
# locals
from ..core.hsg_reaction import HSGReaction

# NOTE: logger
logger = logging.getLogger(__name__)


@measure_time
def dH_rxn_298(
    reaction: Reaction,
    model_source: ModelSource,
    **kwargs
) -> Optional[CustomProp]:
    """
    Calculate the standard enthalpy of reaction at 298.15 K using HSG properties.

    Parameters
    ----------
    reaction : Reaction
        The Reaction object representing the chemical reaction.
    model_source : ModelSource
        The ModelSource object representing the data source for the components.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[CustomProp]
        The standard enthalpy of reaction at 298.15 K in J/mol, or None if calculation fails.
    """
    try:
        # NOTE: calculate standard enthalpy of reaction at 298.15 K
        temperature = Temperature(value=298.15, unit='K')

        dH_rxn_298 = dH_rxn_STD(
            reaction=reaction,
            temperature=temperature,
            model_source=model_source,
            **kwargs
        )

        return dH_rxn_298
    except Exception as e:
        logger.error(
            f"Error calculating standard enthalpy of reaction at 298.15 K: {e}")
        return None


@measure_time
def dH_rxn_STD(
    reaction: Reaction,
    temperature: Temperature,
    model_source: ModelSource,
    **kwargs
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
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[CustomProp]
        The standard enthalpy of reaction in J/mol, or None if calculation fails.

    Notes
    -----
    The standard enthalpy of reaction is calculated for the reaction as written,
    using each species in its specified phase.

    For a reaction:
        aA(x) + bB(x) <=> cC(x) + dD(x)

    the standard reaction enthalpy is:
        ΔH°_rxn(T) = Σ(ν_i * H_i°(T, phase_i))

    where:
        - ν_i is the stoichiometric coefficient
        - H_i°(T, phase_i) is the standard molar enthalpy of species i in its
        stated phase at temperature T

    For ideal gas enthalpy (H_IG) is calculated as:
        H_IG(T) = ΔH_f(298.15) + ∫(Cp_IG dT) from 298.15 K to T

    If phase-specific enthalpy data are unavailable, they may be approximated
    from ideal-gas data:

        H_LIQ(T) ≈ H_IG(T) - ΔH_vap(T)
        H_SOL(T) ≈ H_IG(T) - ΔH_sub(T) or H_IG(T) - ΔH_vap(T) - ΔH_fusion(T)

    Therefore, for a liquid reactant such as CH4(l), vaporization correction is
    required when only ideal-gas enthalpy data are available.

    This method returns the reaction enthalpy for the reaction as written, not
    the reaction enthalpy on a pure ideal-gas reference basis.
    """
    try:
        # SECTION: Prepare source
        # ! component key is set to 'Name-State'
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


@measure_time
def dG_rxn_STD(
    reaction: Reaction,
    temperature: Temperature,
    model_source: ModelSource,
    **kwargs
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
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[CustomProp]
        The standard Gibbs free energy of reaction in J/mol, or None if calculation fails.

    Notes
    -----
    The standard Gibbs free energy of reaction is calculated using ideal-gas
    standard-state Gibbs energies for all species.

    For a reaction:
        aA(x) + bB(x) <=> cC(x) + dD(x)

    the standard Gibbs free energy of reaction is:
        ΔG°_rxn(T) = Σ(ν_i * G_i°,IG(T))

    where:
        - ν_i is the stoichiometric coefficient
        - G_i°,IG(T) is the standard molar Gibbs free energy of species i in
        the ideal-gas standard state at temperature T

    This quantity is used to calculate the equilibrium constant:
        K_eq = exp(-ΔG°_rxn / (R * T))

    Notes on phase information
    --------------------------
    The physical phases written in the reaction are not used directly in the
    Gibbs summation for ΔG°_rxn. All species are evaluated in the ideal-gas
    standard state.

    If only condensed-phase Gibbs data are available for a species, they must
    first be converted to the ideal-gas standard state before use.
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
        dG_rxn_std = hsg_reaction.calc_standard_rxn_gibbs_free_energy(
            temperature=temperature,
        )

        return dG_rxn_std
    except Exception as e:
        logger.error(
            f"Error calculating standard Gibbs free energy of reaction: {e}")
        return None


@measure_time
def Keq_STD(
    reaction: Reaction,
    temperature: Temperature,
    model_source: ModelSource,
    component_key: ComponentKey = 'Name-Formula',
    **kwargs
) -> Optional[CustomProp]:
    """
    Calculate the standard equilibrium constant of reaction using HSG properties using Van't Hoff equation as:
        Keq_std = exp(-ΔG_rxn_std / (R * T))

    Parameters
    ----------
    reaction : Reaction
        The Reaction object representing the chemical reaction.
    temperature : Temperature
        The Temperature object representing the temperature at which to calculate the equilibrium constant of reaction.
    model_source : ModelSource
        The ModelSource object representing the data source for the components.
    component_key : ComponentKey, optional
        The ComponentKey to identify components, by default 'Name-Formula'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[CustomProp]
        The standard equilibrium constant of reaction, or None if calculation fails.
    """
    try:
        # SECTION: Prepare source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # NOTE: initialize HSGReaction
        hsg_reaction = HSGReaction(
            reaction=reaction,
            source=Source_,
        )

        # >> calculate standard equilibrium constant of reaction
        Keq_std = hsg_reaction.calc_equilibrium_constant_STD(
            temperature=temperature,
        )

        return Keq_std
    except Exception as e:
        logger.error(
            f"Error calculating standard equilibrium constant of reaction: {e}")
        return None


@measure_time
def Keq_VH(
    reaction: Reaction,
    temperature: Temperature,
    model_source: ModelSource,
    component_key: ComponentKey = 'Name-Formula',
    **kwargs
) -> Optional[CustomProp]:
    """
    Calculates the equilibrium constant of a reaction at a given temperature using Van't Hoff equation as:

            ln(K_T) = ln(K_ref) + (1/R) ∫(ΔH°(T) / T²) dT from T_ref to T

        where:
        - K_T is the equilibrium constant at temperature T
        - K_ref is the equilibrium constant at reference temperature T_ref
        - ΔH°(T) is the standard enthalpy change of the reaction at temperature T
        - R is the universal gas constant

    Parameters
    ----------
    reaction : Reaction
        The Reaction object representing the chemical reaction.
    temperature : Temperature
        The Temperature object representing the temperature at which to calculate the equilibrium constant of reaction.
    model_source : ModelSource
        The ModelSource object representing the data source for the components.
    component_key : ComponentKey, optional
        The ComponentKey to identify components, by default 'Name-Formula'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[CustomProp]
        The equilibrium constant of reaction at the given temperature, or None if calculation fails.
    """
    try:
        # SECTION: Prepare source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # NOTE: initialize HSGReaction
        hsg_reaction = HSGReaction(
            reaction=reaction,
            source=Source_,
        )

        # >> calculate equilibrium constant of reaction using Van't Hoff equation
        Keq_T = hsg_reaction.calc_equilibrium_constant_vh(
            temperature=temperature,
        )

        return Keq_T
    except Exception as e:
        logger.error(
            f"Error calculating equilibrium constant of reaction using Van't Hoff equation: {e}")
        return None


@measure_time
def Keq(
    dG_rxn_STD: CustomProp,
    temperature: Temperature,
    **kwargs
) -> Optional[CustomProp]:
    """
    Calculate the equilibrium constant of reaction at a given temperature as:
        Keq = exp(-ΔG_rxn_std / (R * T))

    Parameters
    ----------
    dG_rxn_STD : CustomProp
        The standard Gibbs free energy of reaction.
    temperature : Temperature
        The Temperature object representing the temperature at which to calculate the equilibrium constant of reaction.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[CustomProp]
        The equilibrium constant of reaction at the given temperature, or None if calculation fails.
    """
    try:
        # SECTION: input validation
        # NOTE: dGiFrEn_std
        # ! [J/mol]
        if dG_rxn_STD.unit != 'J/mol':
            dGiFrEn_std_value = pycuc.convert_from_to(
                value=dG_rxn_STD.value,
                from_unit=dG_rxn_STD.unit,
                to_unit='J/mol'
            )
        else:
            dGiFrEn_std_value = dG_rxn_STD.value

        # NOTE: temperature
        # ! [K]
        if temperature.unit != 'K':
            temperature_value = pycuc.convert_from_to(
                value=temperature.value,
                from_unit=temperature.unit,
                to_unit='K'
            )
        else:
            temperature_value = temperature.value

        # SECTION: set equation
        Keq_eq = HSGReaction.Keq

        # >> calculate equilibrium constant of reaction using Van't Hoff equation
        return Keq_eq(
            gibbs_energy_of_reaction_std=dGiFrEn_std_value,
            temperature=temperature_value,
        )
    except Exception as e:
        logger.error(
            f"Error calculating equilibrium constant of reaction: {e}")
        return None


@measure_time
def Keq_VH_Shortcut(
    Keq_std: CustomProp,
    dEn_rxn_std: CustomProp,
    temperature: Temperature,
    **kwargs
) -> Optional[CustomProp]:
    """
    Shortcut for Van't Hoff equation to calculate equilibrium constant at different temperatures as:
        Keq = Keq_std * exp( (-ΔH_rxn_std / R) * (1/T - 1/T_ref) )

    where, Keq_std is the equilibrium constant at standard conditions (T_ref), and ΔH_rxn_std is the enthalpy of reaction at standard conditions.
    Tref is set to 298.15 K.

    Parameters
    ----------
    Keq_std : CustomProp
        The standard equilibrium constant of reaction.
    dEn_rxn_std : CustomProp
        The standard enthalpy of reaction.
    temperature : Temperature
        The Temperature object representing the temperature at which to calculate the equilibrium constant of reaction.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[CustomProp]
        The equilibrium constant of reaction at the given temperature, or None if calculation fails.
    """
    try:
        # SECTION: input validation
        # NOTE: Keq_std
        # ! [dimensionless]
        if Keq_std.unit != 'dimensionless':
            logger.error(
                f"Error: Keq_std must be dimensionless, got {Keq_std.unit} instead.")
            return None
        else:
            Keq_std_value = Keq_std.value

        # NOTE: dEn_rxn_std
        # ! [J/mol]
        if dEn_rxn_std.unit != 'J/mol':
            dEn_rxn_std_value = pycuc.convert_from_to(
                value=dEn_rxn_std.value,
                from_unit=dEn_rxn_std.unit,
                to_unit='J/mol'
            )
        else:
            dEn_rxn_std_value = dEn_rxn_std.value

        # NOTE: temperature
        # ! [K]
        if temperature.unit != 'K':
            temperature_value = pycuc.convert_from_to(
                value=temperature.value,
                from_unit=temperature.unit,
                to_unit='K'
            )
        else:
            temperature_value = temperature.value

        # SECTION: vh shortcut equation
        vh_shortcut = HSGReaction.vh_shortcut

        # >> calculate equilibrium constant of reaction using Van't Hoff shortcut equation
        return vh_shortcut(
            enthalpy_of_reaction_std=dEn_rxn_std_value,
            equilibrium_constant_std=Keq_std_value,
            temperature=temperature_value,
        )
    except Exception as e:
        logger.error(
            f"Error calculating equilibrium constant of reaction using Van't Hoff shortcut equation: {e}")
        return None
