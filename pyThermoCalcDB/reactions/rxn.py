# import libs
import logging
from typing import List, Optional, Dict, Any
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Temperature, CustomProp
import pycuc
from math import exp
# locals
from ..configs.constants import (
    R_J_molK,
    T_ref_K,
)
from ..utils.conversions import _to_J__mol
from .reactions import Keq, Keq_VH_Shortcut

# NOTE: logger setup
LOGGER = logging.getLogger(__name__)


class RXN:
    """
    Reaction (RXN) class for handling chemical reactions and their thermodynamic properties as:
    - Standard enthalpy change (dH_rxn_STD)
    - Standard Gibbs free energy change (dG_rxn_STD)
    - Standard entropy change (dS_rxn_STD)
    - Equilibrium constant (Keq)
    - Temperature-dependent equilibrium constant (Keq_T)

    Notes
    -----
    Key Reaction attributes:
    - reaction : str
        The reaction equation string (e.g., "CO2(g) + 3H2(g) => CH3OH(g) + H2O(g)")
    - reactants : List[Dict]
        List of reactant dictionaries with coefficient, molecule, state, and molecule_state
    - products : List[Dict]
        List of product dictionaries with coefficient, molecule, state, and molecule_state
    - reaction_coefficients : float
        Overall stoichiometric coefficients
    - reaction_stoichiometry : Dict[str, float]
        Stoichiometric coefficients for each component (negative for reactants, positive for products)
    - state_count : Dict[str, int]
        Count of components in each state (g, l, aq, s)
    - reaction_phase : str
        Primary phase of the reaction (e.g., 'gas', 'liquid')
    - reaction_state : Dict[str, str]
        State of each component in the reaction
    - reactants_names : List[str]
        List of reactant component IDs
    - products_names : List[str]
        List of product component IDs
    - all_components : List[str]
        List of all component IDs involved in the reaction
    - available_components : List[Component]
        List of available Component objects for the reaction
    - component_checker: bool
        Flag indicating if all components are available
    - map_components: Dict[str, Component]
        Mapping of component IDs to Component objects
    """
    # SECTION: Attributes
    # NOTE: universal gas constant [J/mol.K]
    __R: float = R_J_molK
    # NOTE: reference temperature [K]
    __T_Ref_K: float = T_ref_K

    def __init__(
            self,
            reaction: Reaction,
    ) -> None:
        # NOTE: reaction
        self.reaction: Reaction = reaction

        # SECTION: set components
        self.components = reaction.available_components
        # checker
        if reaction.component_checker is False:
            raise ValueError(
                "Some components in the reaction are not available in the provided components list."
            )

        # NOTE: map components
        self.map_components = reaction.map_components

        # NOTE: all components
        self.all_components: List[str] = reaction.all_components

    def dH_rxn_STD(
        self,
        H_i_IG: Dict[str, CustomProp],
        **kwargs
    ) -> Optional[CustomProp]:
        """
        Calculate the standard enthalpy change of the reaction as:

            dH_rxn = Σ(ν_i * H_IG(T))

        Parameters
        ----------
        H_i_IG : Dict[str, CustomProp]
            A dictionary of standard enthalpies of formation for each component.

        Returns
        -------
        Optional[CustomProp]
            The standard enthalpy change of the reaction.
        """
        try:
            # NOTE: initialize reaction enthalpy
            reaction_enthalpy = 0.0

            # SECTION: calculate reaction enthalpy
            for component_id, coeff in self.reaction.reaction_stoichiometry.items():
                if component_id in H_i_IG:
                    # NOTE: Get component enthalpy value and unit
                    H_component_value = H_i_IG[component_id].value
                    H_component_unit = H_i_IG[component_id].unit

                    # NOTE: Convert to J/mol if necessary
                    if H_component_unit != "J/mol":
                        H_component_value = _to_J__mol(
                            value=H_component_value,
                            from_unit=H_component_unit,
                        )

                    # Accumulate reaction enthalpy
                    reaction_enthalpy += coeff * H_component_value
                else:
                    raise KeyError(
                        f"Component ID '{component_id}' not found in H_i_IG dictionary.")

            return CustomProp(
                value=reaction_enthalpy,
                unit="J/mol",
            )
        except Exception as e:
            LOGGER.error(f"Error in dH_STD calculation: {e}")
            return None

    def dG_rxn_STD(
        self,
        G_i_IG: Dict[str, CustomProp],
        **kwargs
    ) -> Optional[CustomProp]:
        """
        Calculate the standard Gibbs free energy change of the reaction as:

            dG_rxn = Σ(ν_i * G_IG(T))

        Parameters
        ----------
        G_i_IG : Dict[str, CustomProp]
            A dictionary of standard Gibbs free energies of formation for each component.

        Returns
        -------
        Optional[CustomProp]
            The standard Gibbs free energy change of the reaction.
        """
        try:
            # NOTE: initialize reaction Gibbs free energy
            reaction_gibbs = 0.0

            # SECTION: calculate reaction Gibbs free energy
            for component_id, coeff in self.reaction.reaction_stoichiometry.items():
                if component_id in G_i_IG:
                    # NOTE: Get component Gibbs value and unit
                    G_component_value = G_i_IG[component_id].value
                    G_component_unit = G_i_IG[component_id].unit

                    # NOTE: Convert to J/mol if necessary
                    if G_component_unit != "J/mol":
                        G_component_value = _to_J__mol(
                            value=G_component_value,
                            from_unit=G_component_unit,
                        )

                    # Accumulate reaction Gibbs free energy
                    reaction_gibbs += coeff * G_component_value
                else:
                    raise KeyError(
                        f"Component ID '{component_id}' not found in G_i_IG dictionary.")

            return CustomProp(
                value=reaction_gibbs,
                unit="J/mol",
            )
        except Exception as e:
            LOGGER.error(f"Error in dG_STD calculation: {e}")
            return None

    def dS_rxn_STD(
        self,
        S_i_IG: Dict[str, CustomProp],
        **kwargs
    ) -> Optional[CustomProp]:
        """
        Calculate the standard entropy change of the reaction as:

            dS_rxn = Σ(ν_i * S_IG(T))

        Parameters
        ----------
        S_i_IG : Dict[str, CustomProp]
            A dictionary of standard entropies for each component.

        Returns
        -------
        Optional[CustomProp]
            The standard entropy change of the reaction.
        """
        try:
            # NOTE: initialize reaction entropy
            reaction_entropy = 0.0

            # SECTION: calculate reaction entropy
            for component_id, coeff in self.reaction.reaction_stoichiometry.items():
                if component_id in S_i_IG:
                    # NOTE: Get component entropy value and unit
                    S_component_value = S_i_IG[component_id].value
                    S_component_unit = S_i_IG[component_id].unit

                    # NOTE: Convert to J/mol.K if necessary
                    if S_component_unit != "J/mol.K":
                        S_component_value = pycuc.convert_from_to(
                            value=S_component_value,
                            from_unit=S_component_unit,
                            to_unit="J/mol.K",
                        )

                    # Accumulate reaction entropy
                    reaction_entropy += coeff * S_component_value
                else:
                    raise KeyError(
                        f"Component ID '{component_id}' not found in S_i_IG dictionary.")

            return CustomProp(
                value=reaction_entropy,
                unit="J/mol.K",
            )
        except Exception as e:
            LOGGER.error(f"Error in dS_STD calculation: {e}")
            return None

    def Keq(
        self,
        dG_rxn_STD: CustomProp,
        temperature: Temperature,
        **kwargs
    ) -> Optional[CustomProp]:
        """
        Calculate the equilibrium constant of the reaction as:

            Keq = exp(-dG_rxn_STD / (R * T))

        where:
        - R is the universal gas constant [J/mol.K]
        - T is the temperature in Kelvin [K]
        - dG_rxn_STD is the standard Gibbs free energy change of the reaction [J/mol]

        If dG_rxn_STD is not in J/mol, and if temperature is not in K, it will be converted.

        Parameters
        ----------
        dG_rxn_STD : CustomProp
            The standard Gibbs free energy change of the reaction.
        temperature : Temperature
            The temperature at which to calculate the equilibrium constant.

        Returns
        -------
        Optional[CustomProp]
            The equilibrium constant of the reaction.
        """
        try:
            # Calculate Keq
            return Keq(
                dGiFrEn_std=dG_rxn_STD,
                temperature=temperature
            )
        except Exception as e:
            LOGGER.error(f"Error in Keq calculation: {e}")
            return None

    def Keq_vh(
        self,
        Keq_STD: CustomProp,
        dH_rxn_STD: CustomProp,
        temperature: Temperature,
        **kwargs
    ):
        pass

    def Keq_vh_shortcut(
        self,
        Keq_std: CustomProp,
        dH_rxn_STD: CustomProp,
        temperature: Temperature,
    ):
        """
        Shortcut for Van't Hoff equation to calculate equilibrium constant at different temperatures as:
        Keq = Keq_std * exp( (-ΔH_rxn_std / R) * (1/T - 1/T_ref) )

        where, Keq_std is the equilibrium constant at standard conditions (T_ref), and ΔH_rxn_std is the enthalpy of reaction at standard conditions.
        Tref is set to 298.15 K.

        Parameters
        ----------
        Keq_std : CustomProp
            The standard equilibrium constant of reaction.
        dH_rxn_STD : CustomProp
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
            # >> calculate equilibrium constant of reaction using Van't Hoff shortcut equation
            return Keq_VH_Shortcut(
                Keq_std=Keq_std,
                dEn_rxn_std=dH_rxn_STD,
                temperature=temperature,
            )
        except Exception as e:
            LOGGER.error(f"Error in Keq_T calculation: {e}")
            return None
