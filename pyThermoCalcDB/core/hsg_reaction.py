# import libs
import logging
from typing import Optional, List, Dict, Literal, cast, Any
from pythermodb_settings.models import Temperature, Pressure, Component, CustomProp
from pythermodb_settings.utils import set_component_id
import pycuc
from pyThermoLinkDB.thermo import Source
from pyThermoLinkDB.models.component_models import ComponentEquationSource
from pyreactlab_core.models.reaction import Reaction
# locals
from .hsg_properties import HSGProperties
from ..utils.component_tools import map_component_state
# locals

# NOTE: logger setup
logger = logging.getLogger(__name__)


class HSGReaction:
    """
    HSGReaction class to handle chemical reactions using HSG properties. It is used to calculate thermodynamic properties of reactions based on the HSG method as:
    - Enthalpy of reaction
    - Gibbs free energy of reaction

    Attributes
    ----------
    components : List[Component]
        A list of Component objects representing the mixture components.
    reaction : Reaction
        The Reaction object representing the chemical reaction.
    source : Source
        The Source object representing the data source for the components.
    component_key : Literal[...], optional
        The key to identify components, by default 'Name-State'.
    component_ids : List[str]
        A list of component IDs derived from the components.
    mole_fractions : List[float]
        A list of mole fractions for each component in the mixture.

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
    """

    def __init__(
        self,
        components: List[Component],
        reaction: Reaction,
        source: Source,
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
    ):
        """
        Initialize the HSGMixture with a list of components and optional custom properties.

        Parameters
        ----------
        components : List[Component]
            A list of Component objects representing the mixture components.
        reaction : Reaction
            The Reaction object representing the chemical reaction.
        source : Source
            The Source object representing the data source for the components.
        component_key : Literal[...], optional
            The key to identify components, by default 'Name-State'.
        """
        # NOTE: set attributes
        self.components = components
        self.reaction = reaction
        self.source = source
        self.component_key = component_key

        # SECTION: set component IDs
        self.component_ids = [
            set_component_id(
                component=component,
                component_key=self.component_key
            )
            for component in self.components
        ]

        # NOTE: set component phase states
        self.component_states = [
            comp.state for comp in self.components
        ]

        # NOTE: mole fraction components
        self.mole_fractions = [
            c.mole_fraction for c in self.components
        ]

        # SECTION: hsg properties instances
        self.hsg_properties = self._components_hsg_properties()

    def _components_hsg_properties(self) -> Dict[str, HSGProperties]:
        """
        Get the HSGProperties instances for all components in the mixture.

        Returns
        -------
        Dict[str, HSGProperties]
            A dictionary mapping component IDs to their corresponding HSGProperties instances.
        """
        try:
            # NOTE: initialize hsg properties dict
            hsg_properties = {}

            # SECTION: build hsg properties for components
            for component_id, component in zip(self.component_ids, self.components):
                if component_id not in hsg_properties:
                    hsg_properties[component_id] = HSGProperties(
                        component=component,
                        source=self.source,
                        component_key=cast(Literal[
                            'Name-State',
                            'Formula-State',
                            'Name',
                            'Formula',
                            'Name-Formula-State',
                            'Formula-Name-State'
                        ], self.component_key)
                    )

            # >> hsg properties
            return hsg_properties
        except Exception as e:
            logger.error(
                f"Error in getting HSG properties for components: {e}")
            return {}

    def calc_standard_reaction_enthalpy(
        self,
        temperature: Temperature,
    ) -> Optional[Dict[str, Any]]:
        """
        Calculate the standard enthalpy of reaction at the specified temperature and pressure using HSG properties.

            dH_rxn = Σ(ν_i * H_IG(T))

        Parameters
        ----------
        temperature : Temperature
            The temperature at which to calculate the reaction enthalpy.
        rxn_departure_enthalpy : Optional[CustomProp], optional
            The departure enthalpy of the reaction, by default None.
        rxn_excess_enthalpy : Optional[CustomProp], optional
            The excess enthalpy of the reaction, by default None.

        Returns
        -------
        Optional[Dict[str, Any]]
            A dictionary containing the reaction enthalpy and related information, or None if calculation fails.

        Notes
        -----
        - The reaction enthalpy is calculated using the standard enthalpies of formation of the reactants and products, adjusted for departure and excess enthalpies if provided.
        - The calculation may involve ideal and non-ideal contributions depending on the phase.
        - The result is provided in J/mol.
        - R is the universal gas constant (8.3145 J/mol.K).
        - Global reference is ideal gas state at 298.15 K and 1 atm, and all elements are in their standard states (embedded in dEnFo).

        Equations
        ---------
        For a reaction:
            aA + bB => cC + dD
        - The reaction enthalpy (ΔH_rxn) is calculated as:
            dH_rxn = Σ(ν_i * H_IG(T))

        - ideal gas enthalpy (H_IG) is calculated as:
            H_IG(T) = ΔH_f(T) + ∫(Cp_IG dT) from 298.15 K to T

        - liquid enthalpy (H_LIQ) is calculated as:
            H_LIQ(T) = H_IG(T) - ΔH_vap(T)

        - solid enthalpy (H_SOL) is calculated as (not implemented):
            H_SOL(T) = H_IG(T) - ΔH_vap(T) - ΔH_fusion(T)
        """
        try:
            # NOTE: initialize reaction enthalpy
            reaction_enthalpy = 0.0

            # SECTION: calculate reaction enthalpy
            for component_id, coeff in self.reaction.reaction_stoichiometry.items():
                # >> get hsg properties
                hsg_prop = self.hsg_properties.get(component_id, None)
                if hsg_prop is None:
                    logger.error(
                        f"HSG properties not found for component ID: {component_id}"
                    )
                    return None

                # NOTE: get component state
                component_state = self.reaction.reaction_state.get(
                    component_id,
                    None
                )
                if component_state is None:
                    logger.warning(
                        f"Component state not found for component ID: {component_id}. Using default state from component."
                    )
                    return None

                # Cast phase to Literal['IG', 'LIQ', 'SOL']
                phase = cast(Literal['IG', 'LIQ', 'SOL'], "IG")

                # >> calculate component enthalpy
                # NOTE: calculate phase enthalpy
                # ! [J/mol]
                component_enthalpy = hsg_prop.calc_phase_enthalpy(
                    temperature=temperature,
                    phase=phase,
                )
                # >> check component enthalpy
                if component_enthalpy is None:
                    raise ValueError(
                        f"Failed to calculate enthalpy for component ID: {component_id}"
                    )

                # >> accumulate reaction enthalpy
                reaction_enthalpy += coeff * component_enthalpy.value

            # >> return reaction enthalpy
            return {
                'value': reaction_enthalpy,  # J/mol
                'unit': 'J/mol'
            }
        except Exception as e:
            logger.error(f"Error in calculating reaction enthalpy: {e}")
            return None

    def calc_phase_reaction_enthalpy(
        self,
        temperature: Temperature,
        rxn_departure_enthalpy: Optional[CustomProp] = None,
        rxn_excess_enthalpy: Optional[CustomProp] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Calculate the phase (actual) enthalpy of reaction at the specified temperature and pressure using HSG properties.

            dH_rxn = Σ(ν_i * H_IG(T)) + rxn_departure_enthalpy + rxn_excess_enthalpy

        Parameters
        ----------
        temperature : Temperature
            The temperature at which to calculate the reaction enthalpy.
        rxn_departure_enthalpy : Optional[CustomProp], optional
            The departure enthalpy of the reaction, by default None.
        rxn_excess_enthalpy : Optional[CustomProp], optional
            The excess enthalpy of the reaction, by default None.

        Returns
        -------
        Optional[Dict[str, Any]]
            A dictionary containing the reaction enthalpy and related information, or None if calculation fails.

        Notes
        -----
        - The reaction enthalpy is calculated using the standard enthalpies of formation of the reactants and products, adjusted for departure and excess enthalpies if provided.
        - departure_enthalpy is added directly to the reaction enthalpy.
        - excess_enthalpy is also added directly to the reaction enthalpy.
        - The calculation may involve ideal and non-ideal contributions depending on the phase.
        - The result is provided in J/mol.
        - R is the universal gas constant (8.3145 J/mol.K).
        - Global reference is ideal gas state at 298.15 K and 1 atm, and all elements are in their standard states (embedded in dEnFo).

        Equations
        ---------
        For a reaction:
            aA + bB => cC + dD
        - The reaction enthalpy (ΔH_rxn) is calculated as:
            dH_rxn = Σ(ν_i * H_IG(T)) + Σ(ν_i * departure_enthalpy) + Σ(ν_i * excess_enthalpy)

        - ideal gas enthalpy (H_IG) is calculated as:
            H_IG(T) = ΔH_f(T) + ∫(Cp_IG dT) from 298.15 K to T

        - liquid enthalpy (H_LIQ) is calculated as:
            H_LIQ(T) = H_IG(T) - ΔH_vap(T)

        - solid enthalpy (H_SOL) is calculated as (not implemented):
            H_SOL(T) = H_IG(T) - ΔH_vap(T) - ΔH_fusion(T)
        """
        try:
            # NOTE: initialize reaction enthalpy
            reaction_enthalpy = 0.0

            # SECTION: calculate reaction enthalpy
            for component_id, coeff in self.reaction.reaction_stoichiometry.items():
                # >> get hsg properties
                hsg_prop = self.hsg_properties.get(component_id, None)
                if hsg_prop is None:
                    logger.error(
                        f"HSG properties not found for component ID: {component_id}"
                    )
                    return None

                # NOTE: get component state
                component_state = self.reaction.reaction_state.get(
                    component_id,
                    None
                )
                if component_state is None:
                    logger.warning(
                        f"Component state not found for component ID: {component_id}. Using default state from component."
                    )
                    return None

                # ! >> map component state
                phase_str = map_component_state(
                    component=hsg_prop.component
                )
                # Cast phase to Literal['IG', 'LIQ', 'SOL']
                phase = cast(Literal['IG', 'LIQ', 'SOL'], phase_str)

                # >> calculate component enthalpy
                # NOTE: calculate phase enthalpy
                # ! [J/mol]
                component_enthalpy = hsg_prop.calc_phase_enthalpy(
                    temperature=temperature,
                    phase=phase,
                )
                # >> check component enthalpy
                if component_enthalpy is None:
                    raise ValueError(
                        f"Failed to calculate enthalpy for component ID: {component_id}"
                    )

                # >> accumulate reaction enthalpy
                reaction_enthalpy += coeff * component_enthalpy.value

            # SECTION: add departure and excess enthalpy contributions if provided
            # NOTE: departure enthalpy
            if rxn_departure_enthalpy is not None:
                # get departure enthalpy value
                val_ = rxn_departure_enthalpy.value
                unit_ = rxn_departure_enthalpy.unit

                # ! [J/mol]
                if unit_ != 'J/mol':
                    # >> convert to J/mol
                    val_ = pycuc.convert_from_to(
                        value=val_,
                        from_unit=unit_,
                        to_unit='J/mol'
                    )

                reaction_enthalpy += val_

            # NOTE: excess enthalpy
            if rxn_excess_enthalpy is not None:
                # get excess enthalpy value
                val_ = rxn_excess_enthalpy.value
                unit_ = rxn_excess_enthalpy.unit

                # ! [J/mol]
                if unit_ != 'J/mol':
                    # >> convert to J/mol
                    val_ = pycuc.convert_from_to(
                        value=val_,
                        from_unit=unit_,
                        to_unit='J/mol'
                    )

                reaction_enthalpy += val_

            # >> return reaction enthalpy
            return {
                'value': reaction_enthalpy,  # J/mol
                'unit': 'J/mol'
            }
        except Exception as e:
            logger.error(f"Error in calculating reaction enthalpy: {e}")
            return None
