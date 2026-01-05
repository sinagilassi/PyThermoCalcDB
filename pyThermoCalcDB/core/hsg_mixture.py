# import libs
import logging
from typing import Optional, List, Dict, Literal, cast, Any
from pythermodb_settings.models import Temperature, Pressure, Component, CustomProp
from pythermodb_settings.utils import set_component_id
import pycuc
from pyThermoLinkDB.thermo import Source
from pyThermoLinkDB.models.component_models import ComponentEquationSource
# locals
from .hsg_properties import HSGProperties

# NOTE: set up logger
logger = logging.getLogger(__name__)


class HSGMixture:
    """
    Class to represent a mixture of components and calculate its thermodynamic properties using HSG method.

    This class provides thermodynamic property calculations for multi-component mixtures using the HSG
    (Heat capacity, entropy, and Gibbs energy) methodology. It aggregates component properties and
    enables mixture-level calculations with support for custom departure and excess contributions.

    Attributes
    ----------
    components : List[Component]
        A list of Component objects representing the mixture components.
    source : Source
        The Source object representing the data source for component properties.
    component_key : Literal['Name-State', 'Formula-State', 'Name', 'Formula', 'Name-Formula-State', 'Formula-Name-State']
        The key used to uniquely identify components in the source data (default: 'Name-State').
    component_ids : List[str]
        List of unique component identifiers generated from components using the specified component_key.
    mole_fractions : List[float]
        Mole fractions of each component in the mixture (normalized).
    hsg_properties : Dict[str, HSGProperties]
        Dictionary mapping component IDs to their corresponding HSGProperties instances for thermodynamic calculations.

    Methods
    -------
    _components_hsg_properties()
        Private method to initialize HSGProperties instances for all mixture components.
    calc_mixture_enthalpy(temperature, phase, departure_enthalpy, excess_enthalpy)
        Calculate the mixture enthalpy (J/mol) at a given temperature and phase with optional contributions.
    """

    def __init__(
        self,
        components: List[Component],
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
        source : Source
            The Source object representing the data source for the components.
        component_key : Literal[...], optional
            The key to identify components, by default 'Name-State'.
        """
        # NOTE: set attributes
        self.components = components
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

    def calc_mixture_enthalpy(
            self,
            temperature: Temperature,
            phase: Literal['IG', 'LIQ'] = 'IG',
            departure_enthalpy: Optional[CustomProp] = None,
            excess_enthalpy: Optional[CustomProp] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Calculate the mixture enthalpy (J/mol) at a given temperature (K) and pressure for specified phase.

        Parameters
        ----------
        temperature : Temperature
            The temperature of the mixture.
        phase : Literal['IG', 'LIQ', 'SOL'], optional
            The phase of the mixture ('IG' for ideal gas, 'LIQ' for liquid), by default 'IG'.
        departure_enthalpy : Optional[CustomProp], optional
            The departure enthalpy contribution, by default None.
        excess_enthalpy : Optional[CustomProp], optional
            The excess enthalpy contribution, by default None.

        Returns
        -------
        Optional[Dict[str, Any]]
            A dictionary containing the mixture enthalpy value and unit, or None if calculation fails.

        Notes
        -----
        - This method calculates the mixture enthalpy based on the specified phase and contributions.
        - The calculation may involve ideal and non-ideal contributions depending on the phase.
        - The result is provided in J/mol.
        - R is the universal gas constant (8.3145 J/mol.K).
        - Global reference is ideal gas state at 298.15 K and 1 atm, and all elements are in their standard states (embedded in dEnFo).

        Equations
        ---------
        - The mixture enthalpy is calculated using the formula:
            H_mix = Σ(x_i * H_i) + H_departure + H_excess

        - ideal gas enthalpy (H_IG) is calculated as:
            H_IG(T) = ΔH_f(T) + ∫(Cp_IG dT) from 298.15 K to T

        - liquid enthalpy (H_LIQ) is calculated as:
            H_LIQ(T) = H_IG(T) - ΔH_vap(T)

        - solid enthalpy (H_SOL) is calculated as (not implemented):
            H_SOL(T) = H_IG(T) - ΔH_vap(T) - ΔH_fusion(T)
        """
        try:
            # SECTION: calculate mixture enthalpy
            mixture_enthalpy_value = 0.0

            # NOTE: loop through components
            for component_id, mole_fraction in zip(
                self.component_ids,
                self.mole_fractions
            ):
                # ! get component hsg properties
                hsg_prop = self.hsg_properties.get(component_id)
                if hsg_prop is None:
                    logger.error(
                        f"HSG properties not found for component ID: {component_id}")
                    return None

                # NOTE: calculate phase enthalpy for component
                # ! [J/mol]
                component_enthalpy = hsg_prop.calc_phase_enthalpy(
                    temperature=temperature,
                    phase=phase
                )

                # >> check if calculation was successful
                if component_enthalpy is None:
                    logger.error(
                        f"Failed to calculate enthalpy for component ID: {component_id}")
                    return None

                # NOTE: accumulate mixture enthalpy
                # ! [J/mol]
                mixture_enthalpy_value += mole_fraction * component_enthalpy.value

            # SECTION: add departure and excess enthalpy contributions if provided
            # NOTE: departure enthalpy
            if departure_enthalpy is not None:
                # get departure enthalpy value
                val_ = departure_enthalpy.value
                unit_ = departure_enthalpy.unit

                # ! [J/mol]
                if unit_ != 'J/mol':
                    # >> convert to J/mol
                    val_ = pycuc.convert_from_to(
                        value=val_,
                        from_unit=unit_,
                        to_unit='J/mol'
                    )

                mixture_enthalpy_value += val_

            # NOTE: excess enthalpy
            if excess_enthalpy is not None:
                # get excess enthalpy value
                val_ = excess_enthalpy.value
                unit_ = excess_enthalpy.unit

                # ! [J/mol]
                if unit_ != 'J/mol':
                    # >> convert to J/mol
                    val_ = pycuc.convert_from_to(
                        value=val_,
                        from_unit=unit_,
                        to_unit='J/mol'
                    )

                mixture_enthalpy_value += val_

            return {
                'value': mixture_enthalpy_value,
                'unit': 'J/mol'
            }
        except Exception as e:
            logger.error(
                f"Error calculating mixture enthalpy: {e}")
            return None
