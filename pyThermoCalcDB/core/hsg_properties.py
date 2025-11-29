# import libs
import logging
from typing import Dict, Any, Literal
from pythermodb_settings.models import Component, Temperature
from pythermodb_settings.utils import set_component_id
# local
from ..thermo import Source

# NOTE: Logger
logger = logging.getLogger(__name__)


class HSGProperties:
    """
    HSG Properties class to calculate enthalpy and entropy of formation, and Gibbs free energy of formation for a given compound.

    Attributes
    ----------
    component : Component
        The chemical component for which HSG properties are to be calculated.
    source : Source
        The source containing data for calculations.
    """

    def __init__(
        self,
        component: Component,
        source: Source,
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
    ) -> None:
        """
        Initialize HSGProperties with a component and source.

        Parameters
        ----------
        component : Component
            The chemical component for which HSG properties are to be calculated, it consists of the following attributes:
                - name: str
                - formula: str
                - state: str
                - mole_fraction: float, optional
        source : Source
            The source containing data for calculations.
        component_key : Literal
            The key to identify the component in the source data. Defaults to 'Name-State'.
        """
        # NOTE: component
        self.component = component
        # NOTE: source
        self.source = source
        # NOTE: component key
        self.component_key = component_key

        # >> set component id
        self.component_id = set_component_id(
            component=self.component,
            component_key=self.component_key
        )

        # SECTION: extract component source
        self.component_source = self.source.get_component_data(
            component_id=self.component_id,
            components=[self.component],
            component_key=self.component_key
        )

    def _get_formation_data(self, prop_name: str):
        '''
        Retrieve formation data for the specified property.

        Parameters
        ----------
        prop_name : str
            The name of the property for which formation data is to be retrieved.

        Returns
        -------
        Any
            The formation data for the specified property, or None if an error occurs.

        '''
        try:
            res = self.source.data_extractor(
                component_id=self.component_id,
                prop_name=prop_name,
            )

            return res
        except Exception as e:
            logger.exception(
                f"Error retrieving formation data for {prop_name}: {e}")
            return None

    def calc_enthalpy_of_formation(self, temperature: Temperature):
        pass

    def calc_entropy_of_formation(self, temperature: Temperature):
        pass

    def calc_gibbs_free_energy_of_formation(self, temperature: Temperature):
        pass
