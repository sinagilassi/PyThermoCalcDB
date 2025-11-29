# import libs
from typing import Dict, Any
from pythermodb_settings.models import Component, Temperature
# local
from ..thermo import Source


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
        """
        # NOTE: component
        self.component = component
        # NOTE: source
        self.source = source

    def calc_enthalpy_of_formation(self, temperature: Temperature):
        pass

    def calc_entropy_of_formation(self, temperature: Temperature):
        pass

    def calc_gibbs_free_energy_of_formation(self, temperature: Temperature):
        pass
