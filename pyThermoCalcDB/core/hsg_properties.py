# import libs
import logging
from typing import Dict, Any, Literal, Optional
from pythermodb_settings.models import Component, Temperature
from pythermodb_settings.utils import set_component_id
import pycuc
# local
from ..thermo import Source
from ..configs.thermo_props import (
    EnFo_IG_UNIT,
    Ent_STD_UNIT,
    GiEnFo_IG_UNIT,
    EnFo_IG_SYMBOL,
    Ent_STD_SYMBOL,
    GiEnFo_IG_SYMBOL,
    Cp_IG_SYMBOL
)
from ..utils.math_tools import integrate_function

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
    # NOTE: attributes
    T_ref = 298.15  # reference temperature in K
    R = 8.3145  # universal gas constant in J/mol·K

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

    def _get_formation_data(self, prop_name: str) -> Optional[Dict[str, Any]]:
        '''
        Retrieve formation data for the specified property.

        Parameters
        ----------
        prop_name : str
            The name of the property for which formation data is to be retrieved.

        Returns
        -------
        Optional[Dict[str, Any]]
            A dictionary containing the formation data if available, otherwise None.
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
        '''
        Calculate the enthalpy of formation at a given temperature.

        Parameters
        ----------
        temperature : Temperature
            The temperature at which to calculate the enthalpy of formation.

        Returns
        -------
        dict
            A dictionary containing the calculated enthalpy of formation.
        '''
        try:
            # SECTION: get formation data
            # >> enthalpy of formation (EnFo_IG) unit ?
            formation_data = self._get_formation_data(EnFo_IG_SYMBOL)

            # >> check formation data
            if formation_data is None:
                logger.warning(
                    f"No formation data available for enthalpy of formation.")
                return None

            # NOTE: unit conversion if necessary
            # TODO: convert to [J/mol]
            EnFo_val = formation_data['value']
            EnFo_unit = formation_data['unit']
            # ! to [J/mol]
            unit_ = f"{EnFo_unit} => {EnFo_IG_UNIT}"
            EnFo = pycuc.to(
                EnFo_val,
                unit_
            )

            # SECTION: heat capacity calculation
            # NOTE: build equation
            Cp_eq_src = self.source.eq_builder(
                components=[self.component],
                prop_name=Cp_IG_SYMBOL,
                component_key=self.component_key  # type: ignore
            )

            # >> check Cp equation
            if Cp_eq_src is None:
                logger.warning(
                    f"No heat capacity equation available for component {self.component_id}.")
                return None

            # >> for component
            component_Cp_eq_src = Cp_eq_src.get(self.component_id)

            if component_Cp_eq_src is None:
                logger.warning(
                    f"No heat capacity equation available for component {self.component_id}.")
                return None

            # NOTE: extract Cp equation
            # ! args
            Cp_eq_args = component_Cp_eq_src['args']
            # ! arg symbols
            Cp_eq_arg_symbols = component_Cp_eq_src['arg_symbols']
            # ! returns
            Cp_eq_returns = component_Cp_eq_src['returns']
            # >> get return units
            returns_outer_key, returns_inner = next(
                iter(Cp_eq_returns.items()))
            Cp_eq_return_unit = returns_inner['unit']

            # ! return symbols
            Cp_eq_return_symbols = component_Cp_eq_src['return_symbols']
            # ! equation
            Cp_eq = component_Cp_eq_src['value']

            # NOTE: calculate heat capacity at given temperature
            delta_Cp = integrate_function(
                func=Cp_eq,
                var_symbol='T',
                vars=Cp_eq_args,
                a=self.T_ref,
                b=temperature.value,
            )

            # >> unit conversion if necessary
            unit_ = f"{Cp_eq_return_unit} => {EnFo_IG_UNIT}"
            delta_Cp_value = pycuc.to(
                delta_Cp,
                unit_
            )

            # SECTION: calculate enthalpy of formation at temperature
            # ! J/mol
            EnFo_T = EnFo + delta_Cp_value

            # SECTION: prepare result
            result = {
                'value': EnFo_T,
                'unit': EnFo_IG_UNIT,
                'symbol': EnFo_IG_SYMBOL
            }
            return result
        except Exception as e:
            logger.exception(
                f"Error calculating enthalpy of formation: {e}")
            return None

    def calc_gibbs_free_energy_of_formation(self, temperature: Temperature):
        pass
