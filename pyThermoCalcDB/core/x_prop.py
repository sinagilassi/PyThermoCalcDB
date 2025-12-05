# import libs
import logging
from typing import Dict, Literal, Any, Optional
from pythermodb_settings.models import Component
from pyThermoDB.core import TableEquation
from pythermodb_settings.utils import set_component_id
from pyThermoDB.models import EquationResult
import pycuc
from math import pow
# local
from ..thermo import Source
from ..models import ComponentEquationSource, CalcResult


# NOTE: Logger
logger = logging.getLogger(__name__)


class XProp:
    def __init__(
        self,
        prop_name: str,
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
        prop_name : str
            The ID of the equation to be used for calculations, e.g., 'VaPr', 'Cp_IG'.
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
        # NOTE: equation id
        self.prop_name = prop_name
        # NOTE: component
        self.component = component
        # NOTE: source
        self.source = source
        # NOTE: component key
        self.component_key = component_key

        # SECTION: set component id
        self.component_id = set_component_id(
            component=self.component,
            component_key=self.component_key
        )

        # SECTION: get vapor pressure equation source
        self.component_equation: ComponentEquationSource = self._get_equation_from_source()

        # NOTE: equation
        # ! equation
        self._eq: TableEquation = self.component_equation.value
        # ! equation args/returns/symbols
        self._args: Dict[
            str,
            Any
        ] = self.component_equation.args or {}
        self._arg_symbols: Dict[
            str,
            Any
        ] = self.component_equation.arg_symbols or {}
        self._returns: Dict[
            str,
            Any
        ] = self.component_equation.returns or {}
        self._return_symbols: Dict[
            str,
            Any
        ] = self.component_equation.return_symbols or {}

        # NOTE: units and symbols for return
        # >> get return units
        returns_outer_key, returns_inner = next(
            iter(self._returns.items()))
        self._return_unit: str = returns_inner['unit']
        self._return_symbol: str = returns_inner['symbol']

    @property
    def equation(self) -> TableEquation:
        """
        Get the equation object.

        Returns
        -------
        TableEquation
            The equation object.
        """
        return self._eq

    @property
    def body(self) -> str:
        """
        Get the body of the equation as a string.

        Returns
        -------
        str
            The body of the equation.
        """
        return self._eq.body

    @property
    def args(self) -> Dict[str, Any]:
        """
        Get the arguments required for the equation.

        Returns
        -------
        Dict[str, Any]
            A dictionary of argument names and their details.
        """
        return self._args

    @property
    def arg_symbols(self) -> Dict[str, Any]:
        """
        Get the symbols for the equation arguments.

        Returns
        -------
        Dict[str, Any]
            A dictionary of argument names and their symbols.
        """
        return self._arg_symbols

    @property
    def returns(self) -> Dict[str, Any]:
        """
        Get the return values of the equation.

        Returns
        -------
        Dict[str, Any]
            A dictionary of return value names and their details.
        """
        return self._returns

    @property
    def return_symbol(self) -> str:
        """
        Get the symbol for the equation return values.

        Returns
        -------
        str
            The symbol for the return value.
        """
        return self._return_symbol

    def _get_equation_from_source(self) -> ComponentEquationSource:
        """
        Retrieve the component equation for the component from the source.

        Returns
        -------
        ComponentEquationSource
            The component equation source for the specified component and equation ID.

        Raises
        ------
        ValueError
            If no equation is found for the specified component ID.
        """
        try:
            # SECTION: get equations for component
            equations = self.source.eq_builder(
                components=[self.component],
                prop_name=self.prop_name,
                component_key=self.component_key  # type: ignore
            )

            if not equations:
                logger.warning(
                    f'No {self.prop_name} equation found for component ID: {self.component_id}'
                )
                raise ValueError(
                    f'No {self.prop_name} equation found for component ID: {self.component_id}'
                )

            # SECTION: select for component
            eq: ComponentEquationSource | None = equations.get(
                self.component_id
            )

            # >> check
            if eq is None:
                logger.warning(
                    f'No {self.prop_name} equation found for component ID: {self.component_id}'
                )
                raise ValueError(
                    f'No {self.prop_name} equation found for component ID: {self.component_id}')

            return eq
        except Exception as e:
            logger.error(
                f'Error retrieving {self.prop_name} equation for component ID: {self.component_id} - {e}'
            )
            raise e

    def _calc(
        self,
        **input_args
    ) -> Optional[CalcResult]:
        """
        Calculate with args.

        Parameters
        ----------
        **input_args : Dict[str, Any]
            Input arguments for the equation.

        Returns
        -------
        Optional[CalcResult]
            The calculation results, encapsulated in CalcResult as:
            - value: float
            - unit: str
            - symbol: str
        """
        try:
            # SECTION: variables
            if not input_args:
                logger.error(
                    f"No input arguments provided for calculation."
                )
                return None

            # SECTION: calculate
            _result: EquationResult = self._eq.cal(
                **{**self._args, **input_args}
            )
            # >> value
            _value = float(_result['value']) if isinstance(
                _result['value'], (float, int)) else None

            # >> check
            if _value is None:
                logger.error(
                    f"Calculation resulted in None value."
                )
                return None

            # res
            res = {
                'value': _value,
                'unit': self._return_unit,
                'symbol': self._return_symbol
            }

            # >> set
            res = CalcResult(**res)

            return res
        except Exception as e:
            logger.error(
                f"Error calculating {e}")
            return None
