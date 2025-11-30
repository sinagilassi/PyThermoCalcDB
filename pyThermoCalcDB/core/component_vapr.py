# import libs
import logging
from typing import Optional, Dict, List, Literal, Any, Tuple
from pythermodb_settings.models import Component
from pyThermoDB.core import TableEquation
from pythermodb_settings.utils import set_component_id
from pyThermoDB.models import EquationResult
from scipy.optimize import root_scalar
import pycuc
# local
from ..thermo import Source
from ..configs.thermo_props import VaPr_SYMBOL, VaPr_UNIT
from ..models import ComponentEquationSource

# NOTE: logger
logger = logging.getLogger(__name__)


class ComponentVaporPressure:
    """
    Component vapor pressure class used to calculate the following properties:
    - Saturation pressure VaPr(T)
    - Enthalpy of vaporization EnVap(T)
    - Temperature at given vapor pressure TeVaPr(P)
    """

    # NOTE: attributes
    T_ref = 298.15  # reference temperature in K
    R = 8.314462618  # universal gas constant in J/mol.K

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

        # SECTION: set component id
        self.component_id = set_component_id(
            component=self.component,
            component_key=self.component_key
        )

        # SECTION: get vapor pressure equation source
        self.vapor_pressure_equation: ComponentEquationSource = self._get_vapor_pressure_equation()

        # NOTE: equation
        # ! vapor pressure equation
        self.VaPr_eq: TableEquation = self.vapor_pressure_equation.value
        # ! vapor pressure equation args/returns/symbols
        self.VaPr_args: Dict[
            str,
            Any
        ] = self.vapor_pressure_equation.args or {}
        self.VaPr_arg_symbols: Dict[
            str,
            Any
        ] = self.vapor_pressure_equation.arg_symbols or {}
        self.VaPr_returns: Dict[
            str,
            Any
        ] = self.vapor_pressure_equation.returns or {}
        self.VaPr_return_symbols: Dict[
            str,
            Any
        ] = self.vapor_pressure_equation.return_symbols or {}

        # >> get return units
        returns_outer_key, returns_inner = next(
            iter(self.VaPr_returns.items()))
        self.VaPr_return_unit: str = returns_inner['unit']

    def _get_vapor_pressure_equation(self) -> ComponentEquationSource:
        """
        Retrieve the vapor pressure equation for the component from the source.

        Returns
        -------
        ComponentEquationSource
            The vapor pressure equation source for the component.
        """
        try:
            # SECTION: get equations for component
            equations = self.source.eq_builder(
                components=[self.component],
                prop_name=VaPr_SYMBOL,
                component_key=self.component_key  # type: ignore
            )

            if not equations:
                logger.warning(
                    f'No vapor pressure equation found for component ID: {self.component_id}'
                )
                raise ValueError(
                    f'No vapor pressure equation found for component ID: {self.component_id}'
                )

            # SECTION: select for component
            eq: ComponentEquationSource | None = equations.get(
                self.component_id)

            # >> check
            if eq is None:
                logger.warning(
                    f'No vapor pressure equation found for component ID: {self.component_id}'
                )
                raise ValueError(
                    f'No vapor pressure equation found for component ID: {self.component_id}')

            return eq
        except Exception as e:
            logger.error(
                f'Error retrieving vapor pressure equation for component ID: {self.component_id} - {e}'
            )
            raise e

    def calc_dPsat__dT(
            self,
            temperature: float,
            h: Optional[float] = None
    ) -> float:
        """
        Calculate the derivative of saturation pressure with respect to temperature.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin.
        h : Optional[float]
            Small temperature increment for numerical differentiation. If None, a default value is used.

        Returns
        -------
        float
            The derivative dPsat/dT in Pa/K.
        """
        try:
            # SECTION: validate temperature
            if temperature <= 0.0:
                raise ValueError(
                    "Temperature must be greater than zero Kelvin.")

            # SECTION: numerical derivative
            # NOTE: central difference
            delta_T = 1e-5 if h is None else h

            # NOTE: calculate Psat at T + delta_T and T - delta_T
            T_plus = temperature + delta_T
            T_minus = temperature - delta_T

            # >> calculate Psat values
            VaPr_plus: EquationResult = self.VaPr_eq.cal(
                **{**self.VaPr_args, **{"T": T_plus}}
            )
            # >> value
            VaPr_plus_value = float(VaPr_plus['value']) if isinstance(
                VaPr_plus['value'], (float, int)) else None
            # >> check
            if VaPr_plus_value is None:
                raise ValueError(
                    f"Failed to calculate vapor pressure at T={T_plus} K."
                )

            VaPr_minus: EquationResult = self.VaPr_eq.cal(
                **{**self.VaPr_args, **{"T": T_minus}}
            )
            # >> value
            VaPr_minus_value = float(VaPr_minus['value']) if isinstance(
                VaPr_minus['value'], (float, int)) else None
            # >> check
            if VaPr_minus_value is None:
                raise ValueError(
                    f"Failed to calculate vapor pressure at T={T_minus} K."
                )

            # NOTE: unit conversion to Pa
            VaPr_plus_value = pycuc.to(
                VaPr_plus_value,
                f"{self.VaPr_return_unit} => Pa"
            )
            VaPr_minus_value = pycuc.to(
                VaPr_minus_value,
                f"{self.VaPr_return_unit} => Pa"
            )

            res = (VaPr_plus_value - VaPr_minus_value) / (2 * delta_T)

            return res
        except Exception as e:
            logger.error(
                f"Error calculating dPsat/dT at T={temperature} K: {e}")
            raise e

    def calc_EnVap_Clapeyron(
            self,
            temperature: float,
            h: Optional[float] = None
    ) -> float:
        """
        Calculate the enthalpy of vaporization using the Clausius–Clapeyron relation.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin.
        h : Optional[float]
            Small temperature increment for numerical differentiation. If None, a default value is used.

        Returns
        -------
        float
            Enthalpy of vaporization in J/mol.
        """
        try:
            # SECTION: calculate Psat at T
            VaPr_result: EquationResult = self.VaPr_eq.cal(
                **{**self.VaPr_args, **{"T": temperature}}
            )
            # >> value
            VaPr_value = float(VaPr_result['value']) if isinstance(
                VaPr_result['value'], (float, int)) else None
            # >> check
            if VaPr_value is None:
                raise ValueError(
                    f"Failed to calculate vapor pressure at T={temperature} K."
                )

            # NOTE: unit conversion to Pa
            VaPr_value = pycuc.to(
                VaPr_value,
                f"{self.VaPr_return_unit} => Pa"
            )

            # SECTION: calculate dPsat/dT
            dPsat_dT = self.calc_dPsat__dT(temperature, h=h)

            # SECTION: Clausius–Clapeyron relation
            R = self.R  # J/mol.K

            EnVap = (R * temperature**2 / VaPr_value) * dPsat_dT

            return EnVap
        except Exception as e:
            logger.error(
                f"Error calculating EnVap via Clausius–Clapeyron at T={temperature} K: {e}")
            raise e

    def calc_TeVaPr(
            self,
            pressure: float,
            temperature_guess: Optional[float] = None,
            T_bracket: Optional[Tuple[float, float]] = None,
            method: str = "auto",
            tol: float = 1e-6,
            max_iter: int = 50,
            h: Optional[float] = None
    ) -> float:
        """
        Calculate the temperature at a given vapor pressure using the bisection method.

        Parameters
        ----------
        pressure : float
            Vapor pressure in Pa.
        temperature_guess : float
            Initial guess for temperature in Kelvin.
        T_bracket : Tuple[float, float]
            Bracket (min, max) for temperature in Kelvin.
        method : str
            Method to use for root finding. Defaults to "auto".
            - "auto": use brentq if bracket is provided, else Newton with derivative
            - "brentq": robust bracketing method (requires T_bracket)
            - "bisect": slower but very robust (requires T_bracket)
            - "newton": uses derivative dPsat/dT and T_guess
        tol : float
            Tolerance for convergence.
        max_iter : int
            Maximum number of iterations.
        h : Optional[float]
            Small temperature increment for numerical differentiation in dPsat/dT. If None, a default value is used.

        Returns
        -------
        float
            Temperature in Kelvin corresponding to the given vapor pressure.

        Notes
        -----
        - The temperature is found by solving Psat(T) - P = 0 using root scalar method.
        - Bracketing interval [K] such that Psat(T_low) - P and Psat(T_high) - P have opposite signs. Used for robust bracketing methods.
        """
        try:
            # SECTION: define function to find root
            def func(T: float) -> float:
                VaPr_result: EquationResult = self.VaPr_eq.cal(
                    **{**self.VaPr_args, **{"T": T}}
                )
                # >> value
                VaPr_value = float(VaPr_result['value']) if isinstance(
                    VaPr_result['value'], (float, int)) else None
                # >> check
                if VaPr_value is None:
                    raise ValueError(
                        f"Failed to calculate vapor pressure at T={T} K."
                    )

                # NOTE: unit conversion to Pa
                VaPr_value = pycuc.to(
                    VaPr_value,
                    f"{self.VaPr_return_unit} => Pa"
                )

                return VaPr_value - pressure

            # SECTION: define derivative function
            def func_prime(T: float) -> float:
                return self.calc_dPsat__dT(T)

            # SECTION: select method
            method = method.lower()

            # NOTE: auto-selection
            if method == "auto":
                if T_bracket is not None:
                    method = "brentq"
                else:
                    method = "newton"

            # NOTE: root finding
            if method == "newton":
                if temperature_guess is None:
                    raise ValueError(
                        "temperature_guess must be provided for method 'newton'."
                    )

                # ! root finding using Newton's method
                sol = root_scalar(
                    func,
                    fprime=func_prime,
                    method='newton',
                    x0=temperature_guess,
                    xtol=tol,
                    maxiter=max_iter
                )
            elif method in ("brentq", "bisect"):
                # >> check bracket
                if T_bracket is None:
                    raise ValueError(
                        f"T_bracket must be provided for method '{method}'."
                    )

                # ! root finding using bracketing method
                sol = root_scalar(
                    func,
                    bracket=T_bracket,
                    method=method,
                    xtol=tol,
                    maxiter=max_iter
                )

            else:
                raise ValueError(
                    "method must be 'auto', 'newton', 'brentq', or 'bisect'."
                )

            if not sol.converged:
                raise RuntimeError(
                    f"Root finding did not converge for pressure={pressure} Pa."
                )

            return sol.root
        except Exception as e:
            logger.error(
                f"Error calculating TeVaPr for P={pressure} Pa: {e}"
            )
            raise e
