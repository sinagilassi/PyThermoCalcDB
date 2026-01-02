# import libs
import logging
from typing import Literal, Optional, List
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, Temperature, Pressure
from pyThermoLinkDB.thermo import Source
# local
from ..utils.tools import measure_time
from ..core.hsg_properties import HSGProperties
from ..models.component_ref import (
    ComponentGibbsEnergyOfFormation,
    ComponentEnthalpyOfFormation,
    ComponentEnthalpyChange,
    ComponentEntropyChange
)

# NOTE: Logger
logger = logging.getLogger(__name__)


@measure_time
def calc_enthalpy_change(
        component: Component,
        model_source: ModelSource,
        temperature_initial: Temperature,
        temperature_final: Temperature,
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
        **kwargs
) -> Optional[ComponentEnthalpyChange]:
    """
    Calculate the enthalpy change (J/mol) between two temperatures (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the enthalpy change.
    model_source : ModelSource
        The source model containing necessary data.
    temperature_initial : Temperature
        The initial temperature.
    temperature_final : Temperature
        The final temperature.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentEnthalpyChange]
        The enthalpy change value, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the enthalpy change.
    - The enthalpy change is calculated by integrating the heat capacity equation between the two specified temperatures.
    - The result is provided in J/mol.
    - Reference temperature is set to 298.15 K.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(temperature_initial, Temperature):
            logger.error("Invalid initial temperature provided.")
            return None

        if not isinstance(temperature_final, Temperature):
            logger.error("Invalid final temperature provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate enthalpy change
        delta_H = hsg_props.calc_enthalpy_change(
            T1=temperature_initial,
            T2=temperature_final
        )

        return delta_H
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy change for component '{component.name}': {e}")
        return None


@measure_time
def calc_enthalpy_of_formation_at_temperature(
    component: Component,
    model_source: ModelSource,
    temperature: Temperature,
    component_key: Literal[
        'Name-State',
        'Formula-State',
        'Name',
        'Formula',
        'Name-Formula-State',
        'Formula-Name-State'
    ] = 'Name-State',
    **kwargs
) -> Optional[ComponentEnthalpyOfFormation]:
    """
    Calculate the enthalpy of formation (J/mol) at a given temperature (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the enthalpy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperature : Temperature
        The temperature at which to calculate the enthalpy of formation.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentEnthalpyOfFormation]
        A ComponentEnthalpyOfFormation object containing the enthalpy of formation value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the enthalpy of formation.
    - The enthalpy of formation symbol must be consistent with the defined unit in the configuration which is EnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All enthalpy results are provided in J/mol.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(temperature, Temperature):
            logger.error("Invalid temperature provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        EnFo_result = hsg_props.calc_enthalpy_of_formation(
            temperature=temperature
        )

        return EnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy of formation for component '{component.name}': {e}")
        return None


@measure_time
def calc_gibbs_energy_of_formation_at_temperature(
        component: Component,
        model_source: ModelSource,
        temperature: Temperature,
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
        **kwargs
) -> Optional[ComponentGibbsEnergyOfFormation]:
    """
    Calculate the Gibbs free energy of formation (J/mol) at a given temperature (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the Gibbs free energy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperature : Temperature
        The temperature at which to calculate the Gibbs free energy of formation.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentGibbsEnergyOfFormation]
        A ComponentGibbsEnergyOfFormation object containing the Gibbs free energy of formation value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the Gibbs free energy of formation.
    - The Gibbs free energy of formation symbol must be consistent with the defined unit in the configuration which is GiEnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All Gibbs energy results are provided in J/mol.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(temperature, Temperature):
            logger.error("Invalid temperature provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        GiEnFo_result = hsg_props.calc_gibbs_free_energy_of_formation(
            temperature=temperature
        )

        return GiEnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating Gibbs free energy of formation for component '{component.name}': {e}")
        return None


@measure_time
def calc_enthalpy_of_formation_range(
        component: Component,
        model_source: ModelSource,
        temperatures: list[Temperature],
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
        **kwargs
) -> Optional[list[ComponentEnthalpyOfFormation]]:
    """
    Calculate the enthalpy of formation (J/mol) over a range of temperatures (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the enthalpy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperatures : list[Temperature]
        The list of temperatures at which to calculate the enthalpy of formation.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[list[ComponentEnthalpyOfFormation]]
        A list of ComponentEnthalpyOfFormation objects containing the enthalpy of formation values, units, and symbols, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the enthalpy of formation.
    - The enthalpy of formation symbol must be consistent with the defined unit in the configuration which is EnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All enthalpy results are provided in J/mol.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not all(isinstance(temp, Temperature) for temp in temperatures):
            logger.error("Invalid temperatures provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        EnFo_results = hsg_props.calc_enthalpy_of_formation_range(
            temperatures=temperatures
        )

        return EnFo_results
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy of formation range for component '{component.name}': {e}")
        return None


@measure_time
def calc_gibbs_energy_of_formation_range(
        component: Component,
        model_source: ModelSource,
        temperatures: list[Temperature],
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
        **kwargs
) -> Optional[list[ComponentGibbsEnergyOfFormation]]:
    """
    Calculate the Gibbs free energy of formation (J/mol) over a range of temperatures (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the Gibbs free energy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperatures : list[Temperature]
        The list of temperatures at which to calculate the Gibbs free energy of formation.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[list[ComponentGibbsEnergyOfFormation]]
        A list of ComponentGibbsEnergyOfFormation objects containing the Gibbs free energy of formation values, units, and symbols, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the Gibbs free energy of formation.
    - The Gibbs free energy of formation symbol must be consistent with the defined unit in the configuration which is GiEnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All Gibbs energy results are provided in J/mol.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not all(isinstance(temp, Temperature) for temp in temperatures):
            logger.error("Invalid temperatures provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        GiEnFo_results = hsg_props.calc_gibbs_free_energy_of_formation_range(
            temperatures=temperatures
        )

        return GiEnFo_results
    except Exception as e:
        logger.error(
            f"Error calculating Gibbs free energy of formation range for component '{component.name}': {e}")
        return None


@measure_time
def calc_mixture_enthalpy(
        component: List[Component],
        model_source: ModelSource,
        temperature: Temperature,
        phases: List[Literal['Liquid', 'Vapor', 'Vapor-Liquid']],
        include_departure: bool = False,
        include_excess: bool = False,
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
        verbose: bool = False,
        **kwargs
):
    pass


@measure_time
def calc_entropy_change(
        component: Component,
        model_source: ModelSource,
        temperature_initial: Temperature,
        temperature_final: Temperature,
        pressure_initial: Pressure,
        pressure_final: Pressure,
        phase: Literal['IG', 'LIQ'],
        component_key: Literal[
            'Name-State',
            'Formula-State',
            'Name',
            'Formula',
            'Name-Formula-State',
            'Formula-Name-State'
        ] = 'Name-State',
        **kwargs
) -> Optional[ComponentEntropyChange]:
    """
    Calculate the entropy change (J/mol.K) between two temperatures (K) and pressures in any unit for a component.

    The entropy change is calculated by integrating the heat capacity equation between the two specified temperatures and accounting for pressure changes based on the phase of the component. The result is provided in J/mol.K.

    The calculation considers the following:
    - For ideal gases, the pressure change contribution is included using the relation ΔS = nR ln(P2/P1).
    - For liquids, the pressure change contribution is typically negligible and may be omitted.

    The expression for entropy change is given by:
        ΔS = ∫(Cp/T) dT + ln(P2/P1) * R  (for ideal gases)
        ΔS = ∫(Cp/T) dT                    (for liquids)

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the entropy change.
    model_source : ModelSource
        The source model containing necessary data.
    temperature_initial : Temperature
        The initial temperature.
    temperature_final : Temperature
        The final temperature.
    pressure_initial : Pressure
        The initial pressure.
    pressure_final : Pressure
        The final pressure.
    phase : Literal['IG', 'LIQ'], optional
        The phase of the component ('IG' for ideal gas, 'LIQ' for liquid).
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentEntropyChange]
        The entropy change value, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGProperties class and uses it to compute the entropy change.
    - The entropy change is calculated by integrating the heat capacity equation between the two specified temperatures and accounting for pressure changes.
    - The result is provided in J/mol.K.
    - R is the universal gas constant (8.3145 J/mol.K).
    - Pressures can be provided in any unit as they will be converted internally.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(temperature_initial, Temperature):
            logger.error("Invalid initial temperature provided.")
            return None

        if not isinstance(temperature_final, Temperature):
            logger.error("Invalid final temperature provided.")
            return None

        if not isinstance(pressure_initial, Pressure):
            logger.error("Invalid initial pressure provided.")
            return None

        if not isinstance(pressure_final, Pressure):
            logger.error("Invalid final pressure provided.")
            return None

        if phase not in ['IG', 'LIQ']:
            logger.error("Invalid phase provided. Must be 'IG' or 'LIQ'.")
            return None

        # SECTION: Prepare source
        Source_ = Source(model_source=model_source)

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate entropy change
        delta_S = hsg_props.calc_entropy_change(
            T1=temperature_initial,
            T2=temperature_final,
            P1=pressure_initial,
            P2=pressure_final,
            phase=phase
        )

        return delta_S
    except Exception as e:
        logger.error(
            f"Error calculating entropy change for component '{component.name}': {e}")
        return None
