# import libs
import logging
from typing import Literal, Optional, List, cast
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import (
    Component,
    Temperature,
    Pressure,
    CustomProp,
    ComponentKey,
    CustomProperty
)
from pyThermoLinkDB.thermo import Source
import pycuc
# local
from ..configs.constants import T_298_K
from ..utils.tools import measure_time
from ..utils.component_tools import map_state_to_phase
from ..core.hsg_properties import HSGProperties
from ..core.hsg_mixture import HSGMixture
from ..models.component_ref import (
    ComponentGibbsFreeEnergy,
    ComponentEnthalpy,
    ComponentEnthalpyChange,
    ComponentEntropyChange,
    MixtureEnthalpyResult
)
from ..configs.thermo_props import (
    EnFo_IG_SYMBOL,
    EnFo_LIQ_SYMBOL,
    EnFo_SOL_SYMBOL,
    GiEnFo_IG_SYMBOL,
    GiEnFo_LIQ_SYMBOL,
    GiEnFo_SOL_SYMBOL,
    Cp_IG_SYMBOL,
    Cp_LIQ_SYMBOL,
    Cp_SOL_SYMBOL,
    EnVap_SYMBOL,
    EnSub_SYMBOL,
)

# NOTE: Logger
logger = logging.getLogger(__name__)


@measure_time
def build_hsg_properties(
    component: Component,
    model_source: ModelSource,
    component_key: ComponentKey = 'Name-Formula',
    **kwargs
) -> Optional[HSGProperties]:
    """
    Build and return an instance of the HSGProperties class for a given component and model source.

    Parameters
    ----------
    component : Component
        The chemical component for which to build the HSGProperties instance.
    model_source : ModelSource
        The source model containing necessary data.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-Formula'.
    **kwargs
        Additional keyword arguments.

    Returns
    -------
    Optional[HSGProperties]
        An instance of the HSGProperties class, or None if initialization fails.

    Notes
    -----
    - The function prepares the source using the provided model source and component key, then initializes and returns an instance of the HSGProperties class.
    - This function can be used as a helper to create HSGProperties instances for various calculations.
    """
    try:
        # SECTION: Input validation
        if not isinstance(component, Component):
            logger.error("Invalid component provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        return hsg_props
    except Exception as e:
        logger.error(
            f"Error building HSGProperties for component '{component.name}': {e}")
        return None

# SECTION: enthalpy change calculations
# ! enthalpy change between two temperatures


@measure_time
def calc_dEn(
        component: Component,
        model_source: ModelSource,
        temperature_initial: Temperature,
        temperature_final: Temperature,
        component_key: ComponentKey = 'Name-Formula',
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
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

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

# ! enthalpy change between two temperatures using hsg properties


def calc_dEn_hsg(
        hsg_props: HSGProperties,
        temperature_initial: Temperature,
        temperature_final: Temperature,
        **kwargs
) -> Optional[ComponentEnthalpyChange]:
    """
    Calculate the enthalpy change (J/mol) between two temperatures (K) for a component using an existing HSGProperties instance.

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    temperature_initial : Temperature
        The initial temperature.
    temperature_final : Temperature
        The final temperature.
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
    - The function uses the provided HSGProperties instance to compute the enthalpy change by integrating the heat capacity equation between the two specified temperatures.
    - The result is provided in J/mol.
    - Reference temperature is set to 298.15 K.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
            return None

        # NOTE: calculate enthalpy change
        delta_H = hsg_props.calc_enthalpy_change(
            T1=temperature_initial,
            T2=temperature_final
        )

        return delta_H
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy change using HSGProperties: {e}")
        return None

# SECTION: enthalpy calculations
# ! enthalpy at temperature


@measure_time
def calc_En(
    component: Component,
    model_source: ModelSource,
    temperature: Temperature,
    component_key: ComponentKey = 'Name-Formula',
    **kwargs
) -> Optional[ComponentEnthalpy]:
    """
    Calculate the enthalpy (J/mol) at a given temperature (K) for a component.

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
    Optional[ComponentEnthalpy]
        A ComponentEnthalpy object containing the enthalpy of formation value, unit, and symbol, or None if calculation fails.

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
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        EnFo_result = hsg_props.calc_enthalpy(
            temperature=temperature
        )

        return EnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy of formation for component '{component.name}': {e}")
        return None

# ! enthalpy at temperature using hsg properties


def calc_En_hsg(
    hsg_props: HSGProperties,
    temperature: Temperature,
    **kwargs
) -> Optional[ComponentEnthalpy]:
    """
    Calculate the enthalpy (J/mol) at a given temperature (K) for a component using an existing HSGProperties instance.

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    temperature : Temperature
        The temperature at which to calculate the enthalpy of formation.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentEnthalpy]
        A ComponentEnthalpy object containing the enthalpy of formation value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function uses the provided HSGProperties instance to compute the enthalpy of formation at the specified temperature.
    - The enthalpy of formation symbol must be consistent with the defined unit in the configuration which is EnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All enthalpy results are provided in J/mol.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
            return None

        # NOTE: calculate
        EnFo_result = hsg_props.calc_enthalpy(
            temperature=temperature
        )

        return EnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy of formation using HSGProperties: {e}"
        )
        return None


# SECTION: reference enthalpy calculations (phase inferred from component state)
# ! reference enthalpy at temperature

@measure_time
def calc_En_IG_ref(
    component: Component,
    model_source: ModelSource,
    temperature: Temperature,
    component_key: ComponentKey = 'Name-Formula',
    **kwargs
) -> Optional[ComponentEnthalpy]:
    """
    Calculate the reference enthalpy (J/mol) at a given temperature (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the reference enthalpy.
    model_source : ModelSource
        The source model containing necessary data.
    temperature : Temperature
        The temperature at which to calculate the reference enthalpy.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentEnthalpy]
        A ComponentEnthalpy object containing the reference enthalpy value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function initializes HSGProperties and computes phase-specific reference enthalpy via `calc_reference_enthalpy`.
    - The phase is inferred from `component.state` and mapped to one of `IG`, `LIQ`, or `SOL`.
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

        # SECTION: Component phase
        component_state = component.state.upper()
        phase = map_state_to_phase(
            state=component_state
        )

        # SECTION: Prepare source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        EnFo_IG_result = hsg_props.calc_reference_enthalpy(
            temperature=temperature,
            phase=cast(Literal["IG", "LIQ", "SOL"], phase)
        )

        return EnFo_IG_result
    except Exception as e:
        logger.error(
            f"Error calculating reference enthalpy for component '{component.name}': {e}")
        return None

# ! reference enthalpy at temperature using hsg properties


def calc_En_IG_ref_hsg(
    hsg_props: HSGProperties,
    temperature: Temperature,
    **kwargs
) -> Optional[ComponentEnthalpy]:
    """
    Calculate the reference enthalpy (J/mol) at a given temperature (K) using an existing HSGProperties instance.

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    temperature : Temperature
        The temperature at which to calculate the reference enthalpy.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentEnthalpy]
        A ComponentEnthalpy object containing the reference enthalpy value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function uses the provided HSGProperties instance to compute phase-specific reference enthalpy via `calc_reference_enthalpy`.
    - The phase is inferred from `hsg_props.component.state` and mapped to one of `IG`, `LIQ`, or `SOL`.
    - Reference temperature is set to 298.15 K.
    - All enthalpy results are provided in J/mol.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
            return None

        if not isinstance(temperature, Temperature):
            logger.error("Invalid temperature provided.")
            return None

        # SECTION: Component phase
        component_state = hsg_props.component.state.upper()
        phase = map_state_to_phase(
            state=component_state
        )

        # NOTE: calculate
        EnFo_IG_result = hsg_props.calc_reference_enthalpy(
            temperature=temperature,
            phase=cast(Literal["IG", "LIQ", "SOL"], phase)
        )

        return EnFo_IG_result
    except Exception as e:
        logger.error(
            f"Error calculating reference enthalpy using HSGProperties: {e}"
        )
        return None

# ! calculate enthalpy at temperature


@measure_time
def calc_En_IG_ref_hsg_plus(
    hsg_props: HSGProperties,
    EnFo_IG: CustomProperty,
    temperature: Temperature,
    temperature_ref: Temperature = T_298_K,
    output_unit: Optional[str] = None,
    **kwargs
) -> Optional[CustomProperty]:
    """
    Calculate the enthalpy (J/mol) at a given temperature (K) using an existing HSGProperties instance and a reference enthalpy.

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    EnFo_IG : CustomProperty
        The reference enthalpy value at 298.15 K.
    temperature : Temperature
        The temperature at which to calculate the enthalpy.
    temperature_ref : Temperature
        The reference temperature corresponding to the reference enthalpy (typically 298.15 K), default is T_298_K.
    output_unit : Optional[str]
        The desired output unit for the enthalpy result (e.g., 'J/mol', 'kJ/mol'). If None, the unit of the reference enthalpy will be used.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'silent'.

    Returns
    -------
    Optional[CustomProperty]
        A CustomProp object containing the enthalpy value at the specified temperature, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function calculates the enthalpy at the specified temperature by adding the reference enthalpy (EnFo_IG) to the enthalpy change from 298.15 K to the specified temperature.
    - All enthalpy results are provided in J/mol.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
            return None

        if not isinstance(EnFo_IG, CustomProperty):
            logger.error("Invalid reference enthalpy provided.")
            return None

        if not isinstance(temperature, Temperature):
            logger.error("Invalid temperature provided.")
            return None

        if not isinstance(temperature_ref, Temperature):
            logger.error("Invalid reference temperature provided.")
            return None

        # NOTE: enthalpy at reference temperature
        # ! J/mol
        if EnFo_IG.unit != "J/mol":
            val_ = pycuc.convert_from_to(
                value=EnFo_IG.value,
                from_unit=EnFo_IG.unit,
                to_unit="J/mol"
            )

            # upd
            EnFo_IG = CustomProperty(
                value=val_,
                unit="J/mol",
                symbol=EnFo_IG.symbol,
            )

        # >> check symbol
        if EnFo_IG.symbol != EnFo_IG_SYMBOL:
            logger.warning(
                f"Reference enthalpy symbol '{EnFo_IG.symbol}' does not match expected {EnFo_IG_SYMBOL}."
            )

        # NOTE: calculate enthalpy change from 298.15 K to specified temperature
        # ! J/mol
        delta_H = hsg_props.calc_enthalpy_change(
            T1=temperature_ref,
            T2=temperature,
        )

        # >> check
        if delta_H is None:
            logger.error("Failed to calculate enthalpy change.")
            return None

        # >> check unit
        if delta_H.unit != "J/mol":
            val_ = pycuc.convert_from_to(
                value=delta_H.value,
                from_unit=delta_H.unit,
                to_unit="J/mol"
            )

            # upd
            delta_H = CustomProp(
                value=val_,
                unit="J/mol"
            )

        # NOTE: calculate total enthalpy at specified temperature
        # ! J/mol
        En_result = CustomProperty(
            value=float(EnFo_IG.value + delta_H.value),
            unit=EnFo_IG.unit,
            symbol=EnFo_IG.symbol,
        )

        # NOTE: convert to output unit if specified
        if (
            output_unit is not None and
            En_result.unit != output_unit
        ):
            val_converted = pycuc.convert_from_to(
                value=En_result.value,
                from_unit=En_result.unit,
                to_unit=output_unit
            )

            # upd
            En_result = CustomProperty(
                value=val_converted,
                unit=output_unit,
                symbol=En_result.symbol,
            )

        return En_result
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy using HSGProperties and reference enthalpy: {e}"
        )
        return None

# SECTION: Gibbs free energy calculations
# ! Gibbs free energy at temperature


@measure_time
def calc_GiFrEn(
        component: Component,
        model_source: ModelSource,
        temperature: Temperature,
        phase: Literal['IG', 'LIQ', 'SOL'] = 'IG',
        component_key: ComponentKey = 'Name-Formula',
        **kwargs
) -> Optional[ComponentGibbsFreeEnergy]:
    """
    Calculate the Gibbs free energy (J/mol) at a given temperature (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the Gibbs free energy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperature : Temperature
        The temperature at which to calculate the Gibbs free energy of formation.
    phase : Literal['IG', 'LIQ', 'SOL'], optional
        The phase of the component ('IG' for ideal gas, 'LIQ' for liquid, 'SOL' for solid). Default is 'IG'.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-State'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentGibbsFreeEnergy]
        A ComponentGibbsFreeEnergy object containing the Gibbs free energy of formation value, unit, and symbol, or None if calculation fails.

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
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        GiEnFo_result = hsg_props.calc_gibbs_free_energy(
            temperature=temperature,
            phase=phase
        )

        return GiEnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating Gibbs free energy of formation for component '{component.name}': {e}")
        return None

# ! Gibbs free energy at temperature using hsg properties


def calc_GiFrEn_hsg(
        hsg_props: HSGProperties,
        temperature: Temperature,
        phase: Literal['IG', 'LIQ', 'SOL'] = 'IG',
        **kwargs
) -> Optional[ComponentGibbsFreeEnergy]:
    """
    Calculate the Gibbs free energy (J/mol) at a given temperature (K) for a component using an existing HSGProperties instance.

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    temperature : Temperature
        The temperature at which to calculate the Gibbs free energy of formation.
    phase : Literal['IG', 'LIQ', 'SOL'], optional
        The phase of the component ('IG' for ideal gas, 'LIQ' for liquid, 'SOL' for solid). Default is 'IG'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[ComponentGibbsFreeEnergy]
        A ComponentGibbsFreeEnergy object containing the Gibbs free energy of formation value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function uses the provided HSGProperties instance to compute the Gibbs free energy of formation at the specified temperature and phase.
    - The Gibbs free energy of formation symbol must be consistent with the defined unit in the configuration which is GiEnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All Gibbs energy results are provided in J/mol.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
            return None

        # NOTE: calculate
        GiEnFo_result = hsg_props.calc_gibbs_free_energy(
            temperature=temperature,
            phase=phase
        )

        return GiEnFo_result
    except Exception as e:
        logger.error(
            f"Error calculating Gibbs free energy of formation using HSGProperties: {e}")
        return None

# SECTION: enthalpy calculations over a range of temperatures
# ! enthalpy at temperature over a range of temperatures


@measure_time
def calc_En_range(
        component: Component,
        model_source: ModelSource,
        temperatures: list[Temperature],
        component_key: ComponentKey = 'Name-Formula',
        **kwargs
) -> Optional[list[ComponentEnthalpy]]:
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
    Optional[list[ComponentEnthalpy]]
        A list of ComponentEnthalpy objects containing the enthalpy of formation values, units, and symbols, or None if calculation fails.

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
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        EnFo_results = hsg_props.calc_enthalpy_range(
            temperatures=temperatures
        )

        return EnFo_results
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy of formation range for component '{component.name}': {e}")
        return None

# ! Gibbs free energy at temperature over a range of temperatures using hsg properties


def calc_En_range_hsg(
        hsg_props: HSGProperties,
        temperatures: list[Temperature],
        **kwargs
) -> Optional[list[ComponentEnthalpy]]:
    """
    Calculate the enthalpy of formation (J/mol) over a range of temperatures (K) for a component using an existing HSGProperties instance.

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    temperatures : list[Temperature]
        The list of temperatures at which to calculate the enthalpy of formation.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[list[ComponentEnthalpy]]
        A list of ComponentEnthalpy objects containing the enthalpy of formation values, units, and symbols, or None if calculation fails.

    Notes
    -----
    - The function uses the provided HSGProperties instance to compute the enthalpy of formation at the specified range of temperatures.
    - The enthalpy of formation symbol must be consistent with the defined unit in the configuration which is EnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All enthalpy results are provided in J/mol.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
            return None

        if not all(isinstance(temp, Temperature) for temp in temperatures):
            logger.error("Invalid temperatures provided.")
            return None

        # NOTE: calculate
        EnFo_results = hsg_props.calc_enthalpy_range(
            temperatures=temperatures
        )

        return EnFo_results
    except Exception as e:
        logger.error(
            f"Error calculating enthalpy of formation range using HSGProperties: {e}")
        return None

# SECTION: Gibbs free energy calculations over a range of temperatures
# ! Gibbs free energy at temperature over a range of temperatures


@measure_time
def calc_GiFrEn_range(
        component: Component,
        model_source: ModelSource,
        temperatures: list[Temperature],
        component_key: ComponentKey = 'Name-Formula',
        **kwargs
) -> Optional[list[ComponentGibbsFreeEnergy]]:
    """
    Calculate the Gibbs free energy (J/mol) over a range of temperatures (K) for a component.

    Parameters
    ----------
    component : Component
        The chemical component for which to calculate the Gibbs free energy of formation.
    model_source : ModelSource
        The source model containing necessary data.
    temperatures : list[Temperature]
        The list of temperatures at which to calculate the Gibbs free energy of formation.
    component_key : Literal[..., optional]
        The key to identify the component, by default 'Name-Formula'.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[list[ComponentGibbsFreeEnergy]]
        A list of ComponentGibbsFreeEnergy objects containing the Gibbs free energy of formation values, units, and symbols, or None if calculation fails.

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
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # SECTION: Initialize HSGProperties
        hsg_props = HSGProperties(
            component=component,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate
        GiEnFo_results = hsg_props.calc_gibbs_free_energy_range(
            temperatures=temperatures
        )

        return GiEnFo_results
    except Exception as e:
        logger.error(
            f"Error calculating Gibbs free energy of formation range for component '{component.name}': {e}")
        return None

# ! Gibbs free energy at temperature over a range of temperatures using hsg properties


def calc_GiFrEn_range_hsg(
        hsg_props: HSGProperties,
        temperatures: list[Temperature],
        **kwargs
) -> Optional[list[ComponentGibbsFreeEnergy]]:
    """
    Calculate the Gibbs free energy (J/mol) over a range of temperatures (K) for a component using an existing HSGProperties instance.

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    temperatures : list[Temperature]
        The list of temperatures at which to calculate the Gibbs free energy of formation.
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[list[ComponentGibbsFreeEnergy]]
        A list of ComponentGibbsFreeEnergy objects containing the Gibbs free energy of formation values, units, and symbols, or None if calculation fails.

    Notes
    -----
    - The function uses the provided HSGProperties instance to compute the Gibbs free energy of formation at the specified range of temperatures.
    - The Gibbs free energy of formation symbol must be consistent with the defined unit in the configuration which is GiEnFo_IG.
    - Reference temperature is set to 298.15 K.
    - All Gibbs energy results are provided in J/mol.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
            return None

        if not all(isinstance(temp, Temperature) for temp in temperatures):
            logger.error("Invalid temperatures provided.")
            return None

        # NOTE: calculate
        GiEnFo_results = hsg_props.calc_gibbs_free_energy_range(
            temperatures=temperatures
        )

        return GiEnFo_results
    except Exception as e:
        logger.error(
            f"Error calculating Gibbs free energy of formation range using HSGProperties: {e}")
        return None

# SECTION: entropy change calculations
# ! entropy change between two temperatures and pressures


@measure_time
def calc_dEnt(
        component: Component,
        model_source: ModelSource,
        temperature_initial: Temperature,
        temperature_final: Temperature,
        pressure_initial: Pressure,
        pressure_final: Pressure,
        phase: Literal['IG', 'LIQ', 'SOL'],
        component_key: ComponentKey = 'Name-Formula',
        **kwargs
) -> Optional[ComponentEntropyChange]:
    """
    Calculate the entropy change (J/mol.K) between two temperatures (K) and pressures in any unit for a component.

    The entropy change is calculated by integrating the heat capacity equation between the two specified temperatures and accounting for pressure changes based on the phase of the component. The result is provided in J/mol.K.

    The calculation considers the following:
    - For ideal gases, the pressure change contribution is included using the relation ΔS = nR ln(P2/P1).
    - For liquids, the pressure change contribution is typically negligible and may be omitted.
    - For solids, the pressure change contribution is typically negligible and may be omitted.

    The expression for entropy change is given by:
        ΔS = ∫(Cp/T) dT + ln(P2/P1) * R  (for ideal gases)
        ΔS = ∫(Cp/T) dT                    (for liquids)
        ΔS = ∫(Cp/T) dT                    (for solids)

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
    phase : Literal['IG', 'LIQ', 'SOL'], optional
        The phase of the component ('IG' for ideal gas, 'LIQ' for liquid, 'SOL' for solid).
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

        if phase not in ['IG', 'LIQ', 'SOL']:
            logger.error(
                "Invalid phase provided. Must be 'IG', 'LIQ', or 'SOL'.")
            return None

        # SECTION: Prepare source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

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

# ! entropy change between two temperatures and pressures using hsg properties


def calc_dEnt_hsg(
        hsg_props: HSGProperties,
        temperature_initial: Temperature,
        temperature_final: Temperature,
        pressure_initial: Pressure,
        pressure_final: Pressure,
        phase: Literal['IG', 'LIQ', 'SOL'],
        **kwargs
) -> Optional[ComponentEntropyChange]:
    """
    Calculate the entropy change (J/mol.K) between two temperatures (K) and pressures in any unit for a component using an existing HSGProperties instance.

    The entropy change is calculated by integrating the heat capacity equation between the two specified temperatures and accounting for pressure changes based on the phase of the component. The result is provided in J/mol.K.

    The calculation considers the following:
    - For ideal gases, the pressure change contribution is included using the relation ΔS = nR ln(P2/P1).
    - For liquids, the pressure change contribution is typically negligible and may be omitted.
    - For solids, the pressure change contribution is typically negligible and may be omitted.

    The expression for entropy change is given by:
        ΔS = ∫(Cp/T) dT + ln(P2/P1) * R  (for ideal gases)
        ΔS = ∫(Cp/T) dT                    (for liquids)
        ΔS = ∫(Cp/T) dT                    (for solids)

    Parameters
    ----------
    hsg_props : HSGProperties
        An instance of the HSGProperties class for the component of interest.
    temperature_initial : Temperature
        The initial temperature.
    temperature_final : Temperature
        The final temperature.
    pressure_initial : Pressure
        The initial pressure.
    pressure_final : Pressure
        The final pressure.
    phase : Literal['IG', 'LIQ', 'SOL'], optional
        The phase of the component ('IG' for ideal gas, 'LIQ' for liquid, 'SOL' for solid).
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
    - The function uses the provided HSGProperties instance to compute the entropy change.
    - The entropy change is calculated by integrating the heat capacity equation between the two specified temperatures and accounting for pressure changes.
    - The result is provided in J/mol.K.
    - R is the universal gas constant (8.3145 J/mol.K).
    - Pressures can be provided in any unit as they will be converted internally.
    """
    try:
        if not isinstance(hsg_props, HSGProperties):
            logger.error("Invalid HSGProperties instance provided.")
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

        if phase not in ['IG', 'LIQ', 'SOL']:
            logger.error(
                "Invalid phase provided. Must be 'IG', 'LIQ', or 'SOL'.")
            return None

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
            f"Error calculating entropy change using HSGProperties: {e}")
        return None


@measure_time
def calc_En_mix(
        components: List[Component],
        model_source: ModelSource,
        temperature: Temperature,
        pressure: Pressure,
        reference: Literal['IG', 'None'] = 'IG',
        departure_enthalpy: Optional[CustomProp] = None,
        excess_enthalpy: Optional[CustomProp] = None,
        component_key: ComponentKey = 'Name-Formula',
        output_unit: Optional[str] = None,
        **kwargs
) -> Optional[MixtureEnthalpyResult]:
    """
    Calculate the mixture enthalpy at a given temperature (K) and return the result is provided in J/mol.

    The mixture enthalpy is calculated based on the individual component enthalpies and their mole fractions. The function initializes the HSGMixture class and uses it to compute the mixture enthalpy.

    Parameters
    ----------
    components : List[Component]
        The list of chemical components in the mixture.
    model_source : ModelSource
        The source model containing necessary data.
    temperature : Temperature
        The temperature of the mixture.
    pressure : Pressure
        The pressure of the mixture.
    reference : Literal['IG', 'None']
        The reference state for enthalpy calculation.
    departure_enthalpy : Optional[CustomProp], optional
        Custom departure enthalpy property, by default None.
    excess_enthalpy : Optional[CustomProp], optional
        Custom excess enthalpy property, by default None.
    component_key : Literal[..., optional]
        The key to identify the components, by default 'Name-Formula'.
    output_unit : Optional[str], optional
        The desired output unit for mixture enthalpy, by default None (J/mol).
    **kwargs
        Additional keyword arguments.
        - mode : Literal['silent', 'log', 'attach'], optional
            Mode for time measurement logging. Default is 'log'.

    Returns
    -------
    Optional[MixtureEnthalpyResult]
        A MixtureEnthalpyResult object containing the mixture enthalpy value, unit, and symbol, or None if calculation fails.

    Notes
    -----
    - The function initializes the HSGMixture class and uses it to compute the mixture enthalpy.
    - The mixture enthalpy is calculated based on the individual component enthalpies and their mole fractions.
    - The result is provided in J/mol.
    - The calculation may involve ideal and non-ideal contributions depending on the phase.
    - Only 'IG' and 'LIQ' phases are currently supported.
    - Flash calculations are not performed; the phase must be specified.
    """
    try:
        # SECTION: Input validation
        if not all(isinstance(comp, Component) for comp in components):
            logger.error("Invalid components provided.")
            return None

        if not isinstance(model_source, ModelSource):
            logger.error("Invalid model_source provided.")
            return None

        if not isinstance(temperature, Temperature):
            logger.error("Invalid temperature provided.")
            return None

        # SECTION: Prepare source
        Source_ = Source(
            model_source=model_source,
            component_key=component_key
        )

        # SECTION: Initialize HSGMixture
        hsg_mixture = HSGMixture(
            components=components,
            source=Source_,
            component_key=component_key
        )

        # NOTE: calculate mixture enthalpy
        mixture_enthalpy = hsg_mixture.calc_mixture_enthalpy(
            temperature=temperature,
            reference=reference,
            departure_enthalpy=departure_enthalpy,
            excess_enthalpy=excess_enthalpy,
        )

        # >> check
        if mixture_enthalpy is None:
            logger.error("Mixture enthalpy calculation failed.")
            return None

        # set
        mixture_enthalpy_value = mixture_enthalpy['value']
        mixture_enthalpy_unit = mixture_enthalpy['unit']

        # SECTION: convert unit if needed
        if (
            output_unit is not None and
            mixture_enthalpy_unit != output_unit
        ):
            # >> convert
            converted_value = pycuc.convert_from_to(
                value=mixture_enthalpy_value,
                from_unit=mixture_enthalpy_unit,
                to_unit=output_unit
            )
            # >> update
            mixture_enthalpy_value = converted_value
            mixture_enthalpy_unit = output_unit

        # NOTE: prepare result
        res = {
            'temperature': temperature,
            'pressure': pressure,
            'reference': str(reference),
            'value': mixture_enthalpy_value,
            'unit': mixture_enthalpy_unit,
            'symbol': 'En_mix'
        }

        # set
        res = MixtureEnthalpyResult(**res)
        return res
    except Exception as e:
        logger.error(
            f"Error calculating mixture enthalpy: {e}")
        return None
