# import libs
from typing import Any, Dict, Callable
from pydantic import BaseModel, Field, ConfigDict
from pyThermoDB.core import TableEquation
from pythermodb_settings.models import Temperature

# NOTE: Equation Builder Result


class ComponentEquation(BaseModel):
    '''
    Component Equation Result Model

    Attributes
    ----------
    equation : TableEquation
        The equation value.
    body : str
        The body of the equation.
    args : Dict[str, Any]
        The equation arguments.
    arg_symbols : Dict[str, Any]
        The equation argument symbols.
    returns : Dict[str, Any]
        The equation returns.
    return_symbols : Dict[str, Any]
        The equation return symbols.
    '''
    model_config = ConfigDict(arbitrary_types_allowed=True)

    equation: TableEquation = Field(
        ...,
        description="The equation value."
    )
    body: str = Field(
        ...,
        description="The body of the equation."
    )
    args: Dict[str, Any] = Field(
        default_factory=dict,
        description="The equation arguments."
    )
    arg_symbols: Dict[str, Any] = Field(
        default_factory=dict,
        description="The equation argument symbols."
    )
    returns: Dict[str, Any] = Field(
        default_factory=dict,
        description="The equation returns."
    )
    return_symbol: str = Field(
        default_factory=str,
        description="The equation return symbols."
    )


class ComponentEquationSource(BaseModel):
    '''
    Equation Builder Result Model

    Attributes
    ----------
    value : TableEquation
        The equation value.
    args : Dict[str, Any]
        The equation arguments.
    arg_symbols : Dict[str, Any]
        The equation argument symbols.
    returns : Dict[str, Any]
        The equation returns.
    return_symbols : Dict[str, Any]
        The equation return symbols.
    '''
    model_config = ConfigDict(arbitrary_types_allowed=True)

    value: TableEquation = Field(
        ...,
        description="The equation value."
    )
    args: Dict[str, Any] = Field(
        default_factory=dict,
        description="The equation arguments."
    )
    arg_symbols: Dict[str, Any] = Field(
        default_factory=dict,
        description="The equation argument symbols."
    )
    returns: Dict[str, Any] = Field(
        default_factory=dict,
        description="The equation returns."
    )
    return_symbols: Dict[str, Any] = Field(
        default_factory=dict,
        description="The equation return symbols."
    )


# NOTE: enthalpy of formation result model
class ComponentEnthalpyOfFormation(BaseModel):
    '''
    Component Enthalpy of Formation Result Model

    Attributes
    ----------
    temperature : Temperature
        The temperature at which the enthalpy of formation is calculated.
    value : float
        The enthalpy of formation value.
    unit : str
        The enthalpy of formation unit.
    symbol : str
        The enthalpy of formation symbol.
    '''
    model_config = ConfigDict(arbitrary_types_allowed=True)

    temperature: Temperature = Field(
        ...,
        description="The temperature at which the enthalpy of formation is calculated."
    )

    value: float = Field(
        ...,
        description="The enthalpy of formation value."
    )
    unit: str = Field(
        ...,
        description="The enthalpy of formation unit."
    )
    symbol: str = Field(
        ...,
        description="The enthalpy of formation symbol."
    )


# NOTE: gibbs energy of formation result model


class ComponentGibbsEnergyOfFormation(BaseModel):
    '''
    Component Gibbs Energy of Formation Result Model

    Attributes
    ----------
    temperature : Temperature
        The temperature at which the Gibbs energy of formation is calculated.
    value : float
        The Gibbs energy of formation value.
    unit : str
        The Gibbs energy of formation unit.
    symbol : str
        The Gibbs energy of formation symbol.
    '''
    model_config = ConfigDict(arbitrary_types_allowed=True)

    temperature: Temperature = Field(
        ...,
        description="The temperature at which the Gibbs energy of formation is calculated."
    )

    value: float = Field(
        ...,
        description="The Gibbs energy of formation value."
    )
    unit: str = Field(
        ...,
        description="The Gibbs energy of formation unit."
    )
    symbol: str = Field(
        ...,
        description="The Gibbs energy of formation symbol."
    )
