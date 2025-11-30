# import libs
from typing import Any, Dict, Optional, Literal
from pydantic import BaseModel, Field, ConfigDict
from pyThermoDB.core import TableEquation

# NOTE: Equation Builder Result


class ComponentEquationSource(BaseModel):
    '''
    Equation Builder Result Model
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
