# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Optional
import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import Component, ComponentKey, AnnotatedValue
from pythermodb_settings.utils import (
    to_annotated_value,
)
from pythermodb_settings.decorators import calculation_info
from .core.fractions import (
    _calc_fractions,
    _calc_fractions_from_mapping,
    _calc_component_fractions,
)

# NOTE: logger set
logger = logging.getLogger(__name__)


# ======================================================================
# *** Public annotated API
# ======================================================================

# ! ::: annotated for numpy array

@calculation_info(
    name="fraction",
    description="Calculate the fraction of each value relative to the total.",
    equation="fraction_i = value_i / sum(values)",
    inputs={
        "values": "Component values to normalize into fractions."
    },
    outputs={
        "fraction": "Fractions of each value relative to the total."
    },
    aliases=(
        "composition fraction",
        "normalized fraction",
    ),
    notes=(
        "Fractions are unitless.",
        "All input values must be non-negative and their total must be positive.",
    ),
    tags=(
        "fraction",
        "array_like",
        "numpy",
        "component_wise",
        "numeric",
        "unitless",
    )
)
def _fraction_annotated(
        values: Sequence[float | int] | NDArray[np.number],
        *,
        name: str = "fraction",
        description: str = "Calculate the fraction of each value relative to the total.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.floating]]:
    """Calculate annotated fractions from array-like inputs.

    Parameters
    ----------
    values : Sequence[float | int] | NDArray[np.number]
        Values to normalize. May be a Python sequence or a NumPy array.
    name : str, optional
        The name stored in the annotated result. Defaults to ``"fraction"``.
    description : str, optional
        The description stored in the annotated result.
    unit : str, optional
        The unit stored in the annotated result. Fractions are unitless by
        default.
    symbol : str, optional
        The symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[NDArray[np.floating]]
        Fractions with the same shape as ``values``.
    """
    return to_annotated_value(
        _calc_fractions(values=values),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_fractions",
    )

# ! ::: annotated for sequence


@calculation_info(
    name="fraction",
    description="Calculate the fraction of each value relative to the total.",
    equation="fraction_i = value_i / sum(values)",
    inputs={
        "values": "Sequence of component values to normalize into fractions."
    },
    outputs={
        "fraction": "Fractions of each value relative to the total."
    },
    aliases=(
        "composition fraction",
        "normalized fraction",
    ),
    notes=(
        "Fractions are unitless.",
        "All input values must be non-negative and their total must be positive.",
    ),
    tags=(
        "fraction",
        "sequence",
        "component_wise",
        "numeric",
        "unitless",
    )
)
def _fraction_1_annotated(
        values: Sequence[float | int],
        *,
        name: str = "fraction",
        description: str = "Calculate the fraction of each value relative to the total.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Calculate annotated fractions from a sequence of values."""
    return to_annotated_value(
        _calc_fractions(values=values).tolist(),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_fractions",
    )

# ! ::: annotated for mapping


@calculation_info(
    name="fraction",
    description="Calculate the fraction of each keyed value relative to the total.",
    equation="fraction_i = value_i / sum(values)",
    inputs={
        "values": "Mapping of component identifiers to values to normalize into fractions."
    },
    outputs={
        "fraction": "Mapping of component identifiers to fractions."
    },
    aliases=(
        "composition fraction",
        "normalized fraction",
    ),
    notes=(
        "Fractions are unitless.",
        "All input values must be non-negative and their total must be positive.",
    ),
    tags=(
        "numeric",
        "fraction",
        "mapping",
        "component_aware",
        "keyed",
        "unitless",
    )
)
def _fraction_2_annotated(
        values: Mapping[str, float | int],
        *,
        name: str = "fraction",
        description: str = "Calculate the fraction of each value relative to the total.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated fractions from a value mapping."""
    return to_annotated_value(
        _calc_fractions_from_mapping(values=values),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_fractions_from_mapping",
    )

# ! ::: annotated for component mapping


@calculation_info(
    name="component_fraction",
    description="Calculate component fractions with optional component-key remapping and ordering.",
    equation="fraction_i = value_i / sum(values)",
    inputs={
        "values": "Mapping of component identifiers to values to normalize into fractions.",
        "components": "Optional component definitions used to resolve and order component identifiers.",
        "component_key": "Optional component key used for identifier matching."
    },
    outputs={
        "component_fraction": "Mapping of resolved component identifiers to fractions."
    },
    aliases=(
        "composition fraction",
        "normalized fraction",
    ),
    notes=(
        "Fractions are unitless.",
        "All input values must be non-negative and their total must be positive.",
    ),
    tags=(
        "component_fraction",
        "fraction",
        "mapping",
        "component_wise",
        "component_key",
        "component_ordering",
        "unitless",
    )
)
def _fraction_3_annotated(
        values: Mapping[str, float | int],
        components: Optional[Sequence[Component]] = None,
        component_key: Optional[ComponentKey] = None,
        case_sensitive: bool = True,
        sort_by_components_order: bool = True,
        *,
        name: str = "fraction",
        description: str = "Calculate the fraction of each value relative to the total.",
        unit: str | None = None,
        symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated component fractions from a value mapping."""
    return to_annotated_value(
        _calc_component_fractions(
            values=values,
            components=list(components) if components is not None else None,
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_component_fractions",
    )


# ======================================================================
# *** Aliases
# ======================================================================
# >> fractions
calc_fractions = _fraction_annotated

# >> fractions from sequence
calc_fractions_from_sequence = _fraction_1_annotated

# >> fractions from mapping
calc_fractions_from_mapping = _fraction_2_annotated

# >> component fractions
calc_component_fractions = _fraction_3_annotated


# all
__all__ = [
    "calc_fractions",
    "calc_fractions_from_sequence",
    "calc_fractions_from_mapping",
    "calc_component_fractions",
]
