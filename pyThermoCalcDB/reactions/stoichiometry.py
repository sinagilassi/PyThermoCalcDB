"""Stoichiometric reaction-extent public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue
from pythermodb_settings.utils import to_annotated_value

# locals
from .core.stoichiometry import (
    _calc_component_moles_from_reaction_extent,
    _calc_component_moles_from_reaction_extents,
)


# SECTION: Public wrappers
def calc_component_moles_from_reaction_extent(
    initial_moles,
    stoichiometric_coefficients,
    reaction_extent,
    *,
    name: str = "component_moles",
    description: str = "Calculate component moles from a single reaction extent.",
    unit: str | None = "mol",
    symbol: str | None = "N",
) -> AnnotatedValue[list[float]]:
    """Calculate annotated component moles from one reaction extent."""
    result = _calc_component_moles_from_reaction_extent(
        initial_moles,
        stoichiometric_coefficients,
        reaction_extent,
    )
    return to_annotated_value(
        result.tolist(),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_component_moles_from_reaction_extent",
    )


def calc_component_moles_from_reaction_extents(
    initial_moles,
    stoichiometric_matrix,
    reaction_extents,
    *,
    name: str = "component_moles",
    description: str = "Calculate component moles from multiple reaction extents.",
    unit: str | None = "mol",
    symbol: str | None = "N",
) -> AnnotatedValue[list[float]]:
    """Calculate annotated component moles from multiple reaction extents."""
    result = _calc_component_moles_from_reaction_extents(
        initial_moles,
        stoichiometric_matrix,
        reaction_extents,
    )
    return to_annotated_value(
        result.tolist(),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_component_moles_from_reaction_extents",
    )


__all__ = [
    "calc_component_moles_from_reaction_extent",
    "calc_component_moles_from_reaction_extents",
]
