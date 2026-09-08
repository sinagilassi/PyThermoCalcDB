# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional, cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import (
    AnnotatedValue,
    Component,
    ComponentKey,
    CustomProp,
)
from pythermodb_settings.utils import (
    config_components_values,
    to_annotated_value,
)

# locals
from ..utils.conversions import _to_mass, _to_moles, _to_units

# NOTE: logger setup
logger = logging.getLogger(__name__)


# ======================================================================
# *** Helper functions
# ======================================================================
def _validate_moles_and_mass(
    moles: NDArray[np.float64],
    mass: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Validate the shapes and values of moles and solvent mass arrays.

    Parameters
    ----------
    moles : NDArray[np.number]
        Array of component moles.
    mass : NDArray[np.number]
        Array of solvent masses.

    Returns
    -------
    NDArray[np.number]
        Validated mass array. A 1-D per-state mass for 2-D moles is returned
        as ``(n_states, 1)`` so NumPy broadcasts across components.

    Raises
    ------
    ValueError
        If the shapes of moles and mass are incompatible or if any mass is
        not finite and greater than zero.
    """
    if moles.ndim not in (1, 2):
        raise ValueError("component_moles must be a 1-D or 2-D array.")

    if mass.ndim > 2:
        raise ValueError(
            "solvent_mass must be a scalar, 1-D array, or 2-D array.")

    if moles.ndim == 1:
        if mass.ndim not in (0, 1):
            raise ValueError(
                "For 1-D component_moles, solvent_mass must be a scalar "
                "or have the same shape as component_moles."
            )
        if mass.ndim == 1 and mass.shape != moles.shape:
            raise ValueError(
                "For 1-D component_moles, solvent_mass must be a scalar "
                "or have the same shape as component_moles."
            )

    elif moles.ndim == 2:
        if mass.ndim == 0:
            pass
        elif mass.ndim == 1:
            if mass.shape[0] != moles.shape[0]:
                raise ValueError(
                    "For 2-D component_moles, 1-D solvent_mass must have "
                    "one entry per state."
                )
            mass = mass[:, None]

        elif mass.ndim == 2:
            if mass.shape not in (moles.shape, (moles.shape[0], 1)):
                raise ValueError(
                    "For 2-D component_moles, solvent_mass must be a scalar, "
                    "have shape (n_states,), shape (n_states, 1), or the "
                    "same shape as component_moles."
                )

    if not np.all(np.isfinite(mass)):
        raise ValueError("solvent_mass must contain finite values.")

    if np.any(mass <= 0):
        raise ValueError("solvent_mass must be greater than zero.")

    return mass


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================

# ! ::: Molality from array-like inputs
def _calc_molalities(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solvent_mass: float | int | NDArray[np.number],
) -> NDArray[np.float64]:
    """
    Calculate molalities using NumPy vectorization.

    Parameters
    ----------
    component_moles : Sequence[float | int] | NDArray[np.number]
        Component mole amounts. May be a Python sequence or a 1-D/2-D NumPy
        array.
    solvent_mass : float | int | NDArray[np.number]
        Solvent mass. Must be a scalar or have a compatible shape with
        ``component_moles``.

    Returns
    -------
    NDArray[np.floating]
        Molalities with the same shape as ``component_moles``.
    """
    # set
    moles: NDArray[np.float64] = np.asarray(component_moles, dtype=np.float64)
    mass: NDArray[np.float64] = np.asarray(solvent_mass, dtype=np.float64)

    # validate
    mass = _validate_moles_and_mass(moles, mass)

    return cast(NDArray[np.float64], moles / mass)


# ! ::: Molality from sequence
def _calc_molalities_from_sequence(
    component_moles: Sequence[float | int],
    solvent_mass: float,
) -> list[float]:
    """
    Calculate the molality of each component in a solution given the component
    moles and the solvent mass.

    Parameters
    ----------
    component_moles : Sequence[float | int]
        A sequence of moles for each component.
    solvent_mass : float
        The solvent mass.

    Returns
    -------
    list[float]
        A list of molality values for each component.
    """
    # calc
    return _calc_molalities(
        component_moles,
        solvent_mass
    ).tolist()


# ! ::: Molality from mapping
def _calc_molalities_from_mapping(
    component_moles: Mapping[str, float | int],
    solvent_mass: float,
) -> dict[str, float]:
    """
    Calculate the molality of each component in a solution given the component
    moles and the solvent mass.

    Parameters
    ----------
    component_moles : Mapping[str, float | int]
        A mapping of component names to their respective moles.
    solvent_mass : float
        The solvent mass.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molality
        values.
    """
    # calc
    molalities_ = _calc_molalities(
        component_moles=list(component_moles.values()),
        solvent_mass=solvent_mass
    )

    # to dict
    component_molality_dict = dict(
        zip(component_moles.keys(), molalities_.tolist())
    )

    return component_molality_dict


# ! ::: Molality with solvent mass as CustomProp
def _calc_molalities_from_props(
    component_moles: Mapping[str, CustomProp],
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
) -> dict[str, float]:
    """
    Calculate the molality of each component in a solution given the component
    moles and the solvent mass as a CustomProp.

    Parameters
    ----------
    component_moles : Mapping[str, CustomProp]
        A mapping of component names to their respective mole amounts.
    solvent_mass : CustomProp
        The solvent mass as a CustomProp object.
    output_unit : str, optional
        The unit for the output molality values. Defaults to 'mol/kg'.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molality
        values.

    Notes
    -----
    - The solvent mass is expected to be provided as a CustomProp object. If
    the output_unit is not specified, it defaults to mol/kg.
    - Component mole values are expected to be CustomProp objects so their
    units can be converted to the mole unit from output_unit.
    """
    # SECTION: set default units for moles and mass
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    mass_unit = units_[1]

    # NOTE: component moles
    # ! convert component moles to the specified unit if necessary
    component_moles_dict: dict[str, float] = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
    )

    # NOTE: solvent mass unit should match output unit denominator
    # ! convert solvent mass to the specified unit
    solvent_mass_scalar = _to_mass(
        solvent_mass=solvent_mass,
        output_unit=mass_unit
    )

    # SECTION: calculate molality for each component
    return _calc_molalities_from_mapping(
        component_moles=component_moles_dict,
        solvent_mass=solvent_mass_scalar,
    )


# ! ::: Molality with component ID mapping and sorting
def _calc_component_molalities_from_props(
    component_moles: Mapping[str, CustomProp],
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """
    Calculate the molality of each component in a solution given the component
    moles and the solvent mass as a CustomProp.

    Parameters
    ----------
    component_moles : Mapping[str, CustomProp]
        A mapping of component names to CustomProp objects representing the
        moles.
    solvent_mass : CustomProp
        The solvent mass as a CustomProp object.
    output_unit : str, optional
        The unit for the output molality values. Defaults to 'mol/kg'.
    components : Optional[Sequence[Component]], optional
        A sequence of Component objects to map the component moles to, by
        default None.
    component_key : Optional[ComponentKey], optional
        The key to use for mapping component moles to components, by default
        None.
    case_sensitive : bool, optional
        Whether the component mapping should be case sensitive, by default
        True.
    sort_by_components_order : bool, optional
        Whether to sort the component molalities by the order of components,
        by default True.

    Returns
    -------
    dict[str, float]
        A dictionary mapping component names to their respective molality
        values.
    """
    # SECTION: Unit validation
    units_ = _to_units(output_unit)
    # >> set
    mole_unit = units_[0]
    mass_unit = units_[1]

    # SECTION: convert component moles to moles if they are CustomProp objects
    # ! component moles
    component_moles_dict = _to_moles(
        component_moles=component_moles,
        output_unit=mole_unit
    )

    # ! mass
    solvent_mass_scalar = _to_mass(
        solvent_mass=solvent_mass,
        output_unit=mass_unit
    )

    # SECTION: get component values
    if component_key is not None:
        if not components:
            logger.error(
                "Component key is provided but components list is empty."
            )
            raise ValueError(
                "Component key is provided but components list is empty."
            )

        component_values: tuple[
            dict[str, Any],
            list[Any]
        ] | None = config_components_values(
            values=component_moles_dict,
            components=list(components),
            component_key=component_key,
            case_sensitive=case_sensitive,
            sort_by_components_order=sort_by_components_order
        )
        # >> check
        if component_values is None:
            logger.error(
                "Failed to configure component values."
            )
            raise ValueError(
                "Failed to configure component values."
            )

        # unpack
        component_values_dict, _ = component_values
    else:
        component_values_dict = component_moles_dict

    # SECTION: calculate molality for each component
    return _calc_molalities_from_mapping(
        component_moles=component_values_dict,
        solvent_mass=solvent_mass_scalar,
    )


# ======================================================================
# *** Public annotated API
# ======================================================================

# ::: annotated for numpy array
@calculation_info(
    name="molality",
    description="Calculate the molality of each component in a solution.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solvent_mass": "Mass of the solvent."
    },
    outputs={
        "molality": "Molality of each component in the solution."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solvent_mass are already expressed on that molality basis.",
    ),
    tags=(
        "molality",
        "array_like",
        "numpy",
        "component_wise",
        "numeric",
    )
)
def _molality_annotated(
    component_moles: Sequence[float | int] | NDArray[np.number],
    solvent_mass: float | int | NDArray[np.number],
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.floating]]:
    """Calculate annotated molality values from array-like inputs."""
    return to_annotated_value(
        _calc_molalities(
            component_moles=component_moles,
            solvent_mass=solvent_mass
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molalities",
    )


# ::: annotated for sequence
@calculation_info(
    name="molality",
    description="Calculate the molality of each component in a solution.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Component moles in the solution.",
        "solvent_mass": "Mass of the solvent."
    },
    outputs={
        "molality": "Molality of each component in the solution."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solvent_mass are already expressed on that molality basis.",
    ),
    tags=(
        "molality",
        "sequence",
        "component_wise",
        "numeric",
    )
)
def _molality_1_annotated(
    component_moles: Sequence[float],
    solvent_mass: float,
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None
) -> AnnotatedValue[list[float]]:
    """Calculate annotated molality values from a sequence of component moles."""
    return to_annotated_value(
        _calc_molalities_from_sequence(
            component_moles=component_moles,
            solvent_mass=solvent_mass
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molalities_from_sequence"
    )


# ::: annotated for mapping
@calculation_info(
    name="molality",
    description="Calculate the molality of each keyed component in a solution.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Mapping of component identifiers to component moles in the solution.",
        "solvent_mass": "Mass of the solvent."
    },
    outputs={
        "molality": "Mapping of component identifiers to molality values."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "Numeric inputs do not carry unit metadata, so the annotated result unit is not defined by default.",
        "Pass unit only when component_moles and solvent_mass are already expressed on that molality basis.",
    ),
    tags=(
        "numeric",
        "molality",
        "mapping",
        "component_aware",
        "keyed",
    )
)
def _molality_2_annotated(
    component_moles: Mapping[str, float | int],
    solvent_mass: float,
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molality values from a component-mole mapping."""
    return to_annotated_value(
        _calc_molalities_from_mapping(
            component_moles=component_moles,
            solvent_mass=solvent_mass
        ),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molalities_from_mapping"
    )


# ::: annotated for mapping with custom properties
@calculation_info(
    name="molality",
    description="Calculate keyed molality values from unit-aware component moles and solvent mass.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Mapping of component identifiers to unit-aware component mole amounts.",
        "solvent_mass": "Unit-aware mass of the solvent.",
        "output_unit": "Molality unit used to normalize component moles and solvent mass."
    },
    outputs={
        "molality": "Mapping of component identifiers to molality values in output_unit."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "The output_unit must be a ratio such as mol/kg with amount in the numerator and mass in the denominator.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "molality",
        "mapping",
        "unit_aware",
        "component_aware",
        "unit_conversion",
    )
)
def _molality_3_annotated(
    component_moles: Mapping[str, CustomProp],
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molality values from unit-aware component moles."""
    # SECTION: set default unit for output if not provided
    if unit is None:
        unit = output_unit

    # check unit & output unit consistency
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    return to_annotated_value(
        _calc_molalities_from_props(
            component_moles=component_moles,
            solvent_mass=solvent_mass,
            output_unit=output_unit
        ),
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_molalities_from_props"
    )


# ::: annotated for component mapping with custom properties
@calculation_info(
    name="component_molality",
    description="Calculate unit-aware component molalities with optional component-key remapping and ordering.",
    equation="molality = component_moles / solvent_mass",
    inputs={
        "component_moles": "Mapping of component identifiers to unit-aware component mole amounts.",
        "solvent_mass": "Unit-aware mass of the solvent.",
        "output_unit": "Molality unit used to normalize component moles and solvent mass.",
        "components": "Optional component definitions used to resolve and order component identifiers.",
        "component_key": "Optional component key used for identifier matching."
    },
    outputs={
        "component_molality": "Mapping of resolved component identifiers to molality values in output_unit."
    },
    aliases=(
        "molal concentration",
        "amount concentration by solvent mass",
    ),
    notes=(
        "The output_unit must be a ratio such as mol/kg with amount in the numerator and mass in the denominator.",
        "The annotated result unit is output_unit; an explicitly supplied unit must match output_unit.",
    ),
    tags=(
        "component_molality",
        "molality",
        "mapping",
        "component_wise",
        "unit_aware",
        "unit_conversion",
        "component_key",
        "component_ordering",
    )
)
def _molality_4_annotated(
    component_moles: Mapping[str, CustomProp],
    solvent_mass: CustomProp,
    output_unit: str = 'mol/kg',
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "molality",
    description: str = "Calculate the molality of each component in a solution.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Calculate annotated molality values for mapped components with units."""
    # SECTION: set default unit for output if not provided
    if unit is None:
        unit = output_unit

    # check unit & output unit consistency
    if unit != output_unit:
        raise ValueError(
            f"Mismatch between unit ({unit}) and output_unit ({output_unit})"
        )

    res = _calc_component_molalities_from_props(
        component_moles=component_moles,
        solvent_mass=solvent_mass,
        output_unit=output_unit,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order
    )

    # convert the result to an annotated value
    return to_annotated_value(
        value=res,
        name=name,
        description=description,
        unit=output_unit,
        symbol=symbol,
        implementation="_calc_component_molalities_from_props"
    )


# ======================================================================
# *** Aliases
# ======================================================================
# >> molalities
calc_molalities = _molality_annotated

# >> molalities from sequence
calc_molalities_from_sequence = _molality_1_annotated

# >> molalities from mapping
calc_molalities_from_mapping = _molality_2_annotated

# >> molalities with units
calc_molalities_from_props = _molality_3_annotated

# >> component molalities
calc_component_molalities_from_props = _molality_4_annotated

# all
__all__ = [
    "_calc_molalities",
    "calc_molalities",
    "_calc_molalities_from_sequence",
    "calc_molalities_from_sequence",
    "_calc_molalities_from_mapping",
    "calc_molalities_from_mapping",
    "_calc_molalities_from_props",
    "calc_molalities_from_props",
    "_calc_component_molalities_from_props",
    "calc_component_molalities_from_props",
]
