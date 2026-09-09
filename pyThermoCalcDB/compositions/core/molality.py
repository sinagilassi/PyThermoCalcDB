# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Any, Optional, cast

import numpy as np
from numpy.typing import NDArray
from pythermodb_settings.models import (
    Component,
    ComponentKey,
    CustomProp,
)
from pythermodb_settings.utils import (
    config_components_values,
)

# locals
from ...utils.conversions import _to_mass, _to_moles, _to_units

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


# all
__all__ = [
    "_calc_molalities",
    "_calc_molalities_from_sequence",
    "_calc_molalities_from_mapping",
    "_calc_molalities_from_props",
    "_calc_component_molalities_from_props",
]
