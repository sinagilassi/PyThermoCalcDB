"""Composition basis conversion functions.

These helpers perform deterministic arithmetic conversions between common
composition bases. Inputs may be plain numeric values or ``CustomProp`` values.
Mapping inputs are normalized through ``pythermodb_settings.utils.to_amounts``.

Unit convention
---------------
Arguments named ``output_*_unit`` are optional and define the unit used to
normalize the matching input before calculation. When an ``output_*_unit`` is
``None``, the input value is used as-is. When it is provided, numeric values are
assumed to already use that unit and ``CustomProp`` values are converted to that
unit with ``pycuc.convert_from_to``. The low-level ``_calc_*`` functions return
NumPy arrays or floats. The ``_calc_*_from_sequence`` adapters return lists,
and the ``_calc_*_from_mapping`` and ``_calc_*_from_props`` adapters return
dictionaries keyed by component.
"""

from collections.abc import Mapping, Sequence
from typing import Optional, List, cast
import numpy as np
from numpy.typing import NDArray
# >> pythermodb-settings
from pythermodb_settings.models import CustomProp, Component, ComponentKey
from pythermodb_settings.models.units import UnitConversionFn
from pythermodb_settings.utils.quantity import to_dict
# locals
from ...utils.conversions import (
    _configure_component_values,
    _resolve_unit_conversion_fn,
    _validate_same_keys,
    _pos,
)


# ======================================================================
# *** Internal deterministic calculations
# ======================================================================


def _as_float_array(
    values: float | int | Sequence[float | int] | NDArray[np.number],
    identifier: str,
) -> NDArray[np.float64]:
    """Normalize a numeric scalar, sequence, or NumPy array to float64.

    Parameters
    ----------
    values : float | int | Sequence[float | int] | NDArray[np.number]
        Numeric input to normalize.
    identifier : str
        Name used in validation error messages.

    Returns
    -------
    NDArray[np.float64]
        Finite scalar, 1-D, or 2-D values as a float64 NumPy array.
    """
    array_values: NDArray[np.float64] = np.asarray(values, dtype=np.float64)

    if array_values.ndim > 2:
        raise ValueError(
            f"{identifier} must be a scalar, 1-D array, or 2-D array.")

    if array_values.size == 0:
        raise ValueError(f"{identifier} must not be empty.")

    if not np.all(np.isfinite(array_values)):
        raise ValueError(f"{identifier} must contain finite values.")

    return array_values


def _validate_same_array_shape(
    left_values: NDArray[np.float64],
    right_values: NDArray[np.float64],
    left_identifier: str,
    right_identifier: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Broadcast two numeric arrays for component-wise calculations.

    Parameters
    ----------
    left_values, right_values : NDArray[np.float64]
        Numeric arrays that must be broadcast-compatible.
    left_identifier, right_identifier : str
        Names used in validation error messages.

    Returns
    -------
    tuple[NDArray[np.float64], NDArray[np.float64]]
        Broadcast arrays with matching shapes.
    """
    try:
        left_values, right_values = np.broadcast_arrays(
            left_values, right_values)
    except ValueError as exc:
        raise ValueError(
            f"{left_identifier} and {right_identifier} must have "
            "broadcast-compatible shapes."
        ) from exc

    return (
        cast(NDArray[np.float64], left_values),
        cast(NDArray[np.float64], right_values),
    )


def _calc_mole_fraction_to_mass_fraction(
    mole_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
) -> NDArray[np.float64]:
    """Calculate mass fractions from mole fractions.

    Parameters
    ----------
    mole_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mole fractions for one composition or multiple composition rows.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mole_fractions``.

    Returns
    -------
    NDArray[np.float64]
        Mass fractions normalized along the last axis.
    """
    # SECTION: Convert inputs to numeric arrays
    x = _as_float_array(mole_fractions, "mole_fractions")
    mw = _as_float_array(molecular_weights, "molecular_weights")

    # SECTION: Validate inputs
    x, mw = _validate_same_array_shape(
        x,
        mw,
        "mole_fractions",
        "molecular_weights",
    )
    if np.any((x < 0.0) | (x > 1.0)):
        raise ValueError("mole_fractions must be between zero and one.")
    if np.any(mw <= 0.0):
        raise ValueError("molecular_weights must be positive.")

    # SECTION: Calculate mass fractions
    weighted_values = x * mw
    denominator = (
        weighted_values
        if weighted_values.ndim == 0
        else np.sum(weighted_values, axis=-1, keepdims=True)
    )
    if np.any(denominator <= 0.0):
        raise ValueError("The weighted molecular-weight sum must be positive.")

    return cast(
        NDArray[np.float64],
        np.asarray(weighted_values / denominator, dtype=np.float64),
    )


def _calc_mass_fraction_to_mole_fraction(
    mass_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
) -> NDArray[np.float64]:
    """Calculate mole fractions from mass fractions.

    Parameters
    ----------
    mass_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mass fractions for one composition or multiple composition rows.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mass_fractions``.

    Returns
    -------
    NDArray[np.float64]
        Mole fractions normalized along the last axis.
    """
    # SECTION: Convert inputs to numeric arrays
    w = _as_float_array(mass_fractions, "mass_fractions")
    mw = _as_float_array(molecular_weights, "molecular_weights")

    # SECTION: Validate inputs
    w, mw = _validate_same_array_shape(
        w,
        mw,
        "mass_fractions",
        "molecular_weights",
    )
    if np.any((w < 0.0) | (w > 1.0)):
        raise ValueError("mass_fractions must be between zero and one.")
    if np.any(mw <= 0.0):
        raise ValueError("molecular_weights must be positive.")

    # SECTION: Calculate mole fractions
    mole_values = w / mw
    denominator = (
        mole_values
        if mole_values.ndim == 0
        else np.sum(mole_values, axis=-1, keepdims=True)
    )
    if np.any(denominator <= 0.0):
        raise ValueError(
            "The reciprocal molecular-weight sum must be positive.")

    return cast(
        NDArray[np.float64],
        np.asarray(mole_values / denominator, dtype=np.float64),
    )


def _calc_molarities_to_molalities(
    molarities: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    solution_density: float | int | Sequence[float | int] | NDArray[np.number],
) -> NDArray[np.float64]:
    """Calculate molalities from molarities for one or more solutes.

    Parameters
    ----------
    molarities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molarities for one composition or multiple composition rows.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``molarities``.
    solution_density : float | int | Sequence[float | int] | NDArray[np.number]
        Solution density as a scalar or broadcast-compatible array.

    Returns
    -------
    NDArray[np.float64]
        Molalities with the same solute shape as ``molarities``.
    """
    # SECTION: Convert inputs to numeric arrays
    c = _as_float_array(molarities, "molarities")
    mw = _as_float_array(molecular_weights, "molecular_weights")
    rho = _as_float_array(solution_density, "solution_density")

    # SECTION: Validate inputs
    c, mw = _validate_same_array_shape(
        c,
        mw,
        "molarities",
        "molecular_weights",
    )
    if np.any(c < 0.0):
        raise ValueError("molarities must be non-negative.")
    if np.any(mw <= 0.0):
        raise ValueError("molecular_weights must be positive.")
    if np.any(rho <= 0.0):
        raise ValueError("solution_density must be positive.")

    # SECTION: Calculate molalities
    dissolved_mass = c * mw
    solute_mass = (
        dissolved_mass
        if dissolved_mass.ndim == 0
        else np.sum(dissolved_mass, axis=-1, keepdims=True)
    )
    if dissolved_mass.ndim == 2 and rho.ndim == 1:
        if rho.shape[0] != dissolved_mass.shape[0]:
            raise ValueError(
                "For 2-D molarities, 1-D solution_density must have one "
                "entry per state."
            )
        rho = rho[:, None]

    try:
        rho, solute_mass = np.broadcast_arrays(rho, solute_mass)
    except ValueError as exc:
        raise ValueError(
            "solution_density must be broadcast-compatible with the "
            "summed solute mass."
        ) from exc

    solvent_mass = rho - solute_mass
    if np.any(solvent_mass <= 0.0):
        raise ValueError(
            "solution_density - sum(molarity*molecular_weight) must be positive."
        )

    return cast(
        NDArray[np.float64],
        np.asarray(c / solvent_mass, dtype=np.float64),
    )


def _calc_molality_to_mole_fraction(
    molalities: float | int | Sequence[float | int] | NDArray[np.number],
    solvent_molecular_weight: float | int,
) -> NDArray[np.float64]:
    """Calculate solute and solvent mole fractions from molalities.

    Parameters
    ----------
    molalities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molalities for one composition or multiple composition rows.
    solvent_molecular_weight : float | int
        Molecular weight of the solvent.

    Returns
    -------
    NDArray[np.float64]
        Mole fractions with the solvent appended as the final component.
    """
    # SECTION: Convert inputs to numeric arrays
    b = _as_float_array(molalities, "molalities")

    # SECTION: Validate inputs
    if np.any(b < 0.0):
        raise ValueError("molalities must be non-negative.")
    if (
        not np.isfinite(solvent_molecular_weight)
        or solvent_molecular_weight <= 0.0
    ):
        raise ValueError("solvent_molecular_weight must be positive.")

    # SECTION: Calculate mole fractions
    solvent_moles = 1.0 / solvent_molecular_weight
    total = (
        solvent_moles + b
        if b.ndim == 0
        else solvent_moles + np.sum(b, axis=-1, keepdims=True)
    )
    solute_fractions = b / total

    if solute_fractions.ndim == 0:
        return cast(
            NDArray[np.float64],
            np.asarray([float(solute_fractions),
                       solvent_moles / float(total)]),
        )

    solvent_fraction = np.full(
        (*solute_fractions.shape[:-1], 1),
        solvent_moles,
        dtype=np.float64,
    ) / total

    return cast(
        NDArray[np.float64],
        np.concatenate((solute_fractions, solvent_fraction), axis=-1),
    )


def _calc_molarity_to_molality(
    molarity: float | int,
    molecular_weight: float | int,
    solution_density: float | int,
) -> float:
    """Calculate single-solute molality from numeric inputs.

    Parameters
    ----------
    molarity : float | int
        Solute molarity.
    molecular_weight : float | int
        Solute molecular weight.
    solution_density : float | int
        Solution density.

    Returns
    -------
    float
        Solute molality.
    """
    if molarity <= 0.0:
        raise ValueError("molarity must be positive.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")
    if solution_density <= 0.0:
        raise ValueError("solution_density must be positive.")

    solvent_mass = solution_density - molarity * molecular_weight
    if solvent_mass <= 0.0:
        raise ValueError(
            "solution_density - molarity*molecular_weight must be positive."
        )
    return molarity / solvent_mass


def _calc_molality_to_molarity(
    molality: float | int,
    molecular_weight: float | int,
    solution_density: float | int,
) -> float:
    """Calculate single-solute molarity from numeric inputs.

    Parameters
    ----------
    molality : float | int
        Solute molality.
    molecular_weight : float | int
        Solute molecular weight.
    solution_density : float | int
        Solution density.

    Returns
    -------
    float
        Solute molarity.
    """
    if molality <= 0.0:
        raise ValueError("molality must be positive.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")
    if solution_density <= 0.0:
        raise ValueError("solution_density must be positive.")

    return molality * solution_density / (1.0 + molality * molecular_weight)


def _calc_mole_fraction_to_molality(
    solute_mole_fraction: float | int,
    solvent_mole_fraction: float | int,
    solvent_molecular_weight: float | int,
) -> float:
    """Calculate molality from numeric mole fractions.

    Parameters
    ----------
    solute_mole_fraction : float | int
        Mole fraction of the solute.
    solvent_mole_fraction : float | int
        Mole fraction of the solvent.
    solvent_molecular_weight : float | int
        Molecular weight of the solvent.

    Returns
    -------
    float
        Solute molality.
    """
    if solute_mole_fraction < 0.0 or solute_mole_fraction > 1.0:
        raise ValueError("solute_mole_fraction must be between zero and one.")
    if solvent_mole_fraction <= 0.0 or solvent_mole_fraction > 1.0:
        raise ValueError(
            "solvent_mole_fraction must be greater than zero and no greater than one."
        )
    if solvent_molecular_weight <= 0.0:
        raise ValueError("solvent_molecular_weight must be positive.")

    return solute_mole_fraction / (
        solvent_mole_fraction * solvent_molecular_weight
    )


def _calc_molarity_to_mass_fraction(
    molarity: float | int,
    molecular_weight: float | int,
    solution_density: float | int,
) -> float:
    """Calculate mass fraction from molarity.

    Parameters
    ----------
    molarity : float | int
        Solute molarity.
    molecular_weight : float | int
        Solute molecular weight.
    solution_density : float | int
        Solution density.

    Returns
    -------
    float
        Solute mass fraction.
    """
    if molarity <= 0.0:
        raise ValueError("molarity must be positive.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")
    if solution_density <= 0.0:
        raise ValueError("solution_density must be positive.")

    mass_fraction = molarity * molecular_weight / solution_density
    if mass_fraction > 1.0:
        raise ValueError("Calculated mass fraction is greater than one.")
    return mass_fraction


def _calc_mass_fraction_to_molarity(
    mass_fraction: float | int,
    solution_density: float | int,
    molecular_weight: float | int,
) -> float:
    """Calculate molarity from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Solute mass fraction.
    solution_density : float | int
        Solution density.
    molecular_weight : float | int
        Solute molecular weight.

    Returns
    -------
    float
        Solute molarity.
    """
    if mass_fraction < 0.0 or mass_fraction > 1.0:
        raise ValueError("mass_fraction must be between zero and one.")
    if solution_density <= 0.0:
        raise ValueError("solution_density must be positive.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")

    return mass_fraction * solution_density / molecular_weight


def _calc_molality_to_mass_fraction(
    molality: float | int,
    molecular_weight: float | int,
) -> float:
    """Calculate mass fraction from molality.

    Parameters
    ----------
    molality : float | int
        Solute molality.
    molecular_weight : float | int
        Solute molecular weight.

    Returns
    -------
    float
        Solute mass fraction.
    """
    if molality <= 0.0:
        raise ValueError("molality must be positive.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")

    solute_mass = molality * molecular_weight
    return solute_mass / (1.0 + solute_mass)


def _calc_mass_fraction_to_molality(
    mass_fraction: float | int,
    molecular_weight: float | int,
) -> float:
    """Calculate molality from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Solute mass fraction.
    molecular_weight : float | int
        Solute molecular weight.

    Returns
    -------
    float
        Solute molality.
    """
    if mass_fraction < 0.0 or mass_fraction >= 1.0:
        raise ValueError(
            "mass_fraction must be at least zero and less than one.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")

    return mass_fraction / (molecular_weight * (1.0 - mass_fraction))


def _calc_molarity_to_mass_concentration(
    molarity: float | int,
    molecular_weight: float | int,
) -> float:
    """Calculate mass concentration from molarity.

    Parameters
    ----------
    molarity : float | int
        Component molarity.
    molecular_weight : float | int
        Component molecular weight.

    Returns
    -------
    float
        Component mass concentration.
    """
    if molarity <= 0.0:
        raise ValueError("molarity must be positive.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")

    return molarity * molecular_weight


def _calc_mass_concentration_to_molarity(
    mass_concentration: float | int,
    molecular_weight: float | int,
) -> float:
    """Calculate molarity from mass concentration.

    Parameters
    ----------
    mass_concentration : float | int
        Component mass concentration.
    molecular_weight : float | int
        Component molecular weight.

    Returns
    -------
    float
        Component molarity.
    """
    if mass_concentration <= 0.0:
        raise ValueError("mass_concentration must be positive.")
    if molecular_weight <= 0.0:
        raise ValueError("molecular_weight must be positive.")

    return mass_concentration / molecular_weight


def _calc_mass_fraction_to_weight_percent(mass_fraction: float | int) -> float:
    """Calculate weight percent from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Mass fraction on the interval [0, 1].

    Returns
    -------
    float
        Weight percent on the interval [0, 100].
    """
    if mass_fraction < 0.0 or mass_fraction > 1.0:
        raise ValueError("mass_fraction must be between zero and one.")
    return 100.0 * mass_fraction


def _calc_weight_percent_to_mass_fraction(weight_percent: float | int) -> float:
    """Calculate mass fraction from weight percent.

    Parameters
    ----------
    weight_percent : float | int
        Weight percent on the interval [0, 100].

    Returns
    -------
    float
        Mass fraction on the interval [0, 1].
    """
    if weight_percent < 0.0 or weight_percent > 100.0:
        raise ValueError("weight_percent must be between zero and 100.")
    return weight_percent / 100.0


def _calc_mole_fraction_to_mole_percent(mole_fraction: float | int) -> float:
    """Calculate mole percent from mole fraction.

    Parameters
    ----------
    mole_fraction : float | int
        Mole fraction on the interval [0, 1].

    Returns
    -------
    float
        Mole percent on the interval [0, 100].
    """
    if mole_fraction < 0.0 or mole_fraction > 1.0:
        raise ValueError("mole_fraction must be between zero and one.")
    return 100.0 * mole_fraction


def _calc_mole_percent_to_mole_fraction(mole_percent: float | int) -> float:
    """Calculate mole fraction from mole percent.

    Parameters
    ----------
    mole_percent : float | int
        Mole percent on the interval [0, 100].

    Returns
    -------
    float
        Mole fraction on the interval [0, 1].
    """
    if mole_percent < 0.0 or mole_percent > 100.0:
        raise ValueError("mole_percent must be between zero and 100.")
    return mole_percent / 100.0


def _calc_mass_fraction_to_ppm(mass_fraction: float | int) -> float:
    """Calculate mass-based parts per million from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Mass fraction on the interval [0, 1].

    Returns
    -------
    float
        Mass-based parts per million.
    """
    return _calc_mass_fraction_to_weight_percent(mass_fraction) * 10000.0


def _calc_ppm_mass_to_mass_fraction(ppm: float | int) -> float:
    """Calculate mass fraction from mass-based parts per million.

    Parameters
    ----------
    ppm : float | int
        Mass-based parts per million.

    Returns
    -------
    float
        Mass fraction.
    """
    if ppm < 0.0:
        raise ValueError("ppm must be non-negative.")
    return ppm * 1e-6


def _calc_mole_fraction_to_ppm(mole_fraction: float | int) -> float:
    """Calculate mole-based parts per million from mole fraction.

    Parameters
    ----------
    mole_fraction : float | int
        Mole fraction on the interval [0, 1].

    Returns
    -------
    float
        Mole-based parts per million.
    """
    return _calc_mole_fraction_to_mole_percent(mole_fraction) * 10000.0


def _calc_ppm_mole_to_mole_fraction(ppm: float | int) -> float:
    """Calculate mole fraction from mole-based parts per million.

    Parameters
    ----------
    ppm : float | int
        Mole-based parts per million.

    Returns
    -------
    float
        Mole fraction.
    """
    if ppm < 0.0:
        raise ValueError("ppm must be non-negative.")
    return ppm * 1e-6


def _calc_mass_fraction_to_ppb(mass_fraction: float | int) -> float:
    """Calculate mass-based parts per billion from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Mass fraction on the interval [0, 1].

    Returns
    -------
    float
        Mass-based parts per billion.
    """
    return _calc_mass_fraction_to_weight_percent(mass_fraction) * 10000000.0


def _calc_ppb_mass_to_mass_fraction(ppb: float | int) -> float:
    """Calculate mass fraction from mass-based parts per billion.

    Parameters
    ----------
    ppb : float | int
        Mass-based parts per billion.

    Returns
    -------
    float
        Mass fraction.
    """
    if ppb < 0.0:
        raise ValueError("ppb must be non-negative.")
    return ppb * 1e-9


def _calc_mole_fraction_to_ppb(mole_fraction: float | int) -> float:
    """Calculate mole-based parts per billion from mole fraction.

    Parameters
    ----------
    mole_fraction : float | int
        Mole fraction on the interval [0, 1].

    Returns
    -------
    float
        Mole-based parts per billion.
    """
    return _calc_mole_fraction_to_mole_percent(mole_fraction) * 10000000.0


def _calc_ppb_mole_to_mole_fraction(ppb: float | int) -> float:
    """Calculate mole fraction from mole-based parts per billion.

    Parameters
    ----------
    ppb : float | int
        Mole-based parts per billion.

    Returns
    -------
    float
        Mole fraction.
    """
    if ppb < 0.0:
        raise ValueError("ppb must be non-negative.")
    return ppb * 1e-9


# SECTION: Mole fraction and mass fraction conversions

# ! ::: Mole fraction to mass fraction from sequence


def _calc_mole_fraction_to_mass_fraction_from_sequence(
    mole_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
) -> list[float]:
    """Convert array-like mole fractions to a mass-fraction list.

    Parameters
    ----------
    mole_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mole fractions for one composition or multiple composition rows.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mole_fractions``.

    Returns
    -------
    list[float]
        Mass fractions converted from the NumPy calculation result.
    """
    return _calc_mole_fraction_to_mass_fraction(
        mole_fractions=mole_fractions,
        molecular_weights=molecular_weights,
    ).tolist()


# ! ::: Mole fraction to mass fraction from mapping


def _calc_mole_fraction_to_mass_fraction_from_mapping(
    mole_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping mole fractions to mass fractions.

    Parameters
    ----------
    mole_fractions : Mapping[str, float | int]
        Mole fractions keyed by component.
    molecular_weights : Mapping[str, float | int]
        Molecular weights keyed by component.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Mass fractions keyed by component.
    """
    x = _configure_component_values(
        dict(mole_fractions),
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "mole_fractions",
    )
    mw = _configure_component_values(
        dict(molecular_weights),
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molecular_weights",
    )
    _validate_same_keys(x, mw)

    mass_fractions = _calc_mole_fraction_to_mass_fraction(
        list(x.values()),
        [mw[key] for key in x],
    )
    return dict(zip(x.keys(), mass_fractions.tolist()))


# ! ::: Mole fraction to mass fraction from props


def _calc_mole_fraction_to_mass_fraction_from_props(
    mole_fractions: Mapping[str, CustomProp],
    molecular_weights: Mapping[str, CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping mole fractions to mass fractions using unit-aware weights.

    Parameters
    ----------
    mole_fractions : Mapping[str, CustomProp]
        Mole fractions keyed by component.
    molecular_weights : Mapping[str, CustomProp]
        Unit-aware molecular weights keyed by component.
    output_molecular_weight_unit : str | None, optional
        Unit used to normalize molecular weights.
    unit_conversion_fn : UnitConversionFn | None, optional
        Function used to convert ``CustomProp`` values.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Mass fractions keyed by component.
    """
    # SECTION: Normalize unit-aware molecular weights
    n = to_dict(
        mole_fractions,
        None,
        None
    )
    mw = to_dict(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )
    return _calc_mole_fraction_to_mass_fraction_from_mapping(
        mole_fractions=n,
        molecular_weights=mw,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )


# ! ::: Mass fraction to mole fraction from sequence


def _calc_mass_fraction_to_mole_fraction_from_sequence(
    mass_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
) -> list[float]:
    """Convert array-like mass fractions to a mole-fraction list.

    Parameters
    ----------
    mass_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mass fractions for one composition or multiple composition rows.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mass_fractions``.

    Returns
    -------
    list[float]
        Mole fractions converted from the NumPy calculation result.
    """
    return _calc_mass_fraction_to_mole_fraction(
        mass_fractions=mass_fractions,
        molecular_weights=molecular_weights,
    ).tolist()


# ! ::: Mass fraction to mole fraction from mapping


def _calc_mass_fraction_to_mole_fraction_from_mapping(
    mass_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping mass fractions to mole fractions.

    Parameters
    ----------
    mass_fractions : Mapping[str, float | int]
        Mass fractions keyed by component.
    molecular_weights : Mapping[str, float | int]
        Molecular weights keyed by component.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Mole fractions keyed by component.
    """
    w = _configure_component_values(
        dict(mass_fractions),
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "mass_fractions",
    )
    mw = _configure_component_values(
        dict(molecular_weights),
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molecular_weights",
    )
    _validate_same_keys(w, mw)

    mole_fractions = _calc_mass_fraction_to_mole_fraction(
        list(w.values()),
        [mw[key] for key in w],
    )
    return dict(zip(w.keys(), mole_fractions.tolist()))


# ! ::: Mass fraction to mole fraction from props


def _calc_mass_fraction_to_mole_fraction_from_props(
    mass_fractions: Mapping[str, CustomProp],
    molecular_weights: Mapping[str, CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping mass fractions to mole fractions using unit-aware weights.

    Parameters
    ----------
    mass_fractions : Mapping[str, CustomProp]
        Mass fractions keyed by component.
    molecular_weights : Mapping[str, CustomProp]
        Unit-aware molecular weights keyed by component.
    output_molecular_weight_unit : str | None, optional
        Unit used to normalize molecular weights.
    unit_conversion_fn : UnitConversionFn | None, optional
        Function used to convert ``CustomProp`` values.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Mole fractions keyed by component.
    """
    # SECTION: Normalize unit-aware molecular weights
    m = to_dict(
        mass_fractions,
        None,
        None,
    )
    mw = to_dict(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )
    return _calc_mass_fraction_to_mole_fraction_from_mapping(
        mass_fractions=m,
        molecular_weights=mw,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )


# SECTION: Molarity and molality conversions


def _calc_molarities_to_molalities_from_sequence(
    molarities: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    solution_density: float | int | Sequence[float | int] | NDArray[np.number],
) -> list[float]:
    """Convert array-like molarities to a molality list.

    Parameters
    ----------
    molarities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molarities for one composition or multiple composition rows.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``molarities``.
    solution_density : float | int | Sequence[float | int] | NDArray[np.number]
        Solution density as a scalar or broadcast-compatible array.

    Returns
    -------
    list[float]
        Molalities converted from the NumPy calculation result.
    """
    return _calc_molarities_to_molalities(
        molarities=molarities,
        molecular_weights=molecular_weights,
        solution_density=solution_density,
    ).tolist()


def _calc_molarities_to_molalities_from_mapping(
    molarities: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    solution_density: float | int,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping molarities to molalities.

    Parameters
    ----------
    molarities : Mapping[str, float | int]
        Solute molarities keyed by component.
    molecular_weights : Mapping[str, float | int]
        Molecular weights keyed by component.
    solution_density : float | int
        Solution density.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Molalities keyed by component.
    """
    c = _configure_component_values(
        dict(molarities),
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molarities",
    )
    mw = _configure_component_values(
        dict(molecular_weights),
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molecular_weights",
    )
    _validate_same_keys(c, mw)

    molalities = _calc_molarities_to_molalities(
        list(c.values()),
        [mw[key] for key in c],
        solution_density,
    )
    return dict(zip(c.keys(), molalities.tolist()))


def _calc_molarities_to_molalities_from_props(
    molarities: Mapping[str, CustomProp],
    molecular_weights: Mapping[str, CustomProp],
    solution_density: CustomProp,
    output_molarity_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    output_solution_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping molarities to molalities using unit-aware inputs.

    Parameters
    ----------
    molarities : Mapping[str, CustomProp]
        Unit-aware molarities keyed by component.
    molecular_weights : Mapping[str, CustomProp]
        Unit-aware molecular weights keyed by component.
    solution_density : CustomProp
        Unit-aware solution density.
    output_molarity_unit : str | None, optional
        Unit used to normalize molarities.
    output_molecular_weight_unit : str | None, optional
        Unit used to normalize molecular weights.
    output_solution_density_unit : str | None, optional
        Unit used to normalize solution density.
    unit_conversion_fn : UnitConversionFn | None, optional
        Function used to convert ``CustomProp`` values.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Molalities keyed by component.
    """
    c = to_dict(
        molarities,
        output_molarity_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )
    mw = to_dict(
        molecular_weights,
        output_molecular_weight_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )
    rho = _pos(
        solution_density,
        "solution_density",
        output_solution_density_unit,
        unit_conversion_fn,
    )
    return _calc_molarities_to_molalities_from_mapping(
        molarities=c,
        molecular_weights=mw,
        solution_density=rho,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )


# SECTION: Molality and mole fraction conversions

# ! ::: Molality to mole fraction from sequence


def _calc_molality_to_mole_fraction_from_sequence(
    molalities: float | int | Sequence[float | int] | NDArray[np.number],
    solvent_molecular_weight: float | int,
) -> list[float]:
    """Convert array-like molalities to a mole-fraction list.

    Parameters
    ----------
    molalities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molalities for one composition or multiple composition rows.
    solvent_molecular_weight : float | int
        Molecular weight of the solvent.

    Returns
    -------
    list[float]
        Solute mole fractions with solvent mole fraction appended last.
    """
    return _calc_molality_to_mole_fraction(
        molalities=molalities,
        solvent_molecular_weight=solvent_molecular_weight,
    ).tolist()


# ! ::: Molality to mole fraction from mapping


def _calc_molality_to_mole_fraction_from_mapping(
    molalities: Mapping[str, float | int],
    solvent_molecular_weight: float | int,
    solvent_key: str = "solvent",
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping molalities to mole fractions.

    Parameters
    ----------
    molalities : Mapping[str, float | int]
        Solute molalities keyed by component.
    solvent_molecular_weight : float | int
        Molecular weight of the solvent.
    solvent_key : str, optional
        Key used for the solvent mole fraction in the returned mapping.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Solute and solvent mole fractions keyed by component.
    """
    b = _configure_component_values(
        dict(molalities),
        components,
        component_key,
        case_sensitive,
        sort_by_components_order,
        "molalities",
    )
    mole_fractions = _calc_molality_to_mole_fraction(
        list(b.values()),
        solvent_molecular_weight,
    )
    res = dict(zip(b.keys(), mole_fractions[:-1].tolist()))
    res[solvent_key] = float(mole_fractions[-1])
    return res


# ! ::: Molality to mole fraction from props


def _calc_molality_to_mole_fraction_from_props(
    molalities: Mapping[str, CustomProp],
    solvent_molecular_weight: CustomProp,
    solvent_key: str = "solvent",
    output_molality_unit: str | None = None,
    output_solvent_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[List[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
) -> dict[str, float]:
    """Convert mapping molalities to mole fractions using unit-aware inputs.

    Parameters
    ----------
    molalities : Mapping[str, CustomProp]
        Unit-aware solute molalities keyed by component.
    solvent_molecular_weight : CustomProp
        Unit-aware molecular weight of the solvent.
    solvent_key : str, optional
        Key used for the solvent mole fraction in the returned mapping.
    output_molality_unit : str | None, optional
        Unit used to normalize molalities.
    output_solvent_molecular_weight_unit : str | None, optional
        Unit used to normalize solvent molecular weight.
    unit_conversion_fn : UnitConversionFn | None, optional
        Function used to convert ``CustomProp`` values.
    components : Optional[List[Component]], optional
        Component ordering and lookup metadata.
    component_key : Optional[ComponentKey], optional
        Component attribute used for lookup when ``components`` is provided.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether output follows the supplied component order.

    Returns
    -------
    dict[str, float]
        Solute and solvent mole fractions keyed by component.
    """
    b = to_dict(
        molalities,
        output_molality_unit,
        unit_conversion_fn=_resolve_unit_conversion_fn(unit_conversion_fn),
    )
    mw = _pos(
        solvent_molecular_weight,
        "solvent_molecular_weight",
        output_solvent_molecular_weight_unit,
        unit_conversion_fn,
    )
    return _calc_molality_to_mole_fraction_from_mapping(
        molalities=b,
        solvent_molecular_weight=mw,
        solvent_key=solvent_key,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )


# SECTION: Public exports
__all__ = [
    # ! mole fraction to mass fraction
    "_calc_mole_fraction_to_mass_fraction",
    "_calc_mole_fraction_to_mass_fraction_from_sequence",
    "_calc_mole_fraction_to_mass_fraction_from_mapping",
    "_calc_mole_fraction_to_mass_fraction_from_props",
    # ! mass fraction to mole fraction
    "_calc_mass_fraction_to_mole_fraction",
    "_calc_mass_fraction_to_mole_fraction_from_sequence",
    "_calc_mass_fraction_to_mole_fraction_from_mapping",
    "_calc_mass_fraction_to_mole_fraction_from_props",
    # ! molarities to molalities
    "_calc_molarities_to_molalities",
    "_calc_molarities_to_molalities_from_sequence",
    "_calc_molarities_to_molalities_from_mapping",
    "_calc_molarities_to_molalities_from_props",
    # ! molality conversions
    "_calc_molality_to_mole_fraction",
    "_calc_molality_to_mole_fraction_from_sequence",
    "_calc_molality_to_mole_fraction_from_mapping",
    "_calc_molality_to_mole_fraction_from_props",
    "_calc_molarity_to_molality",
    "_calc_molality_to_molarity",
    "_calc_molality_to_mass_fraction",
    "_calc_molarity_to_mass_concentration",
    "_calc_mole_fraction_to_molality",
    "_calc_mass_fraction_to_molality",
    # ! molarity conversions
    "_calc_molarity_to_mass_fraction",
    "_calc_mass_fraction_to_molarity",
    "_calc_mass_concentration_to_molarity",
    # ! general conversions
    "_calc_mass_fraction_to_weight_percent",
    "_calc_weight_percent_to_mass_fraction",
    "_calc_mole_fraction_to_mole_percent",
    "_calc_mole_percent_to_mole_fraction",
    "_calc_mass_fraction_to_ppm",
    "_calc_ppm_mass_to_mass_fraction",
    "_calc_mole_fraction_to_ppm",
    "_calc_ppm_mole_to_mole_fraction",
    "_calc_mass_fraction_to_ppb",
    "_calc_ppb_mass_to_mass_fraction",
    "_calc_mole_fraction_to_ppb",
    "_calc_ppb_mole_to_mole_fraction",
]
