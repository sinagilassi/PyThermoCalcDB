"""Annotated public APIs for composition-basis conversions."""

# import libs
from collections.abc import Mapping, Sequence
from typing import Optional

import numpy as np
from numpy.typing import NDArray

# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp
from pythermodb_settings.models.units import UnitConversionFn

# locals
from ..utils.tools import to_annotated_value
from .core.conversions import (
    _calc_mass_concentration_to_molarity,
    _calc_mass_fraction_to_molality,
    _calc_mass_fraction_to_molarity,
    _calc_mass_fraction_to_mole_fraction,
    _calc_mass_fraction_to_mole_fraction_from_mapping,
    _calc_mass_fraction_to_mole_fraction_from_props,
    _calc_mass_fraction_to_mole_fraction_from_sequence,
    _calc_mass_fraction_to_ppb,
    _calc_mass_fraction_to_ppm,
    _calc_mass_fraction_to_weight_percent,
    _calc_molarities_to_molalities,
    _calc_molarities_to_molalities_from_mapping,
    _calc_molarities_to_molalities_from_props,
    _calc_molarities_to_molalities_from_sequence,
    _calc_molarity_to_mass_concentration,
    _calc_molarity_to_mass_fraction,
    _calc_molarity_to_molality,
    _calc_molality_to_mass_fraction,
    _calc_molality_to_molarity,
    _calc_molality_to_mole_fraction,
    _calc_molality_to_mole_fraction_from_mapping,
    _calc_molality_to_mole_fraction_from_props,
    _calc_molality_to_mole_fraction_from_sequence,
    _calc_mole_fraction_to_mass_fraction,
    _calc_mole_fraction_to_mass_fraction_from_mapping,
    _calc_mole_fraction_to_mass_fraction_from_props,
    _calc_mole_fraction_to_mass_fraction_from_sequence,
    _calc_mole_fraction_to_molality,
    _calc_mole_fraction_to_mole_percent,
    _calc_mole_fraction_to_ppb,
    _calc_mole_fraction_to_ppm,
    _calc_mole_percent_to_mole_fraction,
    _calc_ppb_mass_to_mass_fraction,
    _calc_ppb_mole_to_mole_fraction,
    _calc_ppm_mass_to_mass_fraction,
    _calc_ppm_mole_to_mole_fraction,
    _calc_weight_percent_to_mass_fraction,
)


def _annotate(
    value: object,
    *,
    name: str,
    description: str,
    unit: str | None,
    symbol: str | None,
    implementation: str,
) -> AnnotatedValue[object]:
    """Build an annotated conversion result.

    Parameters
    ----------
    value : object
        Calculated conversion result.
    name : str
        Result name.
    description : str
        Result description.
    unit : str | None
        Result unit metadata.
    symbol : str | None
        Result symbol metadata.
    implementation : str
        Core implementation name stored in result metadata.

    Returns
    -------
    AnnotatedValue[object]
        Annotated calculation result.
    """
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation=implementation,
    )


# ======================================================================
# *** Public annotated API
# ======================================================================


@calculation_info(
    name="mole_fraction_to_mass_fraction",
    description="Convert mole fractions to mass fractions.",
    equation="w_i = x_i*M_i / sum_j(x_j*M_j)",
    inputs={
        "mole_fractions": "Mole fractions.",
        "molecular_weights": "Molecular weights.",
    },
    outputs={"mass_fractions": "Mass fractions."},
    tags=("conversion", "mole_fraction", "mass_fraction", "array_like", "numpy"),
)
def calc_mole_fraction_to_mass_fraction(
    mole_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "mass_fractions",
    description: str = "Mass fractions converted from mole fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.float64]]:
    """Return annotated mass fractions from numeric array-like inputs.

    Parameters
    ----------
    mole_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mole fractions to convert.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mole_fractions``.

    Returns
    -------
    AnnotatedValue[NDArray[np.float64]]
        Annotated mass fractions.
    """
    value = _calc_mole_fraction_to_mass_fraction(mole_fractions, molecular_weights)
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mole_fraction_to_mass_fraction",
    )


@calculation_info(
    name="mole_fraction_to_mass_fraction",
    description="Convert sequence mole fractions to mass fractions.",
    equation="w_i = x_i*M_i / sum_j(x_j*M_j)",
    inputs={
        "mole_fractions": "Sequence or NumPy array of mole fractions.",
        "molecular_weights": "Sequence or NumPy array of molecular weights.",
    },
    outputs={"mass_fractions": "List of mass fractions."},
    tags=("conversion", "mole_fraction", "mass_fraction", "sequence", "numeric"),
)
def calc_mole_fraction_to_mass_fraction_from_sequence(
    mole_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "mass_fractions",
    description: str = "Mass fractions converted from mole fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Return annotated mass fractions from array-like inputs as a list.

    Parameters
    ----------
    mole_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mole fractions to convert.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mole_fractions``.

    Returns
    -------
    AnnotatedValue[list[float]]
        Annotated mass-fraction list.
    """
    value = _calc_mole_fraction_to_mass_fraction_from_sequence(
        mole_fractions=mole_fractions,
        molecular_weights=molecular_weights,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mole_fraction_to_mass_fraction_from_sequence",
    )


@calculation_info(
    name="mole_fraction_to_mass_fraction",
    description="Convert mapping mole fractions to mass fractions.",
    equation="w_i = x_i*M_i / sum_j(x_j*M_j)",
    inputs={
        "mole_fractions": "Mapping of component keys to mole fractions.",
        "molecular_weights": "Mapping of component keys to molecular weights.",
    },
    outputs={"mass_fractions": "Mapping of component keys to mass fractions."},
    tags=("conversion", "mole_fraction", "mass_fraction", "mapping", "numeric"),
)
def calc_mole_fraction_to_mass_fraction_from_mapping(
    mole_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "mass_fractions",
    description: str = "Mass fractions converted from keyed mole fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated mass fractions from numeric mappings.

    Parameters
    ----------
    mole_fractions : Mapping[str, float | int]
        Mole fractions keyed by component.
    molecular_weights : Mapping[str, float | int]
        Molecular weights keyed by component.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated mass fractions keyed by component.
    """
    value = _calc_mole_fraction_to_mass_fraction_from_mapping(
        mole_fractions=mole_fractions,
        molecular_weights=molecular_weights,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mole_fraction_to_mass_fraction_from_mapping",
    )


@calculation_info(
    name="mole_fraction_to_mass_fraction",
    description="Convert mole fractions to mass fractions using unit-aware weights.",
    equation="w_i = x_i*M_i / sum_j(x_j*M_j)",
    inputs={
        "mole_fractions": "Mapping of component keys to mole fractions.",
        "molecular_weights": "Mapping of component keys to unit-aware molecular weights.",
    },
    outputs={"mass_fractions": "Mapping of component keys to mass fractions."},
    tags=("conversion", "mole_fraction", "mass_fraction", "mapping", "unit_aware"),
)
def calc_mole_fraction_to_mass_fraction_from_props(
    mole_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "mass_fractions",
    description: str = "Mass fractions converted from keyed mole fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated mass fractions using unit-aware molecular weights.

    Parameters
    ----------
    mole_fractions : Mapping[str, float | int]
        Mole fractions keyed by component.
    molecular_weights : Mapping[str, CustomProp]
        Unit-aware molecular weights keyed by component.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated mass fractions keyed by component.
    """
    value = _calc_mole_fraction_to_mass_fraction_from_props(
        mole_fractions=mole_fractions,
        molecular_weights=molecular_weights,
        output_molecular_weight_unit=output_molecular_weight_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mole_fraction_to_mass_fraction_from_props",
    )


@calculation_info(
    name="mass_fraction_to_mole_fraction",
    description="Convert mass fractions to mole fractions.",
    equation="x_i = (w_i/M_i) / sum_j(w_j/M_j)",
    inputs={
        "mass_fractions": "Mass fractions.",
        "molecular_weights": "Molecular weights.",
    },
    outputs={"mole_fractions": "Mole fractions."},
    tags=("conversion", "mass_fraction", "mole_fraction", "array_like", "numpy"),
)
def calc_mass_fraction_to_mole_fraction(
    mass_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from mass fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.float64]]:
    """Return annotated mole fractions from numeric array-like inputs.

    Parameters
    ----------
    mass_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mass fractions to convert.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mass_fractions``.

    Returns
    -------
    AnnotatedValue[NDArray[np.float64]]
        Annotated mole fractions.
    """
    value = _calc_mass_fraction_to_mole_fraction(mass_fractions, molecular_weights)
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mass_fraction_to_mole_fraction",
    )


@calculation_info(
    name="mass_fraction_to_mole_fraction",
    description="Convert sequence mass fractions to mole fractions.",
    equation="x_i = (w_i/M_i) / sum_j(w_j/M_j)",
    inputs={
        "mass_fractions": "Sequence or NumPy array of mass fractions.",
        "molecular_weights": "Sequence or NumPy array of molecular weights.",
    },
    outputs={"mole_fractions": "List of mole fractions."},
    tags=("conversion", "mass_fraction", "mole_fraction", "sequence", "numeric"),
)
def calc_mass_fraction_to_mole_fraction_from_sequence(
    mass_fractions: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from mass fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Return annotated mole fractions from array-like inputs as a list.

    Parameters
    ----------
    mass_fractions : float | int | Sequence[float | int] | NDArray[np.number]
        Mass fractions to convert.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``mass_fractions``.

    Returns
    -------
    AnnotatedValue[list[float]]
        Annotated mole-fraction list.
    """
    value = _calc_mass_fraction_to_mole_fraction_from_sequence(
        mass_fractions=mass_fractions,
        molecular_weights=molecular_weights,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mass_fraction_to_mole_fraction_from_sequence",
    )


@calculation_info(
    name="mass_fraction_to_mole_fraction",
    description="Convert mapping mass fractions to mole fractions.",
    equation="x_i = (w_i/M_i) / sum_j(w_j/M_j)",
    inputs={
        "mass_fractions": "Mapping of component keys to mass fractions.",
        "molecular_weights": "Mapping of component keys to molecular weights.",
    },
    outputs={"mole_fractions": "Mapping of component keys to mole fractions."},
    tags=("conversion", "mass_fraction", "mole_fraction", "mapping", "numeric"),
)
def calc_mass_fraction_to_mole_fraction_from_mapping(
    mass_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from keyed mass fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated mole fractions from numeric mappings.

    Parameters
    ----------
    mass_fractions : Mapping[str, float | int]
        Mass fractions keyed by component.
    molecular_weights : Mapping[str, float | int]
        Molecular weights keyed by component.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated mole fractions keyed by component.
    """
    value = _calc_mass_fraction_to_mole_fraction_from_mapping(
        mass_fractions=mass_fractions,
        molecular_weights=molecular_weights,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mass_fraction_to_mole_fraction_from_mapping",
    )


@calculation_info(
    name="mass_fraction_to_mole_fraction",
    description="Convert mass fractions to mole fractions using unit-aware weights.",
    equation="x_i = (w_i/M_i) / sum_j(w_j/M_j)",
    inputs={
        "mass_fractions": "Mapping of component keys to mass fractions.",
        "molecular_weights": "Mapping of component keys to unit-aware molecular weights.",
    },
    outputs={"mole_fractions": "Mapping of component keys to mole fractions."},
    tags=("conversion", "mass_fraction", "mole_fraction", "mapping", "unit_aware"),
)
def calc_mass_fraction_to_mole_fraction_from_props(
    mass_fractions: Mapping[str, float | int],
    molecular_weights: Mapping[str, CustomProp],
    output_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from keyed mass fractions.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated mole fractions using unit-aware molecular weights.

    Parameters
    ----------
    mass_fractions : Mapping[str, float | int]
        Mass fractions keyed by component.
    molecular_weights : Mapping[str, CustomProp]
        Unit-aware molecular weights keyed by component.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated mole fractions keyed by component.
    """
    value = _calc_mass_fraction_to_mole_fraction_from_props(
        mass_fractions=mass_fractions,
        molecular_weights=molecular_weights,
        output_molecular_weight_unit=output_molecular_weight_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_mass_fraction_to_mole_fraction_from_props",
    )


@calculation_info(
    name="molarities_to_molalities",
    description="Convert multisolute molarities to molalities.",
    equation="b_i = C_i / (rho - sum_j(C_j*M_j))",
    inputs={
        "molarities": "Solute molarities.",
        "molecular_weights": "Solute molecular weights.",
        "solution_density": "Solution density.",
    },
    outputs={"molalities": "Solute molalities."},
    tags=("conversion", "molarity", "molality", "array_like", "numpy"),
)
def calc_molarities_to_molalities(
    molarities: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    solution_density: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "molalities",
    description: str = "Molalities converted from molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.float64]]:
    """Return annotated molalities from numeric array-like inputs.

    Parameters
    ----------
    molarities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molarities to convert.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``molarities``.
    solution_density : float | int | Sequence[float | int] | NDArray[np.number]
        Solution density.

    Returns
    -------
    AnnotatedValue[NDArray[np.float64]]
        Annotated molalities.
    """
    value = _calc_molarities_to_molalities(
        molarities=molarities,
        molecular_weights=molecular_weights,
        solution_density=solution_density,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molarities_to_molalities",
    )


@calculation_info(
    name="molarities_to_molalities",
    description="Convert sequence molarities to molalities.",
    equation="b_i = C_i / (rho - sum_j(C_j*M_j))",
    inputs={
        "molarities": "Sequence or NumPy array of solute molarities.",
        "molecular_weights": "Sequence or NumPy array of molecular weights.",
        "solution_density": "Solution density.",
    },
    outputs={"molalities": "List of solute molalities."},
    tags=("conversion", "molarity", "molality", "sequence", "numeric"),
)
def calc_molarities_to_molalities_from_sequence(
    molarities: float | int | Sequence[float | int] | NDArray[np.number],
    molecular_weights: float | int | Sequence[float | int] | NDArray[np.number],
    solution_density: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "molalities",
    description: str = "Molalities converted from molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Return annotated molalities from array-like inputs as a list.

    Parameters
    ----------
    molarities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molarities to convert.
    molecular_weights : float | int | Sequence[float | int] | NDArray[np.number]
        Molecular weights paired with ``molarities``.
    solution_density : float | int | Sequence[float | int] | NDArray[np.number]
        Solution density.

    Returns
    -------
    AnnotatedValue[list[float]]
        Annotated molality list.
    """
    value = _calc_molarities_to_molalities_from_sequence(
        molarities=molarities,
        molecular_weights=molecular_weights,
        solution_density=solution_density,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molarities_to_molalities_from_sequence",
    )


@calculation_info(
    name="molarities_to_molalities",
    description="Convert mapping molarities to molalities.",
    equation="b_i = C_i / (rho - sum_j(C_j*M_j))",
    inputs={
        "molarities": "Mapping of component keys to molarities.",
        "molecular_weights": "Mapping of component keys to molecular weights.",
        "solution_density": "Solution density.",
    },
    outputs={"molalities": "Mapping of component keys to molalities."},
    tags=("conversion", "molarity", "molality", "mapping", "numeric"),
)
def calc_molarities_to_molalities_from_mapping(
    molarities: Mapping[str, float | int],
    molecular_weights: Mapping[str, float | int],
    solution_density: float | int,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "molalities",
    description: str = "Molalities converted from keyed molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated molalities from numeric mappings.

    Parameters
    ----------
    molarities : Mapping[str, float | int]
        Solute molarities keyed by component.
    molecular_weights : Mapping[str, float | int]
        Molecular weights keyed by component.
    solution_density : float | int
        Solution density.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated molalities keyed by component.
    """
    value = _calc_molarities_to_molalities_from_mapping(
        molarities=molarities,
        molecular_weights=molecular_weights,
        solution_density=solution_density,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molarities_to_molalities_from_mapping",
    )


@calculation_info(
    name="molarities_to_molalities",
    description="Convert unit-aware mapping molarities to molalities.",
    equation="b_i = C_i / (rho - sum_j(C_j*M_j))",
    inputs={
        "molarities": "Mapping of component keys to unit-aware molarities.",
        "molecular_weights": "Mapping of component keys to unit-aware molecular weights.",
        "solution_density": "Unit-aware solution density.",
    },
    outputs={"molalities": "Mapping of component keys to molalities."},
    tags=("conversion", "molarity", "molality", "mapping", "unit_aware"),
)
def calc_molarities_to_molalities_from_props(
    molarities: Mapping[str, CustomProp],
    molecular_weights: Mapping[str, CustomProp],
    solution_density: CustomProp,
    output_molarity_unit: str | None = None,
    output_molecular_weight_unit: str | None = None,
    output_solution_density_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "molalities",
    description: str = "Molalities converted from keyed unit-aware molarities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated molalities from unit-aware mappings.

    Parameters
    ----------
    molarities : Mapping[str, CustomProp]
        Unit-aware molarities keyed by component.
    molecular_weights : Mapping[str, CustomProp]
        Unit-aware molecular weights keyed by component.
    solution_density : CustomProp
        Unit-aware solution density.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated molalities keyed by component.
    """
    value = _calc_molarities_to_molalities_from_props(
        molarities=molarities,
        molecular_weights=molecular_weights,
        solution_density=solution_density,
        output_molarity_unit=output_molarity_unit,
        output_molecular_weight_unit=output_molecular_weight_unit,
        output_solution_density_unit=output_solution_density_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molarities_to_molalities_from_props",
    )


@calculation_info(
    name="molality_to_mole_fraction",
    description="Convert molalities to solute and solvent mole fractions.",
    equation="x_i = b_i / (sum_j(b_j) + 1/M_s)",
    inputs={
        "molalities": "Solute molalities.",
        "solvent_molecular_weight": "Solvent molecular weight.",
    },
    outputs={"mole_fractions": "Solute mole fractions with solvent included."},
    tags=("conversion", "molality", "mole_fraction", "array_like", "numpy"),
)
def calc_molality_to_mole_fraction(
    molalities: float | int | Sequence[float | int] | NDArray[np.number],
    solvent_molecular_weight: float | int,
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[NDArray[np.float64]]:
    """Return annotated mole fractions from numeric molality inputs.

    Parameters
    ----------
    molalities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molalities to convert.
    solvent_molecular_weight : float | int
        Solvent molecular weight.

    Returns
    -------
    AnnotatedValue[NDArray[np.float64]]
        Annotated mole fractions with solvent appended last.
    """
    value = _calc_molality_to_mole_fraction(
        molalities=molalities,
        solvent_molecular_weight=solvent_molecular_weight,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molality_to_mole_fraction",
    )


@calculation_info(
    name="molality_to_mole_fraction",
    description="Convert sequence molalities to mole fractions.",
    equation="x_i = b_i / (sum_j(b_j) + 1/M_s)",
    inputs={
        "molalities": "Sequence or NumPy array of solute molalities.",
        "solvent_molecular_weight": "Solvent molecular weight.",
    },
    outputs={"mole_fractions": "List of solute mole fractions with solvent last."},
    tags=("conversion", "molality", "mole_fraction", "sequence", "numeric"),
)
def calc_molality_to_mole_fraction_from_sequence(
    molalities: float | int | Sequence[float | int] | NDArray[np.number],
    solvent_molecular_weight: float | int,
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[list[float]]:
    """Return annotated mole fractions from array-like molalities as a list.

    Parameters
    ----------
    molalities : float | int | Sequence[float | int] | NDArray[np.number]
        Solute molalities to convert.
    solvent_molecular_weight : float | int
        Solvent molecular weight.

    Returns
    -------
    AnnotatedValue[list[float]]
        Annotated mole-fraction list with solvent appended last.
    """
    value = _calc_molality_to_mole_fraction_from_sequence(
        molalities=molalities,
        solvent_molecular_weight=solvent_molecular_weight,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molality_to_mole_fraction_from_sequence",
    )


@calculation_info(
    name="molality_to_mole_fraction",
    description="Convert mapping molalities to mole fractions.",
    equation="x_i = b_i / (sum_j(b_j) + 1/M_s)",
    inputs={
        "molalities": "Mapping of component keys to molalities.",
        "solvent_molecular_weight": "Solvent molecular weight.",
    },
    outputs={"mole_fractions": "Mapping of component keys to mole fractions."},
    tags=("conversion", "molality", "mole_fraction", "mapping", "numeric"),
)
def calc_molality_to_mole_fraction_from_mapping(
    molalities: Mapping[str, float | int],
    solvent_molecular_weight: float | int,
    solvent_key: str = "solvent",
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from keyed molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated mole fractions from numeric molality mappings.

    Parameters
    ----------
    molalities : Mapping[str, float | int]
        Solute molalities keyed by component.
    solvent_molecular_weight : float | int
        Solvent molecular weight.
    solvent_key : str, optional
        Output key used for the solvent.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated mole fractions keyed by component.
    """
    value = _calc_molality_to_mole_fraction_from_mapping(
        molalities=molalities,
        solvent_molecular_weight=solvent_molecular_weight,
        solvent_key=solvent_key,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molality_to_mole_fraction_from_mapping",
    )


@calculation_info(
    name="molality_to_mole_fraction",
    description="Convert unit-aware mapping molalities to mole fractions.",
    equation="x_i = b_i / (sum_j(b_j) + 1/M_s)",
    inputs={
        "molalities": "Mapping of component keys to unit-aware molalities.",
        "solvent_molecular_weight": "Unit-aware solvent molecular weight.",
    },
    outputs={"mole_fractions": "Mapping of component keys to mole fractions."},
    tags=("conversion", "molality", "mole_fraction", "mapping", "unit_aware"),
)
def calc_molality_to_mole_fraction_from_props(
    molalities: Mapping[str, CustomProp],
    solvent_molecular_weight: CustomProp,
    solvent_key: str = "solvent",
    output_molality_unit: str | None = None,
    output_solvent_molecular_weight_unit: str | None = None,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[list[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "mole_fractions",
    description: str = "Mole fractions converted from keyed unit-aware molalities.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[dict[str, float]]:
    """Return annotated mole fractions from unit-aware molality mappings.

    Parameters
    ----------
    molalities : Mapping[str, CustomProp]
        Unit-aware solute molalities keyed by component.
    solvent_molecular_weight : CustomProp
        Unit-aware solvent molecular weight.
    solvent_key : str, optional
        Output key used for the solvent.

    Returns
    -------
    AnnotatedValue[dict[str, float]]
        Annotated mole fractions keyed by component.
    """
    value = _calc_molality_to_mole_fraction_from_props(
        molalities=molalities,
        solvent_molecular_weight=solvent_molecular_weight,
        solvent_key=solvent_key,
        output_molality_unit=output_molality_unit,
        output_solvent_molecular_weight_unit=output_solvent_molecular_weight_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )
    return _annotate(
        value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_molality_to_mole_fraction_from_props",
    )


@calculation_info(
    name="molarity_to_molality",
    description="Convert single-solute molarity to molality.",
    equation="b = C / (rho - C*M)",
    inputs={"molarity": "Solute molarity.", "molecular_weight": "Molecular weight.", "solution_density": "Solution density."},
    outputs={"molality": "Solute molality."},
    tags=("conversion", "molarity", "molality", "scalar", "numeric"),
)
def calc_molarity_to_molality(
    molarity: float | int,
    molecular_weight: float | int,
    solution_density: float | int,
    *,
    name: str = "molality",
    description: str = "Molality converted from molarity.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality from numeric molarity.

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
    AnnotatedValue[float]
        Annotated solute molality.
    """
    value = _calc_molarity_to_molality(molarity, molecular_weight, solution_density)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_molarity_to_molality")


@calculation_info(
    name="molality_to_molarity",
    description="Convert single-solute molality to molarity.",
    equation="C = b*rho / (1 + b*M)",
    inputs={"molality": "Solute molality.", "molecular_weight": "Molecular weight.", "solution_density": "Solution density."},
    outputs={"molarity": "Solute molarity."},
    tags=("conversion", "molality", "molarity", "scalar", "numeric"),
)
def calc_molality_to_molarity(
    molality: float | int,
    molecular_weight: float | int,
    solution_density: float | int,
    *,
    name: str = "molarity",
    description: str = "Molarity converted from molality.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity from numeric molality.

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
    AnnotatedValue[float]
        Annotated solute molarity.
    """
    value = _calc_molality_to_molarity(molality, molecular_weight, solution_density)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_molality_to_molarity")


@calculation_info(
    name="mole_fraction_to_molality",
    description="Convert solute mole fraction to molality.",
    equation="b_i = x_i / (x_s*M_s)",
    inputs={"solute_mole_fraction": "Solute mole fraction.", "solvent_mole_fraction": "Solvent mole fraction.", "solvent_molecular_weight": "Solvent molecular weight."},
    outputs={"molality": "Solute molality."},
    tags=("conversion", "mole_fraction", "molality", "scalar", "numeric"),
)
def calc_mole_fraction_to_molality(
    solute_mole_fraction: float | int,
    solvent_mole_fraction: float | int,
    solvent_molecular_weight: float | int,
    *,
    name: str = "molality",
    description: str = "Molality converted from mole fraction.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality from numeric mole fractions.

    Parameters
    ----------
    solute_mole_fraction : float | int
        Solute mole fraction.
    solvent_mole_fraction : float | int
        Solvent mole fraction.
    solvent_molecular_weight : float | int
        Solvent molecular weight.

    Returns
    -------
    AnnotatedValue[float]
        Annotated solute molality.
    """
    value = _calc_mole_fraction_to_molality(
        solute_mole_fraction,
        solvent_mole_fraction,
        solvent_molecular_weight,
    )
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mole_fraction_to_molality")


@calculation_info(
    name="molarity_to_mass_fraction",
    description="Convert molarity to mass fraction.",
    equation="w_i = C_i*M_i / rho",
    inputs={"molarity": "Solute molarity.", "molecular_weight": "Molecular weight.", "solution_density": "Solution density."},
    outputs={"mass_fraction": "Solute mass fraction."},
    tags=("conversion", "molarity", "mass_fraction", "scalar", "numeric"),
)
def calc_molarity_to_mass_fraction(
    molarity: float | int,
    molecular_weight: float | int,
    solution_density: float | int,
    *,
    name: str = "mass_fraction",
    description: str = "Mass fraction converted from molarity.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass fraction from numeric molarity.

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
    AnnotatedValue[float]
        Annotated mass fraction.
    """
    value = _calc_molarity_to_mass_fraction(molarity, molecular_weight, solution_density)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_molarity_to_mass_fraction")


@calculation_info(
    name="mass_fraction_to_molarity",
    description="Convert mass fraction to molarity.",
    equation="C_i = w_i*rho / M_i",
    inputs={"mass_fraction": "Solute mass fraction.", "solution_density": "Solution density.", "molecular_weight": "Molecular weight."},
    outputs={"molarity": "Solute molarity."},
    tags=("conversion", "mass_fraction", "molarity", "scalar", "numeric"),
)
def calc_mass_fraction_to_molarity(
    mass_fraction: float | int,
    solution_density: float | int,
    molecular_weight: float | int,
    *,
    name: str = "molarity",
    description: str = "Molarity converted from mass fraction.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity from numeric mass fraction.

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
    AnnotatedValue[float]
        Annotated molarity.
    """
    value = _calc_mass_fraction_to_molarity(mass_fraction, solution_density, molecular_weight)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mass_fraction_to_molarity")


@calculation_info(
    name="molality_to_mass_fraction",
    description="Convert molality to mass fraction.",
    equation="w_i = b_i*M_i / (1 + b_i*M_i)",
    inputs={"molality": "Solute molality.", "molecular_weight": "Molecular weight."},
    outputs={"mass_fraction": "Solute mass fraction."},
    tags=("conversion", "molality", "mass_fraction", "scalar", "numeric"),
)
def calc_molality_to_mass_fraction(
    molality: float | int,
    molecular_weight: float | int,
    *,
    name: str = "mass_fraction",
    description: str = "Mass fraction converted from molality.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass fraction from numeric molality.

    Parameters
    ----------
    molality : float | int
        Solute molality.
    molecular_weight : float | int
        Solute molecular weight.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mass fraction.
    """
    value = _calc_molality_to_mass_fraction(molality, molecular_weight)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_molality_to_mass_fraction")


@calculation_info(
    name="mass_fraction_to_molality",
    description="Convert mass fraction to molality.",
    equation="b_i = w_i / (M_i*(1 - w_i))",
    inputs={"mass_fraction": "Solute mass fraction.", "molecular_weight": "Molecular weight."},
    outputs={"molality": "Solute molality."},
    tags=("conversion", "mass_fraction", "molality", "scalar", "numeric"),
)
def calc_mass_fraction_to_molality(
    mass_fraction: float | int,
    molecular_weight: float | int,
    *,
    name: str = "molality",
    description: str = "Molality converted from mass fraction.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molality from numeric mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Solute mass fraction.
    molecular_weight : float | int
        Solute molecular weight.

    Returns
    -------
    AnnotatedValue[float]
        Annotated molality.
    """
    value = _calc_mass_fraction_to_molality(mass_fraction, molecular_weight)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mass_fraction_to_molality")


@calculation_info(
    name="molarity_to_mass_concentration",
    description="Convert molarity to mass concentration.",
    equation="c_m,i = C_i*M_i",
    inputs={"molarity": "Component molarity.", "molecular_weight": "Molecular weight."},
    outputs={"mass_concentration": "Component mass concentration."},
    tags=("conversion", "molarity", "mass_concentration", "scalar", "numeric"),
)
def calc_molarity_to_mass_concentration(
    molarity: float | int,
    molecular_weight: float | int,
    *,
    name: str = "mass_concentration",
    description: str = "Mass concentration converted from molarity.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass concentration from numeric molarity.

    Parameters
    ----------
    molarity : float | int
        Component molarity.
    molecular_weight : float | int
        Component molecular weight.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mass concentration.
    """
    value = _calc_molarity_to_mass_concentration(molarity, molecular_weight)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_molarity_to_mass_concentration")


@calculation_info(
    name="mass_concentration_to_molarity",
    description="Convert mass concentration to molarity.",
    equation="C_i = c_m,i / M_i",
    inputs={"mass_concentration": "Component mass concentration.", "molecular_weight": "Molecular weight."},
    outputs={"molarity": "Component molarity."},
    tags=("conversion", "mass_concentration", "molarity", "scalar", "numeric"),
)
def calc_mass_concentration_to_molarity(
    mass_concentration: float | int,
    molecular_weight: float | int,
    *,
    name: str = "molarity",
    description: str = "Molarity converted from mass concentration.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated molarity from numeric mass concentration.

    Parameters
    ----------
    mass_concentration : float | int
        Component mass concentration.
    molecular_weight : float | int
        Component molecular weight.

    Returns
    -------
    AnnotatedValue[float]
        Annotated molarity.
    """
    value = _calc_mass_concentration_to_molarity(mass_concentration, molecular_weight)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mass_concentration_to_molarity")


@calculation_info(
    name="mass_fraction_to_weight_percent",
    description="Convert mass fraction to weight percent.",
    equation="wt_percent = 100*w",
    inputs={"mass_fraction": "Mass fraction."},
    outputs={"weight_percent": "Weight percent."},
    tags=("conversion", "mass_fraction", "weight_percent", "scalar", "numeric"),
)
def calc_mass_fraction_to_weight_percent(
    mass_fraction: float | int,
    *,
    name: str = "weight_percent",
    description: str = "Weight percent converted from mass fraction.",
    unit: str | None = "%",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated weight percent from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Mass fraction to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated weight percent.
    """
    value = _calc_mass_fraction_to_weight_percent(mass_fraction)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mass_fraction_to_weight_percent")


@calculation_info(
    name="weight_percent_to_mass_fraction",
    description="Convert weight percent to mass fraction.",
    equation="w = wt_percent / 100",
    inputs={"weight_percent": "Weight percent."},
    outputs={"mass_fraction": "Mass fraction."},
    tags=("conversion", "weight_percent", "mass_fraction", "scalar", "numeric"),
)
def calc_weight_percent_to_mass_fraction(
    weight_percent: float | int,
    *,
    name: str = "mass_fraction",
    description: str = "Mass fraction converted from weight percent.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass fraction from weight percent.

    Parameters
    ----------
    weight_percent : float | int
        Weight percent to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mass fraction.
    """
    value = _calc_weight_percent_to_mass_fraction(weight_percent)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_weight_percent_to_mass_fraction")


@calculation_info(
    name="mole_fraction_to_mole_percent",
    description="Convert mole fraction to mole percent.",
    equation="mol_percent = 100*x",
    inputs={"mole_fraction": "Mole fraction."},
    outputs={"mole_percent": "Mole percent."},
    tags=("conversion", "mole_fraction", "mole_percent", "scalar", "numeric"),
)
def calc_mole_fraction_to_mole_percent(
    mole_fraction: float | int,
    *,
    name: str = "mole_percent",
    description: str = "Mole percent converted from mole fraction.",
    unit: str | None = "%",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mole percent from mole fraction.

    Parameters
    ----------
    mole_fraction : float | int
        Mole fraction to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mole percent.
    """
    value = _calc_mole_fraction_to_mole_percent(mole_fraction)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mole_fraction_to_mole_percent")


@calculation_info(
    name="mole_percent_to_mole_fraction",
    description="Convert mole percent to mole fraction.",
    equation="x = mol_percent / 100",
    inputs={"mole_percent": "Mole percent."},
    outputs={"mole_fraction": "Mole fraction."},
    tags=("conversion", "mole_percent", "mole_fraction", "scalar", "numeric"),
)
def calc_mole_percent_to_mole_fraction(
    mole_percent: float | int,
    *,
    name: str = "mole_fraction",
    description: str = "Mole fraction converted from mole percent.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mole fraction from mole percent.

    Parameters
    ----------
    mole_percent : float | int
        Mole percent to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mole fraction.
    """
    value = _calc_mole_percent_to_mole_fraction(mole_percent)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mole_percent_to_mole_fraction")


@calculation_info(
    name="mass_fraction_to_ppm",
    description="Convert mass fraction to mass-based ppm.",
    equation="ppm = 1e6*w",
    inputs={"mass_fraction": "Mass fraction."},
    outputs={"ppm": "Mass-based parts per million."},
    tags=("conversion", "mass_fraction", "ppm", "scalar", "numeric"),
)
def calc_mass_fraction_to_ppm(
    mass_fraction: float | int,
    *,
    name: str = "ppm",
    description: str = "Mass-based ppm converted from mass fraction.",
    unit: str | None = "ppm",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass-based ppm from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Mass fraction to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mass-based ppm.
    """
    value = _calc_mass_fraction_to_ppm(mass_fraction)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mass_fraction_to_ppm")


@calculation_info(
    name="ppm_mass_to_mass_fraction",
    description="Convert mass-based ppm to mass fraction.",
    equation="w = ppm*1e-6",
    inputs={"ppm": "Mass-based parts per million."},
    outputs={"mass_fraction": "Mass fraction."},
    tags=("conversion", "ppm", "mass_fraction", "scalar", "numeric"),
)
def calc_ppm_mass_to_mass_fraction(
    ppm: float | int,
    *,
    name: str = "mass_fraction",
    description: str = "Mass fraction converted from mass-based ppm.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass fraction from mass-based ppm.

    Parameters
    ----------
    ppm : float | int
        Mass-based parts per million.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mass fraction.
    """
    value = _calc_ppm_mass_to_mass_fraction(ppm)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_ppm_mass_to_mass_fraction")


@calculation_info(
    name="mole_fraction_to_ppm",
    description="Convert mole fraction to mole-based ppm.",
    equation="ppm = 1e6*x",
    inputs={"mole_fraction": "Mole fraction."},
    outputs={"ppm": "Mole-based parts per million."},
    tags=("conversion", "mole_fraction", "ppm", "scalar", "numeric"),
)
def calc_mole_fraction_to_ppm(
    mole_fraction: float | int,
    *,
    name: str = "ppm",
    description: str = "Mole-based ppm converted from mole fraction.",
    unit: str | None = "ppm",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mole-based ppm from mole fraction.

    Parameters
    ----------
    mole_fraction : float | int
        Mole fraction to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mole-based ppm.
    """
    value = _calc_mole_fraction_to_ppm(mole_fraction)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mole_fraction_to_ppm")


@calculation_info(
    name="ppm_mole_to_mole_fraction",
    description="Convert mole-based ppm to mole fraction.",
    equation="x = ppm*1e-6",
    inputs={"ppm": "Mole-based parts per million."},
    outputs={"mole_fraction": "Mole fraction."},
    tags=("conversion", "ppm", "mole_fraction", "scalar", "numeric"),
)
def calc_ppm_mole_to_mole_fraction(
    ppm: float | int,
    *,
    name: str = "mole_fraction",
    description: str = "Mole fraction converted from mole-based ppm.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mole fraction from mole-based ppm.

    Parameters
    ----------
    ppm : float | int
        Mole-based parts per million.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mole fraction.
    """
    value = _calc_ppm_mole_to_mole_fraction(ppm)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_ppm_mole_to_mole_fraction")


@calculation_info(
    name="mass_fraction_to_ppb",
    description="Convert mass fraction to mass-based ppb.",
    equation="ppb = 1e9*w",
    inputs={"mass_fraction": "Mass fraction."},
    outputs={"ppb": "Mass-based parts per billion."},
    tags=("conversion", "mass_fraction", "ppb", "scalar", "numeric"),
)
def calc_mass_fraction_to_ppb(
    mass_fraction: float | int,
    *,
    name: str = "ppb",
    description: str = "Mass-based ppb converted from mass fraction.",
    unit: str | None = "ppb",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass-based ppb from mass fraction.

    Parameters
    ----------
    mass_fraction : float | int
        Mass fraction to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mass-based ppb.
    """
    value = _calc_mass_fraction_to_ppb(mass_fraction)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mass_fraction_to_ppb")


@calculation_info(
    name="ppb_mass_to_mass_fraction",
    description="Convert mass-based ppb to mass fraction.",
    equation="w = ppb*1e-9",
    inputs={"ppb": "Mass-based parts per billion."},
    outputs={"mass_fraction": "Mass fraction."},
    tags=("conversion", "ppb", "mass_fraction", "scalar", "numeric"),
)
def calc_ppb_mass_to_mass_fraction(
    ppb: float | int,
    *,
    name: str = "mass_fraction",
    description: str = "Mass fraction converted from mass-based ppb.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mass fraction from mass-based ppb.

    Parameters
    ----------
    ppb : float | int
        Mass-based parts per billion.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mass fraction.
    """
    value = _calc_ppb_mass_to_mass_fraction(ppb)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_ppb_mass_to_mass_fraction")


@calculation_info(
    name="mole_fraction_to_ppb",
    description="Convert mole fraction to mole-based ppb.",
    equation="ppb = 1e9*x",
    inputs={"mole_fraction": "Mole fraction."},
    outputs={"ppb": "Mole-based parts per billion."},
    tags=("conversion", "mole_fraction", "ppb", "scalar", "numeric"),
)
def calc_mole_fraction_to_ppb(
    mole_fraction: float | int,
    *,
    name: str = "ppb",
    description: str = "Mole-based ppb converted from mole fraction.",
    unit: str | None = "ppb",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mole-based ppb from mole fraction.

    Parameters
    ----------
    mole_fraction : float | int
        Mole fraction to convert.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mole-based ppb.
    """
    value = _calc_mole_fraction_to_ppb(mole_fraction)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_mole_fraction_to_ppb")


@calculation_info(
    name="ppb_mole_to_mole_fraction",
    description="Convert mole-based ppb to mole fraction.",
    equation="x = ppb*1e-9",
    inputs={"ppb": "Mole-based parts per billion."},
    outputs={"mole_fraction": "Mole fraction."},
    tags=("conversion", "ppb", "mole_fraction", "scalar", "numeric"),
)
def calc_ppb_mole_to_mole_fraction(
    ppb: float | int,
    *,
    name: str = "mole_fraction",
    description: str = "Mole fraction converted from mole-based ppb.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """Return annotated mole fraction from mole-based ppb.

    Parameters
    ----------
    ppb : float | int
        Mole-based parts per billion.

    Returns
    -------
    AnnotatedValue[float]
        Annotated mole fraction.
    """
    value = _calc_ppb_mole_to_mole_fraction(ppb)
    return _annotate(value, name=name, description=description, unit=unit, symbol=symbol, implementation="_calc_ppb_mole_to_mole_fraction")


# SECTION: Public exports
__all__ = [
    "calc_mole_fraction_to_mass_fraction",
    "calc_mole_fraction_to_mass_fraction_from_sequence",
    "calc_mole_fraction_to_mass_fraction_from_mapping",
    "calc_mole_fraction_to_mass_fraction_from_props",
    "calc_mass_fraction_to_mole_fraction",
    "calc_mass_fraction_to_mole_fraction_from_sequence",
    "calc_mass_fraction_to_mole_fraction_from_mapping",
    "calc_mass_fraction_to_mole_fraction_from_props",
    "calc_molarities_to_molalities",
    "calc_molarities_to_molalities_from_sequence",
    "calc_molarities_to_molalities_from_mapping",
    "calc_molarities_to_molalities_from_props",
    "calc_molality_to_mole_fraction",
    "calc_molality_to_mole_fraction_from_sequence",
    "calc_molality_to_mole_fraction_from_mapping",
    "calc_molality_to_mole_fraction_from_props",
    "calc_molarity_to_molality",
    "calc_molality_to_molarity",
    "calc_mole_fraction_to_molality",
    "calc_molarity_to_mass_fraction",
    "calc_mass_fraction_to_molarity",
    "calc_molality_to_mass_fraction",
    "calc_mass_fraction_to_molality",
    "calc_molarity_to_mass_concentration",
    "calc_mass_concentration_to_molarity",
    "calc_mass_fraction_to_weight_percent",
    "calc_weight_percent_to_mass_fraction",
    "calc_mole_fraction_to_mole_percent",
    "calc_mole_percent_to_mole_fraction",
    "calc_mass_fraction_to_ppm",
    "calc_ppm_mass_to_mass_fraction",
    "calc_mole_fraction_to_ppm",
    "calc_ppm_mole_to_mole_fraction",
    "calc_mass_fraction_to_ppb",
    "calc_ppb_mass_to_mass_fraction",
    "calc_mole_fraction_to_ppb",
    "calc_ppb_mole_to_mole_fraction",
]
