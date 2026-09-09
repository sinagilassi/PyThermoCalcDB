"""Charge-balance calculations."""

# import libs
import logging
from collections.abc import Mapping, Sequence
from typing import Optional

import numpy as np
from numpy.typing import NDArray

# >> pythermodb-settings
from pythermodb_settings.decorators import calculation_info
from pythermodb_settings.models import AnnotatedValue, Component, ComponentKey, CustomProp
from pythermodb_settings.models.units import UnitConversionFn

# locals
from ..utils.conversions import (
    _resolve_result_unit,
)
from ..utils.tools import to_annotated_value
from .core.charge_balance import (
    _calc_charge_balance,
    _calc_charge_balance_from_sequence,
    _calc_charge_balance_from_props,
    _calc_charge_balance_from_mapping,
    _calc_charge_balance_from_props,
    _check_electroneutrality,
    _check_electroneutrality_from_sequence,
    _check_electroneutrality_from_mapping,
    _check_electroneutrality_from_props
)

# NOTE: logger setup
logger = logging.getLogger(__name__)


# =====================================================================
# *** Public annotated API
# =====================================================================
# ::: annotated for numpy array


@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "array_like",
        "numpy",
        "numeric",
    ),
)
def calc_charge_balance(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float | NDArray[np.float64]]:
    """
    Return annotated charge-balance residual from numeric inputs.

    Parameters
    ----------
    concentrations : float | int | Sequence[float | int] | NDArray[np.number]
        Species concentrations on a common concentration basis.
    charges : float | int | Sequence[float | int] | NDArray[np.number]
        Species charges aligned with ``concentrations``.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    unit : str, optional
        Unit stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[float | NDArray[np.float64]]
        Annotated charge-balance residual.
    """
    # SECTION: Calculate value
    value = _calc_charge_balance(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_charge_balance",
    )


# ::: annotated for sequence

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from sequence inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "sequence",
        "numeric",
    ),
)
def calc_charge_balance_from_sequence(
    concentrations: Sequence[float | int],
    charges: Sequence[float | int],
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """
    Return annotated charge-balance residual from numeric sequences.

    Parameters
    ----------
    concentrations : Sequence[float | int]
        Species concentrations in sequence order.
    charges : Sequence[float | int]
        Species charges aligned with ``concentrations``.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    unit : str, optional
        Unit stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[float]
        Annotated charge-balance residual for the sequence.
    """
    # SECTION: Calculate value
    value = _calc_charge_balance_from_sequence(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_charge_balance_from_sequence",
    )


# ::: annotated for mapping

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from mapping inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Keyed species concentrations.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "mapping",
        "numeric",
    ),
)
def calc_charge_balance_from_mapping(
    concentrations: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    unit: str | None = None,
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """
    Return annotated charge-balance residual from numeric mappings.

    Parameters
    ----------
    concentrations : Mapping[str, float | int]
        Species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Species charges keyed by the same identifiers as ``concentrations``.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    unit : str, optional
        Unit stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[float]
        Annotated charge-balance residual for the keyed species set.
    """
    # SECTION: Calculate value
    value = _calc_charge_balance_from_mapping(
        concentrations=concentrations,
        charges=charges,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_charge_balance_from_mapping",
    )


# ::: annotated for CustomProp mapping

@calculation_info(
    name="charge_balance",
    description="Calculate charge-balance residual from unit-aware mapping inputs.",
    equation="sum_i(z_i * c_i)",
    inputs={
        "concentrations": "Keyed unit-aware species concentrations.",
        "charges": "Keyed species charges.",
    },
    outputs={
        "charge_balance": "Charge-balance residual."
    },
    tags=(
        "charge_balance",
        "mapping",
        "unit_aware",
    ),
)
def calc_charge_balance_from_props(
    concentrations: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_concentration_unit: str,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "charge_balance",
    description: str = "Charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[float]:
    """
    Return annotated charge-balance residual from unit-aware mappings.

    Parameters
    ----------
    concentrations : Mapping[str, CustomProp]
        Unit-aware species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Numeric species charges keyed by component identifier.
    output_concentration_unit : str
        Unit used to normalize concentration values before calculation.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to the configured conversion helper.
    components : Sequence[Component], optional
        Component metadata used for key remapping and ordering.
    component_key : ComponentKey, optional
        Component identifier format used for mapping keys.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether values should follow the order of ``components``.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[float]
        Annotated charge-balance residual on the normalized concentration basis.
    """
    # SECTION: Calculate value
    value = _calc_charge_balance_from_props(
        concentrations=concentrations,
        charges=charges,
        output_concentration_unit=output_concentration_unit,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )

    # SECTION: Build annotated value
    return to_annotated_value(
        value=value,
        name=name,
        description=description,
        unit=_resolve_result_unit(
            "concentrations",
            concentrations,
            output_concentration_unit,
        ),
        symbol=symbol,
        implementation="_calc_charge_balance_from_props",
    )


# ! ::: Electroneutrality check

@calculation_info(
    name="electroneutrality",
    description="Check electroneutrality from numeric inputs.",
    equation="abs(sum_i(z_i * c_i)) <= tolerance",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
        "tolerance": "Maximum allowed absolute charge-balance residual.",
    },
    outputs={
        "electroneutrality": "Whether the charge-balance residual is within tolerance."
    },
    tags=(
        "electroneutrality",
        "array_like",
        "numpy",
        "numeric",
        "boolean",
    ),
)
def check_electroneutrality(
    concentrations: float | int | Sequence[float | int] | NDArray[np.number],
    charges: float | int | Sequence[float | int] | NDArray[np.number],
    tolerance: float = 1e-12,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool | NDArray[np.bool_]]:
    """
    Return annotated electroneutrality check from numeric inputs.

    Parameters
    ----------
    concentrations : float | int | Sequence[float | int] | NDArray[np.number]
        Species concentrations on a common concentration basis.
    charges : float | int | Sequence[float | int] | NDArray[np.number]
        Species charges aligned with ``concentrations``.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[bool | NDArray[np.bool_]]
        Annotated boolean result of the electroneutrality check.
    """
    # SECTION: Check electroneutrality
    result = _check_electroneutrality(
        concentrations=concentrations,
        charges=charges,
        tolerance=tolerance,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=result,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# ! ::: Sequence check

@calculation_info(
    name="electroneutrality",
    description="Check electroneutrality from sequence inputs.",
    equation="abs(sum_i(z_i * c_i)) <= tolerance",
    inputs={
        "concentrations": "Species concentrations.",
        "charges": "Species charges.",
        "tolerance": "Maximum allowed absolute charge-balance residual.",
    },
    outputs={
        "electroneutrality": "Whether the charge-balance residual is within tolerance."
    },
    tags=(
        "electroneutrality",
        "sequence",
        "numeric",
        "boolean",
    ),
)
def check_electroneutrality_from_sequence(
    concentrations: Sequence[float | int],
    charges: Sequence[float | int],
    tolerance: float = 1e-12,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool]:
    """
    Return annotated electroneutrality check from numeric sequences.

    Parameters
    ----------
    concentrations : Sequence[float | int]
        Species concentrations in sequence order.
    charges : Sequence[float | int]
        Species charges aligned with ``concentrations``.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[bool]
        Annotated boolean result for the sequence check.
    """
    # SECTION: Check electroneutrality
    result = _check_electroneutrality_from_sequence(
        concentrations=concentrations,
        charges=charges,
        tolerance=tolerance,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=result,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# ! ::: Mapping check

@calculation_info(
    name="electroneutrality",
    description="Check electroneutrality from mapping inputs.",
    equation="abs(sum_i(z_i * c_i)) <= tolerance",
    inputs={
        "concentrations": "Keyed species concentrations.",
        "charges": "Keyed species charges.",
        "tolerance": "Maximum allowed absolute charge-balance residual.",
    },
    outputs={
        "electroneutrality": "Whether the charge-balance residual is within tolerance."
    },
    tags=(
        "electroneutrality",
        "mapping",
        "numeric",
        "boolean",
    ),
)
def check_electroneutrality_from_mapping(
    concentrations: Mapping[str, float | int],
    charges: Mapping[str, float | int],
    tolerance: float = 1e-12,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool]:
    """
    Return annotated electroneutrality check from numeric mappings.

    Parameters
    ----------
    concentrations : Mapping[str, float | int]
        Species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Species charges keyed by the same identifiers as ``concentrations``.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[bool]
        Annotated boolean result for the keyed check.
    """
    # SECTION: Check electroneutrality
    result = _check_electroneutrality_from_mapping(
        concentrations=concentrations,
        charges=charges,
        tolerance=tolerance,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=result,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# ! ::: Unit-aware mapping check

@calculation_info(
    name="electroneutrality",
    description="Check electroneutrality from unit-aware mapping inputs.",
    equation="abs(sum_i(z_i * c_i)) <= tolerance",
    inputs={
        "concentrations": "Keyed unit-aware species concentrations.",
        "charges": "Keyed species charges.",
        "output_concentration_unit": "Unit used to normalize concentrations.",
        "tolerance": "Maximum allowed absolute charge-balance residual.",
    },
    outputs={
        "electroneutrality": "Whether the charge-balance residual is within tolerance."
    },
    tags=(
        "electroneutrality",
        "mapping",
        "unit_aware",
        "boolean",
    ),
)
def check_electroneutrality_from_props(
    concentrations: Mapping[str, CustomProp],
    charges: Mapping[str, float | int],
    output_concentration_unit: str,
    tolerance: float = 1e-12,
    unit_conversion_fn: UnitConversionFn | None = None,
    components: Optional[Sequence[Component]] = None,
    component_key: Optional[ComponentKey] = None,
    case_sensitive: bool = True,
    sort_by_components_order: bool = True,
    *,
    name: str = "electroneutrality",
    description: str = "Electroneutrality check based on the charge-balance residual.",
    symbol: str | None = None,
) -> AnnotatedValue[bool]:
    """
    Return annotated electroneutrality check from unit-aware mappings.

    Parameters
    ----------
    concentrations : Mapping[str, CustomProp]
        Unit-aware species concentrations keyed by component identifier.
    charges : Mapping[str, float | int]
        Numeric species charges keyed by component identifier.
    output_concentration_unit : str
        Unit used to normalize concentration values before checking.
    tolerance : float, optional
        Maximum absolute charge-balance residual allowed.
    unit_conversion_fn : UnitConversionFn, optional
        Unit conversion function. Defaults to the configured conversion helper.
    components : Sequence[Component], optional
        Component metadata used for key remapping and ordering.
    component_key : ComponentKey, optional
        Component identifier format used for mapping keys.
    case_sensitive : bool, optional
        Whether component matching is case-sensitive.
    sort_by_components_order : bool, optional
        Whether values should follow the order of ``components``.
    name : str, optional
        Name stored in the annotated result.
    description : str, optional
        Description stored in the annotated result.
    symbol : str, optional
        Symbol stored in the annotated result.

    Returns
    -------
    AnnotatedValue[bool]
        Annotated boolean result for the unit-aware check.
    """
    # SECTION: Check electroneutrality
    result = _check_electroneutrality_from_props(
        concentrations=concentrations,
        charges=charges,
        output_concentration_unit=output_concentration_unit,
        tolerance=tolerance,
        unit_conversion_fn=unit_conversion_fn,
        components=components,
        component_key=component_key,
        case_sensitive=case_sensitive,
        sort_by_components_order=sort_by_components_order,
    )

    # SECTION: Build annotated boolean value
    return to_annotated_value(
        value=result,
        name=name,
        description=description,
        unit=None,
        symbol=symbol,
    )


# =====================================================================
# *** Aliases
# =====================================================================

calc_mapping_charge_balance = calc_charge_balance_from_mapping
calc_sequence_charge_balance = calc_charge_balance_from_sequence


# SECTION: Public exports
__all__ = [
    # public
    "calc_charge_balance",
    "calc_charge_balance_from_sequence",
    "calc_charge_balance_from_mapping",
    "calc_charge_balance_from_props",
    "check_electroneutrality",
    "check_electroneutrality_from_sequence",
    "check_electroneutrality_from_mapping",
    "check_electroneutrality_from_props",

    # compatibility
    "calc_mapping_charge_balance",
    "calc_sequence_charge_balance",
]
