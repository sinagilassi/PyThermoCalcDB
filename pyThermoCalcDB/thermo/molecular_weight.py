# import libs
import logging
from typing import Optional
from pyreactlab_core.core import parse_elemental_composition, parse_ionic_charge
from pythermodb_settings.models import CustomProp
# ! locals
from ..configs.atomic_weights import ATOMIC_WEIGHTS, ELECTRON_MOLAR_MASS

# NOTE: logger
logger = logging.getLogger(__name__)

# ! ::: Calculate Molecular Weight


def calc_molecular_weight(
    formula: str,
    include_electron_mass: bool = False,
    decimal_digits: Optional[int] = None
) -> Optional[CustomProp]:
    """
    Calculate the molecular weight of a chemical species in `g/mol`.

    The elemental composition is obtained using
    ``parse_elemental_composition`` from PyReactLab-Core.

    For ionic species, electron-mass correction can optionally
    be included.

    Parameters
    ----------
    formula : str
        Chemical formula.

        Examples
        --------
        H2O
        CO2
        Fe(OH)3
        Ca3(PO4)2
        CuSO4*5H2O
        Fe{3+}
        SO4{2-}
        e{-}

    include_electron_mass : bool, default=False
        If True, correct the molecular weight according to the
        ionic charge.

        For charge ``z``:

            MW_ion = MW_neutral - z * M_e

        where ``M_e`` is the electron molar mass.

        Therefore:

        - cations lose electron mass,
        - anions gain electron mass,
        - neutral species are unchanged.
    decimal_digits : int | None, optional
        Number of decimal digits to round the result to. If None, the result is not rounded.

    Returns
    -------
    CustomProp | None
        Molecular weight as a CustomProp instance, or None if calculation fails.

    Raises
    ------
    ValueError
        If an element is not available in ``ATOMIC_WEIGHTS``.

    Examples
    --------
    >>> calculate_molecular_weight("H2O")
    18.015

    >>> calculate_molecular_weight("SO4{2-}")
    96.056

    >>> calculate_molecular_weight(
    ...     "SO4{2-}",
    ...     include_electron_mass=True,
    ... )
    96.05709715981813

    >>> calculate_molecular_weight(
    ...     "e{-}",
    ...     include_electron_mass=True,
    ... )
    0.000548579909065
    """
    try:
        # SECTION: elemental composition
        composition = parse_elemental_composition(formula)
        # >> check
        if composition is None or not isinstance(composition, dict):
            logger.warning(
                f"Failed to parse elemental composition for formula '{formula}'."
            )
            return None

        molecular_weight = 0.0

        # SECTION: atomic contribution
        for element, count in composition.items():

            if element not in ATOMIC_WEIGHTS:
                logger.warning(
                    f"Atomic weight is not available for element '{element}'."
                )
                return None

            molecular_weight += (
                ATOMIC_WEIGHTS[element] * count
            )

        # SECTION: electron-mass correction
        if include_electron_mass:
            charge = parse_ionic_charge(formula)
            # >> check
            if charge is None:
                logger.warning(
                    f"Failed to parse ionic charge for formula '{formula}'."
                )
                return None

            molecular_weight -= (
                charge * ELECTRON_MOLAR_MASS
            )

        # NOTE: round the result to the specified number of decimal digits
        if decimal_digits is not None:
            molecular_weight = round(molecular_weight, decimal_digits)

        # res
        res = CustomProp(
            value=molecular_weight,
            unit="g/mol",
        )

        return res
    except Exception as e:
        logger.warning(
            f"An error occurred while calculating molecular weight for formula '{formula}': {e}"
        )
        return None


# SECTION: alias
calc_MW = calc_molecular_weight
