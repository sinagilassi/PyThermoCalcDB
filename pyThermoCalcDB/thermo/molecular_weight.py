# import libs
import logging
from typing import Optional
from pythermodb_settings.models import CustomProp
# ! locals
from .core.molecular_weight import _calc_molecular_weight

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
        # SECTION: calc
        molecular_weight = _calc_molecular_weight(
            formula=formula,
            include_electron_mass=include_electron_mass,
            decimal_digits=decimal_digits,
        )

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
