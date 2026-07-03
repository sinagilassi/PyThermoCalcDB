# import libs
import logging
import math
from typing import Optional, Sequence, Literal
from pythermodb_settings.models import (
    CustomProp,
)
# locals

# NOTE: set up logger
logger = logging.getLogger(__name__)


# SECTION: Gas mixture viscosity calculation
def gas_mixture_viscosity(
    mole_fractions: Sequence[float],
    viscosities: Sequence[CustomProp],
    molecular_weights: Sequence[CustomProp],
    mode: Literal["wilke", "linear"] = "wilke",
    normalize: bool = True,
) -> Optional[CustomProp]:
    """
    Calculate gas-mixture viscosity.

    Parameters
    ----------
    mole_fractions : sequence of float
        Mole fractions of components.
    viscosities : sequence of CustomProp
        Pure-component gas viscosities. All values must have the same unit.
    molecular_weights : sequence of CustomProp
        Molecular weights. All values must have the same unit.
    mode : {"wilke", "linear"}, optional
        "wilke"  : Wilke mixing rule for gas mixtures.
        "linear" : simple mole-fraction average.
    normalize : bool, optional
        If True, mole fractions are normalized to sum to 1.

    Returns
    -------
    Optional[CustomProp]
        Mixture viscosity using the same unit as the input viscosities, or None
        if an error occurs.

    Notes
    -----
    Wilke equation:

        mu_mix = sum_i y_i mu_i / sum_j y_j phi_ij

    where:

        phi_ij =
        [1 + (mu_i/mu_j)^0.5 * (MW_j/MW_i)^0.25]^2
        / [8 * (1 + MW_i/MW_j)]^0.5

    For i = j, phi_ij = 1.
    """
    try:
        # NOTE: collect inputs as lists for length checks and repeated access
        y = list(mole_fractions)
        viscosity_props = list(viscosities)
        molecular_weight_props = list(molecular_weights)

        n = len(y)

        # NOTE: validate list sizes and mole fractions
        if not (len(viscosity_props) == len(molecular_weight_props) == n):
            logger.warning(
                "mole_fractions, viscosities, and molecular_weights must have the same length."
            )
            return None

        if n == 0:
            logger.warning("Input lists cannot be empty.")
            return None

        if any(value < 0 for value in y):
            logger.warning("Mole fractions cannot be negative.")
            return None

        # NOTE: validate CustomProp units before using values
        viscosity_units = {prop.unit.strip() for prop in viscosity_props}
        if len(viscosity_units) != 1:
            logger.warning("All viscosities must have the same unit.")
            return None

        molecular_weight_units = {
            prop.unit.strip() for prop in molecular_weight_props}
        if len(molecular_weight_units) != 1:
            logger.warning("All molecular weights must have the same unit.")
            return None

        # NOTE: extract numeric values from CustomProp inputs
        viscosity_unit = viscosity_props[0].unit.strip()
        mu = [float(prop.value) for prop in viscosity_props]
        mw = [float(prop.value) for prop in molecular_weight_props]

        # NOTE: validate physical property values
        if any(value <= 0 for value in mu):
            logger.warning("Viscosities must be positive.")
            return None

        if any(value <= 0 for value in mw):
            logger.warning("Molecular weights must be positive.")
            return None

        y_sum = sum(y)

        if y_sum <= 0:
            logger.warning("Sum of mole fractions must be positive.")
            return None

        # NOTE: normalize or validate mole-fraction sum
        if normalize:
            y = [value / y_sum for value in y]
        else:
            if abs(y_sum - 1.0) > 1e-8:
                logger.warning(
                    "Mole fractions must sum to 1, or use normalize=True.")
                return None

        # NOTE: calculate mixture viscosity with the selected mixing rule
        if mode == "linear":
            mu_mix = sum(y_i * mu_i for y_i, mu_i in zip(y, mu))

        elif mode == "wilke":
            mu_mix = 0.0

            for i in range(n):
                denominator = 0.0

                for j in range(n):
                    if i == j:
                        phi_ij = 1.0
                    else:
                        phi_ij = (
                            (
                                1.0
                                + math.sqrt(mu[i] / mu[j])
                                * (mw[j] / mw[i]) ** 0.25
                            )
                            ** 2
                        ) / math.sqrt(8.0 * (1.0 + mw[i] / mw[j]))

                    denominator += y[j] * phi_ij

                mu_mix += y[i] * mu[i] / denominator

        else:
            logger.warning("mode must be either 'wilke' or 'linear'.")
            return None

        # NOTE: return result as CustomProp using the input viscosity unit
        return CustomProp(
            value=float(mu_mix),
            unit=viscosity_unit,
        )
    except Exception as e:
        logger.error(f"Error in gas mixture viscosity calculation: {e}")
        return None

# SECTION: Liquid mixture viscosity calculation


def liquid_mixture_viscosity(
    mole_fractions: Sequence[float],
    viscosities: Sequence[CustomProp],
    mode: Literal["log", "linear"] = "log",
    normalize: bool = True,
) -> Optional[CustomProp]:
    """
    Calculate ideal/simple liquid-mixture viscosity.

    Parameters
    ----------
    mole_fractions : sequence of float
        Liquid mole fractions of components.
    viscosities : sequence of CustomProp
        Pure-component liquid viscosities. All values must have the same unit.
    mode : {"log", "linear"}, optional
        "log":
            Logarithmic mixing rule:

                ln(mu_mix) = sum_i x_i ln(mu_i)

            Recommended for ideal or near-ideal liquid mixtures.

        "linear":
            Simple mole-fraction average:

                mu_mix = sum_i x_i mu_i

            Rough approximation.

    normalize : bool, optional
        If True, mole fractions are normalized to sum to 1.
        If False, mole fractions must already sum to 1.

    Returns
    -------
    Optional[CustomProp]
        Liquid-mixture viscosity using the same unit as the input viscosities,
        or None if an error occurs.
    """
    try:
        # NOTE: collect inputs as lists for length checks and repeated access
        x = list(mole_fractions)
        viscosity_props = list(viscosities)

        n = len(x)

        # NOTE: validate list sizes and mole fractions
        if len(viscosity_props) != n:
            logger.warning(
                "mole_fractions and viscosities must have the same length.")
            return None

        if n == 0:
            logger.warning("Input lists cannot be empty.")
            return None

        if any(value < 0 for value in x):
            logger.warning("Mole fractions cannot be negative.")
            return None

        # NOTE: validate CustomProp units before using values
        viscosity_units = {prop.unit.strip() for prop in viscosity_props}
        if len(viscosity_units) != 1:
            logger.warning("All viscosities must have the same unit.")
            return None

        # NOTE: extract numeric values from CustomProp inputs
        viscosity_unit = viscosity_props[0].unit.strip()
        mu = [float(prop.value) for prop in viscosity_props]

        # NOTE: validate physical property values
        if any(value <= 0 for value in mu):
            logger.warning("Viscosities must be positive.")
            return None

        x_sum = sum(x)

        if x_sum <= 0:
            logger.warning("Sum of mole fractions must be positive.")
            return None

        # NOTE: normalize or validate mole-fraction sum
        if normalize:
            x = [value / x_sum for value in x]
        else:
            if abs(x_sum - 1.0) > 1e-8:
                logger.warning(
                    "Mole fractions must sum to 1, or use normalize=True.")
                return None

        # NOTE: calculate mixture viscosity with the selected mixing rule
        if mode == "log":
            ln_mu_mix = sum(x_i * math.log(mu_i)
                            for x_i, mu_i in zip(x, mu))
            mu_mix = math.exp(ln_mu_mix)

        elif mode == "linear":
            mu_mix = sum(x_i * mu_i for x_i, mu_i in zip(x, mu))

        else:
            logger.warning("mode must be either 'log' or 'linear'.")
            return None

        # NOTE: return result as CustomProp using the input viscosity unit
        return CustomProp(
            value=float(mu_mix),
            unit=viscosity_unit,
        )
    except Exception as e:
        logger.error(f"Error in liquid mixture viscosity calculation: {e}")
        return None
