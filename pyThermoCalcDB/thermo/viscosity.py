# import libs
import logging
import math
from typing import Dict, Optional, Any, Sequence, Literal, cast
from pythermodb_settings.models import (
    Temperature,
    Pressure,
    CustomProp,
)
# locals

# NOTE: set up logger
logger = logging.getLogger(__name__)


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
    y : sequence of float
        Mole fractions of components.
    mu : sequence of float
        Pure-component gas viscosities [Pa.s].
    mw : sequence of float
        Molecular weights [g/mol] or [kg/kmol].
        Only ratios are used, so either unit is fine.
    mode : {"wilke", "linear"}, optional
        "wilke"  : Wilke mixing rule for gas mixtures.
        "linear" : simple mole-fraction average.
    normalize : bool, optional
        If True, mole fractions are normalized to sum to 1.

    Returns
    -------
    float
        Mixture viscosity [Pa.s].

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

    y = list(y)
    mu = list(mu)
    mw = list(mw)

    n = len(y)

    if not (len(mu) == len(mw) == n):
        raise ValueError("y, mu, and mw must have the same length.")

    if n == 0:
        raise ValueError("Input lists cannot be empty.")

    if any(value < 0 for value in y):
        raise ValueError("Mole fractions cannot be negative.")

    if any(value <= 0 for value in mu):
        raise ValueError("Viscosities must be positive.")

    if any(value <= 0 for value in mw):
        raise ValueError("Molecular weights must be positive.")

    y_sum = sum(y)

    if y_sum <= 0:
        raise ValueError("Sum of mole fractions must be positive.")

    if normalize:
        y = [value / y_sum for value in y]
    else:
        if abs(y_sum - 1.0) > 1e-8:
            raise ValueError(
                "Mole fractions must sum to 1, or use normalize=True.")

    if mode == "linear":
        return sum(y_i * mu_i for y_i, mu_i in zip(y, mu))

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

        return mu_mix

    else:
        raise ValueError("mode must be either 'wilke' or 'linear'.")
