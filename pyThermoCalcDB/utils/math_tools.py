from typing import Callable, TypedDict, Literal
import numpy as np
from scipy.integrate import quad, simps


class QuantityResult(TypedDict):
    value: float
    unit: str
    symbol: str
    prop_name: str


def integrate_quantity_function(
    func: Callable[..., QuantityResult],
    var: str,
    a: float,
    b: float,
    *,
    integral_mode: Literal["quad", "trapz", "simps"] = "quad",
    n_points: int = 1000,
    **kwargs
) -> QuantityResult:
    """
    Integrate a quantity-returning function over [a, b].

    Parameters
    ----------
    func : Callable[..., QuantityResult]
        Function that returns a dict with keys:
        value, unit, symbol, prop_name.
        The integration is applied only to the `value` field.
    var : str
        Name of the integration variable (e.g. 'T', 'x').
    a : float
        Lower integration limit.
    b : float
        Upper integration limit.
    integral_mode : {'quad', 'trapz', 'simps'}, optional
        Numerical integration method:
        - 'quad'  : scipy.integrate.quad (adaptive, default)
        - 'trapz' : numpy.trapz on an evenly spaced grid
        - 'simps' : scipy.integrate.simps on an evenly spaced grid
    n_points : int, optional
        Number of grid points for 'trapz' and 'simps'.
    **kwargs :
        Extra parameters forwarded to `func`.

    Returns
    -------
    QuantityResult
        Same structure as func(...) but with `value` replaced by
        the definite integral over [a, b].
    """

    # ---------- Choose integration core ----------
    if integral_mode == "quad":
        # quad expects a float -> float function
        def integrand(x: float) -> float:
            result = func(x, **kwargs)
            return result["value"]

        integral_value, _err = quad(integrand, a, b)

    else:
        # Grid-based methods: build x-grid and evaluate func
        x = np.linspace(a, b, n_points)
        y = np.empty_like(x)

        for i, xi in enumerate(x):
            result = func(xi, **kwargs)
            y[i] = result["value"]

        if integral_mode == "trapz":
            integral_value = np.trapz(y, x)
        elif integral_mode == "simps":
            integral_value = simps(y, x)
        else:
            raise ValueError(f"Unknown integral_mode: {integral_mode}")

    # ---------- Build output metadata ----------
    # Call once (e.g. at midpoint) to copy metadata
    sample = func((a + b) / 2.0, **kwargs)

    return QuantityResult(
        value=integral_value,
        # You can adjust this to your unit system
        unit=f"{sample['unit']} * {var}",
        symbol=f"∫ {sample['symbol']} d{var}",
        prop_name=f"integral_of_{sample['prop_name']}",
    )
