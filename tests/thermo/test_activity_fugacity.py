import math

import numpy as np
import pytest

from pythermocalcdb.thermo import (
    calc_activity_from_concentration,
    calc_activity_from_mole_fraction,
    calc_effective_concentration,
    calc_liquid_fugacity_coefficient,
    calc_liquid_partial_fugacity,
    calc_poynting_factor_from_integral,
    calc_poynting_factor_incompressible,
)
from pythermocalcdb.thermo.core import (
    _calc_activity_from_mole_fraction,
    _calc_poynting_factor_incompressible,
)


def test_activity_from_mole_fraction_scalar_and_vector():
    assert calc_activity_from_mole_fraction(0.4, 1.2) == pytest.approx(0.48)
    result = _calc_activity_from_mole_fraction([0.2, 0.8], [1.5, 0.75])
    assert isinstance(result, np.ndarray)
    assert np.allclose(result, [0.3, 0.6])


def test_activity_from_concentration_and_effective_concentration():
    assert calc_activity_from_concentration(2.0, 0.5, 4.0) == pytest.approx(0.25)
    assert calc_effective_concentration(0.25, 4.0) == pytest.approx(1.0)


def test_activity_rejects_invalid_domain():
    with pytest.raises(ValueError):
        calc_activity_from_mole_fraction(-0.1, 1.0)
    with pytest.raises(ValueError):
        calc_activity_from_concentration(1.0, 0.0, 1.0)


def test_poynting_factor_helpers():
    expected = math.exp(1.0e-4 * (2.0e5 - 1.0e5) / (8.314462618 * 300.0))
    assert calc_poynting_factor_incompressible(1.0e-4, 2.0e5, 1.0e5, 300.0) == pytest.approx(expected)
    assert calc_poynting_factor_from_integral(10.0, 300.0) == pytest.approx(
        math.exp(10.0 / (8.314462618 * 300.0))
    )
    vector = _calc_poynting_factor_incompressible([1.0e-4, 2.0e-4], 2.0e5, 1.0e5, 300.0)
    assert isinstance(vector, np.ndarray)
    assert vector.shape == (2,)


def test_liquid_fugacity_transformations():
    coeff = calc_liquid_fugacity_coefficient(1.2, 0.95, 5.0e4, 1.0e5, 1.1)
    assert coeff == pytest.approx(1.2 * 0.95 * 0.5 * 1.1)
    fugacity = calc_liquid_partial_fugacity(0.25, 1.2, 0.95, 5.0e4, 1.1)
    assert fugacity == pytest.approx(0.25 * 1.2 * 0.95 * 5.0e4 * 1.1)


def test_fugacity_rejects_invalid_domain():
    with pytest.raises(ValueError):
        calc_poynting_factor_incompressible(1.0e-4, 1.0e5, 1.0e5, 0.0)
    with pytest.raises(ValueError):
        calc_liquid_fugacity_coefficient(1.0, 1.0, 1.0e5, 0.0, 1.0)
