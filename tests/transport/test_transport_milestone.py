import math

import numpy as np
import pytest

from pythermocalcdb.transport import (
    calc_fuller_schettler_giddings_diffusivity,
    calc_ideal_liquid_thermal_conductivity,
    calc_ideal_vapor_liquid_surface_tension,
    calc_lennard_jones_pair_diameter,
    calc_lennard_jones_pair_energy,
    calc_neufeld_diffusion_collision_integral,
    calc_reduced_collision_temperature,
    calc_stiel_thodos_gas_thermal_conductivity,
    calc_wilke_chang_diffusivity,
    calc_wilke_lee_diffusivity,
    calc_winterfeld_vapor_liquid_surface_tension,
)
from pythermocalcdb.transport.core import (
    _calc_fuller_schettler_giddings_diffusivity,
    _calc_ideal_vapor_liquid_surface_tension,
)


def test_collision_helpers():
    assert calc_lennard_jones_pair_diameter(3.0, 5.0) == pytest.approx(4.0)
    assert calc_lennard_jones_pair_energy(100.0, 400.0) == pytest.approx(200.0)
    assert calc_reduced_collision_temperature(300.0, 150.0) == pytest.approx(2.0)
    omega = calc_neufeld_diffusion_collision_integral(2.0)
    expected = (
        1.06036 / 2.0**0.15610
        + 0.19300 * math.exp(-0.47635 * 2.0)
        + 1.03587 * math.exp(-1.52996 * 2.0)
        + 1.76474 * math.exp(-3.89411 * 2.0)
    )
    assert omega == pytest.approx(expected)


def test_fuller_diffusivity_scalar_and_vector():
    result = calc_fuller_schettler_giddings_diffusivity(
        298.15,
        101325.0,
        0.0280134,
        0.0319988,
        18.5,
        16.3,
    )
    assert result > 0.0
    vector = _calc_fuller_schettler_giddings_diffusivity(
        [298.15, 320.0],
        101325.0,
        0.0280134,
        0.0319988,
        18.5,
        16.3,
    )
    assert isinstance(vector, np.ndarray)
    assert vector[1] > vector[0]


def test_wilke_lee_and_wilke_chang_diffusivity_are_positive():
    gas = calc_wilke_lee_diffusivity(
        298.15,
        101325.0,
        0.0280134,
        0.0319988,
        3.798,
        3.467,
        71.4,
        106.7,
    )
    liquid = calc_wilke_chang_diffusivity(
        298.15,
        8.9e-4,
        0.018015,
        7.5e-5,
        solvent_association_factor=2.6,
    )
    assert gas > 0.0
    assert liquid > 0.0


def test_thermal_conductivity_helpers():
    value = calc_stiel_thodos_gas_thermal_conductivity(1.0e-5, 0.028, 29.0)
    expected = (1.0e-5 / 0.028) * (1.15 * 29.0 + 0.88 * 8.314462618)
    assert value == pytest.approx(expected)
    assert calc_ideal_liquid_thermal_conductivity([0.25, 0.75], [0.1, 0.2]) == pytest.approx(0.175)


def test_surface_tension_helpers():
    ideal = calc_ideal_vapor_liquid_surface_tension([0.25, 0.75], [0.02, 0.03])
    assert ideal == pytest.approx(0.0275)
    core = _calc_ideal_vapor_liquid_surface_tension([[0.5, 0.5], [0.25, 0.75]], [[0.02, 0.03], [0.02, 0.03]])
    assert np.allclose(core, [0.025, 0.0275])
    winterfeld = calc_winterfeld_vapor_liquid_surface_tension(
        [0.5, 0.5],
        [0.02, 0.03],
        [10000.0, 10000.0],
    )
    expected = 0.25 * 0.02 + 0.5 * math.sqrt(0.02 * 0.03) + 0.25 * 0.03
    assert winterfeld == pytest.approx(expected)


def test_transport_rejects_invalid_domain():
    with pytest.raises(ValueError):
        calc_neufeld_diffusion_collision_integral(0.0)
    with pytest.raises(ValueError):
        calc_ideal_liquid_thermal_conductivity([0.2, 0.2], [0.1, 0.2])
    with pytest.raises(ValueError):
        calc_ideal_vapor_liquid_surface_tension([1.0], [0.0])
