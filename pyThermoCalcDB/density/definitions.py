"""Density definition public wrappers."""

# import libs
from pythermodb_settings.models import AnnotatedValue, ScalarValue
from pythermodb_settings.utils import to_annotated_value

# locals
from ..utils.conversions import _pos
from .core.definitions import (
    _calc_density_from_molar_volume,
    _calc_gas_density_from_compressibility,
    _calc_multiphase_density_from_volume_fractions,
)


# SECTION: Public wrappers
def calc_density_from_molar_volume(
    molar_mass: ScalarValue,
    molar_volume: ScalarValue,
    *,
    name: str = "density",
    description: str = "Calculate mass density from molar mass and molar volume.",
    unit: str | None = "kg/m3",
    symbol: str | None = "rho",
) -> AnnotatedValue[float]:
    """Calculate annotated density ``rho = M / Vm``."""
    return to_annotated_value(
        float(_calc_density_from_molar_volume(
            _pos(molar_mass, "molar_mass"),
            _pos(molar_volume, "molar_volume"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_density_from_molar_volume",
    )


def calc_gas_density_from_compressibility(
    molar_mass: ScalarValue,
    pressure: ScalarValue,
    temperature: ScalarValue,
    compressibility_factor: ScalarValue,
    gas_constant: ScalarValue = 8.31446261815324,
    *,
    name: str = "gas_density",
    description: str = "Calculate gas density from compressibility factor.",
    unit: str | None = "kg/m3",
    symbol: str | None = "rho",
) -> AnnotatedValue[float]:
    """Calculate annotated gas density ``rho = M*P/(Z*R*T)``."""
    return to_annotated_value(
        float(_calc_gas_density_from_compressibility(
            _pos(molar_mass, "molar_mass"),
            _pos(pressure, "pressure"),
            _pos(temperature, "temperature"),
            _pos(compressibility_factor, "compressibility_factor"),
            _pos(gas_constant, "gas_constant"),
        )),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_gas_density_from_compressibility",
    )


def calc_multiphase_density_from_volume_fractions(
    phase_densities,
    phase_volume_fractions,
    *,
    name: str = "multiphase_density",
    description: str = "Calculate multiphase density from phase volume fractions.",
    unit: str | None = "kg/m3",
    symbol: str | None = "rho_m",
) -> AnnotatedValue[float]:
    """Calculate annotated multiphase density from phase volume fractions."""
    return to_annotated_value(
        float(_calc_multiphase_density_from_volume_fractions(phase_densities, phase_volume_fractions)),
        name=name,
        description=description,
        unit=unit,
        symbol=symbol,
        implementation="_calc_multiphase_density_from_volume_fractions",
    )


__all__ = [
    "calc_density_from_molar_volume",
    "calc_gas_density_from_compressibility",
    "calc_multiphase_density_from_volume_fractions",
]
