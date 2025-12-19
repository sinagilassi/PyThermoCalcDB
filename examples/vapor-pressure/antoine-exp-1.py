# import libs
from pyThermoCalcDB.thermo.vapor_pressure import antoine
from pythermodb_settings.models import Temperature
from rich import print

# NOTE: Example usage of Antoine vapor pressure calculation
# antoine constants for water (example values)
A = 8.07131
B = 1730.63
C = 233.426

# temperature at which to calculate vapor pressure
temperature = Temperature(value=100, unit='C')  # 100 °C

# calculate vapor pressure using Antoine equation
result = antoine(
    A=A,
    B=B,
    C=C,
    temperature=temperature,
    output_unit='kPa', base='log10',
    message="Calculated vapor pressure for water at 100 °C"
)
if result:
    print("Antoine Vapor Pressure Calculation Result:")
    print(result.model_dump())
