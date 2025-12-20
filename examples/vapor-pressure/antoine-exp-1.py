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

# temperature range for validity check (optional)
T_min = Temperature(value=0, unit='C')
T_max = Temperature(value=374, unit='C')
temperature_range = (T_min, T_max)

# calculate vapor pressure using Antoine equation
result = antoine(
    A=A,
    B=B,
    C=C,
    temperature=temperature,
    temperature_range=temperature_range,
    output_unit='kPa', base='log10',
    message="Calculated vapor pressure for water at 100 °C"
)
if result:
    print("Antoine Vapor Pressure Calculation Result:")
    print(result)
