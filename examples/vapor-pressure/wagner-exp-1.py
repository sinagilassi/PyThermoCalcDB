# import libs
from pyThermoCalcDB.thermo.vapor_pressure import wagner
from pythermodb_settings.models import Temperature, Pressure
from rich import print

# NOTE: Example usage of Wagner vapor pressure calculation
# Wagner constants for Propane (from provided data)
A = -6.76368
B = 1.55481
C = -1.5872
D = -2.024

# temperature at which to calculate vapor pressure
# Example temperature for Propane
temperature = Temperature(value=300, unit='K')
# critical temperature and pressure for Propane
critical_temperature = Temperature(value=369.85, unit='K')  # K
critical_pressure = Pressure(value=4247, unit='kPa')  # 42.47 bar = 4247 kPa

# temperature range for validity check (optional)
T_min = Temperature(value=200, unit='K')
T_max = Temperature(value=265, unit='K')
temperature_range = (T_min, T_max)

# calculate vapor pressure using Wagner equation
result = wagner(
    A=A,
    B=B,
    C=C,
    D=D,
    temperature=temperature,
    critical_temperature=critical_temperature,
    critical_pressure=critical_pressure,
    output_unit='kPa',
    message="Calculated vapor pressure for Propane"
)
if result:
    print("Wagner Vapor Pressure Calculation Result:")
    print(result)
