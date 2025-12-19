# import libs
from pyThermoCalcDB.thermo.heat_capacity import Cp_IG_polynomial, Cp_IG_NASA_polynomial
from pythermodb_settings.models import Temperature
from rich import print

# NOTE: Example coefficients for a hypothetical substance
A = 25.0  # J/mol.K
B = 1.0e-2  # J/mol.K^2
C = 2.0e-5  # J/mol.K^3
D = -1.0e-8  # J/mol.K^4
E = 0.0  # J·K/mol

# NOTE: Define the temperature at which to calculate heat capacity
temperature = Temperature(value=300.0, unit="K")

# temperature range for validity check (optional)
T_min = Temperature(value=200.0, unit="K")
T_max = Temperature(value=1000.0, unit="K")
temperature_range = (T_min, T_max)

# NOTE: Calculate the ideal gas heat capacity using the polynomial equation
result = Cp_IG_polynomial(
    A=A,
    B=B,
    C=C,
    D=D,
    E=E,
    temperature=temperature,
    temperature_range=temperature_range,
    output_unit="J/mol.K",
    message="Calculation successful"
)
# Print the result
print(result)


# NOTE: Example NASA polynomial coefficients for a hypothetical substance
a1 = 30.09200
a2 = 6.832514
a3 = 6.793435
a4 = -2.534480
a5 = 0.082139

# NOTE: Calculate the ideal gas heat capacity using the NASA polynomial equation
nasa_result = Cp_IG_NASA_polynomial(
    a1=a1,
    a2=a2,
    a3=a3,
    a4=a4,
    a5=a5,
    temperature=temperature,
    temperature_range=temperature_range,
    output_unit="J/mol.K",
    message="NASA calculation successful"
)
# Print the NASA result
print(nasa_result)
