# import libs
from pythermodb_settings.models import Temperature
from pyThermoCalcDB.thermo.enthalpy import (
    En_IG_NASA9_polynomial,
    En_IG_shomate,
    En_IG_NASA9_polynomial_ranges,
    En_IG_NASA7_polynomial_ranges
)
from rich import print

# NOTE: Example NASA polynomial coefficients for a hypothetical substance
a1 = 30.09200
a2 = 6.832514
a3 = 6.793435
a4 = -2.534480
a5 = 0.082139
a6 = -250.8810
a7 = 223.3967
b1 = -242.7400
b2 = 0.0  # not used in enthalpy calculation

# NOTE: Define the temperature at which to calculate enthalpy
temperature = Temperature(value=300.0, unit="K")
# temperature range for validity check (optional)
T_min = Temperature(value=200.0, unit="K")
T_max = Temperature(value=1000.0, unit="K")
temperature_range = (T_min, T_max)

# NOTE: Calculate the ideal gas enthalpy using the NASA polynomial equation
result = En_IG_NASA9_polynomial(
    a1=a1,
    a2=a2,
    a3=a3,
    a4=a4,
    a5=a5,
    a6=a6,
    a7=a7,
    b1=b1,
    b2=b2,
    temperature=temperature,
    temperature_range=temperature_range,
    output_unit="J/mol",
    message="NASA enthalpy calculation successful"
)
# Print the result
print(result)


# NOTE: Example Shomate coefficients for a hypothetical substance
A = 25.56759
B = 6.096130
C = 4.054656
D = -2.671301
E = 0.131021
F = -30.327060
G = 223.3967

# NOTE: Calculate the ideal gas enthalpy using the Shomate equation
result_shomate = En_IG_shomate(
    A=A,
    B=B,
    C=C,
    D=D,
    E=E,
    F=F,
    G=G,
    temperature=temperature,
    temperature_range=temperature_range,
    output_unit="kJ/mol",
    message="Shomate enthalpy calculation successful"
)
# Print the result
print(result_shomate)

# SECTION: calculate enthalpy at different temperature ranges
# NOTE: Define temperature ranges for NASA 9-coefficient polynomial
temperatures = [
    Temperature(value=250.0, unit="K"),
    Temperature(value=500.0, unit="K"),
    Temperature(value=800.0, unit="K"),
    Temperature(value=1000.0, unit="K"),
]

result = En_IG_NASA9_polynomial_ranges(
    a1=a1,
    a2=a2,
    a3=a3,
    a4=a4,
    a5=a5,
    a6=a6,
    a7=a7,
    b1=b1,
    b2=b2,
    temperatures=temperatures,
)
# Print the results
print(result)
