# import libs
from pythermodb_settings.models import Temperature
from pyThermoCalcDB.thermo.enthalpy import (
    En_IG_NASA9_polynomial,
    En_IG_shomate,
    En_IG_NASA9_polynomial_ranges,
    En_IG_NASA7_polynomial_ranges,
    En_IG_NASA9_polynomial_range,
    En_IG_NASA7_polynomial_range,
)
from rich import print

# NOTE: Example NASA polynomial coefficients for a substance
# ! CO2 coefficients (NASA9)
#  49437.8364, -626.429208, 5.30181336, 0.002503601, -2.12e-07, -7.69e-10, 2.85e-13, -45281.8986, -7.0487901,
a1 = 49437.8364
a2 = -626.429208
a3 = 5.30181336
a4 = 0.002503601
a5 = -2.12e-07
a6 = -7.69e-10
a7 = 2.85e-13
b1 = -45281.8986
b2 = -7.0487901

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


result = En_IG_NASA9_polynomial_range(
    a1=a1,
    a2=a2,
    a3=a3,
    a4=a4,
    a5=a5,
    a6=a6,
    a7=a7,
    b1=b1,
    b2=b2,
    T_low=Temperature(value=300.0, unit="K"),
    T_high=Temperature(value=800.0, unit="K"),
)
# Print the result
print(result)
