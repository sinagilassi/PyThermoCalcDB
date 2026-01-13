# import libs
from rich import print
from pythermodb_settings.models import Temperature
from pyThermoCalcDB.thermo.entropy import (
    S_IG_NASA9_polynomial,
    S_IG_NASA9_polynomial_ranges,
    dS_IG_NASA7_polynomial,
    S_IG_NASA7_polynomial_ranges,
)

# NOTE: Example NASA polynomial coefficients for a hypothetical substance
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
# temperature range
temperatures = [
    Temperature(value=250.0, unit="K"),
    Temperature(value=500.0, unit="K"),
    Temperature(value=800.0, unit="K"),
    Temperature(value=1000.0, unit="K"),
]

# NOTE: calculate the ideal gas Entropy using the NASA polynomial equation
result = S_IG_NASA9_polynomial(
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
    output_unit="J/mol.K",
    message="NASA entropy calculation successful"
)
# Print the result
print(result)

# NOTE: calculate the ideal gas Entropy over a range of temperatures
result_range = S_IG_NASA9_polynomial_ranges(
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
    output_unit="J/mol.K",
    message="NASA entropy range calculation successful"
)
# Print the result
print(result_range)
