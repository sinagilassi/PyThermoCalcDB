# import libs
from pyThermoCalcDB.thermo.density import rackett
from pythermodb_settings.models import (
    Temperature,
    Pressure,
    CustomProp,
)
from rich import print

# NOTE: Define input parameters for the Rackett density calculation
temperature = Temperature(value=350.0, unit="K")
critical_temperature = Temperature(value=500.0, unit="K")
critical_pressure = Pressure(value=40.0, unit="bar")

molecular_weight = CustomProp(value=44.01, unit="g/mol")  # e.g., CO2
critical_compressibility = CustomProp(value=0.274, unit="dimensionless")

# NOTE: Perform the Rackett density calculation
result = rackett(
    temperature=temperature,
    critical_temperature=critical_temperature,
    critical_pressure=critical_pressure,
    molecular_weight=molecular_weight,
    critical_compressibility=critical_compressibility,
    message="Rackett density calculation successful"
)
# Print the result
print(result)
