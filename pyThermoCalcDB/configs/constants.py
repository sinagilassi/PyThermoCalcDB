# import libs
from pythermodb_settings.models import Temperature

# SECTION: PyThermoDBLink/PyThermoDB
DATASOURCE = "datasource"
EQUATIONSOURCE = "equationsource"


# SECTION: constants
R_J_molK = 8.314462618  # universal gas constant in J/mol.K
T_ref_K = 298.15  # reference temperature in K
P_ref_Pa = 101325.0  # reference pressure in Pa

# NOTE: reference temperature
T_298_K: Temperature = Temperature(
    value=298.15,
    unit="K"
)
