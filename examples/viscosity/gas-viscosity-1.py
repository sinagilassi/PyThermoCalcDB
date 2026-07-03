# import libs
from pythermocalcdb.thermo.viscosity import gas_mixture_viscosity
from pythermodb_settings.models import CustomProp
from rich import print


# NOTE: Define gas-mixture composition
mole_fractions = [0.50, 0.30, 0.20]

# NOTE: Define pure-component gas viscosities with the same unit
viscosities = [
    CustomProp(value=1.10e-5, unit="Pa.s"),  # e.g., methane
    CustomProp(value=1.76e-5, unit="Pa.s"),  # e.g., nitrogen
    CustomProp(value=2.05e-5, unit="Pa.s"),  # e.g., oxygen
]

# NOTE: Define molecular weights with the same unit
molecular_weights = [
    CustomProp(value=16.04, unit="g/mol"),  # methane
    CustomProp(value=28.01, unit="g/mol"),  # nitrogen
    CustomProp(value=32.00, unit="g/mol"),  # oxygen
]

# NOTE: Perform Wilke gas-mixture viscosity calculation
result = gas_mixture_viscosity(
    mole_fractions=mole_fractions,
    viscosities=viscosities,
    molecular_weights=molecular_weights,
    mode="wilke",
)

# NOTE: Print the result
print(result)
