# import libs
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from pythermocalcdb.thermo.viscosity import liquid_mixture_viscosity
from pythermodb_settings.models import CustomProp
from rich import print


# NOTE: Define liquid-mixture composition
mole_fractions = [0.40, 0.60]

# NOTE: Define pure-component liquid viscosities with the same unit
viscosities = [
    CustomProp(value=0.544, unit="cP"),  # e.g., methanol
    CustomProp(value=1.074, unit="cP"),  # e.g., ethanol
]

# NOTE: Perform logarithmic liquid-mixture viscosity calculation
result = liquid_mixture_viscosity(
    mole_fractions=mole_fractions,
    viscosities=viscosities,
    mode="log",
)

# NOTE: Print the result
print(result)
