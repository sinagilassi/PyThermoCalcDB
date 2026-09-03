
# import libs
from rich import print
from pythermodb_settings.models import CustomProp
from pythermocalcdb.compositions import calc_ionic_strength_molarity


# CustomProp values are normalized to the requested output concentration unit.
molarities = {
    "Na+": CustomProp(value=0.10, unit="mol/L"),
    "Cl-": CustomProp(value=0.10, unit="mol/L"),
}
charges = {"Na+": 1.0, "Cl-": -1.0}

ionic_strength = calc_ionic_strength_molarity(
    molarities,
    charges,
    output_molarity_unit="mol/L",
)

print(f"Molarity-based ionic strength: {ionic_strength:.3f} mol/L")
assert ionic_strength == 0.1
