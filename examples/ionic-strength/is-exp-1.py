# import libs
from rich import print
from pythermodb_settings.models import Component
from pythermocalcdb.compositions import (
    _calc_ionic_strength_molality,
    calc_ionic_strength_molality,
    calc_ionic_strength_molality_2
)


# Ionic strength on a molality basis:
# I = 0.5 * (m_Na+ z_Na+^2 + m_Cl- z_Cl-^2 + m_Ca2+ z_Ca2+^2)
molalities = {
    "Na{+}": 0.10,
    "Cl{-}": 0.10,
    "Ca{2+}": 0.05,
}

# list
molality_list = list(molalities.values())

charges = {
    "Na{+}": 1.0,
    "Cl{-}": -1.0,
    "Ca{2+}": 2.0,
}

# Components
sodium = Component(
    name="Sodium",
    formula="Na{+}",
    state="aq"
)
chloride = Component(
    name="Chloride",
    formula="Cl{-}",
    state="aq"
)

calcium = Component(
    name="Calcium",
    formula="Ca{2+}",
    state="aq"
)
# component list
component_list = [sodium, chloride, calcium]

# NOTE: return float
ionic_strength = _calc_ionic_strength_molality(molalities, charges)

print(f"Molality-based ionic strength (float): {ionic_strength} mol/kg")

# NOTE: return AnnotatedValue
ionic_strength = calc_ionic_strength_molality(molalities, charges)

print(f"Molality-based ionic strength: {ionic_strength} mol/kg")

# NOTE: return AnnotatedValue using the alternative molality calculation method
ionic_strength_2 = calc_ionic_strength_molality_2(
    molality_list,
    component_list
)

print(
    f"Molality-based ionic strength (method 2): {ionic_strength_2} mol/kg"
)
