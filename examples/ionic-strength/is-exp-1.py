# import libs
from rich import print
from pythermocalcdb.compositions import calc_ionic_strength_molality


# Ionic strength on a molality basis:
# I = 0.5 * (m_Na+ z_Na+^2 + m_Cl- z_Cl-^2 + m_Ca2+ z_Ca2+^2)
molalities = {
    "Na+": 0.10,
    "Cl-": 0.10,
    "Ca2+": 0.05,
}
charges = {
    "Na+": 1.0,
    "Cl-": -1.0,
    "Ca2+": 2.0,
}

ionic_strength = calc_ionic_strength_molality(molalities, charges)

print(f"Molality-based ionic strength: {ionic_strength:.3f} mol/kg")
assert ionic_strength == 0.2
