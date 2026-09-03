"""Basic composition conversion examples using sequence inputs."""
from rich import print
from pythermocalcdb.compositions.conversions import (
    mass_concentration_to_molarity,
    mass_fraction_to_molality,
    mass_fraction_to_mole_fraction,
    mass_fraction_to_ppm,
    mass_fraction_to_weight_percent,
    molality_to_molarity,
    molality_to_mass_fraction,
    molarity_to_mass_concentration,
    molarity_to_mass_fraction,
    molarity_to_molality,
    mole_fraction_to_mass_fraction,
    mole_fraction_to_molality,
    mole_fraction_to_mole_percent,
    ppm_mass_to_mass_fraction,
    weight_percent_to_mass_fraction,
)


# NOTE: Sequence inputs return lists in the same order as the input sequence.
mole_fractions = [0.5, 0.5]
molecular_weights = [0.018015, 0.04607]  # kg/mol
mass_fractions = mole_fraction_to_mass_fraction(
    mole_fractions,
    molecular_weights,
)
print("Mole fraction to mass fraction:")
print(mass_fractions)

mole_fractions_round_trip = mass_fraction_to_mole_fraction(
    mass_fractions,
    molecular_weights,
)
print("Mass fraction to mole fraction:")
print(mole_fractions_round_trip)

# NOTE: Scalar conversions use numerically consistent units.
molarity = 1.0  # mol/L
solute_molecular_weight = 0.05844  # kg/mol
solution_density = 1.02  # kg/L

molality = molarity_to_molality(
    molarity,
    solute_molecular_weight,
    solution_density,
)
print("Molarity to molality:")
print(molality)

molarity_round_trip = molality_to_molarity(
    molality,
    solute_molecular_weight,
    solution_density,
)
print("Molality to molarity:")
print(molarity_round_trip)

mass_fraction = molarity_to_mass_fraction(
    molarity,
    solute_molecular_weight,
    solution_density,
)
print("Molarity to mass fraction:")
print(mass_fraction)

print("Mass fraction to molality:")
print(mass_fraction_to_molality(mass_fraction, solute_molecular_weight))

print("Molality to mass fraction:")
print(molality_to_mass_fraction(molality, solute_molecular_weight))

print("Molarity to mass concentration:")
print(molarity_to_mass_concentration(molarity, solute_molecular_weight))

print("Mass concentration to molarity:")
print(mass_concentration_to_molarity(0.05844, solute_molecular_weight))

print("Mass fraction to weight percent:")
print(mass_fraction_to_weight_percent(mass_fraction))

print("Weight percent to mass fraction:")
print(weight_percent_to_mass_fraction(5.0))

print("Mass fraction to ppm:")
print(mass_fraction_to_ppm(1e-4))

print("ppm to mass fraction:")
print(ppm_mass_to_mass_fraction(100.0))

print("Mole fraction to molality:")
print(mole_fraction_to_molality(0.02, 0.98, 0.01801528))

print("Mole fraction to mole percent:")
print(mole_fraction_to_mole_percent(0.25))
