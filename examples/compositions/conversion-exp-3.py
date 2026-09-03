"""Unit-aware composition conversion examples using CustomProp values."""
from rich import print
from pythermodb_settings.models import CustomProp

from pythermocalcdb.compositions.conversions import (
    mass_concentration_to_molarity,
    mass_fraction_to_molarity,
    molality_to_mass_fraction,
    molality_to_molarity,
    molarity_to_mass_concentration,
    molarity_to_mass_fraction,
    molarity_to_molality,
    mole_fraction_to_mass_fraction,
)


# NOTE: CustomProp values are converted only when the matching output unit is set.
molecular_weights = {
    "water": CustomProp(value=18.015, unit="g/mol"),
    "ethanol": CustomProp(value=46.07, unit="g/mol"),
}

mass_fractions = mole_fraction_to_mass_fraction(
    mole_fractions={"water": 0.5, "ethanol": 0.5},
    molecular_weights=molecular_weights,
    output_molecular_weight_unit="kg/mol",
)
print("Mole fraction to mass fraction with molecular weight conversion:")
print(mass_fractions)

# NOTE: Molarity, molecular weight, and density are converted before calculation.
solute_mass_fraction = molarity_to_mass_fraction(
    molarity=CustomProp(value=1.0, unit="mol/L"),
    molecular_weight=CustomProp(value=58.44, unit="g/mol"),
    solution_density=CustomProp(value=1020.0, unit="kg/m^3"),
    output_molarity_unit="mol/L",
    output_molecular_weight_unit="kg/mol",
    output_solution_density_unit="kg/L",
)
print("Molarity to mass fraction with unit conversion:")
print(solute_mass_fraction)

solute_molality = molarity_to_molality(
    molarity=CustomProp(value=1.0, unit="mol/L"),
    molecular_weight=CustomProp(value=58.44, unit="g/mol"),
    solution_density=CustomProp(value=1020.0, unit="kg/m^3"),
    output_molarity_unit="mol/L",
    output_molecular_weight_unit="kg/mol",
    output_solution_density_unit="kg/L",
)
print("Molarity to molality with unit conversion:")
print(solute_molality)

molarity_round_trip = molality_to_molarity(
    molality=solute_molality,
    molecular_weight=CustomProp(value=58.44, unit="g/mol"),
    solution_density=CustomProp(value=1020.0, unit="kg/m^3"),
    output_molecular_weight_unit="kg/mol",
    output_solution_density_unit="kg/L",
)
print("Molality back to molarity with unit conversion:")
print(molarity_round_trip)

print("Mass fraction to molarity with unit conversion:")
print(
    mass_fraction_to_molarity(
        mass_fraction=solute_mass_fraction,
        solution_density=CustomProp(value=1020.0, unit="kg/m^3"),
        molecular_weight=CustomProp(value=58.44, unit="g/mol"),
        output_solution_density_unit="kg/L",
        output_molecular_weight_unit="kg/mol",
    )
)

print("Molality to mass fraction with unit conversion:")
print(
    molality_to_mass_fraction(
        molality=solute_molality,
        molecular_weight=CustomProp(value=58.44, unit="g/mol"),
        output_molecular_weight_unit="kg/mol",
    )
)

mass_concentration = molarity_to_mass_concentration(
    molarity=CustomProp(value=1.0, unit="mol/L"),
    molecular_weight=CustomProp(value=58.44, unit="g/mol"),
    output_molarity_unit="mol/L",
    output_molecular_weight_unit="kg/mol",
)
print("Molarity to mass concentration with unit conversion:")
print(mass_concentration)

print("Mass concentration to molarity with unit conversion:")
print(
    mass_concentration_to_molarity(
        mass_concentration=CustomProp(value=58.44, unit="g/L"),
        molecular_weight=CustomProp(value=58.44, unit="g/mol"),
        output_mass_concentration_unit="g/L",
        output_molecular_weight_unit="g/mol",
    )
)
