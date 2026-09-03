"""Component-keyed composition conversion examples using mapping inputs."""
from rich import print
from pythermodb_settings.models import Component

from pythermocalcdb.compositions.conversions import (
    mass_fraction_to_mole_fraction,
    molality_to_mole_fraction,
    molarities_to_molalities,
    mole_fraction_to_mass_fraction,
)


components = [
    Component(name="carbon dioxide", formula="CO2", state="g"),
    Component(name="methane", formula="CH4", state="g"),
    Component(name="oxygen", formula="O2", state="g"),
]

# NOTE: Inputs are intentionally not in component order. When component_key is
# provided, mapping outputs are remapped and ordered by the components list.
mole_fractions = {
    "oxygen": 0.20,
    "carbon dioxide": 0.30,
    "methane": 0.50,
}
molecular_weights = {
    "oxygen": 0.03200,
    "carbon dioxide": 0.04401,
    "methane": 0.01604,
}

mass_fractions = mole_fraction_to_mass_fraction(
    mole_fractions=mole_fractions,
    molecular_weights=molecular_weights,
    components=components,
    component_key="Formula",
)
print("Mole fraction to mass fraction by formula:")
print(mass_fractions)
# print(list(mass_fractions.values()))

mole_fractions_from_mass = mass_fraction_to_mole_fraction(
    mass_fractions=mass_fractions,
    molecular_weights={
        "CO2": 0.04401,
        "CH4": 0.01604,
        "O2": 0.03200,
    },
    components=components,
    component_key="Formula",
)
print("Mass fraction to mole fraction by formula:")
print(mole_fractions_from_mass)
# print(list(mole_fractions_from_mass.values()))

molarities = {
    "oxygen": 0.10,
    "carbon dioxide": 0.25,
    "methane": 0.15,
}
molalities = molarities_to_molalities(
    molarities=molarities,
    molecular_weights=molecular_weights,
    solution_density=1.10,
    components=components,
    component_key="Formula",
)
print("Molarities to molalities by formula:")
print(molalities)
# print(list(molalities.values()))

mole_fractions_with_solvent = molality_to_mole_fraction(
    molalities=molalities,
    solvent_molecular_weight=0.01801528,
    solvent_key="H2O",
    components=components,
    component_key="Formula",
)
print("Molalities to mole fractions with solvent:")
print(mole_fractions_with_solvent)
# print(list(mole_fractions_with_solvent.values()))

# NOTE: Without component_key, mapping output keeps the input mapping keys and order.
plain_order = molality_to_mole_fraction(
    molalities={"oxygen": 0.10, "carbon dioxide": 0.25, "methane": 0.15},
    solvent_molecular_weight=0.01801528,
)
print("Molalities to mole fractions without component mapping:")
print(plain_order)
