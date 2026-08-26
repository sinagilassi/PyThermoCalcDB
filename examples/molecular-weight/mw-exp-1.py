# import libs
from rich import print
from pythermocalcdb.thermo.molecular_weight import calc_molecular_weight, calc_MW

# NOTE: Compounds

compounds = [
    "Fe(OH)3",
    "Ca3(PO4)2",
    "CuSO4*5H2O",
    "SO4{2-}",
    "e{-}",
    "H2O",
    "C6H12O6",
    "Fe{3+}",
    "H{+}",
    "OH{-}",
]

# NOTE: Calculate and print molecular weights
for compound in compounds:
    mw = calc_MW(
        compound,
        include_electron_mass=True,
        decimal_digits=3
    )
    print(f"Molecular weight of {compound}: {mw}")
