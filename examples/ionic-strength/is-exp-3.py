# import libs
from rich import print
from pythermocalcdb.compositions import calc_charge_balance, check_electroneutrality


concentrations = {"Na+": 0.10, "Cl-": 0.10}
charges = {"Na+": 1.0, "Cl-": -1.0}

charge_balance = calc_charge_balance(concentrations, charges)
neutral = check_electroneutrality(
    concentrations,
    charges,
    tolerance=1e-12,
)

print(f"Charge-balance residual: {charge_balance:.3g} mol/L")
print(f"Electroneutral: {neutral}")
assert charge_balance == 0.0
assert neutral
