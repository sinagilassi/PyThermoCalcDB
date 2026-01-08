# Examples

The repository ships runnable scripts under `examples/` that load ThermoDB pickle files, build a `ModelSource`, and call the helper functions. Make sure the `examples/thermodb` directory is available or point the scripts at your own ThermoDB pickles.

Run any script from the project root:

```bash
cd examples
python exp-3.py
```

All scripts share the same setup pattern:

```python
import os
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource

parent_dir = os.path.dirname(os.path.abspath(__file__))
thermodb_dir = os.path.join(parent_dir, "thermodb")

CO2 = Component(name="carbon dioxide", formula="CO2", state="g")
CO2_src = ComponentThermoDBSource(
    component=CO2,
    source=os.path.join(thermodb_dir, "carbon dioxide-g.pkl"),
)

rules = {
    "ALL": {
        "DATA": {"critical-pressure": "Pc", "critical-temperature": "Tc", "acentric-factor": "AcFa"},
        "EQUATIONS": {"CUSTOM-REF-1::vapor-pressure": "VaPr", "CUSTOM-REF-1::ideal-gas-heat-capacity": "Cp_IG"},
    }
}

model_source: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=[CO2_src],
    rules=rules,
)
```

## Script tour

| Script | Focus | Key helpers |
| --- | --- | --- |
| `examples/exp-1.py` | Sanity check that PyThermoCalcDB imports and reports its version. | n/a |
| `examples/exp-2.py` | Builds a multi-component `ModelSource`, evaluates vapor pressure / Cp equations, and computes enthalpy and Gibbs free energy at a single temperature and across a temperature range. | `thermo.calc_En`, `thermo.calc_GiFrEn`, `thermo.calc_En_range`, `thermo.calc_GiFrEn_range` |
| `examples/exp-3.py` | Vapor-pressure focused sweep for CO2: Psat at several temperatures, enthalpy of vaporization, saturated temperature vs pressure, and dPsat/dT. | `sat.calc_VaPr`, `sat.calc_EnVap`, `sat.calc_T_sat`, `sat.calc_VaPr_sensitivity` |
| `examples/calc-hsg.py` | Computes enthalpy and entropy changes relative to reference T/P and reproduces formation enthalpy/Gibbs values. | `thermo.calc_dEn`, `thermo.calc_dEnt`, `thermo.calc_En`, `thermo.calc_GiFrEn` |
| `examples/calc-hsg-mix.py` | Mixture enthalpy for a CO2/CO ideal-gas blend with mole fractions and optional unit conversion. | `thermo.calc_En_mix` |
| `examples/calc-sat.py` | Saturation properties for ethanol across a temperature sweep; also inspects vapor-pressure equation ranges. | `sat.calc_VaPr`, `sat.calc_EnVap`, `sat.calc_T_sat`, `sat.calc_VaPr_sensitivity` |
| `examples/reaction/exp-1.py` | Standard enthalpy of reaction for methane combustion at 298.15 K. | `reactions.dEn_rxn_STD` |

Notes:
- Scripts import verbose names like `calc_vapor_pressure_at_temperature`; the underlying implementations live in `pyThermoCalcDB.docs.sat` and `pyThermoCalcDB.docs.thermo` as listed in the method reference.
- Temperatures are specified in Kelvin and pressures in Pascal unless noted; conversions are handled internally.
