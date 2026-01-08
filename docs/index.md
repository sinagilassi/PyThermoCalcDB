# PyThermoCalcDB

PyThermoCalcDB is a thermodynamic calculation layer built on PyThermoDB and PyThermoLinkDB. It wraps curated equations and data into quick helpers for enthalpy, Gibbs energy, entropy, vapor pressure, and reaction energetics.

## What you get
- Uses ThermoDB pickle sources through PyThermoLinkDB `ModelSource`.
- Typed inputs from `pythermodb_settings` (`Component`, `Temperature`, `Pressure`, `CustomProp`).
- Helper modules for pure components, mixtures, saturation/phase-change, and reactions.

## Install

```bash
pip install pyThermoCalcDB
```

## Quickstart

1. Collect ThermoDB pickle paths for the components you care about.
2. Build a `ModelSource` with `pyThermoLinkDB.load_and_build_model_source`.
3. Call helpers from `pyThermoCalcDB.docs.thermo` or `pyThermoCalcDB.docs.sat`.

```python
import os
import pyThermoLinkDB as ptdblink
from pyThermoLinkDB.models import ModelSource
from pythermodb_settings.models import Component, ComponentRule, ComponentThermoDBSource, Temperature, Pressure
from pyThermoCalcDB.docs import thermo, sat

CO2 = Component(name="carbon dioxide", formula="CO2", state="g")
thermodb_dir = "/path/to/thermodb"
CO2_src = ComponentThermoDBSource(
    component=CO2,
    source=os.path.join(thermodb_dir, "carbon dioxide-g.pkl"),
)

model_source: ModelSource = ptdblink.load_and_build_model_source(
    thermodb_sources=[CO2_src],
    original_equation_label=False
)

T = Temperature(value=300.0, unit="K")
P = Pressure(value=101325.0, unit="Pa")

enthalpy = thermo.calc_En(
    component=CO2,
    model_source=model_source,
    temperature=T,
)
Pvap = sat.calc_VaPr(
    component=CO2,
    model_source=model_source,
    temperature=T,
)
Tsat = sat.calc_T_sat(
    component=CO2,
    model_source=model_source,
    pressure=P,
    temperature_guess=T,
)
print(enthalpy)
print(Pvap)
print(Tsat)
```

## Helper modules

| Area | Module | Highlights |
| --- | --- | --- |
| Thermodynamic properties | `pyThermoCalcDB.docs.thermo` | enthalpy/Gibbs at a temperature or range, enthalpy and entropy changes, mixture enthalpy with optional departure/excess terms |
| Saturation & phase change | `pyThermoCalcDB.docs.sat` | vapor pressure, enthalpy of vaporization, saturated temperature solving, vapor-pressure sensitivity |
| Reactions | `pyThermoCalcDB.reactions.reactions` | standard reaction enthalpy and Gibbs free energy from stoichiometry |

## Next steps
- See `docs.md` for the full method reference and signatures.
- Skim `examples.md` to pick the script that matches your use case.
