# Method Reference 🧱

All helpers take typed values from `pythermodb_settings` and expect a `ModelSource` built with PyThermoLinkDB. The optional `mode` kwarg supports `silent`, `log`, or `attach` timing/logging behavior.

## Thermodynamic properties (`pyThermoCalcDB.docs.thermo`) 🔥

| Function | What it does | Key inputs | Output (units) |
| --- | --- | --- | --- |
| `calc_dEn` | Enthalpy change (Delta H) between two temperatures (integrates Cp) | `component`, `model_source`, `temperature_initial`, `temperature_final`, `component_key` | `ComponentEnthalpyChange` in J/mol |
| `calc_En` | Enthalpy at a temperature (formation enthalpy plus Cp integral) | `component`, `model_source`, `temperature`, `component_key` | `ComponentEnthalpy` in J/mol |
| `calc_GiFrEn` | Gibbs free energy at a temperature for a phase | `component`, `model_source`, `temperature`, `phase` (`IG`/`LIQ`/`SOL`), `component_key` | `ComponentGibbsFreeEnergy` in J/mol |
| `calc_En_range` | Enthalpy values across a list of temperatures | `component`, `model_source`, `temperatures` | `list[ComponentEnthalpy]` |
| `calc_GiFrEn_range` | Gibbs free energies across a list of temperatures | `component`, `model_source`, `temperatures` | `list[ComponentGibbsFreeEnergy]` |
| `calc_dEnt` | Entropy change (Delta S) between two T/P states | `component`, `model_source`, `temperature_initial`, `temperature_final`, `pressure_initial`, `pressure_final`, `phase` | `ComponentEntropyChange` in J/(mol*K) |
| `calc_En_mix` | Mixture enthalpy with optional departure/excess and unit conversion | `components` (with `mole_fraction`), `model_source`, `temperature`, `pressure`, `reference` (`IG` or `None`), optional `departure_enthalpy` / `excess_enthalpy`, `component_key`, `output_unit` | `MixtureEnthalpyResult` |

Notes 🧠:
- `component_key` defaults to "Name-Formula" for thermo helpers; pass a different key if your ThermoDB uses another identifier.
- Reference temperature is 298.15 K inside the HSG formulations.
- Use `output_unit` in `calc_En_mix` to get convenient units such as `kJ/mol`.

Example 🎯:

```python
from pyThermoCalcDB.docs import thermo
from pythermodb_settings.models import Component, Temperature, Pressure

# model_source built via pyThermoLinkDB
CO2 = Component(name="carbon dioxide", formula="CO2", state="g")
CO = Component(name="carbon monoxide", formula="CO", state="g")
T1 = Temperature(value=298.15, unit="K")
T2 = Temperature(value=350.0, unit="K")

delta_H = thermo.calc_dEn(
    component=CO2,
    model_source=model_source,
    temperature_initial=T1,
    temperature_final=T2,
)

gibbs_range = thermo.calc_GiFrEn_range(
    component=CO2,
    model_source=model_source,
    temperatures=[T1, T2],
    phase="IG",
)

mix_H = thermo.calc_En_mix(
    components=[CO2, CO],
    model_source=model_source,
    temperature=T2,
    pressure=Pressure(value=1.0, unit="atm"),
    reference="IG",
    output_unit="kJ/mol",
)
```

## Saturation and phase-change (`pyThermoCalcDB.docs.sat`) 💧

| Function | What it does | Key inputs | Output (units) |
| --- | --- | --- | --- |
| `calc_VaPr` | Vapor pressure (Psat) at temperature | `component`, `model_source`, `temperature`, `component_key` | `CalcResult` in Pa |
| `calc_VaPr_range` | Vapor pressure across temperatures | `component`, `model_source`, `temperatures`, `component_key` | `list[CalcResult]` |
| `calc_EnVap` | Enthalpy of vaporization via Clapeyron using dPsat/dT | `component`, `model_source`, `temperature`, `component_key` | `CalcResult` (J/mol) |
| `calc_EnVap_range` | Enthalpy of vaporization across temperatures | `component`, `model_source`, `temperatures`, `component_key` | `list[CalcResult]` |
| `calc_T_sat` | Saturated temperature at a target pressure (root finding) | `component`, `model_source`, `pressure`, optional `temperature_guess`, `T_bracket`, `method` (`auto`, `brentq`, `bisect`, `newton`, `least_squares`), `tol`, `max_iter`, `h`, `component_key` | `CalcResult` in K |
| `calc_VaPr_sensitivity` | dPsat/dT at temperature | `component`, `model_source`, `temperature`, `component_key` | `CalcResult` |
| `calc_VaPr_sensitivity_range` | dPsat/dT across temperatures | `component`, `model_source`, `temperatures`, `component_key` | `list[CalcResult]` |

Notes 💡:
- Temperatures should be provided in Kelvin and pressures in Pascal; conversion is handled internally by `pycuc`.
- The model source must contain a vapor-pressure equation (`VaPr`) for the component.
- `calc_T_sat` accepts either a temperature guess or a bracket; `least_squares` is helpful when a bracket is not obvious.

Example 🧪:

```python
from pyThermoCalcDB.docs import sat
from pythermodb_settings.models import Pressure, Temperature

T_grid = [Temperature(value=v, unit="K") for v in (280.0, 300.0, 320.0)]
Pvaps = sat.calc_VaPr_range(component=CO2, model_source=model_source, temperatures=T_grid)

Tsat = sat.calc_T_sat(
    component=CO2,
    model_source=model_source,
    pressure=Pressure(value=1.5e6, unit="Pa"),
    temperature_guess=Temperature(value=290.0, unit="K"),
    method="least_squares",
)
```

## Reaction energetics (`pyThermoCalcDB.reactions.reactions`) 🔁

| Function | What it does | Key inputs | Output (units) |
| --- | --- | --- | --- |
| `dEn_rxn_STD` | Standard enthalpy of reaction at a temperature | `reaction` (`pyreactlab_core.models.reaction.Reaction`), `temperature`, `model_source` | `CustomProp` in J/mol |
| `dGiFrEn_rxn_STD` | Standard Gibbs free energy of reaction at a temperature | `reaction`, `temperature`, `model_source` | `CustomProp` in J/mol |

Example 🔥:

```python
from pyThermoCalcDB.reactions.reactions import dEn_rxn_STD
from pyreactlab_core.models.reaction import Reaction
from pythermodb_settings.models import Component, Temperature

# model_source built via pyThermoLinkDB
reaction = Reaction(
    name="Combustion of Methane",
    reaction="CH4(g) + 2O2(g) => CO2(g) + 2H2O(g)",
    components=[
        Component(name="methane", formula="CH4", state="g"),
        Component(name="oxygen", formula="O2", state="g"),
        Component(name="carbon dioxide", formula="CO2", state="g"),
        Component(name="water", formula="H2O", state="g"),
    ],
)

dH_rxn = dEn_rxn_STD(
    reaction=reaction,
    temperature=Temperature(value=298.15, unit="K"),
    model_source=model_source,
)
```
