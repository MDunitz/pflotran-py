# PFLOTRAN Module References

## Atmospheric CO₂(aq) Constraint

The atmospheric CO₂(aq) boundary condition uses a total concentration constraint (`T`)
equilibrated with CO₂(g).

### Henry's Law Calculation

```
[CO₂(aq)] = K_H × p_CO₂

K_H(T) = Henry's law constant for CO₂ in water [mol/(L·atm)]
p_CO₂  = atmospheric partial pressure of CO₂ [atm]
```

| Parameter | Value | Source |
|-----------|-------|--------|
| K_H (25°C) | 3.4 × 10⁻² mol/(L·atm) | Sander, R. (2023). Compilation of Henry's law constants (v5.0.0). *Atmos. Chem. Phys.*, 23, 10901–12440. https://doi.org/10.5194/acp-23-10901-2023 |
| K_H (20°C) | 3.7 × 10⁻² mol/(L·atm) | ibid. |
| K_H (8°C)  | ~5.3 × 10⁻² mol/(L·atm) | Interpolated from Weiss (1974) temperature dependence |
| p_CO₂ (atmosphere) | 4.2 × 10⁻⁴ atm (~420 ppm) | NOAA Global Monitoring Laboratory, 2024 |

**Expected [CO₂(aq)] at simulation temperature (8°C):**
```
[CO₂(aq)] = 5.3e-2 × 4.2e-4 ≈ 2.2e-5 M
```

**Value in generator:** `1.906e-05 M` — reasonable for slightly warmer conditions or
lower assumed pCO₂. Consistent with O'Meara et al. (2024) SWaMP constraint using
`G CO2(g)` gas equilibrium.

**BUG:** The generator wrote `'1.906-05 T'` (missing `e` in scientific notation).
Fixed to `'1.906e-05 T'`.

### Key Reference
- Weiss, R.F. (1974). Carbon dioxide in water and seawater: the solubility of a
  non-ideal gas. *Marine Chemistry*, 2(3), 203–215.
  https://doi.org/10.1016/0304-4203(74)90015-2

## Microbial Reaction Rate Constants

Rate constants are in PFLOTRAN units: `mol/(L·s)` for MICROBIAL_REACTION without biomass.
See PFLOTRAN docs: https://www.pflotran.org/documentation/user_guide/cards/subsurface/chemistry/microbial_reaction_card.html

### Fe(III) Reduction

| Source | DOM Pool | Rate Constant | Notes |
|--------|----------|--------------|-------|
| O'Meara et al. (2024) SWaMP — DOM1 | DOM1 | 3.0 × 10⁻¹⁰ | Calibrated to GCReW porewater data |
| O'Meara et al. (2024) SWaMP — DOM2 | DOM2 | 1.67 × 10⁻¹⁰ | Calibrated to GCReW porewater data |
| O'Meara SWaMP (commented out) | DOM1 | 5.0 × 10⁻⁹ | Initial estimate, Furukawa et al. (2004) |
| This generator | single pool | **2.25 × 10⁻¹⁰** | Consistent with O'Meara calibrated range |
| Main branch `.in` files | single pool | 2.25 × 10⁻⁸ | **100× higher** — likely from pre-calibration |

**Conclusion:** The generator value (2.25e-10) is consistent with O'Meara's calibrated
SWaMP network. The main branch value (2.25e-08) appears to be an earlier,
uncalibrated estimate. The generator value should be retained.

### Full Rate Constant Table

| Reaction | Generator Value | O'Meara SWaMP Reference | Literature Source |
|----------|----------------|------------------------|-------------------|
| DOM aerobic respiration | 1.80 × 10⁻⁷ | 6.3 × 10⁻⁹ (DOM1) | Thompson et al. (1995) |
| Fermentation | 6.00 × 10⁻⁸ | — | — |
| Fe(III) reduction | 2.25 × 10⁻¹⁰ | 3.0 × 10⁻¹⁰ (DOM1) | Furukawa et al. (2004), calibrated |
| Sulfate reduction | 1.50 × 10⁻⁹ | 1.5 × 10⁻¹³ (DOM1) | Furukawa et al. (2004), calibrated |
| Hydrogenotrophic methanogenesis | 7.20 × 10⁻⁹ | 9.33 × 10⁻¹⁶ (HCO₃⁻) | Thompson et al. (1995), calibrated |
| Acetoclastic methanogenesis | 1.50 × 10⁻⁸ | 8.83 × 10⁻¹⁶ (DOM1) | Thompson et al. (1995), calibrated |
| Methylotrophic methanogenesis | 9.10 × 10⁻⁶ | — | — |
| Nitrification | — | 3.3 × 10⁻¹⁰ | Dettmann (2001) |
| Denitrification | — | 1.5 × 10⁻⁹ (DOM1) | Dettmann (2001) |

**Note:** Generator rates are substantially higher than O'Meara's calibrated rates for most
reactions, particularly sulfate reduction (4 OOM higher) and methanogenesis (7+ OOM higher).
This likely reflects the different domain (salt pond vs. tidal marsh), higher DOM concentrations,
and the absence of tidal flushing in the generator's batch system.

### Half-Saturation Constants

| Species | Generator Value (M) | O'Meara SWaMP Reference (M) | Source |
|---------|-------------------|---------------------------|--------|
| Fe(OH)₃ | 1.0 × 10⁻¹⁰ | 4.7 × 10⁻⁴ | Gao et al. (2010) |
| SO₄²⁻ | 1.0 × 10⁻⁴ | 2.0 × 10⁻⁴ | Gao et al. (2010) |
| DOM1 (fermentation) | 5.0 × 10⁻² | — | — |
| O₂ | 1.0 × 10⁻⁴ | 1.0 × 10⁻⁵ | — |

## Salinity Inhibition

### Cl⁻ Inhibition Threshold
Generator uses 0.2 M Cl⁻ as a single threshold for all pathways (fermentation,
hydrogenotrophic, acetoclastic, methylotrophic methanogenesis).

Seawater Cl⁻ at 1× = 0.536 M → already above threshold at baseline.

### Water Activity (a_w) Inhibition
Custom Fortran reaction sandboxes (AWINHIBIT, AWINHIBITACETATE, AWINHIBITMETHYL)
apply a_w threshold of 0.5. Water activity captures the combined effect of all dissolved
ions (not just Cl⁻), making it a more physically complete inhibition mechanism.

**Potential double-counting concern:** See design note in PR #50 review.

## Literature Cited

- Furukawa, Y., Smith, A.C., Kostka, J.E., Watkins, J. & Alexander, C.R. (2004).
  Quantification of macrobenthic effects on diagenesis using a multicomponent
  inverse model in salt marsh sediments. *Limnol. Oceanogr.*, 49(6), 2058–2072.

- Gao, H., et al. (2010). Aerobic denitrification in permeable Wadden Sea sediments.
  *ISME Journal*, 4, 417–426.

- O'Meara, T.A., Yuan, F., Sulman, B.N., Noyce, G.L., Rich, R., Thornton, P.E.,
  & Megonigal, J.P. (2024). Developing a Redox Network for Coastal Saltmarsh Systems
  in the PFLOTRAN Reaction Model. *J. Geophys. Res.: Biogeosciences*, 129(3),
  e2023JG007633. https://doi.org/10.1029/2023JG007633
  - GitHub: https://github.com/omearata/REDOX-PFLOTRAN
  - Data: https://data.ess-dive.lbl.gov/datasets/doi:10.15485/2294096

- Sander, R. (2023). Compilation of Henry's law constants (version 5.0.0) for water as
  solvent. *Atmos. Chem. Phys.*, 23, 10901–12440.
  https://doi.org/10.5194/acp-23-10901-2023

- Sulman, B.N., et al. (2022). Integrating Tide-Driven Wetland Soil Redox and
  Biogeochemical Interactions Into a Land Surface Model. *J. Adv. Modeling Earth Systems*,
  14(4), e2021MS002916.

- Thompson, S.P., Paerl, H.W. & Go, M.C. (1995). Seasonal patterns of nitrification and
  denitrification in a natural and a restored salt marsh. *Estuaries*, 18(3), 399–408.

- Weiss, R.F. (1974). Carbon dioxide in water and seawater: the solubility of a non-ideal
  gas. *Marine Chemistry*, 2(3), 203–215.
