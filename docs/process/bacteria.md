# Bacteria

**FABM module name:** `CalculateBacteria`  
**Source:** `fortran/bacteria/bacteria.F90`

## Purpose

> TODO (decorate): one or two sentences describing what this block represents in BAMHBI.

### State variables used / updated

**Core:**
- `BAC` (`bac`, mmol C m-3): Bacterial biomass

**Other tracers referenced:**
- `DCL` (`dcl`, mmol C m-3): Labile detritus concentration in carbon
- `DCS` (`dcs`, mmol C m-3): Semi-labile detritus concentration in carbon
- `DIC` (`dic`, mmol C m-3): Dissolved inorganic carbon concentration
- `DNL` (`dnl`, mmol N m-3): Labile detritus concentration in nitrogen
- `DNS` (`dns`, mmol C m-3): Semi-labile detritus concentration in nitrogen
- `DOX` (`dox`, mmol O2 m-3): Dissolved oxygen concentration
- `NHS` (`nhs`, mmol N m-3): Ammonium concentration
- `NOS` (`nos`, mmol N m-3): Nitrate concentration
- `ODU` (`odu`, mmol ODU m-3): Oxygen demand unit concentration
- `PHO` (`pho`, mmol P m-3): Phosphorus

(read from state/tracer arrays in the host model)

### Internal rates (names as in `Databases/Modules.txt`)

- `Denitrificaiton`
- `Excretion`
- `Growth`
- `Iron`
- `Limitation_iron`
- `Limitation_nutrient`
- `Mortality_C`
- `Mortality_N`
- `Remineralization_anoxic_loc`
- `Respiration`
- `Respiration_loc`
- `Uptake_DCL`
- `Uptake_DNL`
- `Uptake_NHS`
- `Uptake_NHS_pot`
- `testratio`
- `tf`

### Diagnostics exposed

- `bacteria_anoxrem`
- `bacteria_oxygenconsumption`
- `Bacteria_Respiration`
- `denitrification`
- `Uptake_DOCL`

### Parameters consumed

- `eff_gr_bac_c` (default `0.17`, -): BAC gross growth efficiency on C (Anderson and Pondhaven, 2003)
- `f_dl_dom` (default `0.7`, -): Labile fraction of PHY- and nonPHY-produced DOM (A&P, 2003)
- `f_solid_odu` (default `0.2`, -): Percentage of solid ODU formation
- `i1_curve` (default `25000.0`, -): Parameter of the curve simulating the iron concentration
- `i2_curve` (default `50.0`, -): Parameter of the curve simulating the iron c
- `iron` (default `10.0`, mmolFe m-3): Concentration of iron in surface water
- `ki_anox_nos` (default `0.0005`, mmolN m-3): Half-sat. constant for NOS inhibition in anoxic remineralization
- `ki_anox_o2` (default `0.0005`, mmolO2 m-3): Half-sat. constant for O2 inhibition in anoxic remineralization
- `ki_denit_o2` (default `0.5`, mmolO2 m-3): Half-sat. constant for O2 inhibition in denitrification
- `ks_denitr_nos` (default `0.3`, mmolN m-3): Half-sat. constant for NOS lim. in denitrif. (Soetaert et al., 1996)
- `ks_dls_bac` (default `25.0`, mmolC m-3): Half-sat. constant for DLC uptake by BAC (Anderson and Pondhaven, 2003)
- `ks_nhs_bac` (default `0.5`, mmolN m-3): Half-sat. constant for NHS uptake by BAC (Anderson and Pondhaven, 2003)
- `ks_odu_iron` (default `100.0`, mmolFe m-3): Half-sat. constant for iron lim. in solid ODU formation
- `ks_oxic_o2` (default `3.0`, mmolO2 m-3): Half-sat. constant for O2 lim. in oxic min. (Soetaert et al., 1996)
- `ks_po4_bac` (default `0.031`, mmolP m-3): Half-saturation constant for PO4 uptake by BAC
- `mo_bac` (default `0.05`, d-1): Bacteria natural mortality (Anderson and Pondhaven, 2003)
- `mumax_bac` (default `0.000154`, d-1): Maximum labile DOC or NHS uptake by BAC (A&P, 2003)
- `q10_bac` (default `2.0`, -): Temperature factor for BAC
- `r_n_c_bac` (default `0.196`, molN molC-1): N:C (Goldman) ratio in BAC (Anderson and Pondhaven, 2003)
- `r_n_c_denit` (default `0.8`, molN molC-1): N:C ratio of denitrification
- `r_o2_c_resp` (default `1.0`, molO2 molC-1): O2:C ratio of respiration process
- `r_odu_c_anox` (default `1.0`, molODU molC-1): ODU:C ratio in anoxic remineralization
- `r_p_n_redfield` (default `0.0625`, molP molN-1): N:P Redfield ratio in PHY

### Routines

- `initialize`
- `do`

#### Routine hat template

```fortran
!-------------------------------------------------------------------------------
! <ROUTINE NAME> — <one-line purpose>
!
! Inputs (read-only)
! - <state/dependencies used> [units]
!
! Updates (write)
! - <which tracers / tendencies / diagnostics> [units]
!
! Notes
! - Conservation: <C/N/P/Si/O2/ODU; exceptions>
! - Numerics: <explicit/implicit; clipping/positivity; dt constraints>
! - Sign conventions: <flux positive upward/downward, etc.>
!-------------------------------------------------------------------------------
```
