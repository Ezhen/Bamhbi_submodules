# Microzooplankton

**FABM module name:** `CalculateMicroZoo`  
**Source:** `fortran/zooplankton/microzoo.F90`

## Purpose

> TODO (decorate): one or two sentences describing what this block represents in BAMHBI.

### State variables used / updated

**Core:**
- `MIC` (`mic`, mmol C m-3): Microzooplakton biomass

**Other tracers referenced:**
- `BAC` (`bac`, mmol C m-3): Bacterial biomass
- `CDI` (`cdi`, mmol C m-3): Diatom biomass in carbon
- `CEM` (`cem`, mmol C m-3): Small flagellate biomass in carbon
- `CFL` (`cfl`, mmol C m-3): Small flagellate biomass in carbon
- `DCL` (`dcl`, mmol C m-3): Labile detritus concentration in carbon
- `DCS` (`dcs`, mmol C m-3): Semi-labile detritus concentration in carbon
- `DIC` (`dic`, mmol C m-3): Dissolved inorganic carbon concentration
- `DNL` (`dnl`, mmol N m-3): Labile detritus concentration in nitrogen
- `DNS` (`dns`, mmol C m-3): Semi-labile detritus concentration in nitrogen
- `DOX` (`dox`, mmol O2 m-3): Dissolved oxygen concentration
- `NDI` (`ndi`, mmol N m-3): Diatom biomass in nitrogen
- `NEM` (`nem`, mmol N m-3): Large flagellate biomass in nitrogen
- `NFL` (`nfl`, mmol N m-3): Large flagellate biomass in nitrogen
- `NHS` (`nhs`, mmol N m-3): Ammonium concentration
- `PHO` (`pho`, mmol P m-3): Phosphorus
- `POC` (`poc`, mmol C m-3): Particulate organic carbon concentration
- `PON` (`pon`, mmol N m-3): Particulate organic nitrogen concentration
- `SID` (`sid`, mmol Si m-3): Detrital silicate concentration

(read from state/tracer arrays in the host model)

### Internal rates (names as in `Databases/Modules.txt`)

- `Egestion_C`
- `Egestion_N`
- `Excretion`
- `Grazing_C`
- `Grazing_N`
- `Growth`
- `Intake_C`
- `Intake_N`
- `Messy_feeding_C`
- `Messy_feeding_N`
- `Mortality_C`
- `Mortality_N`
- `Prey_C`
- `Prey_N`
- `Ratio_N_C`
- `Ratio_N_C_thr`
- `Respiration_C`
- `tf`

### Diagnostics exposed

- `bac_to_ZOO`
- `phy_to_ZOO`
- `POC_to_ZOO`

### Parameters consumed

- `doxsatmort` (default `7.8125`, mmolO2 m-3 (?)): Perc. of sat. where metabolic respiration is 1/2 the one under O2 sat.
- `eff_ass_zoo_c` (default `0.64`, -): ZOO assimilation efficiency on C (Anderson and Pondhaven, 2003)
- `eff_ass_zoo_n` (default `0.77`, -): ZOO assimilation efficiencies on N (Anderson and Pondhaven, 2003)
- `eff_gr_mic_c` (default `0.8`, -): MIC net growth efficiency on C (Anderson and Pondhaven, 2003)
- `eff_mic_bac` (default `0.7`, -): Capture efficiency of MIC on BAC
- `eff_mic_dia` (default `0.0`, -): Capture efficiency of MIC on DI
- `eff_mic_emi` (default `1.0`, -): Capture efficiency of MIC on EM
- `eff_mic_fla` (default `0.0`, -): Capture efficiency of MIC on FL
- `eff_mic_mes` (default `0.0`, -): Capture efficiency of MIC on MES
- `eff_mic_mic` (default `0.0`, -): Capture efficiency of MIC on MIC
- `eff_mic_pom` (default `0.0`, -): Capture efficiency of MIC on POM
- `f_dl_dom` (default `0.7`, -): Labile fraction of PHY- and nonPHY-produced DOM (A&P, 2003)
- `gmax_mic` (default `3.6`, d-1): Maximum grazing rate of MIC (Strom and Morello, 1998)
- `ks_mort_mic` (default `1.0`, mmolC m-3): (?) Mortality half-saturation rate of MIC
- `ks_prey_mic` (default `5.0`, mmolC m-3): Half-saturation constant for MIC grazing (Soetart et al., 2001)
- `mess_prey_mic` (default `0.23`, -): Messy feeding fraction of MIC grazing (Anderson and Pondhaven, 2003)
- `mo_anox_pred` (default `0.25`, d-1): Mortality rate in anoxia
- `moexp_mic` (default `2.0`, -): Order of the non-linearity of mortality rate for MIC
- `momax_mic` (default `0.3`, d-1): Maximum mortality rate of MIC (Anderson and William)
- `q10_zoo` (default `2.0`, -): Temperature factor Soetart et al., 2001
- `r_n_c_bac` (default `0.196`, molN molC-1): N:C (Goldman) ratio in BAC (Anderson and Pondhaven, 2003)
- `r_n_c_mic` (default `0.18`, molN molC-1): N:C molar ratio in MIC (Anderson and Pondhaven, 2003)
- `r_o2_c_resp` (default `1.0`, molO2 molC-1): O2:C ratio of respiration process
- `r_p_n_redfield` (default `0.0625`, molP molN-1): N:P Redfield ratio in PHY
- `r_si_n_dia` (default `0.83`, molSi molN-1): Si:N ratio in DI (Aksnes et al. 1994)
- `w_dia` (default `-1.0`, m d-1): Sinking velocity of DI

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
