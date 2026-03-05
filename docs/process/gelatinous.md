# Gelatinous

**FABM module name:** `CalculateGelatinous`  
**Source:** `fortran/jellyfish/gelatinous.F90`

## Purpose

> TODO (decorate): one or two sentences describing what this block represents in BAMHBI.

### State variables used / updated

**Core:**
- `GEL` (`gel`, mmol C m-3): Gelatinous omnivorous biomass

**Other tracers referenced:**
- `CDI` (`cdi`, mmol C m-3): Diatom biomass in carbon
- `CEM` (`cem`, mmol C m-3): Small flagellate biomass in carbon
- `CFL` (`cfl`, mmol C m-3): Small flagellate biomass in carbon
- `DIC` (`dic`, mmol C m-3): Dissolved inorganic carbon concentration
- `DOX` (`dox`, mmol O2 m-3): Dissolved oxygen concentration
- `MES` (`mes`, mmol C m-3): Mesozooplakton biomass
- `MIC` (`mic`, mmol C m-3): Microzooplakton biomass
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
- `Excretion_C_adj`
- `Excretion_N_adj`
- `Grazing_C`
- `Grazing_N`
- `Mortality_C`
- `Mortality_N`
- `Prey_C`
- `Ratio_N_C`
- `Ratio_N_C_test`
- `Respiration_C`
- `tf`

### Diagnostics exposed

- `TotalRespiration_Gel`

### Parameters consumed

- `doxsatmort` (default `7.8125`, mmolO2 m-3 (?)): Perc. of sat. where metabolic respiration is 1/2 the one under O2 sat.
- `eff_ass_gel_prey` (default `0.75`, -): GEL assimilation efficiency on prey (Lancelot et al., 2002)
- `eff_gel_dia` (default `0.0`, -): Capture efficiency of GEL on DI
- `eff_gel_emi` (default `0.0`, -): Capture efficiency of GEL on EM
- `eff_gel_fla` (default `0.0`, -): Capture efficiency of GEL on FL
- `eff_gel_mes` (default `1.0`, -): Capture efficiency of GEL on MES
- `eff_gel_mic` (default `0.0`, -): Capture efficiency of GEL on MIC
- `eff_gel_pom` (default `0.0`, -): Capture efficiency of GEL on POM
- `eff_gr_gel_c` (default `0.2`, -): Part of the assimil. food used for GEL growth (Lancelot et al., 2002)
- `gmax_gel` (default `0.3`, d-1): Maximum grazing rate of GEL
- `ks_mort_gel` (default `0.0`, mmolC m-3): (?) Mortality half-saturation rate of GEL
- `mo_anox_pred` (default `0.25`, d-1): Mortality rate in anoxia
- `momax_gel` (default `0.009`, d-1): Maximum mortality rate of GEL (Lancelot et al., 2002)
- `q10_gel` (default `3.5`, -): Temperature dependency for GEL (Kremer, 1977)
- `r_n_c_gel` (default `0.25`, molN molC-1): N:C molar ratio in GEL
- `r_o2_c_resp` (default `1.0`, molO2 molC-1): O2:C ratio of respiration process
- `r_p_n_redfield` (default `0.0625`, molP molN-1): N:P Redfield ratio in PHY
- `r_si_n_dia` (default `0.83`, molSi molN-1): Si:N ratio in DI (Aksnes et al. 1994)
- `respb_gel` (default `0.0001`, d-1): Basal respiration rate of GEL
- `t_g_gel` (default `0`, mmolC m-3): Feeding threshold for GEL grazing (Lancelot et al., 2002)

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
