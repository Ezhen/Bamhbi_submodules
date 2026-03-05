# Chemistry

**FABM module name:** `ChemicalDynamics`  
**Source:** `fortran/chem/chemical.F90`

## Purpose

> TODO (decorate): one or two sentences describing what this block represents in BAMHBI.

### State variables used / updated

- `DOX` (`dox`, mmol O2 m-3): Dissolved oxygen concentration
- `NHS` (`nhs`, mmol N m-3): Ammonium concentration
- `NOS` (`nos`, mmol N m-3): Nitrate concentration
- `ODU` (`odu`, mmol ODU m-3): Oxygen demand unit concentration

(read from state/tracer arrays in the host model)

### Internal rates (names as in `Databases/Modules.txt`)

- `Rate_nitrification`
- `Rate_oxidation_NHS_NOS`
- `Rate_oxidation_ODU_NO3`
- `Rate_oxidation_ODU_O2`
- `tf`

### Diagnostics exposed

- `ANAMMOX`
- `Nitrification`
- `Oxidation_by_nitrate`
- `Oxidation_by_oxygen`

### Parameters consumed

- `ki_nhs_o2` (default `8.0`, mmolO2 m-3): Half-sat. constant for O2 inhibition in NHS oxidation by NOS
- `ki_nhs_odu` (default `0.5`, mmolO2 m-3): Half-sat. constant for O2 inhibition in NHS oxidation by ODU
- `ki_odu_o2` (default `5.0`, mmolO2 m-3): Half-sat. constant for O2 inhibition in ODU oxidation by NOS
- `ks_nhs_o2` (default `3.0`, mmolO2 m-3): Half-sat. constant for O2 lim. in NHS oxidation by O2
- `ks_odu_nos` (default `2.0`, mmolN m-3): Half-sat. constant for NOS lim. in ODU oxidation by NOS
- `ks_odu_o2` (default `1.0`, mmolO2 m-3): Half-sat. constant for O2 lim. in ODU oxidation (Soetaert et al., 1996)
- `q10_che` (default `2.0`, -): Temperature factor for chemical processes
- `r_nos_nhs_oxid` (default `0.6`, molNOS molNHS-1): NOS:NHS ratio in NHS oxidation
- `r_nos_odu_oxid` (default `0.8`, molNOS molODU-1): NOS:ODU ratio in ODU oxidation
- `r_o2_nhs_nitr` (default `2.0`, molO2 molNS-1): O2:NHS ratio in NHS oxidation in nitrification
- `r_o2_odu_oxid` (default `1.0`, molO2 molODU-1): O2:ODU ratio in ODU oxidation
- `rox_nhs_nos` (default `0.05`, d-1): Maximum NHS oxidation rate by NOS
- `rox_nhs_o2` (default `0.03`, d-1): Maximum NHS oxidation rate of NHS by O2
- `rox_odu_nos` (default `0.05`, d-1): Maximum ODU oxidation rate by NOS
- `rox_odu_o2` (default `0.1`, d-1): Maximum ODU oxidation rate by O2 (Oguz et al. 2000)

### Routines

- `initialize`
- `do`
- `do_surface`

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
