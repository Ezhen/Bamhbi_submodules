# DOM/POM

**FABM module name:** `CalculateDOM`  
**Source:** `fortran/dom/dom.F90`

## Purpose

> TODO (decorate): one or two sentences describing what this block represents in BAMHBI.

### State variables used / updated

- `DCL` (`dcl`, mmol C m-3): Labile detritus concentration in carbon
- `DCS` (`dcs`, mmol C m-3): Semi-labile detritus concentration in carbon
- `DNL` (`dnl`, mmol N m-3): Labile detritus concentration in nitrogen
- `DNS` (`dns`, mmol C m-3): Semi-labile detritus concentration in nitrogen
- `POC` (`poc`, mmol C m-3): Particulate organic carbon concentration
- `PON` (`pon`, mmol N m-3): Particulate organic nitrogen concentration

### External dependencies (read-only)

- `BAC` (`bac`, mmol C m-3): Bacterial biomass
- `DOX` (`dox`, mmol O2 m-3): Dissolved oxygen concentration

### Internal rates (names as in `Databases/Modules.txt`)

- `Hydrolysis_DNS_DNL`
- `Hydrolysis_DSC_DCL`
- `Hydrolysis_POC_DOC`
- `Hydrolysis_PON_DON`
- `Limitation_hydrolysis_POM`
- `tf`

### Diagnostics exposed

- (none declared in `Databases/Modules.txt`)

### Parameters consumed

- `dr_dom` (default `5.5e-06`, mol d-1): Deposition rate of dissolved detritus
- `dr_pom` (default `5.6e-06`, mol d-1): Deposition rate of particulate detritus
- `f_dl_dom` (default `0.7`, -): Labile fraction of PHY- and nonPHY-produced DOM (A&P, 2003)
- `hmax_dsl` (default `4.0`, d-1): Maximum DSL hydrolysis (Anderson and Pondhaven, 2003)
- `hmax_poc` (default `0.04`, d-1): POC hydrolysis rate (Anderson and Pondhaven, 2003)
- `hmax_pon` (default `0.055`, d-1): PON hydrolysis rate (Anderson and Pondhaven, 2003)
- `k_d` (default `0.03`, m-1): Background light attanuation coefficient
- `ks_dsc_bac` (default `417.0`, mmolC m-3): Half-sat. constant for DSC uptake by BAC (Anderson and Pondhaven, 2003)
- `ks_hydr_o2` (default `2.7`, mmolO2 m-3): Half-saturation constant for oxic hydrolysis rate
- `q10_che` (default `2.0`, -): Temperature factor for chemical processes
- `w_dom` (default `-1.0`, m d-1): Sinking velocity of dissolved detritus
- `w_pom` (default `-2.0`, m d-1): Sinking velocity of particulate detritus

### Routines

- `initialize`
- `do`
- `do_bottom`

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
