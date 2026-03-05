# Flagellates

**FABM module name:** `CalculateFlagellates`  
**Source:** `fortran/phytoplankton/flagellates.F90`

## Purpose

> TODO (decorate): one or two sentences describing what this block represents in BAMHBI.

### State variables used / updated

**Core:**
- `CFL` (`cfl`, mmol C m-3): Small flagellate biomass in carbon
- `NFL` (`nfl`, mmol N m-3): Large flagellate biomass in nitrogen

**Other tracers referenced:**
- `CEM` (`cem`, mmol C m-3): Small flagellate biomass in carbon
- `DCL` (`dcl`, mmol C m-3): Labile detritus concentration in carbon
- `DCS` (`dcs`, mmol C m-3): Semi-labile detritus concentration in carbon
- `DIC` (`dic`, mmol C m-3): Dissolved inorganic carbon concentration
- `DNL` (`dnl`, mmol N m-3): Labile detritus concentration in nitrogen
- `DNS` (`dns`, mmol C m-3): Semi-labile detritus concentration in nitrogen
- `DOX` (`dox`, mmol O2 m-3): Dissolved oxygen concentration
- `NEM` (`nem`, mmol N m-3): Large flagellate biomass in nitrogen
- `NHS` (`nhs`, mmol N m-3): Ammonium concentration
- `NOS` (`nos`, mmol N m-3): Nitrate concentration
- `PHO` (`pho`, mmol P m-3): Phosphorus
- `POC` (`poc`, mmol C m-3): Particulate organic carbon concentration
- `PON` (`pon`, mmol N m-3): Particulate organic nitrogen concentration

(read from state/tracer arrays in the host model)

### Internal rates (names as in `Databases/Modules.txt`)

- `Excretion_ext`
- `Growth`
- `Leakage_DOC`
- `Leakage_DON`
- `Limitation_light`
- `Limitation_nutrient`
- `Mortality`
- `Mortality_C`
- `Mortality_N`
- `Ratio_Chl_C`
- `Ratio_N_C`
- `Respiration_tot`
- `Uptake_C`
- `Uptake_NHS`
- `Uptake_NO3`
- `Uptake_Nit`
- `Uptake_PO4`
- `Uptake_nutrient`
- `tf`

### Diagnostics exposed

- `Carbon_UptakeFlagellates`
- `Nitrogen_Uptake_Flagellates`
- `NPP`
- `PhytoNitrateReduction`
- `TotalRespirationFlagellates`

### Parameters consumed

- `dr_fla` (default `5.8e-06`, mol d-1): Deposition rate of FL
- `exc_extra_doc` (default `0.05`, -): Extra-photosynthetic DOC excretion (Van der Molen et al, 2004)
- `f_dl_dom` (default `0.7`, -): Labile fraction of PHY- and nonPHY-produced DOM (A&P, 2003)
- `f_dl_phy_ex` (default `0.65`, -): Labile fraction phytoxcreted DOC (Anderson and Pondhaven, 2003)
- `f_dl_phy_mo` (default `0.34`, -): DOM fraction of phytoplankton mortality
- `f_leak_phy` (default `0.02`, -): Phytoplankton leakage fraction (Van der Molen et al., 2004)
- `f_pp_resp_fla` (default `0.1`, -): Part of primary production used for respiration by FL
- `k_d` (default `0.03`, m-1): Background light attanuation coefficient
- `ki_nhs_phy` (default `0.5`, mmolN m-3): Inhib. constant of NHS for NOS uptake by PHY (Soetaert et al., 2001)
- `ks_nhs_fla` (default `3.0`, mmolN m-3): Half-saturation constant for NHS uptake by FL
- `ks_nos_fla` (default `3.0`, mmolN m-3): Half-saturation constant for NOS uptake by FL
- `ks_po4_fla` (default `0.2`, mmolP m-3): Half-saturation constant for PO4 uptake by FL
- `mo_fla` (default `0.03`, d-1): Mortality rate of FL (Asknes et al., 1994)
- `mumax_fla` (default `1.0`, d-1): Maximum specific growth rate of FL
- `pi_fla` (default `0.2153`, m2 W-1 d-1): Initial slope of photosynthesis-light curve for FL
- `q10_phy` (default `2.0`, -): Temperature factor (Soetaert et al., 2001)
- `r_o2_c_resp` (default `1.0`, molO2 molC-1): O2:C ratio of respiration process
- `r_o2_nhs_nitr` (default `2.0`, molO2 molNS-1): O2:NHS ratio in NHS oxidation in nitrification
- `r_p_n_redfield` (default `0.0625`, molP molN-1): N:P Redfield ratio in PHY
- `respb_fla` (default `0.009`, d-1): Basal respiration rate of FL (Lancelot et al., 2002)
- `rmax_chl_n_fla` (default `2.0`, g Chla molN-1): Maximum Chl:N ratio in FL (Soetaert et al., 2001)
- `rmax_n_c_fla` (default `0.2`, mol N molC-1): Maximum N:C ratio in FL (Soetaert et al., 2001)
- `rmin_chl_n_fla` (default `1.0`, g Chla molN-1): Minimum Chl:N ratio in FL (Soetaert et al., 2001)
- `rmin_n_c_fla` (default `0.05`, molN molC-1): Minimum N:C ratio in FL (Soetaert et al., 2001)
- `umax_nhs_fla` (default `0.5`, molN molC-1 d-1): Maximal NHS uptake rate by FL
- `umax_nos_fla` (default `0.50`, molN molC-1 d-1): Maximal NOS uptake rate by FL
- `umax_po4_fla` (default `0.03125`, molP molC-1 d-1): Maximal PO4 uptake rate by FL
- `w_fla` (default `-2.0`, m d-1): Sinking velocity of FL

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
