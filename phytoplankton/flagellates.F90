#include "fabm_driver.h" 
 
!#########################################################################################
!
! 1-D ecosystem model - Biological model of Tett
!
! References to the microplankton model:
! Tett, P., 1998. Parameterising a microplankton model.
! Department of Biological Sciences, Napier University,
! Report ISBN 0 902703 60 9, 60 pp.
!
! Sharples Tett (1994). Modeling the effect of physical! variability on the midwater chlorophyll maximum.
! Journal of marine research 52: 219-238
!
! Implementation: Marilaure Gregoire,                 NIOO-CEME
! Translation into FABM : Evgeny Ivanov, ULg / MAST

! Contains the pelagic submodel, as used in Soetaert et al 2001.                                                                                     !
!
!--------------------------------------------------------------------*

   module fabm_ulg_flagellates 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_flagellates 
      type (type_state_variable_id)         :: id_cfl,id_nfl
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_uptake_c_fla,id_uptake_n_fla,id_npp,id_reduction_nitrate_phy,id_respiration_fla
      type (type_diagnostic_variable_id)    :: id_chla

!     Model parameters 
      real(rk)     :: dr_fla, exc_extra_doc, f_dl_dom, f_dl_phy_ex
      real(rk)     :: f_dl_phy_mo, f_leak_phy, f_pp_resp_fla, k_d
      real(rk)     :: ki_nhs_phy, ks_nhs_fla, ks_nos_fla, ks_po4_fla
      real(rk)     :: mo_fla, mumax_fla, pi_fla, q10_phy, r_o2_c_resp
      real(rk)     :: r_o2_nhs_nitr, r_p_n_redfield, respb_fla
      real(rk)     :: rmax_chl_n_fla, rmax_n_c_fla, rmin_chl_n_fla
      real(rk)     :: rmin_n_c_fla, umax_nhs_fla, umax_nos_fla
      real(rk)     :: umax_po4_fla, w_fla

      contains 

      procedure :: initialize 
      procedure :: do 
      procedure :: do_bottom
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

   contains
   ! Initialise the Flagellates model

   subroutine initialize(self,configunit)
   class (type_ulg_flagellates), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

!     Model parameters 
      real(rk)     :: dr_fla, exc_extra_doc, f_dl_dom, f_dl_phy_ex
      real(rk)     :: f_dl_phy_mo, f_leak_phy, f_pp_resp_fla, k_d
      real(rk)     :: ki_nhs_phy, ks_nhs_fla, ks_nos_fla, ks_po4_fla
      real(rk)     :: mo_fla, mumax_fla, pi_fla, q10_phy, r_o2_c_resp
      real(rk)     :: r_o2_nhs_nitr, r_p_n_redfield, respb_fla
      real(rk)     :: rmax_chl_n_fla, rmax_n_c_fla, rmin_chl_n_fla
      real(rk)     :: rmin_n_c_fla, umax_nhs_fla, umax_nos_fla
      real(rk)     :: umax_po4_fla, w_fla

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%dr_fla, 'dr_fla', 'mol d-1', 'Deposition rate of FL', default=5.8e-06_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%exc_extra_doc, 'exc_extra_doc', '-', 'Extra-photosynthetic DOC excretion', default=0.05_rk) 
   call self%get_parameter(self%f_dl_dom, 'f_dl_dom', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%f_dl_phy_ex, 'f_dl_phy_ex', '-', 'Labile fraction phytoxcreted DOC', default=0.65_rk) 
   call self%get_parameter(self%f_dl_phy_mo, 'f_dl_phy_mo', '-', 'DOM fraction of phytoplankton mortality', default=0.34_rk) 
   call self%get_parameter(self%f_leak_phy, 'f_leak_phy', '-', 'Phytoplankton leakage fraction', default=0.02_rk) 
   call self%get_parameter(self%f_pp_resp_fla, 'f_pp_resp_fla', '-', 'Part of primary production used for respiration by FL', default=0.1_rk) 
   call self%get_parameter(self%k_d, 'k_d', 'm-1', 'Background light attanuation coefficient', default=0.03_rk) 
   call self%get_parameter(self%ki_nhs_phy, 'ki_nhs_phy', 'mmolN m-3', 'Inhib. constant of NHS for NOS uptake by PHY', default=0.5_rk) 
   call self%get_parameter(self%ks_nhs_fla, 'ks_nhs_fla', 'mmolN m-3', 'Half-saturation constant for NHS uptake by FL', default=3.0_rk) 
   call self%get_parameter(self%ks_nos_fla, 'ks_nos_fla', 'mmolN m-3', 'Half-saturation constant for NOS uptake by FL', default=3.0_rk) 
   call self%get_parameter(self%ks_po4_fla, 'ks_po4_fla', 'mmolP m-3', 'Half-saturation constant for PO4 uptake by FL', default=0.2_rk) 
   call self%get_parameter(self%mo_fla, 'mo_fla', 'd-1', 'Mortality rate of FL', default=0.03_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%mumax_fla, 'mumax_fla', 'd-1', 'Maximum specific growth rate of FL', default=1.0_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%pi_fla, 'pi_fla', 'm2 W-1 d-1', 'Initial slope of photosynthesis-light curve for FL', default=0.2153_rk, scale_factor=one_pr_day) 
   call self%get_parameter(self%q10_phy, 'q10_phy', '-', 'Temperature factor', default=2.0_rk) 
   call self%get_parameter(self%r_o2_c_resp, 'r_o2_c_resp', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%r_o2_nhs_nitr, 'r_o2_nhs_nitr', 'molO2 molNS-1', 'O2:NHS ratio in NHS oxidation in nitrification', default=2.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%respb_fla, 'respb_fla', 'd-1', 'Basal respiration rate of FL', default=0.009_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rmax_chl_n_fla, 'rmax_chl_n_fla', 'g Chla molN-1', 'Maximum Chl:N ratio in FL', default=2.0_rk) 
   call self%get_parameter(self%rmax_n_c_fla, 'rmax_n_c_fla', 'mol N molC-1', 'Maximum N:C ratio in FL', default=0.2_rk) 
   call self%get_parameter(self%rmin_chl_n_fla, 'rmin_chl_n_fla', 'g Chla molN-1', 'Minimum Chl:N ratio in FL', default=1.0_rk) 
   call self%get_parameter(self%rmin_n_c_fla, 'rmin_n_c_fla', 'molN molC-1', 'Minimum N:C ratio in FL', default=0.05_rk) 
   call self%get_parameter(self%umax_nhs_fla, 'umax_nhs_fla', 'molN molC-1 d-1', 'Maximal NHS uptake rate by FL', default=0.5_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%umax_nos_fla, 'umax_nos_fla', 'molN molC-1 d-1', 'Maximal NOS uptake rate by FL', default=0.50_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%umax_po4_fla, 'umax_po4_fla', 'molP molC-1 d-1', 'Maximal PO4 uptake rate by FL', default=0.03125_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%w_fla, 'w_fla', 'm d-1', 'Sinking velocity of FL', default=-2.0_rk, scale_factor=one_pr_day) 

   ! Register state variables 

   call self%register_state_variable(self%id_cfl, 'CFL', 'mmol C m-3', 'Small flagellate biomass in carbon', minimum=0.0e-7_rk, vertical_movement=self%w_fla) 
   call self%register_state_variable(self%id_nfl, 'NFL', 'mmol N m-3', 'Large flagellate biomass in nitrogen', minimum=0.0e-7_rk, vertical_movement=self%w_fla) 

   call self%register_state_dependency(self%id_dcl, 'DCL', 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'DCS', 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'DIC', 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'DNL', 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'DNS', 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'DOX', 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_nhs, 'NHS', 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nos, 'NOS', 'Nitrate concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'PHO', 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'POC', 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'PON', 'Particulate organic nitrogen concentration', 'mmol N m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 

    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_uptake_c_fla, 'uptake_c_fla', 'mmol C m-3 d-1', & 
      'Carbon uptake by Flagellates', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_uptake_n_fla, 'uptake_n_fla', 'mmol N m-3 d-1', & 
      'Nitrogen uptake by Flagellates', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_npp, 'npp', 'mmol N m-3 d-1', & 
      ' Primary production of nitrogen by all types of phytoplankton', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_reduction_nitrate_phy, 'reduction_nitrate_phy', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_respiration_fla, 'respiration_fla', 'mmol C m-3 d-1', & 
      'Total Respiration of Flagellates', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_chla, 'chla','mg chl a m-3', 'Chlorophyll concentration')

    ! Add to aggregate variables 
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_cfl)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nfl)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_nfl, scale_factor=self%r_p_n_redfield) 
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_chlorophyll',units="mg chl a m-3",aggregate_variable=.true.),self%id_chla,scale_factor=1._rk)

   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_nfl, scale_factor=self%k_d)

   return 

99 call self%fatal_error('Flagellates', 'Error reading namelist ulg_flagellates') 

   end subroutine initialize 


   ! Right hand sides of Flagellates model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_flagellates), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CFL,NFL
      real(rk) ::   Excretion_ext	  ! mmol C d-1, Phytoplankton extra excretion
      real(rk) ::   Growth	  ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   Leakage_DOC	  ! mmol C d-1, Phytoplankton passive leakage rate for carbon
      real(rk) ::   Leakage_DON	  ! mmol N d-1, Phytoplankton passive leakage rate for nitrogen
      real(rk) ::   Limitation_light	  ! -, Light limitation for small flagellates
      real(rk) ::   Limitation_nutrient	  ! -, Nutrient limitation for small flagellates
      real(rk) ::   Mortality	  ! mmol m-3, Phytoplankton mortality rate
      real(rk) ::   Mortality_C	  ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Mortality_N	  ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   Ratio_Chl_C	  ! g Chla mol C-1, Chl/C ratio in small flagellates
      real(rk) ::   Ratio_N_C	  ! mol N mol C-1, N/C ratio in small flagellates
      real(rk) ::   Ratio_Si_C	  ! g Chla mol C-1, Chl/C ratio in small flagellates
      real(rk) ::   Respiration_tot	  ! mmol C m-3, Total phytoplankton respiration (basal & activity)
      real(rk) ::   Uptake_C	  ! mmol C m-3, C assimilation of phytoplankton
      real(rk) ::   Uptake_NHS	  ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   Uptake_NO3	  ! mmol N m-3, Nitrate uptake of phytoplankton
      real(rk) ::   Uptake_Nit	  ! mmol N m-3, Nitrogen uptake of phytoplankton
      real(rk) ::   Uptake_PO4	  ! mmol P m-3, Phosphate uptake by large flagellates
      real(rk) ::   Uptake_nutrient	  ! mmol m-3, Nutrient uptake of phytoplankton
      real(rk) ::   tf	  ! -, Temperature factor

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_cfl,CFL)       ! Large flagellate biomass in carbon
   _GET_(self%id_nfl,NFL)       ! Large flagellate biomass in nitrogen
   _GET_(self%id_dcl,DCL)       ! Labile detritus concentration in carbon
   _GET_(self%id_dcs,DCS)       ! Semi-labile detritus concentration in carbon
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dnl,DNL)       ! Labile detritus concentration in nitrogen
   _GET_(self%id_dns,DNS)       ! Semi-labile detritus concentration in nitrogen
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_nos,NOS)       ! Nitrate concentration
   _GET_(self%id_pho,PHO)       ! Phosphorus
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor(temp,self%q10_phy)
    
   ! Calculate ratios in phytoplankton 
    Ratio_N_C  = Ratio(NFL,CFL)
    Ratio_Chl_C = ratio_chl_c_phyt(Ratio_N_C,self%rmax_n_c_fla,self%rmin_n_c_fla,self%rmin_chl_n_fla,self%rmax_chl_n_fla)
    
   ! Nitrate uptake rate 
    Uptake_NO3 = uptake_nitrate_phyt(Ratio_N_C,tf,self%rmax_n_c_fla,self%umax_nos_fla) * Michaelis(NOS,self%ks_nos_fla) * Inhibition(NHS,self%ki_nhs_phy) * CFL
    
   ! Ammonium uptake rate 
    Uptake_NHS = uptake_nutrient_phyt(Ratio_N_C,(NHS-0.0_rk),tf,self%rmax_n_c_fla,self%umax_nhs_fla,self%ks_nhs_fla) * CFL
    
   ! Phosphate uptake rate 
    Uptake_PO4 = uptake_nutrient_phyt(Ratio_N_C,(PHO-0.0_rk),tf,self%rmax_n_c_fla,self%umax_po4_fla,self%ks_po4_fla) * CFL
    
   ! Potential nitrogen uptake 
    Uptake_Nit = Uptake_NHS + Uptake_NO3
    
   ! Nutrient uptake 
    Uptake_nutrient = min(Uptake_Nit,Uptake_PO4/self%r_p_n_redfield)
                                                                                                                                                                                               
   ! Compute actual N:C and Si:C ratios 
    Ratio_N_C = ratio_adjustment(Ratio_N_C,self%rmax_n_c_fla,self%rmin_n_c_fla)
    Ratio_Si_C = ratio_adjustment(1.0_rk,1.0_rk,0.0_rk)
    
   ! Compute nutrient and light limitation 
    Limitation_nutrient = limitation_by_nutrient(self%rmin_n_c_fla,0.0_rk,Ratio_N_C,Ratio_Si_C)
    Limitation_light = 1.0-exp(-self%pi_fla * par / self%mumax_fla)
    
   ! Compute carbon uptake 
    Uptake_C = self%mumax_fla  * Limitation_nutrient * CFL * tf * Limitation_light
    
   ! Compute respiration 
    Respiration_tot = Uptake_C * self%f_pp_resp_fla + self%respb_fla * CFL * tf
    
   ! Compute growth 
    Growth = Uptake_C - Respiration_tot
    
   ! Compute the extra DOC excretion from the difference of growth in nutrient limited (actual NC ratio) and nutrient saturated (max NC ratio) 
    Excretion_ext = CFL * tf * self%exc_extra_doc * self%mumax_fla * extra_excretion(Limitation_light,self%rmin_n_c_fla,self%rmax_n_c_fla,Ratio_N_C)
    
   ! Compute the leakage 
    Leakage_DOC = self%f_leak_phy * Uptake_C
    Leakage_DON = self%f_leak_phy * abs(Uptake_nutrient)
    
   ! Compute mortality 
    Mortality_C  = mortality_phyt(self%mo_fla,tf) * CFL
    Mortality_N  = mortality_phyt(self%mo_fla,tf) * NFL
    
   ! Carbon/nitrogen in flagellates increases by growth and decreases by leakage and mortality 
   _ADD_SOURCE_(self%id_cfl, Growth - Mortality_C - Leakage_DOC) 
   _ADD_SOURCE_(self%id_nfl, Uptake_nutrient - Mortality_N - Leakage_DON) 
    
   ! IF CN ratio of phytoplankton higher than CNmin, than nitrogen is taken up, unless it gets excreted 
   IF (Uptake_nutrient.gt.0) THEN 
     _ADD_SOURCE_(self%id_nos, -Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit))
     _ADD_SOURCE_(self%id_dox, Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit) * self%r_o2_nhs_nitr)
     _ADD_SOURCE_(self%id_nhs, -Uptake_nutrient * Ratio(Uptake_NHS,Uptake_Nit))
     _ADD_SOURCE_(self%id_pho, -Uptake_nutrient * self%r_p_n_redfield)
   ELSE 
     _ADD_SOURCE_(self%id_nhs, -Uptake_nutrient)
     _ADD_SOURCE_(self%id_pho, -Uptake_nutrient * self%r_p_n_redfield)
   ENDIF 
    
   ! Mortality increases the pool of POM and DOM with proper partitioning and leakage adds to the labile pool, extra-excretion and leakage add to the labile and semi-labile detritus 
   _ADD_SOURCE_(self%id_poc, (1.0 - self%f_dl_phy_mo) * Mortality_C)
   _ADD_SOURCE_(self%id_pon, (1.0 - self%f_dl_phy_mo) * Mortality_N)
   _ADD_SOURCE_(self%id_dcl, self%f_dl_phy_mo * Mortality_C * self%f_dl_dom + Leakage_DOC + self%f_leak_phy * Excretion_ext + self%f_dl_phy_ex * (1.0 - self%f_leak_phy) * Excretion_ext)
   _ADD_SOURCE_(self%id_dnl, self%f_dl_phy_mo * Mortality_N * self%f_dl_dom + Leakage_DON) 
   _ADD_SOURCE_(self%id_dcs, self%f_dl_phy_mo * Mortality_C * (1.0 - self%f_dl_dom) + (1.0 - self%f_dl_phy_ex) * (1.0 - self%f_leak_phy) * Excretion_ext)
   _ADD_SOURCE_(self%id_dns, self%f_dl_phy_mo * Mortality_N * (1.0 - self%f_dl_dom))
   _ADD_SOURCE_(self%id_dox, (Growth + Excretion_ext) * self%r_o2_c_resp)
   _ADD_SOURCE_(self%id_dic, -Growth - Excretion_ext)

   _SET_DIAGNOSTIC_(self%id_uptake_n_fla, Uptake_NHS + Uptake_NO3)
   _SET_DIAGNOSTIC_(self%id_uptake_c_fla, Uptake_C)
   _SET_DIAGNOSTIC_(self%id_respiration_fla, Respiration_tot)
   _SET_DIAGNOSTIC_(self%id_reduction_nitrate_phy, Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit) * self%r_o2_nhs_nitr)
   _SET_DIAGNOSTIC_(self%id_npp, Growth)
   _SET_DIAGNOSTIC_(self%id_chla, Ratio_Chl_C*CFL)

   _LOOP_END_

   end subroutine do


   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_ulg_flagellates), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk),save :: nfl, cfl

   _HORIZONTAL_LOOP_BEGIN_
 
   _GET_(self%id_nfl,nfl) 
   _GET_(self%id_cfl,cfl)    

   _ADD_BOTTOM_FLUX_(self%id_nfl,-self%dr_fla*nfl*nfl)
   _ADD_BOTTOM_FLUX_(self%id_cfl,-self%dr_fla*cfl*cfl)   

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

   end module fabm_ulg_flagellates 
