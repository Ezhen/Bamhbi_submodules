#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DemilianaF90
!
! 1-D ecosystem model - Biological model of Tett
!
! Tett, P., 1998. Parameterising a microplankton model.
! Department of Biological Sciences, Napier University,
! Report ISBN 0 902703 60 9, 60 pp.
!
! Sharples Tett (1994). Modeling the effect of physical! variability on the midwater chlorophyll maximum.
! Journal of marine research 52: 219-238
!
! Implementation: Marilaure Gregoire,             NIOO-CEME
! Translation into FABM: E. Ivanov, ULg / MAST
!
! Contains the pelagic submodel, as used in Soetaert et al., 2001.
!
!--------------------------------------------------------------------*

   module fabm_ulg_emiliana 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_emiliana 
      type (type_state_variable_id)         :: id_cem,id_nem
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_uptake_c_emi,id_uptake_n_emi,id_npp,id_reduction_nitrate_phy,id_respiration_emi
      type (type_diagnostic_variable_id)    :: id_chla

!     Model parameters 
      real(rk)     :: dr_emi, exc_extra_doc, f_dl_dom, f_dl_phy_ex
      real(rk)     :: f_dl_phy_mo, f_leak_phy, f_pp_resp_emi, k_d
      real(rk)     :: ki_nhs_phy, ks_nhs_emi, ks_nos_emi, ks_po4_emi
      real(rk)     :: mo_emi, mumax_emi, pi_emi, q10_phy, r_o2_c_resp
      real(rk)     :: r_o2_nhs_nitr, r_p_n_redfield, respb_emi
      real(rk)     :: rmax_chl_n_emi, rmax_n_c_emi, rmin_chl_n_emi
      real(rk)     :: rmin_n_c_emi, umax_nhs_emi, umax_nos_emi
      real(rk)     :: umax_po4_emi, w_emi

      contains 

      procedure :: initialize 
      procedure :: do 
      procedure :: do_bottom  
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

   contains
   ! Initialise the Emiliana model

   subroutine initialize(self,configunit)
   class (type_ulg_emiliana), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

!     Model parameters 
      real(rk)     :: dr_emi, exc_extra_doc, f_dl_dom, f_dl_phy_ex
      real(rk)     :: f_dl_phy_mo, f_leak_phy, f_pp_resp_emi, k_d
      real(rk)     :: ki_nhs_phy, ks_nhs_emi, ks_nos_emi, ks_po4_emi
      real(rk)     :: mo_emi, mumax_emi, pi_emi, q10_phy, r_o2_c_resp
      real(rk)     :: r_o2_nhs_nitr, r_p_n_redfield, respb_emi
      real(rk)     :: rmax_chl_n_emi, rmax_n_c_emi, rmin_chl_n_emi
      real(rk)     :: rmin_n_c_emi, umax_nhs_emi, umax_nos_emi
      real(rk)     :: umax_po4_emi, w_emi


   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%dr_emi, 'dr_emi', 'mol d-1', 'Deposition rate of EM', default=5.7e-06_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%exc_extra_doc, 'exc_extra_doc', '-', 'Extra-photosynthetic DOC excretion', default=0.05_rk) 
   call self%get_parameter(self%f_dl_dom, 'f_dl_dom', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%f_dl_phy_ex, 'f_dl_phy_ex', '-', 'Labile fraction phytoxcreted DOC', default=0.65_rk) 
   call self%get_parameter(self%f_dl_phy_mo, 'f_dl_phy_mo', '-', 'DOM fraction of phytoplankton mortality', default=0.34_rk) 
   call self%get_parameter(self%f_leak_phy, 'f_leak_phy', '-', 'Phytoplankton leakage fraction', default=0.02_rk) 
   call self%get_parameter(self%f_pp_resp_emi, 'f_pp_resp_emi', '-', 'Part of primary production used for respiration by EM ', default=0.1_rk) 
   call self%get_parameter(self%k_d, 'k_d', 'm-1', 'Background light attanuation coefficient', default=0.03_rk) 
   call self%get_parameter(self%ki_nhs_phy, 'ki_nhs_phy', 'mmolN m-3', 'Inhib. constant of NHS for NOS uptake by PHY', default=0.5_rk) 
   call self%get_parameter(self%ks_nhs_emi, 'ks_nhs_emi', 'mmolN m-3', 'Half-saturation constant for NHS uptake by EM', default=0.05_rk) 
   call self%get_parameter(self%ks_nos_emi, 'ks_nos_emi', 'mmolN m-3', 'Half-saturation constant for NOS uptake by EM', default=0.05_rk) 
   call self%get_parameter(self%ks_po4_emi, 'ks_po4_emi', 'mmolP m-3', 'Half-saturation constant for PO4 uptake by EM', default=0.02_rk) 
   call self%get_parameter(self%mo_emi, 'mo_emi', 'd-1', 'Mortality rate of EM', default=0.03_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%mumax_emi, 'mumax_emi', 'd-1', 'Maximum specific growth rate of EM', default=2.5_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%pi_emi, 'pi_emi', 'm2 W-1 d-1', 'Initial slope of photosynthesis-light curve for EM', default=0.3_rk, scale_factor=one_pr_day) 
   call self%get_parameter(self%q10_phy, 'q10_phy', '-', 'Temperature factor', default=2.0_rk) 
   call self%get_parameter(self%r_o2_c_resp, 'r_o2_c_resp', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%r_o2_nhs_nitr, 'r_o2_nhs_nitr', 'molO2 molNS-1', 'O2:NHS ratio in NHS oxidation in nitrification', default=2.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%respb_emi, 'respb_emi', 'd-1', 'Basal respiration rate of EM', default=0.009_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rmax_chl_n_emi, 'rmax_chl_n_emi', 'g Chla molN-1', 'Maximum Chl:N ratio in EM', default=2.0_rk) 
   call self%get_parameter(self%rmax_n_c_emi, 'rmax_n_c_emi', 'molN molC-1', 'Maximum N:C ratio in EM', default=0.2_rk) 
   call self%get_parameter(self%rmin_chl_n_emi, 'rmin_chl_n_emi', 'g Chla molN-1', 'Minimum Chl:N ratio in EM', default=1.0_rk) 
   call self%get_parameter(self%rmin_n_c_emi, 'rmin_n_c_emi', 'molN molC-1', 'Minimum N:C ratio in EM', default=0.05_rk) 
   call self%get_parameter(self%umax_nhs_emi, 'umax_nhs_emi', 'molN molC-1 d-1', 'Maximal NHS uptake rate by EM', default=1.5_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%umax_nos_emi, 'umax_nos_emi', 'molN molC-1 d-1', 'Maximal NOS uptake rate by EM', default=1.5_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%umax_po4_emi, 'umax_po4_emi', 'molP molC-1 d-1', 'Maximal PO4 uptake rate by EM', default=0.09375_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%w_emi, 'w_emi', 'm d-1', 'Sinking velocity of EM', default=-1.5_rk, scale_factor=one_pr_day) 

   ! Register state variables 

   call self%register_state_variable(self%id_cem, 'CEM', 'mmol C m-3', 'Small flagellate biomass in carbon', minimum=0.0e-7_rk, vertical_movement=self%w_emi) 
   call self%register_state_variable(self%id_nem, 'NEM', 'mmol N m-3', 'Small flagellate biomass in nitrogen', minimum=0.0e-7_rk, vertical_movement=self%w_emi) 

   call self%register_state_dependency(self%id_dcl, 'DCL', 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'DCS', 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'DIC', 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'DNL', 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'DNS', 'Semi-labile detritus concentration in nitrogen', 'mmol N m-3') 
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
   call self%register_diagnostic_variable(self%id_uptake_c_emi, 'uptake_c_emi', 'mmol C m-3 d-1', & 
      'Carbon uptake by Emiliana', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_uptake_n_emi, 'uptake_n_emi', 'mmol N m-3 d-1', & 
      'Nitrogen uptake by Emiliana', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_npp, 'npp', 'mmol N m-3 d-1', & 
      ' Primary production of nitrogen by all types of phytoplankton', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_reduction_nitrate_phy, 'reduction_nitrate_phy', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_respiration_emi, 'respiration_emi', 'mmol C m-3 d-1', & 
      'Total Respiration of Emiliana', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_chla, 'chla','mg chl a m-3', 'Chlorophyll concentration')

    ! Add to aggregate variables 
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_cem)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nem)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_nem, scale_factor=self%r_p_n_redfield)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_chlorophyll',units="mg chl a m-3",aggregate_variable=.true.),self%id_chla,scale_factor=1._rk)

   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_nem, scale_factor=self%k_d)

   return 

99 call self%fatal_error('Emiliana', 'Error reading namelist ulg_emiliana') 

   end subroutine initialize 


   ! Right hand sides of Emiliana model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_emiliana), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CEM,NEM
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
      real(rk) ::   Ratio_Si_C	  ! mol Si mol C-1, Si/C ratio in diatoms
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
   _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
   _GET_(self%id_nem,NEM)       ! Small flagellate biomass in nitrogen
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
    Ratio_N_C  = Ratio(NEM,CEM)
    Ratio_Chl_C = ratio_chl_c_phyt(Ratio_N_C,self%rmax_n_c_emi,self%rmin_n_c_emi,self%rmin_chl_n_emi,self%rmax_chl_n_emi)
    
   ! Nitrate uptake rate 
    Uptake_NO3 = uptake_nitrate_phyt(Ratio_N_C,tf,self%rmax_n_c_emi,self%umax_nos_emi) * Michaelis(NOS,self%ks_nos_emi) * Inhibition(NHS,self%ki_nhs_phy) * CEM
    
   ! Ammonium uptake rate 
    Uptake_NHS = uptake_nutrient_phyt(Ratio_N_C,(NHS-0.0_rk),tf,self%rmax_n_c_emi,self%umax_nhs_emi,self%ks_nhs_emi) * CEM
    
   ! Phosphate uptake rate 
    Uptake_PO4 = uptake_nutrient_phyt(Ratio_N_C,(PHO-0.0_rk),tf,self%rmax_n_c_emi,self%umax_po4_emi,self%ks_po4_emi) * CEM
    
    
   ! Potential nitrogen uptake 
    Uptake_Nit = Uptake_NHS + Uptake_NO3
    
   ! Nutrient uptake 
    Uptake_nutrient = min(Uptake_Nit,Uptake_PO4/self%r_p_n_redfield)
    
   ! Compute actual N:C and Si:C ratios 
    Ratio_N_C = ratio_adjustment(Ratio_N_C,self%rmax_n_c_emi,self%rmin_n_c_emi)
    Ratio_Si_C = ratio_adjustment(1.0_rk,1.0_rk,0.0_rk)
    
    
   ! Compute nutrient and light limitation 
    Limitation_nutrient = limitation_by_nutrient(self%rmin_n_c_emi,0.0_rk,Ratio_N_C,Ratio_Si_C)
    Limitation_light = 1.0-exp(-self%pi_emi * par * 4.56 / self%mumax_emi) ! WattToPhoton = 4.56
    
   ! Compute carbon uptake 
    Uptake_C = self%mumax_emi  * CEM * tf * Limitation_nutrient * Limitation_light 
    
   ! Compute respiration 
    Respiration_tot = Uptake_C * self%f_pp_resp_emi + self%respb_emi * CEM * tf
    
   ! Compute growth 
    Growth = Uptake_C - Respiration_tot
    
   ! Compute the extra DOC excretion from the difference of growth in nutrient limited (actual NC ratio) and nutrient saturated (max NC ratio) 
    Excretion_ext = CEM * tf * self%exc_extra_doc * self%mumax_emi * extra_excretion(Limitation_light,self%rmin_n_c_emi,self%rmax_n_c_emi,Ratio_N_C)        
    
   ! Compute the leakage 
    Leakage_DOC = self%f_leak_phy * Uptake_C
    Leakage_DON = self%f_leak_phy * abs(Uptake_nutrient)
    
   ! Compute mortality 
    Mortality_C  = mortality_phyt(self%mo_emi,tf) * CEM
    Mortality_N  = mortality_phyt(self%mo_emi,tf) * NEM
    
   ! Carbon in Emiliana increases by growth and decreases by leakage and mortality 
   _ADD_SOURCE_(self%id_cem, Growth - Mortality_C - Leakage_DOC) 
   _ADD_SOURCE_(self%id_nem, Uptake_nutrient - Mortality_N - Leakage_DON) 
    
   ! IF CN ratio of phytoplankton higher than CNmin, than nitrogen is taken up, unless it gets excreted 
   IF (Uptake_nutrient.gt.0) THEN 
     _ADD_SOURCE_(self%id_nos, -Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit))
     _ADD_SOURCE_(self%id_dox, Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit) * self%r_o2_nhs_nitr) 
     _ADD_SOURCE_(self%id_nhs, -Uptake_nutrient * Ratio(Uptake_NHS,Uptake_Nit))
     _ADD_SOURCE_(self%id_pho, -Uptake_nutrient * self%r_p_n_redfield)
   ELSE 
     _ADD_SOURCE_(self%id_nhs, -Uptake_nutrient)
     _ADD_SOURCE_(self%id_pho, -Uptake_nutrient * self%r_p_n_redfield)
   END IF 
    
   ! Mortality increases the pool of POM and DOM with proper partitioning and leakage adds to the labile pool, extra-excretion and leakage add to the labile and semi-labile detritus 
   _ADD_SOURCE_(self%id_poc, (1.0 - self%f_dl_phy_mo) * Mortality_C)
   _ADD_SOURCE_(self%id_pon, (1.0 - self%f_dl_phy_mo) * Mortality_N)
   _ADD_SOURCE_(self%id_dcl, self%f_dl_phy_mo * Mortality_C * self%f_dl_dom + Leakage_DOC + self%f_leak_phy * Excretion_ext + self%f_dl_phy_ex * (1.0 - self%f_leak_phy) * Excretion_ext)
   _ADD_SOURCE_(self%id_dnl, self%f_dl_phy_mo * Mortality_N * self%f_dl_dom + Leakage_DON) 
   _ADD_SOURCE_(self%id_dcs, self%f_dl_phy_mo * Mortality_C * (1.0 - self%f_dl_dom) + (1.0 - self%f_dl_phy_ex) * (1.0 - self%f_leak_phy) * Excretion_ext) 
   _ADD_SOURCE_(self%id_dns, self%f_dl_phy_mo * Mortality_N * (1.0 - self%f_dl_dom))
   _ADD_SOURCE_(self%id_dox, (Growth + Excretion_ext) * self%r_o2_c_resp)
   _ADD_SOURCE_(self%id_dic, -Growth - Excretion_ext) 
    
   _SET_DIAGNOSTIC_(self%id_uptake_n_emi, Uptake_NHS + Uptake_NO3)
   _SET_DIAGNOSTIC_(self%id_uptake_c_emi, Uptake_C)
   _SET_DIAGNOSTIC_(self%id_respiration_emi, Respiration_tot)
   _SET_DIAGNOSTIC_(self%id_reduction_nitrate_phy, Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit) * self%r_o2_nhs_nitr)
   _SET_DIAGNOSTIC_(self%id_npp, Growth)
   _SET_DIAGNOSTIC_(self%id_chla, Ratio_Chl_C*CEM)

   _LOOP_END_

   end subroutine do


   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_ulg_emiliana), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk),save :: nem, cem

   _HORIZONTAL_LOOP_BEGIN_
 
   _GET_(self%id_nem,nem) 
   _GET_(self%id_cem,cem)    

   _ADD_BOTTOM_FLUX_(self%id_nem,-self%dr_emi*nem*nem)
   _ADD_BOTTOM_FLUX_(self%id_cem,-self%dr_emi*cem*cem)   

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

   end module fabm_ulg_emiliana 
