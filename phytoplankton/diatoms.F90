#include "fabm_driver.h" 
 
!#########################################################################################
!                              3Ddiatoms.F90
!
! 1-D ecosystem model - Biological model of Tett

! References to the microplankton model:
! Tett, P., 1998. Parameterising a microplankton model.
! Department of Biological Sciences, Napier University,
! Report ISBN 0 902703 60 9, 60 pp.
!
! Sharples Tett (1994). Modeling the effect of physical! variability on the midwater chlorophyll maximum.
! Journal of marine research 52: 219-238
!
! Implementation: Marilaure Gregoire,                 NIOO-CEME
! Translation in FABM: Evgeny Ivanov,   Universite de Liege, MAST
!
! Contains the pelagic submodel, as used in Soetaert et al., 2001.
!######################################################################

   module fabm_ulg_diatoms 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_diatoms 
      type (type_state_variable_id)         :: id_cdi,id_ndi,id_sid,id_sio
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      !type (type_bottom_state_variable_id)  :: id_cflux, id_nflux, id_sflux
      type (type_bottom_dependency_id)      :: id_taub
      type (type_diagnostic_variable_id)    :: id_uptake_c_dia,id_uptake_n_dia,id_npp,id_reduction_nitrate_phy,id_uptake_sio_dia,id_respiration_dia
      type (type_diagnostic_variable_id)    :: id_chla

!     Model parameters 
      real(rk)     :: dr_dia, dr_sid, exc_extra_doc, f_dl_dom
      real(rk)     :: f_dl_phy_ex, f_dl_phy_mo, f_leak_phy, f_pp_resp_dia
      real(rk)     :: hmax_sid, k_d, ki_nhs_phy, ks_nhs_dia, ks_nos_dia
      real(rk)     :: ks_po4_dia, ks_sio_dia, mo_dia, mumax_dia
      real(rk)     :: pi_dia, q10_dia, q10_si_diss, r_o2_c_resp
      real(rk)     :: r_o2_nhs_nitr, r_p_n_redfield, r_si_n_dia
      real(rk)     :: respb_dia, rmax_chl_n_dia, rmax_n_c_dia
      real(rk)     :: rmin_chl_n_dia, rmin_n_c_dia, taucr_dep, umax_nhs_dia
      real(rk)     :: umax_nos_dia, umax_po4_dia, umax_si_dia
      real(rk)     :: w_dia_min, w_dia_max, w_sid

      contains 

      procedure :: initialize 
      procedure :: do 
      procedure :: get_vertical_movement
      procedure :: do_bottom 
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

   contains
   ! Initialise the Diatoms model

   subroutine initialize(self,configunit)
   class (type_ulg_diatoms), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%dr_dia, 'dr_dia', 'mol d-1', 'Deposition rate of FL', default=5.8e-06_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%dr_sid, 'dr_sid', 'mol d-1', 'Deposition rate of silicious detritus', default=5.5e-06_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%exc_extra_doc, 'exc_extra_doc', '-', 'Extra-photosynthetic DOC excretion', default=0.05_rk) 
   call self%get_parameter(self%f_dl_dom, 'f_dl_dom', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%f_dl_phy_ex, 'f_dl_phy_ex', '-', 'Labile fraction phytoxcreted DOC', default=0.65_rk) 
   call self%get_parameter(self%f_dl_phy_mo, 'f_dl_phy_mo', '-', 'DOM fraction of phytoplankton mortality', default=0.34_rk) 
   call self%get_parameter(self%f_leak_phy, 'f_leak_phy', '-', 'Phytoplankton leakage fraction', default=0.02_rk) 
   call self%get_parameter(self%f_pp_resp_dia, 'f_pp_resp_dia', '-', 'Part of primary production used for respiration by DI', default=0.1_rk) 
   call self%get_parameter(self%hmax_sid, 'hmax_sid', 'd-1', 'Rate of dissolution of silicious detritus', default=0.08_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%k_d, 'k_d', 'm-1', 'Background light attanuation coefficient', default=0.03_rk) 
   call self%get_parameter(self%ki_nhs_phy, 'ki_nhs_phy', 'mmolN m-3', 'Inhib. constant of NHS for NOS uptake by PHY', default=0.5_rk) 
   call self%get_parameter(self%ks_nhs_dia, 'ks_nhs_dia', 'mmolN m-3', 'Half-saturation constant for NHS uptake by DI', default=1.0_rk) 
   call self%get_parameter(self%ks_nos_dia, 'ks_nos_dia', 'mmolN m-3', 'Half-saturation constant for NOS uptake by DI', default=1.0_rk) 
   call self%get_parameter(self%ks_po4_dia, 'ks_po4_dia', 'mmolP m-3', 'Half-saturation constant for PO4 uptake by DI', default=0.1_rk) 
   call self%get_parameter(self%ks_sio_dia, 'ks_sio_dia', 'mmolSi m-3', 'Half-saturation constant for SiOs uptake by DI', default=3.5_rk) 
   call self%get_parameter(self%mo_dia, 'mo_dia', 'd-1', 'Mortality rate of DI', default=0.03_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%mumax_dia, 'mumax_dia', 'd-1', 'Maximum specific growth rate of DI', default=3.5_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%pi_dia, 'pi_dia', 'm2 W-1 d-1', 'Initial slope of photosynthesis-light curve for DI', default=0.3312_rk, scale_factor=one_pr_day) 
   call self%get_parameter(self%q10_dia, 'q10_dia', '-', 'Temperature factor for DI', default=1.8_rk) 
   call self%get_parameter(self%q10_si_diss, 'q10_si_diss', '-', 'Temperature factor for chemical processes', default=3.3_rk) 
   call self%get_parameter(self%r_o2_c_resp, 'r_o2_c_resp', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%r_o2_nhs_nitr, 'r_o2_nhs_nitr', 'molO2 molNS-1', 'O2:NHS ratio in NHS oxidation in nitrification', default=2.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%r_si_n_dia, 'r_si_n_dia', 'molSi molN-1', 'Si:N ratio in DI', default=0.83_rk) 
   call self%get_parameter(self%respb_dia, 'respb_dia', 'd-1', 'Basal respiration rate of DI', default=0.009_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rmax_chl_n_dia, 'rmax_chl_n_dia', 'g Chla molN-1', 'Maximum Chl:N ratio in DI', default=2.0_rk) 
   call self%get_parameter(self%rmax_n_c_dia, 'rmax_n_c_dia', 'molN molC-1', 'Maximum N:C ratio in DI', default=0.2_rk) 
   call self%get_parameter(self%rmin_chl_n_dia, 'rmin_chl_n_dia', 'g Chla molN-1', 'Minimum Chl:N ratio in DI', default=1.0_rk) 
   call self%get_parameter(self%rmin_n_c_dia, 'rmin_n_c_dia', 'molN molC-1', 'Minimum N:C ratio in DI', default=0.05_rk) 
   call self%get_parameter(self%taucr_dep, 'taucr_dep', 'N m-2', 'critical shear stress for deposition', default=0.02_rk)
   call self%get_parameter(self%umax_nhs_dia, 'umax_nhs_dia', 'molN molC-1 d-1', 'Maximal NHS uptake rate by DI', default=1.0_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%umax_nos_dia, 'umax_nos_dia', 'molN molC-1 d-1', 'Maximal NOS uptake rate by DI', default=1.0_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%umax_po4_dia, 'umax_po4_dia', 'molP molC-1 d-1', 'Maximal PO4 uptake rate by DI', default=0.0625_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%umax_si_dia, 'umax_si_dia', 'molSi molC-1 d-1', 'Maximal SiOs uptake rate by DI', default=0.5_rk, scale_factor=one_pr_day)
   !call self%get_parameter(self%w_dia, 'w_dia', 'm d-1', 'Sinking velocity of DI', default=-1.0_rk, scale_factor=one_pr_day) 
   call self%get_parameter(self%w_dia_min, 'w_dia_min', 'm d-1', 'Sinking velocity of DI', default=-0.1_rk, scale_factor=one_pr_day) 
   call self%get_parameter(self%w_dia_max, 'w_dia_max', 'm d-1', 'Sinking velocity of DI', default=-1.0_rk, scale_factor=one_pr_day) 
   call self%get_parameter(self%w_sid, 'w_sid', 'm d-1', 'Sinking velocity of silicious detritus', default=-2.0_rk, scale_factor=one_pr_day) 

   ! Register state variables 
   call self%register_state_variable(self%id_cdi, 'CDI', 'mmol C m-3', 'Diatom biomass in carbon', minimum=0.0e-7_rk) !, vertical_movement=self%w_dia) 
   call self%register_state_variable(self%id_ndi, 'NDI', 'mmol N m-3', 'Diatom biomass in nitrogen', minimum=0.0e-7_rk) !, vertical_movement=self%w_dia) 
   call self%register_state_variable(self%id_sid, 'SID', 'mmol Si m-3', 'Detrital silicate concentration', minimum=0.0e-7_rk, vertical_movement=self%w_sid) 
   call self%register_state_variable(self%id_sio, 'SIO', 'mmol Si m-3', 'Silicilic acid concentration', minimum=0.0e-7_rk)

   ! Couple to benthic pools to deposit sinking material in.
   !call self%register_state_dependency(self%id_cflux, 'cflux', '??', 'not a flux but a pool C')
   !call self%register_state_dependency(self%id_nflux, 'nflux', '??', 'not a flux but a pool N')
   !call self%register_state_dependency(self%id_sflux, 'sflux', '??', 'not a flux but a pool Si')

   call self%register_dependency(self%id_taub, standard_variables%bottom_stress)

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
   call self%register_diagnostic_variable(self%id_uptake_c_dia, 'uptake_c_dia', 'mmol C m-3 d-1', & 
      'Carbon uptake by diatoms', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_uptake_n_dia, 'uptake_n_dia', 'mmol N m-3 d-1', & 
      'Nitrogen uptake by diatoms', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_npp, 'npp', 'mmol N m-3 d-1', & 
      ' Primary production of nitrogen by all types of phytoplankton', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_reduction_nitrate_phy, 'reduction_nitrate_phy', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_uptake_sio_dia, 'uptake_sio_dia', 'mmol Si m-3 d-1', & 
      'Uptake of silicates by diatoms', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_respiration_dia, 'respiration_dia', 'mmol C m-3 d-1', & 
      'Total Respiration of Diatoms', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_chla, 'chla','mg chl a m-3', 'Chlorophyll concentration')

    ! Add to aggregate variables 
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_cdi)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ndi)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_ndi, scale_factor=self%r_p_n_redfield)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_ndi, scale_factor=self%r_si_n_dia)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_sio)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_sid)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_chlorophyll',units="mg chl a m-3",aggregate_variable=.true.),self%id_chla,scale_factor=1._rk)

   call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_ndi, scale_factor=self%k_d)

   return 

99 call self%fatal_error('Diatoms', 'Error reading namelist ulg_diatoms') 

   end subroutine initialize 


   ! Right hand sides of Diatoms model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_diatoms), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CDI,NDI,SID,SIO
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
      real(rk) ::   Ratio_max_Si_C	  ! mol Si mol C-1, Maximum Si/C ratio in diatoms
      real(rk) ::   Ratio_min_SiO_C	  ! mol Si mol C-1, Minimum Si/C ratio in diatoms
      real(rk) ::   Respiration_tot	  ! mmol C m-3, Total phytoplankton respiration (basal & activity)
      real(rk) ::   Uptake_C	  ! mmol C m-3, C assimilation of phytoplankton
      real(rk) ::   Uptake_NHS	  ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   Uptake_NO3	  ! mmol N m-3, Nitrate uptake of phytoplankton
      real(rk) ::   Uptake_Nit	  ! mmol N m-3, Nitrogen uptake of phytoplankton
      real(rk) ::   Uptake_PO4	  ! mmol P m-3, Phosphate uptake by large flagellates
      real(rk) ::   Uptake_SiO	  ! mmol Si m-3, Diatoms silicate uptake
      real(rk) ::   Uptake_nutrient	  ! mmol m-3, Nutrient uptake of phytoplankton
      real(rk) ::   tf	  ! -, Temperature factor
      real(rk) ::   tf_SiO_dissolution	  ! -, Silicate dissolution adjusted for temperature

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
   _GET_(self%id_ndi,NDI)       ! Diatom biomass in nitrogen
   _GET_(self%id_sid,SID)       ! Detrital silicate concentration
   _GET_(self%id_sio,SIO)       ! Silicilic acid concentration
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
    
    tf = Q10Factor (temp,self%q10_dia)
    tf_SiO_dissolution = Q10Factor(temp,self%q10_si_diss)
    
   ! Calculate ratios in phytoplankton 
    Ratio_N_C  = Ratio(NDI,CDI)
    Ratio_Chl_C = ratio_chl_c_phyt(Ratio_N_C,self%rmax_n_c_dia,self%rmin_n_c_dia,self%rmin_chl_n_dia,self%rmax_chl_n_dia)
    
   ! Nitrate uptake rate 
    Uptake_NO3 = uptake_nitrate_phyt(Ratio_N_C,tf,self%rmax_n_c_dia,self%umax_nos_dia) * Michaelis(NOS,self%ks_nos_dia) * Inhibition(NHS,self%ki_nhs_phy) * CDI
    
   ! Ammonium uptake rate 
    Uptake_NHS = uptake_nutrient_phyt(Ratio_N_C,(NHS-0.0),tf,self%rmax_n_c_dia,self%umax_nhs_dia,self%ks_nhs_dia) * CDI
    
   ! Phosphate uptake rate 
    Uptake_PO4 = uptake_nutrient_phyt(Ratio_N_C,(PHO-0.0),tf,self%rmax_n_c_dia,self%umax_po4_dia,self%ks_po4_dia) * CDI
    
   ! Compute silicate uptake 
    Ratio_max_Si_C = self%rmax_n_c_dia * self%r_si_n_dia
    Ratio_min_SiO_C = self%rmin_n_c_dia * self%r_si_n_dia
    Ratio_Si_C = Ratio_N_C * self%r_si_n_dia
    Uptake_SiO = uptake_nutrient_phyt(Ratio_Si_C,(SIO-0.0),tf,Ratio_max_Si_C,self%umax_si_dia,self%ks_sio_dia) * CDI
    
   ! Potential nitrogen uptake 
    Uptake_Nit = Uptake_NHS + Uptake_NO3
    
   ! Nutrient uptake 
    Uptake_nutrient = min(Uptake_Nit,Uptake_SiO/self%r_si_n_dia,Uptake_PO4/self%r_p_n_redfield)
    
   ! Compute actual N:C and Si:C ratios 
    Ratio_N_C = ratio_adjustment(Ratio_N_C,self%rmax_n_c_dia,self%rmin_n_c_dia)
    Ratio_Si_C = ratio_adjustment(Ratio_Si_C,Ratio_max_Si_C,Ratio_min_SiO_C)
    
   ! Compute nutrient and light limitation 
    Limitation_nutrient = limitation_by_nutrient(self%rmin_n_c_dia,Ratio_min_SiO_C,Ratio_N_C,Ratio_Si_C)
    Limitation_light = 1.-exp(-self%pi_dia * par * 4.56 / self%mumax_dia) ! WattToPhoton = 4.56
    
   ! Compute carbon uptake 
    Uptake_C = self%mumax_dia  * Limitation_nutrient * Limitation_light * CDI * tf 
    
   ! Compute respiration 
    Respiration_tot = Uptake_C * self%f_pp_resp_dia + self%respb_dia * CDI * tf
    
   ! Compute growth 
    Growth = Uptake_C - Respiration_tot
    
   ! Compute the extra DOC excretion from the difference of growth in nutrient limited (actual NC ratio) and nutrient saturated (max NC ratio)        
    Excretion_ext = CDI * tf * self%exc_extra_doc * self%mumax_dia * extra_excretion(Limitation_light,self%rmin_n_c_dia,self%rmax_n_c_dia,Ratio_N_C)    
    
   ! Compute the leakage 
    Leakage_DOC = self%f_leak_phy * Uptake_C
    Leakage_DON = self%f_leak_phy * abs(Uptake_nutrient)
    
   ! Compute mortality 
    Mortality_C  = mortality_phyt(self%mo_dia,tf) * CDI
    Mortality_N  = mortality_phyt(self%mo_dia,tf) * NDI
    
   !Carbon in Diatoms increases by growth and decreases by leakage and mortality 
   _ADD_SOURCE_(self%id_cdi, Growth - Mortality_C - Leakage_DOC)
   _ADD_SOURCE_(self%id_ndi, Uptake_nutrient - Mortality_N - Leakage_DON) 
    
   ! IF CN ratio of phytoplankton higher than CNmin, than nitrogen is taken up, unless it gets excreted 
   IF (Uptake_nutrient.gt.0) THEN 
     _ADD_SOURCE_(self%id_nos, -Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit)) 
     _ADD_SOURCE_(self%id_dox,  Uptake_nutrient * Ratio(Uptake_NO3,Uptake_Nit)*self%r_o2_nhs_nitr) 
     _ADD_SOURCE_(self%id_nhs, -Uptake_nutrient * Ratio(Uptake_NHS,Uptake_Nit))
     _ADD_SOURCE_(self%id_pho, -Uptake_nutrient * self%r_p_n_redfield)
   ELSE 
     _ADD_SOURCE_(self%id_nhs, -Uptake_nutrient) 
     _ADD_SOURCE_(self%id_pho, -Uptake_nutrient * self%r_p_n_redfield) 
   ENDIF 

    Uptake_SiO = (Uptake_nutrient - Leakage_DON) * self%r_si_n_dia

   ! Mortality increases the pool of POM and DOM with proper partitioning and leakage adds to the labile pool, extra-excretion and leakage add to the labile and semi-labile detritus 
   _ADD_SOURCE_(self%id_sio, tf_SiO_dissolution * self%hmax_sid * SID - Uptake_SiO) 
   _ADD_SOURCE_(self%id_sid, Mortality_N * self%r_si_n_dia - tf_SiO_dissolution * self%hmax_sid * SID) 
   _ADD_SOURCE_(self%id_poc, (1.0 - self%f_dl_phy_mo) * Mortality_C)
   _ADD_SOURCE_(self%id_pon, (1.0 - self%f_dl_phy_mo) * Mortality_N)
   _ADD_SOURCE_(self%id_dcl, self%f_dl_phy_mo * Mortality_C * self%f_dl_dom + Leakage_DOC + self%f_leak_phy*Excretion_ext + self%f_dl_phy_ex * (1.0 - self%f_leak_phy) * Excretion_ext)
   _ADD_SOURCE_(self%id_dnl, self%f_dl_phy_mo * Mortality_N * self%f_dl_dom + Leakage_DON)
   _ADD_SOURCE_(self%id_dcs, self%f_dl_phy_mo * Mortality_C * (1.0 - self%f_dl_dom) + (1.0 - self%f_dl_phy_ex) * (1.0 - self%f_leak_phy) * Excretion_ext) 
   _ADD_SOURCE_(self%id_dns, self%f_dl_phy_mo * Mortality_N * (1.0 - self%f_dl_dom))
   _ADD_SOURCE_(self%id_dox, (Growth + Excretion_ext) * self%r_o2_c_resp)
   _ADD_SOURCE_(self%id_dic, -Growth - Excretion_ext)

   _SET_DIAGNOSTIC_(self%id_uptake_n_dia, Uptake_nutrient)
   _SET_DIAGNOSTIC_(self%id_uptake_sio_dia, Uptake_SiO)
   _SET_DIAGNOSTIC_(self%id_uptake_c_dia, Uptake_C)
   _SET_DIAGNOSTIC_(self%id_respiration_dia, Respiration_tot)
   _SET_DIAGNOSTIC_(self%id_npp, Growth)
   _SET_DIAGNOSTIC_(self%id_reduction_nitrate_phy, Uptake_nutrient * Ratio(Uptake_NO3,(1.0e-10+Uptake_Nit)) * self%r_o2_nhs_nitr)
   _SET_DIAGNOSTIC_(self%id_chla, Ratio_Chl_C*CDI)


   _LOOP_END_

   end subroutine do


   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class (type_ulg_diatoms), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: Ratio_N_C, sink_rate, w_dia
      real(rk) :: cdi, ndi

      _LOOP_BEGIN_

        _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
        _GET_(self%id_ndi,NDI)       ! Diatom biomass in nitrogen

        Ratio_N_C  = Ratio(NDI,CDI)

	! In BAMHBI, sinking rates are positive, and in FABM - negative. Here, we account for that.
        sink_rate = -1.0*(self%w_dia_max-self%w_dia_min)*(self%rmax_n_c_dia+Ratio_N_C)/(self%rmin_n_c_dia-self%rmax_n_c_dia)
        w_dia = min(max(sink_rate,-1.0*self%w_dia_min),-1.0*self%w_dia_max)

        _ADD_VERTICAL_VELOCITY_(self%id_cdi, -w_dia)
        _ADD_VERTICAL_VELOCITY_(self%id_ndi, -w_dia)

      _LOOP_END_

   end subroutine get_vertical_movement


   ! Uncomment the commented lines to get benthic-pelagic coupling

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_ulg_diatoms), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_

   real(rk),save :: cdi, ndi, sid
   real(rk) :: taub, f_resp

   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_cdi,cdi)   
   _GET_(self%id_ndi,ndi) 
   _GET_(self%id_sid,sid)    

   !_GET_BOTTOM_(self%id_taub,  taub)

   !f_resp = taub / self%taucr_dep

   ! if (f_resp .le. 1.0_rk) THEN

   !  f_resp = 1.0_rk - f_resp
   
   ! default if no betnhic-pelagic coupling defined
   _ADD_BOTTOM_FLUX_(self%id_cdi,-self%dr_dia*cdi*cdi)   
   _ADD_BOTTOM_FLUX_(self%id_ndi,-self%dr_dia*ndi*ndi)
   _ADD_BOTTOM_FLUX_(self%id_sid,-self%dr_sid*sid*sid)   

   !  _ADD_BOTTOM_FLUX_(self%id_cdi,-self%dr_dia*cdi*f_resp)
   !  _ADD_BOTTOM_FLUX_(self%id_ndi,-self%dr_dia*ndi*f_resp)
   !  _ADD_BOTTOM_FLUX_(self%id_sid,-self%dr_sid*sid*f_resp)
 
   !  _ADD_BOTTOM_SOURCE_(self%id_cflux, self%dr_dia*cdi*f_resp)
   !  _ADD_BOTTOM_SOURCE_(self%id_nflux, self%dr_dia*ndi*f_resp)
   !  _ADD_BOTTOM_SOURCE_(self%id_sflux, (self%dr_sid*sid + self%r_si_n_dia*self%dr_dia*ndi) * f_resp)

   ! else
   !  _ADD_BOTTOM_SOURCE_(self%id_cflux, 0.0_rk)
   !  _ADD_BOTTOM_SOURCE_(self%id_nflux, 0.0_rk)
   !  _ADD_BOTTOM_SOURCE_(self%id_sflux, 0.0_rk)

   ! endif

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

   end module fabm_ulg_diatoms 
