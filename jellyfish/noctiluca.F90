#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DNOCTILUCAF90
!
! Contains the pelagic submodel, for gelatinous (noctiluca, aurelia, mnemiopsis) based on
! the model of Lancelot et al. (2002) published in ESCS. Parameters have to be adjusted.
! The gelatinous model differs from the classic zooplankton model essentially by the term of grazing.
! The grazing rate increases linearly with the prey concentration as in the litterature.
! Otherwise, the mortality is a non-linear term as fro the zooplankton. This term will have to be thoroughly calibrated since it is the
! closure of the model. The egestion is assumed to be the non-assimilated food which directly goes to POM. Respiration
! is considered to be the part of the assimilated food which is not used for the growth (as for the zooplankton).
! as in Lancelot et al. (2002), each gelatinous group has no predator. Noctiluca feeds on the three phytoplnakton group,
! the microzooplaqnkton, bacteria and the two POMS. Aurelia and Mnemiopsis feed both only on mesozoo.
! A flux of adjustment considered as an extra-respiration or an excretion is computed to maimtain constant the N:C ratio.
! This part has to be improved on the model of the zooplankton
!--------------------------------------------------------------------*

   module fabm_ulg_noctiluca 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_noctiluca 
      type (type_state_variable_id)         :: id_noc
      type (type_state_variable_id)         :: id_cdi,id_cem,id_cfl,id_dic,id_dox,id_mes,id_mic,id_ndi,id_nem,id_nfl,id_nhs,id_pho,id_poc,id_pon,id_sid
      type (type_dependency_id)             :: id_temp 
      type (type_diagnostic_variable_id)    :: id_respiration_gel

!     Model parameters 
      real(rk)     :: doxsatmort, eff_ass_noc_prey, eff_gr_noc_c
      real(rk)     :: eff_noc_dia, eff_noc_emi, eff_noc_fla, eff_noc_mes
      real(rk)     :: eff_noc_mic, eff_noc_pom, gmax_noc, ks_mort_noc
      real(rk)     :: mo_anox_pred, momax_noc, q10_zoo, r_n_c_mes  
      real(rk)     :: r_n_c_mic, r_n_c_noc, r_o2_c_resp, r_p_n_redfield
      real(rk)     :: r_si_n_dia, respb_noc, t_g_noc

      contains 

      procedure :: initialize 
      procedure :: do
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

   contains
   ! Initialise the Noctiluca model

   subroutine initialize(self,configunit)
   class (type_ulg_noctiluca), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

!     Model parameters 
      real(rk)     :: doxsatmort, eff_ass_noc_prey, eff_gr_noc_c
      real(rk)     :: eff_noc_dia, eff_noc_emi, eff_noc_fla, eff_noc_mes
      real(rk)     :: eff_noc_mic, eff_noc_pom, gmax_noc, ks_mort_noc
      real(rk)     :: mo_anox_pred, momax_noc, q10_zoo, r_n_c_mes  
      real(rk)     :: r_n_c_mic, r_n_c_noc, r_o2_c_resp, r_p_n_redfield
      real(rk)     :: r_si_n_dia, respb_noc, t_g_noc

   namelist /ulg_noctiluca/ doxsatmort, 	 & 
                      eff_ass_noc_prey, eff_gr_noc_c, 	 & 
                      eff_noc_dia, eff_noc_emi, eff_noc_fla, 	 & 
                      eff_noc_mes, eff_noc_mic, eff_noc_pom, 	 & 
                      gmax_noc, ks_mort_noc, mo_anox_pred, 	 & 
                      momax_noc, q10_zoo, r_n_c_mes, r_n_c_mic,  & 
                      r_n_c_noc, r_o2_c_resp, r_p_n_redfield,  	 & 
                      r_si_n_dia, respb_noc, t_g_noc

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%doxsatmort, 'doxsatmort', 'mmolO2 m-3', 'Perc. of sat. where metabolic respiration is 1/2 the one under O2 sat.', default=7.8125_rk) 
   call self%get_parameter(self%eff_ass_noc_prey, 'eff_ass_noc_prey', '-', 'NOC assimilation efficiency on prey', default=0.75_rk) 
   call self%get_parameter(self%eff_gr_noc_c, 'eff_gr_noc_c', '-', 'Part of the assimil. food used for GEL growth', default=0.2_rk) 
   call self%get_parameter(self%eff_noc_dia, 'eff_noc_dia', '-', 'Capture efficiency of NOC on DIA', default=1.0_rk) 
   call self%get_parameter(self%eff_noc_emi, 'eff_noc_emi', '-', 'Capture efficiency of NOC on EMI', default=1.0_rk) 
   call self%get_parameter(self%eff_noc_fla, 'eff_noc_fla', '-', 'Capture efficiency of NOC on FLA', default=0.5_rk) 
   call self%get_parameter(self%eff_noc_mes, 'eff_noc_mes', '-', 'Capture efficiency of NOC on MES', default=0.0_rk) 
   call self%get_parameter(self%eff_noc_mic, 'eff_noc_mic', '-', 'Capture efficiency of NOC on MIC', default=1.0_rk) 
   call self%get_parameter(self%eff_noc_pom, 'eff_noc_pom', '-', 'Capture efficiency of NOC on POM', default=1.0_rk) 
   call self%get_parameter(self%gmax_noc, 'gmax_noc', 'd-1', 'Maximum grazing rate of NOC', default=0.06_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%ks_mort_noc, 'ks_mort_noc', 'mmolC m-3', '', default=0.0_rk) 
   call self%get_parameter(self%mo_anox_pred, 'mo_anox_pred', 'd-1', 'Mortality rate in anoxia', default=0.25_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%momax_noc, 'momax_noc', 'd-1', 'Maximum mortality rate of NOC', default=0.06_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%q10_zoo, 'q10_zoo', '-', 'Temperature factor Soetart et al., 2001', default=2.0_rk) 
   call self%get_parameter(self%r_n_c_mes, 'r_n_c_mes', 'molN molC-1', 'N:C molar ratio in MES', default=0.21_rk) 
   call self%get_parameter(self%r_n_c_mic, 'r_n_c_mic', 'molN molC-1', 'N:C molar ratio in MIC', default=0.18_rk) 
   call self%get_parameter(self%r_n_c_noc, 'r_n_c_noc', 'molN molC-1', 'N:C molar ratio in NOC', default=0.21_rk) 
   call self%get_parameter(self%r_o2_c_resp, 'r_o2_c_resp', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%r_si_n_dia, 'r_si_n_dia', 'molSi molN-1', 'Si:N ratio in DI', default=0.83_rk) 
   call self%get_parameter(self%respb_noc, 'respb_noc', 'd-1', 'Basal respiration rate of NOC', default=0.0001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%t_g_noc, 't_g_noc', 'mmolC m-3', 'Feeding threshold for NOC grazing', default=0.0_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_noc, 'NOC', 'mmol C m-3', 'Gelatinous carnivorous biomass', minimum=0.0e-7_rk)

   call self%register_state_dependency(self%id_cdi, 'CDI', 'Diatom biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cem, 'CEM', 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cfl, 'CFL', 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'DIC', 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'DOX', 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_mes, 'MES', 'Mesozooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_mic, 'MIC', 'Microzooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_ndi, 'NDI', 'Diatom biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nem, 'NEM', 'Small flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nfl, 'NFL', 'Large flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nhs, 'NHS', 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'PHO', 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'POC', 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'PON', 'Particulate organic nitrogen concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_sid, 'SID', 'Detrital silicate concentration', 'mmol Si m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 

    ! Add to aggregate variables 
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_noc)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_noc, scale_factor=self%r_n_c_noc)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_noc, scale_factor=self%r_n_c_noc*self%r_p_n_redfield)

    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_respiration_gel, 'respiration_gel', 'mmol C m-3 d-1', & 
      'Total Respiration of Gelatinous', output=output_instantaneous) 
   return 

99 call self%fatal_error('Noctiluca', 'Error reading namelist ulg_noctiluca') 

   end subroutine initialize 


   ! Right hand sides of Noctiluca model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_noctiluca), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  CDI,CEM,CFL,DIC,DOX,MES,MIC,NDI,NEM,NFL,NHS,PHO,POC,PON,SID
      real(rk) ::  temp
      real(rk) ::  NOC
      real(rk) ::   Egestion_C	  ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   Egestion_N	  ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   Excretion_C_adj	  ! mmol N m-3, Potential excretion rate necessary to keep the N/C ratio of Gelationous constant
      real(rk) ::   Excretion_N_adj	  ! mmol C m-3, Potential additional respiration flux to keep the N/C ratio of Gelationous constant
      real(rk) ::   Grazing_C	  ! mmol C m-3, Grazing in carbon by microzooplankaton
      real(rk) ::   Grazing_N	  ! mmol N m-3, Grazing in nitrogen all zooplankaton
      real(rk) ::   Mortality_C	  ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Mortality_N	  ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   Prey_N	  ! mmol C m-3, Flux of ingested preys in nitrogen 
      real(rk) ::   Prey_C	  ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   Ratio_N_C	  ! mol N mol C-1, N/C ratio in small flagellates
      real(rk) ::   Ratio_N_C_test	  ! -, N/C ratio of the zooplankton before adjustment
      real(rk) ::   Respiration_C	  ! mmol C m-3, Zooplankton respiration flux
      real(rk) ::   tf	  ! -, Temperature factor

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_noc,NOC)       ! Gelatinous carnivorous biomass
   _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
   _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
   _GET_(self%id_cfl,CFL)       ! Large flagellate biomass in carbon
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_mes,MES)       ! Mesozooplakton biomass
   _GET_(self%id_mic,MIC)       ! Microzooplakton biomass
   _GET_(self%id_ndi,NDI)       ! Diatom biomass in nitrogen
   _GET_(self%id_nem,NEM)       ! Small flagellate biomass in nitrogen
   _GET_(self%id_nfl,NFL)       ! Large flagellate biomass in nitrogen
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_pho,PHO)       ! Phosphorus
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration
   _GET_(self%id_sid,SID)       ! Detrital silicate concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor (temp,self%q10_zoo)
    
   ! Flux of consummed preys in carbon 
    Prey_C = self%eff_noc_fla*CFL+self%eff_noc_emi*CEM+self%eff_noc_dia*CDI+self%eff_noc_mic*MIC+self%eff_noc_mes*MES+self%eff_noc_pom*POC
    
   ! Flux of consummed preys in nitrogen 
    Prey_N = self%eff_noc_fla*NFL+self%eff_noc_emi*NEM+self%eff_noc_dia*NDI+self%eff_noc_mic*MIC*self%r_n_c_mic + self%eff_noc_mes*MES*self%r_n_c_mes+self%eff_noc_pom*PON
    
   ! Grazing rate in carbon 
    if (Prey_C>self%t_g_noc) then 
      Grazing_C = tf * self%gmax_noc * (Prey_C - self%t_g_noc) * NOC
    else 
      Grazing_C = 0.0
    endif 
    
   ! Grazing rate of nitrogen 
    Ratio_N_C = Ratio(Prey_N,Prey_C)
    Grazing_N = Grazing_C * Ratio_N_C
    
   ! Egestion rate of zooplankton 
    Egestion_C = zoo_egestion(self%eff_ass_noc_prey,Grazing_C)
    Egestion_N = Egestion_C * Ratio_N_C
    
   ! Zooplankton respiration 
    Respiration_C = respiration_jelly(tf,self%eff_ass_noc_prey,self%eff_gr_noc_c,self%respb_noc,NOC,Grazing_C)
    
   ! Zooplankton mortality rate 
    Mortality_C = mortality_consument(self%ks_mort_noc,self%momax_noc,1.0_rk,self%doxsatmort,self%mo_anox_pred,tf,NOC,DOX)
    Mortality_N = Mortality_C * self%r_n_c_noc
    
   ! Computes the N/C fluxes necessary to conserve the N/C ratio 
    Ratio_N_C_test = Ratio(Grazing_C * Ratio_N_C - Egestion_N,Grazing_C - Egestion_C - Respiration_C)
    if (Ratio_N_C_test > self%r_n_c_noc) then 
      Excretion_C_adj = 0.0
      Excretion_N_adj = (Grazing_C * Ratio_N_C - Egestion_N) - (Grazing_C - Egestion_C - Respiration_C) * self%r_n_c_noc
    else 
      Excretion_N_adj = 0.0
    endif 
    
   ! Carbon content increases by intake of preys and decreases by egestion, respiration, mortality and adjustement

   _ADD_SOURCE_(self%id_noc, Grazing_C - Egestion_C - Respiration_C - Mortality_C - Excretion_C_adj) 

   ! Particulate detritus pool if formed from non-assimilated grazed food and dead zooplatkton; and grazing on zoooplankton and on phytoplankton
   ! When eating diatoms, zooplankton, ejects silicate as silicious_detritus; DOX decreases due to respiration of zooplankton 

   _ADD_SOURCE_(self%id_poc, Egestion_C + Mortality_C - Ratio(Grazing_C,Prey_C) * self%eff_noc_pom*POC) 
   _ADD_SOURCE_(self%id_pon, Egestion_N + Mortality_N - Ratio(Grazing_C,Prey_C) * self%eff_noc_pom*PON) 
   _ADD_SOURCE_(self%id_nhs, Excretion_N_adj)
   _ADD_SOURCE_(self%id_pho, Excretion_N_adj * self%r_p_n_redfield) 
   _ADD_SOURCE_(self%id_mic, -Ratio(Grazing_C,Prey_C) * self%eff_noc_mic*MIC)
   _ADD_SOURCE_(self%id_mes, -Ratio(Grazing_C,Prey_C) * self%eff_noc_mes*MES) 
   _ADD_SOURCE_(self%id_cfl, -Ratio(Grazing_C,Prey_C) * self%eff_noc_fla*CFL) 
   _ADD_SOURCE_(self%id_cem, -Ratio(Grazing_C,Prey_C) * self%eff_noc_emi*CEM)
   _ADD_SOURCE_(self%id_cdi, -Ratio(Grazing_C,Prey_C) * self%eff_noc_dia*CDI)
   _ADD_SOURCE_(self%id_nfl, -Ratio(Grazing_C,Prey_C) * self%eff_noc_fla*NFL)
   _ADD_SOURCE_(self%id_nem, -Ratio(Grazing_C,Prey_C) * self%eff_noc_emi*NEM)
   _ADD_SOURCE_(self%id_ndi, -Ratio(Grazing_C,Prey_C) * self%eff_noc_dia*NDI) 
   _ADD_SOURCE_(self%id_sid, Ratio(Grazing_C,Prey_C) * self%eff_noc_dia*NDI*self%r_si_n_dia)
   _ADD_SOURCE_(self%id_dox,-self%r_o2_c_resp*(Respiration_C + Excretion_C_adj)) 
   _ADD_SOURCE_(self%id_dic, Respiration_C + Excretion_C_adj)

   _SET_DIAGNOSTIC_(self%id_respiration_gel, Respiration_C + Excretion_C_adj)
 
   _LOOP_END_

   end subroutine do


   end module fabm_ulg_noctiluca 
