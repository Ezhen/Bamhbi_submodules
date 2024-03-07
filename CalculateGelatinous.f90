#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DGELATINOUS
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

   module fabm_ulg_gelatinous 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_gelatinous 
      type (type_state_variable_id)         :: id_gel
      type (type_state_variable_id)         :: id_cdi,id_cem,id_cfl,id_dic,id_dox,id_mes,id_mic,id_ndi,id_nem,id_nfl,id_nhs,id_pho,id_poc,id_pon,id_sid
      type (type_dependency_id)             :: id_temp 
      type (type_diagnostic_variable_id)    :: id_respiration_gel

!     Model parameters 
      real(rk)     :: doxsatmort, eff_ass_gel_prey, eff_gel_dia
      real(rk)     :: eff_gel_emi, eff_gel_fla, eff_gel_mes, eff_gel_mic
      real(rk)     :: eff_gel_pom, eff_gr_gel_c, gmax_gel, ks_mort_gel
      real(rk)     :: mo_anox_pred, momax_gel, q10_gel, r_n_c_gel
      real(rk)     :: r_o2_c_resp, r_p_n_redfield, r_si_n_dia
      real(rk)     :: respb_gel, t_g_gel

      contains 

      procedure :: initialize 
      procedure :: do 
      procedure :: do_bottom 
      procedure :: check_surface_state 
      procedure :: check_bottom_state 
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

   contains
   ! Initialise the Gelatinous model

   subroutine initialize(self,configunit)
   class (type_ulg_gelatinous), intent(inout), target :: self
   integer,                        intent(in)          :: configunit



   namelist /ulg_Gelatinous/ doxsatmort, 	 & 
                      eff_ass_gel_prey, eff_gel_dia, 	 & 
                      eff_gel_emi, eff_gel_fla, eff_gel_mes, 	 & 
                      eff_gel_mic, eff_gel_pom, eff_gr_gel_c, 	 & 
                      gmax_gel, ks_mort_gel, mo_anox_pred, 	 & 
                      momax_gel, q10_gel, r_n_c_gel, 	 & 
                      r_o2_c_resp, r_p_n_redfield, r_si_n_dia, 	 & 
                      respb_gel, t_g_gel

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%doxsatmort, 'doxsatmort', 'mmolO2 m-3 (?)', 'Perc. of sat. where metabolic respiration is 1/2 the one under O2 sat.', default=7.8125_rk) 
   call self%get_parameter(self%eff_ass_gel_prey, 'eff_ass_gel_prey', '-', 'GEL assimilation efficiency on prey', default=0.75_rk) 
   call self%get_parameter(self%eff_gel_dia, 'eff_gel_dia', '-', 'Capture efficiency of GEL on DI', default=0.0_rk) 
   call self%get_parameter(self%eff_gel_emi, 'eff_gel_emi', '-', 'Capture efficiency of GEL on EM', default=0.0_rk) 
   call self%get_parameter(self%eff_gel_fla, 'eff_gel_fla', '-', 'Capture efficiency of GEL on FL', default=0.0_rk) 
   call self%get_parameter(self%eff_gel_mes, 'eff_gel_mes', '-', 'Capture efficiency of GEL on MES', default=1.0_rk) 
   call self%get_parameter(self%eff_gel_mic, 'eff_gel_mic', '-', 'Capture efficiency of GEL on MIC', default=0.0_rk) 
   call self%get_parameter(self%eff_gel_pom, 'eff_gel_pom', '-', 'Capture efficiency of GEL on POM', default=0.0_rk) 
   call self%get_parameter(self%eff_gr_gel_c, 'eff_gr_gel_c', '-', 'Part of the assimil. food used for GEL growth', default=0.2_rk) 
   call self%get_parameter(self%gmax_gel, 'gmax_gel', 'd-1', 'Maximum grazing rate of GEL', default=0.3_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%ks_mort_gel, 'ks_mort_gel', 'mmolC m-3', '', default=0.0_rk) 
   call self%get_parameter(self%mo_anox_pred, 'mo_anox_pred', 'd-1', 'Mortality rate in anoxia', default=0.25_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%momax_gel, 'momax_gel', 'd-1', 'Maximum mortality rate of GEL', default=0.009_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%q10_gel, 'q10_gel', '-', 'Temperature dependency for GEL', default=3.5_rk) 
   call self%get_parameter(self%r_n_c_gel, 'r_n_c_gel', 'molN molC-1', 'N:C molar ratio in GEL', default=0.25_rk) 
   call self%get_parameter(self%r_o2_c_resp, 'r_o2_c_resp', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%r_si_n_dia, 'r_si_n_dia', 'molSi molN-1', 'Si:N ratio in DI', default=0.83_rk) 
   call self%get_parameter(self%respb_gel, 'respb_gel', 'd-1', 'Basal respiration rate of GEL', default=0.0001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%t_g_gel, 't_g_gel', 'mmolC m-3', 'Feeding threshold for GEL grazing', default=0_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_gel, 'GEL', 'mmol C m-3', 'Gelatinous omnivorous biomass',minimum=0.0e-7_rk)

   call self%register_state_dependency(self%id_cdi, 'Diatom biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cem, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cfl, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_mes, 'Mesozooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_mic, 'Microzooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_nhs, 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_ndi, 'Diatom biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nem, 'Small flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nfl, 'Large flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'Particulate organic nitrogen concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_sid, 'Detrital silicate concentration', 'mmol Si m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_respiration_gel, 'respiration_gel', 'mmol C m-3 d-1', & 
      'Total Respiration of Gelatinous', output=output_instantaneous) 
   return 

99 call self%fatal_error('Gelatinous', 'Error reading namelist ulg_gelatinous') 

   end subroutine initialize 


   ! Right hand sides of Gelatinous model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_gelatinous), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  CDI,CEM,CFL,DIC,DOX,MES,MIC,NDI,NEM,NFL,NHS,PHO,POC,PON,SID
      real(rk) ::  temp
      real(rk) ::  GEL
      real(rk) ::   Egestion_C	  ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   Egestion_N	  ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   Excretion_C_adj	  ! mmol N m-3, Potential excretion rate necessary to keep the N/C ratio of Gelationous constant
      real(rk) ::   Excretion_N_adj	  ! mmol C m-3, Potential additional respiration flux to keep the N/C ratio of Gelationous constant
      real(rk) ::   Grazing_C	  ! mmol C m-3, Grazing in carbon by microzooplankaton
      real(rk) ::   Grazing_N	  ! mmol N m-3, Grazing in nitrogen all zooplankaton
      real(rk) ::   Mortality_C	  ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Mortality_N	  ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   Prey_C	  ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   Ratio_N_C	  ! mol N mol C-1, N/C ratio in small flagellates
      real(rk) ::   Ratio_N_C_test	  ! -, N/C ratio of the zooplankton before adjustment
      real(rk) ::   Respiration_C	  ! mmol C m-3, Zooplankton respiration flux
      real(rk) ::   tf	  ! -, Temperature factor

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_gel,GEL)       ! Gelatinous omnivorous biomass
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
    
    tf = Q10Factor(temp,self%q10_gel)
    
   ! Flux of consummed preys in carbon 
    Prey_C = self%eff_gel_fla*CFL+self%eff_gel_emi*CEM+self%eff_gel_dia*CDI+self%eff_gel_mic*MIC+self%eff_gel_mes*MES+self%eff_gel_pom*POC
    
   ! Flux of consummed preys in nitrogen 
    Prey_N = self%eff_gel_fla*NFL+self%eff_gel_emi*NEM+self%eff_gel_dia*NDI+self%eff_gel_mic*MIC*self%r_n_c_mic + self%eff_gel_mes*MES*self%r_n_c_mes+self%eff_gel_pom*PON
    
   ! Grazing rate in carbon 
    if (Prey_C>self%t_g_gel) then 
      Grazing_C = tf * self%gmax_gel * (Prey_C - self%t_g_gel) * GEL
    else 
      Grazing_C = 0.0
    endif 
    
   ! Grazing rate of nitrogen 
    Ratio_N_C = Ratio(Prey_N,Prey_C)
    Grazing_N = Grazing_C * Ratio_N_C
    
   ! Egestion rate of zooplankton 
    Egestion_C = ZOOEgestion(self%eff_ass_gel_prey,Grazing_C)
    Egestion_N = Egestion_C * Ratio_N_C
    
   ! Zooplankton respiration  
    Respiration_C = GELRespiration(tf,self%eff_ass_gel_prey,self%eff_gr_gel_c,self%respb_gel,GEL,Grazing_C)
    
   ! Zooplankton mortality rate 
    Mortality_C = Mortality_consument(self%ks_mort_gel,self%momax_gel,1.0,self%doxsatmort,self%mo_anox_pred,tf,GEL,DOX)
    Mortality_N = Mortality_C * self%r_n_c_gel
    
   ! Computes the N/C fluxes necessary to conserve the N/C ratio 
    Ratio_N_C_test = (Grazing_C * Ratio_N_C - Egestion_N) / (Grazing_C - Egestion_C-Respiration_C)
    if (Ratio_N_C_test > self%r_n_c_gel) then 
      Excretion_C_adj = 0.0
      Excretion_N_adj = (Grazing_C * Ratio_N_C - Egestion_N) - (Grazing_C - Egestion_C-Respiration_C) * self%r_n_c_gel
    else 
      Excretion_N_adj = 0.0
    endif 
    
   ! Carbon content increases by intake of preys and decreases by egestion, respiration, mortality, adjustement 
   ! Particulate detritus pool if formed from non-assimilated grazed food and dead zooplatkton
   _ADD_SOURCE_(self%id_gel, Grazing_C - Egestion_C - Respiration_C - Mortality_C - Excretion_C_adj) 
   _ADD_SOURCE_(self%id_poc, Egestion_C + Mortality_C - Grazing_C / Prey_C * self%eff_gel_pom*POC) 
   _ADD_SOURCE_(self%id_pon, Egestion_N + Mortality_N - Grazing_C / Prey_C * self%eff_gel_pom*PON) 
   _ADD_SOURCE_(self%id_nhs, Excretion_N_adj) 
   _ADD_SOURCE_(self%id_pho, Excretion_N_adj * self%r_p_n_redfield)
   _ADD_SOURCE_(self%id_cfl, -Grazing_C / Prey_C * self%eff_gel_fla * CFL) 
   _ADD_SOURCE_(self%id_cem, -Grazing_C / Prey_C * self%eff_gel_emi * CEM)
   _ADD_SOURCE_(self%id_cdi, -Grazing_C / Prey_C * self%eff_gel_dia * CDI)
   _ADD_SOURCE_(self%id_nfl, -Grazing_C / Prey_C * self%eff_gel_fla * NFL)
   _ADD_SOURCE_(self%id_nem, -Grazing_C / Prey_C * self%eff_gel_emi * NEM)
   _ADD_SOURCE_(self%id_ndi,- Grazing_C / Prey_C * self%eff_gel_dia * NDI) 
   _ADD_SOURCE_(self%id_mic, -Grazing_C / Prey_C * self%eff_gel_mic * MIC) 
   _ADD_SOURCE_(self%id_mes, -Grazing_C / Prey_C * self%eff_gel_mes * MES)
   _ADD_SOURCE_(self%id_sid, Grazing_C / Prey_C * self%eff_gel_dia * NDI * self%r_si_n_dia)
   _ADD_SOURCE_(self%id_dox, -self%r_o2_c_resp * (Respiration_C + Excretion_C_adj)) 
   _ADD_SOURCE_(self%id_dic, Respiration_C + Excretion_C_adj)

   _SET_DIAGNOSTIC_(self%id_respiration_gel, Respiration_C + Excretion_C_adj)

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_gelatinous 
