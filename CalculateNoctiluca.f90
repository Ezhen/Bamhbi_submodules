#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DGELATINOUS
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

   module fabm_ulg_Noctiluca 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Noctiluca 
      type (type_state_variable_id)         :: id_noc
      type (type_state_variable_id)         :: id_cdi,id_cem,id_cfl,id_dic,id_dox,id_mes,id_mic,id_ndi,id_nem,id_nfl,id_nhs,id_pho,id_poc,id_pon,id_sid
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_TotalRespiration_Gel
      type (type_diagnostic_variable_id)    :: id_TotalRespiration_GelIntegrated

!     Model parameters 
      real(rk)     :: Ass_Eff_Noctiluca, basal_Resp_Noctiluca
      real(rk)     :: Capt_eff_Noctiluca_Diatoms, Capt_eff_Noctiluca_Emiliana
      real(rk)     :: Capt_eff_Noctiluca_Flagellates, Capt_eff_Noctiluca_Mesozoo
      real(rk)     :: Capt_eff_Noctiluca_Microzoo, Capt_eff_Noctiluca_POM
      real(rk)     :: DOXsatmort, efficiency_growth_Noctiluca
      real(rk)     :: expmortNoctiluca, HalfSatMort_Noctiluca
      real(rk)     :: MaxgrazingrateNoctiluca, Mortanoxic, NCrNoctiluca
      real(rk)     :: NLin_Mort_Noctiluca, OCr, PNRedfield, Q10Zoo
      real(rk)     :: SiNrDiatoms, threshold_feeding_Noctiluca

      contains 

      procedure :: initialize 
      procedure :: do 
      procedure :: do_bottom 
      procedure :: check_surface_state 
      procedure :: check_bottom_state 
      procedure :: get_light_extinction 
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

   contains
   ! Initialise the Noctiluca model

   subroutine initialize(self,configunit)
   class (type_ulg_Noctiluca), intent(inout), target :: self
   integer,                        intent(in)          :: configunit


   namelist /ulg_Noctiluca/ Ass_Eff_Noctiluca, 	 & 
                      basal_Resp_Noctiluca, 	 & 
                      Capt_eff_Noctiluca_Diatoms, 	 & 
                      Capt_eff_Noctiluca_Emiliana, 	 & 
                      Capt_eff_Noctiluca_Flagellates, 	 & 
                      Capt_eff_Noctiluca_Mesozoo, 	 & 
                      Capt_eff_Noctiluca_Microzoo, 	 & 
                      Capt_eff_Noctiluca_POM, DOXsatmort, 	 & 
                      efficiency_growth_Noctiluca, 	 & 
                      expmortNoctiluca, HalfSatMort_Noctiluca, 	 & 
                      MaxgrazingrateNoctiluca, Mortanoxic, 	 & 
                      NCrNoctiluca, NLin_Mort_Noctiluca, OCr, 	 & 
                      PNRedfield, Q10Zoo, SiNrDiatoms, 	 & 
                      threshold_feeding_Noctiluca

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%Ass_Eff_Noctiluca, 'Ass_Eff_Noctiluca', '-', 'NOC assimilation efficiency on prey', default=0.75_rk) 
   call self%get_parameter(self%basal_Resp_Noctiluca, 'basal_Resp_Noctiluca', 'd-1', 'Basal respiration rate of NOC', default=0.0001_rk) 
   call self%get_parameter(self%Capt_eff_Noctiluca_Diatoms, 'Capt_eff_Noctiluca_Diatoms', '-', 'Capture efficiency of NOC on DI', default=1.0_rk) 
   call self%get_parameter(self%Capt_eff_Noctiluca_Emiliana, 'Capt_eff_Noctiluca_Emiliana', '-', 'Capture efficiency of NOC on EM', default=1.0_rk) 
   call self%get_parameter(self%Capt_eff_Noctiluca_Flagellates, 'Capt_eff_Noctiluca_Flagellates', '-', 'Capture efficiency of NOC on FL', default=0.5_rk) 
   call self%get_parameter(self%Capt_eff_Noctiluca_Mesozoo, 'Capt_eff_Noctiluca_Mesozoo', '-', 'Capture efficiency of NOC on MES', default=0.0_rk) 
   call self%get_parameter(self%Capt_eff_Noctiluca_Microzoo, 'Capt_eff_Noctiluca_Microzoo', '-', 'Capture efficiency of NOC on MIC', default=1.0_rk) 
   call self%get_parameter(self%Capt_eff_Noctiluca_POM, 'Capt_eff_Noctiluca_POM', '-', 'Capture efficiency of NOC on POM', default=1.0_rk) 
   call self%get_parameter(self%DOXsatmort, 'DOXsatmort', 'mmolO2 m-3', 'Percentage of saturation where metabolic respiration is half the one under oxygen satyrated conditions', default=7.8125_rk) 
   call self%get_parameter(self%efficiency_growth_Noctiluca, 'efficiency_growth_Noctiluca', '-', 'Part of the assimil. food used for GEL growth', default=0.2_rk) 
   call self%get_parameter(self%expmortNoctiluca, 'expmortNoctiluca', 'd-2 (?)', 'Quadratic mortality of NOC', default=2.0_rk) 
   call self%get_parameter(self%HalfSatMort_Noctiluca, 'HalfSatMort_Noctiluca', 'mmolC m-3', '', default=0.0_rk) 
   call self%get_parameter(self%MaxgrazingrateNoctiluca, 'MaxgrazingrateNoctiluca', 'd-1', 'Maximum grazing rate of NOC', default=0.06_rk) 
   call self%get_parameter(self%Mortanoxic, 'Mortanoxic', 'd-1', 'Mortality rate in anoxia', default=0.25_rk) 
   call self%get_parameter(self%NCrNoctiluca, 'NCrNoctiluca', 'molN molC-1', 'N:C molar ratio in NOC', default=0.21_rk) 
   call self%get_parameter(self%NLin_Mort_Noctiluca, 'NLin_Mort_Noctiluca', '(?) d-1', 'Quadratic mortality rate of NOC', default=0.06_rk) 
   call self%get_parameter(self%OCr, 'OCr', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%Q10Zoo, 'Q10Zoo', '-', 'Temperature factor Soetart et al., 2001', default=2.0_rk) 
   call self%get_parameter(self%SiNrDiatoms, 'SiNrDiatoms', 'molSi molN-1', 'Si:N ratio in DI', default=0.83_rk) 
   call self%get_parameter(self%threshold_feeding_Noctiluca, 'threshold_feeding_Noctiluca', '?', 'Feeding threshold for NOC grazing', default=?_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_noc, 'NOC'  & 
         , 'mmol C m-3', 'Gelatinous carnivorous biomass' & 
         minimum=0.0e-7_rk)
   call self%register_state_dependency(self%id_cdi, 'Diatom biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cem, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cfl, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_mes, 'Mesozooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_mic, 'Microzooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_ndi, 'Diatom biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nem, 'Small flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nfl, 'Large flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nhs, 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'Particulate organic nitrogen concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_sid, 'Detrital silicate concentration', 'mmol Si m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_TotalRespiration_Gel, 'TotalRespiration_Gel', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_TotalRespiration_GelIntegrated, 'TotalRespiration_GelIntegrated', '-', & 
      '-', output=output_instantaneous) 

   return 

99 call self%fatal_error('Noctiluca', 'Error reading namelist ulg_Noctiluca') 

   end subroutine initialize 


   ! Right hand sides of Noctiluca model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_Noctiluca), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  CDI,CEM,CFL,DIC,DOX,MES,MIC,NDI,NEM,NFL,NHS,PHO,POC,PON,SID
      real(rk) ::  par,temp
      real(rk) ::  NOC
      real(rk) ::   TotalRespiration_Gel
      real(rk) ::   TotalRespiration_GelIntegrated
      real(rk) ::   C_ZOOAdjust	  ! mmol N m-3, Potential excretion rate necessary to keep the N/C ratio of Gelationous constant
      real(rk) ::   C_ZOOEgest	  ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   C_ZOOMort	  ! mmol C m-3, Zooplankton mortality flux in carbon
      real(rk) ::   C_ZOOResp	  ! flux, Zooplankton respiration
      real(rk) ::   FluxPrey_carbon	  ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   FluxPrey_nitrogen	  ! mmol N m-3, Flux of consummed preys in nitrogen
      real(rk) ::   grazing_carbonNoctiluca	  ! mmol C m-3, Grazing in carbon by Noctiluca
      real(rk) ::   grazing_nitrogen	  ! mmol N m-3, Grazing in nitrogen by gelatinous
      real(rk) ::   N_ZOOAdjust	  ! mmol C m-3, Potential additional respiration flux to keep the N/C ratio of Gelationous constant
      real(rk) ::   N_ZOOEgest	  ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   N_ZOOMort	  ! mmol N m-3, Zooplankton mortality flux in nitrogen
      real(rk) ::   NCrfoodNoctiluca	  ! mmol N mmol C-1, N/C ratio in food of noctiluca
      real(rk) ::   NCrzootest	  ! -, N/C ratio of the zooplankton before adjustment
      real(rk) ::   tf	  ! -, Temperature factor
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_noc,NOC)       ! Gelatinous carnivorous biomass
   _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
   _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
   _GET_(self%id_cfl,CFL)       ! Small flagellate biomass in carbon
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
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor (temp,Q10Zoo)
    
   ! Flux of consummed preys in carbon 
    FluxPrey_carbon = self%Capt_eff_Noctiluca_Flagellates*CFL+self%Capt_eff_Noctiluca_Emiliana*CEM+self%Capt_eff_Noctiluca_Diatoms*CDI+self%Capt_eff_Noctiluca_Microzoo*MIC+self%Capt_eff_Noctiluca_Mesozoo*MES+self%Capt_eff_Noctiluca_POM*POC
    
   ! Flux of consummed preys in nitrogen 
    FluxPrey_nitrogen = self%Capt_eff_Noctiluca_Flagellates*NFL+self%Capt_eff_Noctiluca_Emiliana*NEM+self%Capt_eff_Noctiluca_Diatoms*NDI+self%Capt_eff_Noctiluca_Microzoo*MIC*self%NCrMicroZoo + self%Capt_eff_Noctiluca_Mesozoo*MES*self%NCrMesoZoo+self%Capt_eff_Noctiluca_POM*PON
    
   ! Grazing rate in carbon 
    if (FluxPrey_carbon>threshold_feeding_Noctiluca) then 
      grazing_carbonNoctiluca = tf*self%MaxgrazingrateNoctiluca*(FluxPrey_carbon-self%threshold_feeding_Noctiluca)*NOC
    else 
      grazing_carbonNoctiluca = 0.0
    endif 
    
   ! Grazing rate of nitrogen 
    NCrfoodNoctiluca = FluxPrey_nitrogen/FluxPrey_carbon
    grazing_nitrogen = grazing_carbonNoctiluca*NCrfoodNoctiluca
    
   ! Egestion rate of zooplankton 
    C_ZOOEgest = ZOOEgestion(Ass_Eff_Noctiluca,grazing_carbonNoctiluca)
    N_ZOOEgest = C_ZOOEgest*NCrfoodNoctiluca
    
   ! Zooplankton respiration 
    C_ZOOResp = GELRespiration(tf,Ass_Eff_Noctiluca,efficiency_growth_Noctiluca,basal_Resp_Noctiluca,NOC,grazing_carbonNoctiluca)
    
   ! Zooplankton mortality rate 
    C_ZOOMort = Mortality_consument(HalfSatMort_Noctiluca,NLin_Mort_Noctiluca,DOXsatmort,Mortanoxic,tf,NOC,DOX)
    N_ZOOMort = C_ZOOMort* self%NCrNoctiluca
    
   ! Computes the N/C fluxes necessary to conserve the N/C ratio 
    NCrzootest = (grazing_carbonNoctiluca*NCrfoodNoctiluca-N_ZOOEgest)/(grazing_carbonNoctiluca-C_ZOOEgest-C_ZOOResp)
    if (NCrzootest > NCrNoctiluca) then 
      C_ZOOAdjust = 0
      N_ZOOAdjust = (grazing_carbonNoctiluca*NCrfoodNoctiluca-N_ZOOEgest)-(grazing_carbonNoctiluca-C_ZOOEgest-C_ZOOResp)*self%NCrNoctiluca
    else 
      N_ZOOAdjust = 0
    endif 
    
   ! Carbon content increases by intake of preys and decreases by egestion, respiration, mortality and adjustement

   _ADD_SOURCE_(self%id_noc, grazing_carbonNoctiluca - C_ZOOEgest - C_ZOOResp - C_ZOOMort - C_ZOOAdjust) 

   ! Particulate detritus pool if formed from non-assimilated grazed food and dead zooplatkton; and grazing on zoooplankton and on phytoplankton
   ! When eating diatoms, zooplankton, ejects silicate as silicious_detritus; DOX decreases due to respiration of zooplankton 

   _ADD_SOURCE_(self%id_poc, C_ZOOEgest + C_ZOOMort - grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_POM*POC) 
   _ADD_SOURCE_(self%id_pon, N_ZOOEgest + N_ZOOMort - grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_POM*PON) 
   _ADD_SOURCE_(self%id_nhs, N_ZOOAdjust)
   _ADD_SOURCE_(self%id_pho, N_ZOOAdjust*self%PNRedfield) 
   _ADD_SOURCE_(self%id_mic, -grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Microzoo*MIC)
   _ADD_SOURCE_(self%id_mes, -grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Mesozoo*MES) 
   _ADD_SOURCE_(self%id_cfl, -grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Flagellates*CFL) 
   _ADD_SOURCE_(self%id_cem, -grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Emiliana*CEM)
   _ADD_SOURCE_(self%id_cdi, -grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Diatoms*CDI)
   _ADD_SOURCE_(self%id_nfl, -grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Flagellates*NFL)
   _ADD_SOURCE_(self%id_nem, -grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Emiliana*NEM)
   _ADD_SOURCE_(self%id_ndi,- grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Diatoms*NDI) 
   _ADD_SOURCE_(self%id_sid, grazing_carbonNoctiluca/FluxPrey_carbon*self%Capt_eff_Noctiluca_Diatoms*NDI*self%SiNrDiatoms)
   _ADD_SOURCE_(self%id_dox,-self%OCr*(C_ZOOResp + C_ZOOAdjust)) 
   _ADD_SOURCE_(self%id_dic, C_ZOOResp + C_ZOOAdjust)

   _SET_DIAGNOSTIC_(self%id_TotalRespiration_Gel, C_ZOOResp + C_ZOOAdjust)
 
   _LOOP_END_

   end subroutine do


   end module fabm_ulg_Noctiluca 
