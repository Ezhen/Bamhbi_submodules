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

   module fabm_ulg_Gelatinous 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Gelatinous 
      type (type_state_variable_id)         :: id_gel
      type (type_state_variable_id)         :: id_dic,id_dox,id_mes,id_mic,id_nhs,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_TotalRespiration_Gel
      type (type_diagnostic_variable_id)    :: 

!     Model parameters 
      real(rk)     :: Ass_Eff_Gelatinous, basal_Resp_Gelatinous
      real(rk)     :: Capt_eff_Gelatinous_Diatoms, Capt_eff_Gelatinous_Emiliana
      real(rk)     :: Capt_eff_Gelatinous_Flagellates, Capt_eff_Gelatinous_Mesozoo
      real(rk)     :: Capt_eff_Gelatinous_Microzoo, Capt_eff_Gelatinous_POM
      real(rk)     :: DOXsatmort, efficiency_growth_Gelatinous
      real(rk)     :: expmortGelatinous, HalfSatMort_Gelatinous
      real(rk)     :: MaxgrazingrateGelatinous, Mortanoxic, NCrGelatinous
      real(rk)     :: NLin_Mort_Gelatinous, OCr, PNRedfield, Q10Gelatinous
      real(rk)     :: threshold_feeding_Gelatinous

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
   ! Initialise the Gelatinous model

   subroutine initialize(self,configunit)
   class (type_ulg_Gelatinous), intent(inout), target :: self
   integer,                        intent(in)          :: configunit


   namelist /ulg_Gelatinous/ Ass_Eff_Gelatinous, 	 & 
                      basal_Resp_Gelatinous, 	 & 
                      Capt_eff_Gelatinous_Diatoms, 	 & 
                      Capt_eff_Gelatinous_Emiliana, 	 & 
                      Capt_eff_Gelatinous_Flagellates, 	 & 
                      Capt_eff_Gelatinous_Mesozoo, 	 & 
                      Capt_eff_Gelatinous_Microzoo, 	 & 
                      Capt_eff_Gelatinous_POM, DOXsatmort, 	 & 
                      efficiency_growth_Gelatinous, 	 & 
                      expmortGelatinous, 	 & 
                      HalfSatMort_Gelatinous, 	 & 
                      MaxgrazingrateGelatinous, Mortanoxic, 	 & 
                      NCrGelatinous, NLin_Mort_Gelatinous, 	 & 
                      OCr, PNRedfield, Q10Gelatinous, 	 & 
                      threshold_feeding_Gelatinous

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%Ass_Eff_Gelatinous, 'Ass_Eff_Gelatinous', '-', 'GEL assimilation efficiency on prey', default=0.75_rk) 
   call self%get_parameter(self%basal_Resp_Gelatinous, 'basal_Resp_Gelatinous', 'd-1', 'Basal respiration rate of GEL', default=0.0001_rk) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Diatoms, 'Capt_eff_Gelatinous_Diatoms', '-', 'Capture efficiency of GEL on DI', default=0.0_rk) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Emiliana, 'Capt_eff_Gelatinous_Emiliana', '-', 'Capture efficiency of GEL on EM', default=0.0_rk) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Flagellates, 'Capt_eff_Gelatinous_Flagellates', '-', 'Capture efficiency of GEL on FL', default=0.0_rk) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Mesozoo, 'Capt_eff_Gelatinous_Mesozoo', '-', 'Capture efficiency of GEL on MES', default=1.0_rk) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Microzoo, 'Capt_eff_Gelatinous_Microzoo', '-', 'Capture efficiency of GEL on MIC', default=0.0_rk) 
   call self%get_parameter(self%Capt_eff_Gelatinous_POM, 'Capt_eff_Gelatinous_POM', '-', 'Capture efficiency of GEL on POM', default=0.0_rk) 
   call self%get_parameter(self%DOXsatmort, 'DOXsatmort', 'mmolO2 m-3', 'Percentage of saturation where metabolic respiration is half the one under oxygen satyrated conditions', default=7.8125_rk) 
   call self%get_parameter(self%efficiency_growth_Gelatinous, 'efficiency_growth_Gelatinous', '-', 'Part of the assimil. food used for GEL growth', default=0.2_rk) 
   call self%get_parameter(self%expmortGelatinous, 'expmortGelatinous', 'd-2 (?)', 'Quadratic mortality of GEL', default=2.0_rk) 
   call self%get_parameter(self%HalfSatMort_Gelatinous, 'HalfSatMort_Gelatinous', 'mmolC m-3', '', default=0.0_rk) 
   call self%get_parameter(self%MaxgrazingrateGelatinous, 'MaxgrazingrateGelatinous', 'd-1', 'Maximum grazing rate of GEL', default=0.3_rk) 
   call self%get_parameter(self%Mortanoxic, 'Mortanoxic', 'd-1', 'Mortality rate in anoxia', default=0.25_rk) 
   call self%get_parameter(self%NCrGelatinous, 'NCrGelatinous', 'molN molC-1', 'N:C molar ratio in GEL', default=0.25_rk) 
   call self%get_parameter(self%NLin_Mort_Gelatinous, 'NLin_Mort_Gelatinous', '(?) d-1', 'Quadratic mortality rate of GEL', default=0.009_rk) 
   call self%get_parameter(self%OCr, 'OCr', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%Q10Gelatinous, 'Q10Gelatinous', '-', 'Temperature dependency for GEL', default=3.5_rk) 
   call self%get_parameter(self%threshold_feeding_Gelatinous, 'threshold_feeding_Gelatinous', '? mmolC m-3', 'Feeding threshold for GEL grazing', default=0_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_gel, 'GEL'  & 
         , 'mmol C m-3', 'Gelatinous omnivorous biomass' & 
         minimum=0.0e-7_rk)
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_mes, 'Mesozooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_mic, 'Microzooplakton biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_nhs, 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'Particulate organic nitrogen concentration', 'mmol N m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_TotalRespiration_Gel, 'TotalRespiration_Gel', '-', & 
      '-', output=output_instantaneous) 

   return 

99 call self%fatal_error('Gelatinous', 'Error reading namelist ulg_Gelatinous') 

   end subroutine initialize 


   ! Right hand sides of Gelatinous model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_Gelatinous), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DIC,DOX,MES,MIC,NHS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  GEL
      real(rk) ::   TotalRespiration_Gel
      real(rk) ::   
      real(rk) ::   C_ZOOAdjust	  ! mmol N m-3, Potential excretion rate necessary to keep the N/C ratio of Gelationous constant
      real(rk) ::   C_ZOOEgest	  ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   C_ZOOMort	  ! mmol C m-3, Zooplankton mortality flux in carbon
      real(rk) ::   C_ZOOResp	  ! flux, Zooplankton respiration
      real(rk) ::   FluxPrey_carbon	  ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   FluxPrey_nitrogen	  ! mmol N m-3, Flux of consummed preys in nitrogen
      real(rk) ::   grazing_carbonGelatinous	  ! mmol C m-3, Grazing in carbon by Gelatinous
      real(rk) ::   grazing_nitrogen	  ! mmol N m-3, Grazing in nitrogen by gelatinous
      real(rk) ::   N_ZOOAdjust	  ! mmol C m-3, Potential additional respiration flux to keep the N/C ratio of Gelationous constant
      real(rk) ::   N_ZOOEgest	  ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   N_ZOOMort	  ! mmol N m-3, Zooplankton mortality flux in nitrogen
      real(rk) ::   NCrfoodGelatinous	  ! mmol N mmol C-1, N/C ratio in food of gelatinous
      real(rk) ::   NCrzootest	  ! -, N/C ratio of the zooplankton before adjustment
      real(rk) ::   tf	  ! -, Temperature factor
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_gel,GEL)       ! Gelatinous omnivorous biomass
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_mes,MES)       ! Mesozooplakton biomass
   _GET_(self%id_mic,MIC)       ! Microzooplakton biomass
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_pho,PHO)       ! Phosphorus
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor(temp,Q10Gelatinous)
    
   ! Flux of consummed preys in carbon 
    FluxPrey_carbon = self%Capt_eff_Gelatinous_Flagellates*CFL+self%Capt_eff_Gelatinous_Emiliana*CEM+self%Capt_eff_Gelatinous_Diatoms*CDI+self%Capt_eff_Gelatinous_Microzoo*MIC+self%Capt_eff_Gelatinous_Mesozoo*MES+self%Capt_eff_Gelatinous_POM*POC
    
   ! Flux of consummed preys in nitrogen 
    FluxPrey_nitrogen = self%Capt_eff_Gelatinous_Flagellates*NFL+self%Capt_eff_Gelatinous_Emiliana*NEM+self%Capt_eff_Gelatinous_Diatoms*NDI+self%Capt_eff_Gelatinous_Microzoo*MIC*self%NCrMicroZoo + self%Capt_eff_Gelatinous_Mesozoo*MES*self%NCrMesoZoo+self%Capt_eff_Gelatinous_POM*PON
    
   ! Grazing rate in carbon 
             if (FluxPrey_carbon>threshold_feeding_Gelatinous) then 
    grazing_carbonGelatinous = tf*self%MaxgrazingrateGelatinous*(FluxPrey_carbon-self%threshold_feeding_Gelatinous)*GEL
             else 
    grazing_carbonGelatinous = 0.0
          endif 
    
   ! Grazing rate of nitrogen 
    NCrfoodGelatinous = FluxPrey_nitrogen/FluxPrey_carbon
    grazing_nitrogen = grazing_carbonGelatinous*NCrfoodGelatinous
    
   ! Egestion rate of zooplankton 
    C_ZOOEgest = ZOOEgestion(Ass_Eff_Gelatinous,grazing_carbonGelatinous)
    N_ZOOEgest = C_ZOOEgest*NCrfoodGelatinous
    
   ! Zooplankton respiration  
    C_ZOOResp = GELRespiration(tf,Ass_Eff_Gelatinous,efficiency_growth_Gelatinous,basal_Resp_Gelatinous,trb(I,J,K,GEL),grazing_carbonGelatinous)
    
   ! Zooplankton mortality rate 
    C_ZOOMort = Mortality_consument(HalfSatMort_Gelatinous,NLin_Mort_Gelatinous,DOXsatmort,Mortanoxic,tf,GEL,DOX)
    N_ZOOMort = C_ZOOMort * self%NCrGelatinous
    
   ! Computes the N/C fluxes necessary to conserve the N/C ratio 
    NCrzootest = (grazing_carbonGelatinous*NCrfoodGelatinous-N_ZOOEgest)/(grazing_carbonGelatinous-C_ZOOEgest-C_ZOOResp)
             if (NCrzootest> NCrGelatinous) then 
    C_ZOOAdjust = 0.0
    N_ZOOAdjust = (grazing_carbonGelatinous*NCrfoodGelatinous-N_ZOOEgest)-(grazing_carbonGelatinous-C_ZOOEgest-C_ZOOResp)*self%NCrGelatinous
             else 
    N_ZOOAdjust = 0.0
          endif 
    
   ! Carbon content increases by intake of preys and decreases by egestion, respiration, mortality, adjustement 
   _ADD_SOURCE_(self%id_gel,1.0*( grazing_carbonGelatinous)) 
   _ADD_SOURCE_(self%id_gel,-1.0*( C_ZOOEgest + C_ZOOResp + C_ZOOMort + C_ZOOAdjust)) 
    
   ! Particulate detritus pool if formed from non-assimilated grazed food and dead zooplatkton 
   _ADD_SOURCE_(self%id_poc,1.0*( C_ZOOEgest + C_ZOOMort)) 
   _ADD_SOURCE_(self%id_pon,1.0*( N_ZOOEgest + N_ZOOMort)) 
    
   ! Ammonium is excreted by zooplankton 
   _ADD_SOURCE_(self%id_nhs,1.0*( N_ZOOAdjust)) 
   _ADD_SOURCE_(self%id_pho,1.0*( N_ZOOAdjust*self%PNRedfield)) 
    
   ! Grazing on zoooplankton 
   _ADD_SOURCE_(self%id_mic,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_Microzoo*MIC)) 
   _ADD_SOURCE_(self%id_mes,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_Mesozoo*MES)) 
    
   ! Grazing on detritus 
   _ADD_SOURCE_(self%id_poc,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_POM*POC)) 
   _ADD_SOURCE_(self%id_pon,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_POM*PON)) 
    
   ! DOX decreases due to respiration of zooplankton 
   _ADD_SOURCE_(self%id_dox,-1.0*( (C_ZOOResp + C_ZOOAdjust)*self%OCr)) 
   _ADD_SOURCE_(self%id_dic,1.0*( C_ZOOResp + C_ZOOAdjust)) 
    
#ifdef biodiag1 
          TotalRespiration_Gel=C_ZOOResp + C_ZOOAdjust
#endif 

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_Gelatinous 
