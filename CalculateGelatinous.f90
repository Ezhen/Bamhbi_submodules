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
      type (type_state_variable_id)         :: id_agg,id_dic,id_dox,id_mes,id_mic,id_nhs,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_TotalRespiration_Gel
      type (type_diagnostic_variable_id)    :: 

!     Model parameters 
      real(rk) :: Ass_Eff_Gelatinous
      real(rk) :: basal_Resp_Gelatinous
      real(rk) :: Capt_eff_Gelatinous_Diatoms
      real(rk) :: Capt_eff_Gelatinous_Emiliana
      real(rk) :: Capt_eff_Gelatinous_Flagellates
      real(rk) :: Capt_eff_Gelatinous_Mesozoo
      real(rk) :: Capt_eff_Gelatinous_Microzoo
      real(rk) :: Capt_eff_Gelatinous_POM
      real(rk) :: DOXsatmort
      real(rk) :: efficiency_growth_Gelatinous
      real(rk) :: expmortGelatinous
      real(rk) :: HalfSatMort_Gelatinous
      real(rk) :: MaxgrazingrateGelatinous
      real(rk) :: Mortanoxic
      real(rk) :: NCrGelatinous
      real(rk) :: NLin_Mort_Gelatinous
      real(rk) :: OCr
      real(rk) :: PNRedfield
      real(rk) :: Q10Gelatinous
      real(rk) :: threshold_feeding_Gelatinous

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

   real(rk)     :: Ass_Eff_Gelatinous=0.75
   real(rk)     :: basal_Resp_Gelatinous=0.00/daytosecond
   real(rk)     :: Capt_eff_Gelatinous_Diatoms=0.0
   real(rk)     :: Capt_eff_Gelatinous_Emiliana=0.0
   real(rk)     :: Capt_eff_Gelatinous_Flagellates=0.0
   real(rk)     :: Capt_eff_Gelatinous_Mesozoo=1.0
   real(rk)     :: Capt_eff_Gelatinous_Microzoo=0.0
   real(rk)     :: Capt_eff_Gelatinous_POM=0.0
   real(rk)     :: DOXsatmort=7.8125
   real(rk)     :: efficiency_growth_Gelatinous=0.2
   real(rk)     :: expmortGelatinous=2.0
   real(rk)     :: HalfSatMort_Gelatinous=0.0
   real(rk)     :: MaxgrazingrateGelatinous=?
   real(rk)     :: Mortanoxic=0.25/daytosecond
   real(rk)     :: NCrGelatinous=0.25
   real(rk)     :: NLin_Mort_Gelatinous=0.009/daytosecond
   real(rk)     :: OCr=1.0
   real(rk)     :: PNRedfield=1.0/16.0
   real(rk)     :: Q10Gelatinous=3.5
   real(rk)     :: threshold_feeding_Gelatinous=0./12.

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
   call self%get_parameter(self%Ass_Eff_Gelatinous, 'Ass_Eff_Gelatinous', default=Ass_Eff_Gelatinous) 
   call self%get_parameter(self%basal_Resp_Gelatinous, 'basal_Resp_Gelatinous', default=basal_Resp_Gelatinous) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Diatoms, 'Capt_eff_Gelatinous_Diatoms', default=Capt_eff_Gelatinous_Diatoms) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Emiliana, 'Capt_eff_Gelatinous_Emiliana', default=Capt_eff_Gelatinous_Emiliana) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Flagellates, 'Capt_eff_Gelatinous_Flagellates', default=Capt_eff_Gelatinous_Flagellates) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Mesozoo, 'Capt_eff_Gelatinous_Mesozoo', default=Capt_eff_Gelatinous_Mesozoo) 
   call self%get_parameter(self%Capt_eff_Gelatinous_Microzoo, 'Capt_eff_Gelatinous_Microzoo', default=Capt_eff_Gelatinous_Microzoo) 
   call self%get_parameter(self%Capt_eff_Gelatinous_POM, 'Capt_eff_Gelatinous_POM', default=Capt_eff_Gelatinous_POM) 
   call self%get_parameter(self%DOXsatmort, 'DOXsatmort', default=DOXsatmort) 
   call self%get_parameter(self%efficiency_growth_Gelatinous, 'efficiency_growth_Gelatinous', default=efficiency_growth_Gelatinous) 
   call self%get_parameter(self%expmortGelatinous, 'expmortGelatinous', default=expmortGelatinous) 
   call self%get_parameter(self%HalfSatMort_Gelatinous, 'HalfSatMort_Gelatinous', default=HalfSatMort_Gelatinous) 
   call self%get_parameter(self%MaxgrazingrateGelatinous, 'MaxgrazingrateGelatinous', default=MaxgrazingrateGelatinous) 
   call self%get_parameter(self%Mortanoxic, 'Mortanoxic', default=Mortanoxic) 
   call self%get_parameter(self%NCrGelatinous, 'NCrGelatinous', default=NCrGelatinous) 
   call self%get_parameter(self%NLin_Mort_Gelatinous, 'NLin_Mort_Gelatinous', default=NLin_Mort_Gelatinous) 
   call self%get_parameter(self%OCr, 'OCr', default=OCr) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', default=PNRedfield) 
   call self%get_parameter(self%Q10Gelatinous, 'Q10Gelatinous', default=Q10Gelatinous) 
   call self%get_parameter(self%threshold_feeding_Gelatinous, 'threshold_feeding_Gelatinous', default=threshold_feeding_Gelatinous) 

   ! Register state variables 

   call self%register_state_variable(self%id_gel, 'GEL'  & 
         , 'mmol C m-3', 'Gelatinous omnivorous biomass' & 
         minimum=0.0e-7_rk)
   call self%register_state_dependency(self%id_agg, 'Aggregates', 'm-3') 
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
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  AGG,DIC,DOX,MES,MIC,NHS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  GEL
      real(rk) ::   TotalRespiration_Gel
      real(rk) ::   
      real(rk) ::   C_ZOOAdjust	 + ! mmol N m-3, Potential excretion rate necessary to keep the N/C ratio of Gelationous constant
      real(rk) ::   C_ZOOEgest	 + ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   C_ZOOMort	 + ! mmol C m-3, Zooplankton mortality flux in carbon
      real(rk) ::   C_ZOOResp	 + ! flux, Zooplankton respiration
      real(rk) ::   FluxPrey_carbon	 + ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   grazing_carbonGelatinous	 + ! mmol C m-3, Grazing in carbon by Gelatinous
      real(rk) ::   grazing_nitrogen	 + ! mmol N m-3 d-1, Grazing in nitrogen by gelatinous
      real(rk) ::   N_ZOOAdjust	 + ! mmol C m-3, Potential additional respiration flux to keep the N/C ratio of Gelationous constant
      real(rk) ::   N_ZOOEgest	 + ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   N_ZOOMort	 + ! mmol N m-3, Zooplankton mortality flux in nitrogen
      real(rk) ::   NCrfoodGelatinous	 + ! mmol N mmol C-1, N/C ratio in food of gelatinous
      real(rk) ::   NCrzootest	 + ! -, N/C ratio of the zooplankton before adjustment
      real(rk) ::   tf	 + ! -, Temperature factor
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_agg,AGG)       ! Aggregates
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
   ! TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor(temp,Q10Gelatinous)
   ! ZOOPLANKTON 
   ! Grazing rate of zooplankton (grazing_carbonZoo(I),NCrfoodZoo(I),mmolC/day) 
   CALL GELATINOUS_GRAZING_RATE(tf,MaxgrazingrateGelatinous,threshold_feeding_Gelatinous,Capt_eff_Gelatinous_Flagellates,Capt_eff_Gelatinous_Emiliana,Capt_eff_Gelatinous_Diatoms,Capt_eff_Gelatinous_Microzoo,Capt_eff_Gelatinous_Mesozoo,Capt_eff_Gelatinous_POM,GEL,grazing_carbonGelatinous,NCrfoodGelatinous,FluxPrey_carbon,i,j,k) ! REMOVE (POSSIBLY)
    grazing_nitrogen=grazing_carbonGelatinous*NCrfoodGelatinous
   ! Egestion rate of zooplankton(C_ZOOEgest,N_ZOOEgest,mmol/day) 
             CALL GELATINOUS_EGESTION_RATE(Ass_Eff_Gelatinous,grazing_carbonGelatinous,C_ZOOEgest) ! REMOVE (POSSIBLY)
    N_ZOOEgest = C_ZOOEgest*NCrfoodGelatinous
   ! Zooplankton respiration(C_ZOOResp, mmolC/m3/day) 
   CALL GELATINOUS_RESPIRATION_RATE(tf,Ass_Eff_Gelatinous,efficiency_growth_Gelatinous,basal_Resp_Gelatinous,trb(I,J,K,GEL),grazing_carbonGelatinous,C_ZOOResp) ! REMOVE (POSSIBLY)
   ! Zooplankton mortality rate (C_ZOOMort,N_ZOOMort, /day) 
   ! Attention do not forget to add OXYGEN 
   CALL GELATINOUS_MORTALITY_RATE(HalfSatMort_Gelatinous,NLin_Mort_Gelatinous,expmortGelatinous,DOXsatmort,Mortanoxic,tf,GEL,DOX,C_ZOOMort) ! REMOVE (POSSIBLY)
   ! Mortality in nitrogen units 
    N_ZOOMort  = C_ZOOMort * self%NCrGelatinous
   !Computes the nitrogen/carbon fluxes necessary to conserve the NCrzoo (C_ZOOAdjust,N_ZOOAdjust) 
   !The adjustment is realised through an extra-respiration or an excretion term on the model of the zooplankton in ERSEM 
   ! in this case excretion can be zero. We have to compute the N:C ratio of zooplankton with the computed fluxes 
   
             if (NCrzootest> NCrGelatinous) then ! REMOVE (POSSIBLY)
    C_ZOOAdjust=0
   
             else ! REMOVE (POSSIBLY)
    N_ZOOAdjust=0
   
          endif 
   ! ADJUSTING THE RATE OF CHANGE 
   ! zoooplankton C increases by intake of preys, 
   ! it decreases by egestion,respiration,mortality,predation,adjustement 
   _ADD_SOURCE_(self%id_gel,1.0*( grazing_carbonGelatinous)) 
   _ADD_SOURCE_(self%id_gel,-1.0*( C_ZOOEgest + C_ZOOResp + C_ZOOMort + C_ZOOAdjust)) 
   ! Detritus is formed by the non-assimilated zooplankton grazing, 
   ! when phytoplankton dies, when zooplankton dies 
   _ADD_SOURCE_(self%id_poc,1.0*( C_ZOOEgest+C_ZOOMort)) 
   _ADD_SOURCE_(self%id_pon,1.0*( N_ZOOEgest+N_ZOOMort)) 
   ! The number of aggregates increases with PON 
   ! if diatoms can form aggregate 
             SELECT CASE (SinkingVelocityType) ! REMOVE (POSSIBLY)
             CASE ('aggregation') ! REMOVE (POSSIBLY)
               !dPAGG(i,j,k) = dPAGG(i,j,k) + (N_ZOOEgest+N_ZOOMort)*AGG(i,j,k)/(PON(i,j,k)+PNS(i,j,k)) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,1.0*( (N_ZOOEgest+N_ZOOMort)*AGG/(PON))) 
             END SELECT ! REMOVE (POSSIBLY)
   ! Ammonium is excreyed by zooplankton 
   _ADD_SOURCE_(self%id_nhs,1.0*( N_ZOOAdjust)) 
   _ADD_SOURCE_(self%id_pho,1.0*( N_ZOOAdjust*self%PNRedfield)) 
   ! Grazing on zoooplankton 
   _ADD_SOURCE_(self%id_mic,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_Microzoo*MIC)) 
   _ADD_SOURCE_(self%id_mes,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_Mesozoo*MES)) 
   ! Grazing on detritus 
   _ADD_SOURCE_(self%id_poc,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_POM*POC)) 
   _ADD_SOURCE_(self%id_pon,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_POM*PON)) 
   ! POMNOS decreases by the grazing 
   ! If diatoms can form aggregate 
             SELECT CASE (SinkingVelocityType) ! REMOVE (POSSIBLY)
             CASE ('aggregation') ! REMOVE (POSSIBLY)
               !dDAGG(i,j,k)=dDAGG(i,j,k) + grazing_carbonGelatinous/FluxPrey_carbon*Capt_eff_Gelatinous_POM*PON(i,j,k)*AGG(i,j,k)/(PON(i,j,k)+PNS(i,j,k)) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,-1.0*( grazing_carbonGelatinous/FluxPrey_carbon*self%Capt_eff_Gelatinous_POM*AGG)) 
             END SELECT ! REMOVE (POSSIBLY)
   ! DOX decreases due to respiration of zooplankton 
   _ADD_SOURCE_(self%id_dox,-1.0*( (C_ZOOResp + C_ZOOAdjust)*self%OCr)) 
   _ADD_SOURCE_(self%id_dic,1.0*( (C_ZOOResp + C_ZOOAdjust))) 
   ! diagnostics 
#ifdef biodiag1 
          TotalRespiration_Gel=C_ZOOResp + C_ZOOAdjust
#endif 
           end if ! REMOVE (POSSIBLY)
         END DO ! REMOVE (POSSIBLY)
       END DO ! REMOVE (POSSIBLY)
     end do ! REMOVE (POSSIBLY)
   _LOOP_END_

   end subroutine do

