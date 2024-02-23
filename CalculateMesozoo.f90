#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DmesozooF90
!
! 1-D ecosystem model - Biological model of zooplankton
!
! Anderson, 1992,  Modelling the influence of food C:N ratio, and respiration on growth
! and nitrogen excretion in marine zooplankton and bacteria,
! Journal of Plankton Research, vol. 14, n 12, pp. 1645-1671, 1992
!
! Anderson and Hessen (1995), Carbon and Nitrogen limitation in marine copepods,
! Journal of Plankton Research, vol 17, n 2, pp. 317 - 331.
!
! Anderson amd Pondhaven (2003),Non-redfield carbon and nitrogen cycling in the Sarasso Sea :
! pelagic imbalances and export flux, in press in DSR
!
! Implementation: Marilaure Gregoire,                 NIOO-CEME
! Translation into FABM: Evgeny Ivanov, ULg / MAST
!
!--------------------------------------------------------------------
! Contains the pelagic submodel, for zooplankton (micro- and meso- zooplankton)
! This model is based on the Anderson and Hensen (1995) model described in JPR and also
! in Anderson and Ponhaven described in DSR I. This model is a nitrogen-carbon balanced model.
! It is assumed that the organic matter is composed of nitrogenous (proteins, amino acids)
! and non nitrogenous compounds (carbohydrate, lipids).
! The hypothesis is that zooplankton preferentially used nitrogeneous compounds for growth and carbon compounds for
! respiration (nore energy in non nitrogenous substrate). Based on this hypothesis and also on the
!fact that nitrogen and carbon have different assimilation efficiencies, the model estimates the zooplankton
!growth, respiration and excretion so as to maintain their internal ratio.
!--------------------------------------------------------------------*

   module fabm_ulg_Mesozoo 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Mesozoo 
      type (type_state_variable_id)         :: id_mes
      type (type_state_variable_id)         :: id_agg,id_bac,id_cdi,id_cem,id_cfl,id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_mic,id_ndi,id_nem,id_nfl,id_nhs,id_pho,id_poc,id_pon,id_sid
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_bac_to_ZOO,id_phy_to_ZOO,id_POC_to_ZOO,id_Totalrespiration_ZOO
      type (type_diagnostic_variable_id)    :: id_bac_to_ZOOIntegrated,id_phy_to_ZOOIntegrated,id_POC_to_ZOOIntegrated,id_Totalrespiration_zooIntegrated

!     Model parameters 
      real(rk) :: Ass_Eff_OnCarbon
      real(rk) :: Ass_Eff_OnNitrogen
      real(rk) :: Capt_eff_MesoZoo_BAC
      real(rk) :: Capt_eff_MesoZoo_Diatoms
      real(rk) :: Capt_eff_MesoZoo_Emiliana
      real(rk) :: Capt_eff_MesoZoo_Flagellates
      real(rk) :: Capt_eff_MesoZoo_MesoZoo
      real(rk) :: Capt_eff_MesoZoo_MicroZoo
      real(rk) :: Capt_eff_MesoZoo_POM
      real(rk) :: DOXsatmort
      real(rk) :: efficiency_growth_MesoZoo
      real(rk) :: expmortMesoZoo
      real(rk) :: Half_Saturation_MesoZoo
      real(rk) :: HalfSatMort_MesoZoo
      real(rk) :: labilefraction
      real(rk) :: MaxgrazingrateMesoZoo
      real(rk) :: Messy_feeding_MesoZoo
      real(rk) :: Mortanoxic
      real(rk) :: NCrMesoZoo
      real(rk) :: NLin_Mort_MesoZoo
      real(rk) :: OCr
      real(rk) :: PNRedfield
      real(rk) :: Q10Zoo
      real(rk) :: SiNrDiatoms

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
   ! Initialise the Mesozoo model

   subroutine initialize(self,configunit)
   class (type_ulg_Mesozoo), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   real(rk)     :: Ass_Eff_OnCarbon=0.64
   real(rk)     :: Ass_Eff_OnNitrogen=0.77
   real(rk)     :: Capt_eff_MesoZoo_BAC=0.0
   real(rk)     :: Capt_eff_MesoZoo_Diatoms=1.0
   real(rk)     :: Capt_eff_MesoZoo_Emiliana=0.4
   real(rk)     :: Capt_eff_MesoZoo_Flagellates=0.4
   real(rk)     :: Capt_eff_MesoZoo_MesoZoo=0.0
   real(rk)     :: Capt_eff_MesoZoo_MicroZoo=1.0
   real(rk)     :: Capt_eff_MesoZoo_POM=0.8
   real(rk)     :: DOXsatmort=7.8125
   real(rk)     :: efficiency_growth_MesoZoo=0.8
   real(rk)     :: expmortMesoZoo=2.0
   real(rk)     :: Half_Saturation_MesoZoo=5.0
   real(rk)     :: HalfSatMort_MesoZoo=1.0
   real(rk)     :: labilefraction=0.7
   real(rk)     :: MaxgrazingrateMesoZoo=1.2/daytosecond
   real(rk)     :: Messy_feeding_MesoZoo=0.23
   real(rk)     :: Mortanoxic=0.25/daytosecond
   real(rk)     :: NCrMesoZoo=1./4.67
   real(rk)     :: NLin_Mort_MesoZoo=0.3/daytosecond
   real(rk)     :: OCr=1.0
   real(rk)     :: PNRedfield=1.0/16.0
   real(rk)     :: Q10Zoo=2.0
   real(rk)     :: SiNrDiatoms=5./6.

   namelist /ulg_Mesozoo/ Ass_Eff_OnCarbon, 	 & 
                      Ass_Eff_OnNitrogen, 	 & 
                      Capt_eff_MesoZoo_BAC, 	 & 
                      Capt_eff_MesoZoo_Diatoms, 	 & 
                      Capt_eff_MesoZoo_Emiliana, 	 & 
                      Capt_eff_MesoZoo_Flagellates, 	 & 
                      Capt_eff_MesoZoo_MesoZoo, 	 & 
                      Capt_eff_MesoZoo_MicroZoo, 	 & 
                      Capt_eff_MesoZoo_POM, DOXsatmort, 	 & 
                      efficiency_growth_MesoZoo, 	 & 
                      expmortMesoZoo, Half_Saturation_MesoZoo, 	 & 
                      HalfSatMort_MesoZoo, labilefraction, 	 & 
                      MaxgrazingrateMesoZoo, 	 & 
                      Messy_feeding_MesoZoo, Mortanoxic, 	 & 
                      NCrMesoZoo, NLin_Mort_MesoZoo, OCr, 	 & 
                      PNRedfield, Q10Zoo, SiNrDiatoms

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%Ass_Eff_OnCarbon, 'Ass_Eff_OnCarbon', default=Ass_Eff_OnCarbon) 
   call self%get_parameter(self%Ass_Eff_OnNitrogen, 'Ass_Eff_OnNitrogen', default=Ass_Eff_OnNitrogen) 
   call self%get_parameter(self%Capt_eff_MesoZoo_BAC, 'Capt_eff_MesoZoo_BAC', default=Capt_eff_MesoZoo_BAC) 
   call self%get_parameter(self%Capt_eff_MesoZoo_Diatoms, 'Capt_eff_MesoZoo_Diatoms', default=Capt_eff_MesoZoo_Diatoms) 
   call self%get_parameter(self%Capt_eff_MesoZoo_Emiliana, 'Capt_eff_MesoZoo_Emiliana', default=Capt_eff_MesoZoo_Emiliana) 
   call self%get_parameter(self%Capt_eff_MesoZoo_Flagellates, 'Capt_eff_MesoZoo_Flagellates', default=Capt_eff_MesoZoo_Flagellates) 
   call self%get_parameter(self%Capt_eff_MesoZoo_MesoZoo, 'Capt_eff_MesoZoo_MesoZoo', default=Capt_eff_MesoZoo_MesoZoo) 
   call self%get_parameter(self%Capt_eff_MesoZoo_MicroZoo, 'Capt_eff_MesoZoo_MicroZoo', default=Capt_eff_MesoZoo_MicroZoo) 
   call self%get_parameter(self%Capt_eff_MesoZoo_POM, 'Capt_eff_MesoZoo_POM', default=Capt_eff_MesoZoo_POM) 
   call self%get_parameter(self%DOXsatmort, 'DOXsatmort', default=DOXsatmort) 
   call self%get_parameter(self%efficiency_growth_MesoZoo, 'efficiency_growth_MesoZoo', default=efficiency_growth_MesoZoo) 
   call self%get_parameter(self%expmortMesoZoo, 'expmortMesoZoo', default=expmortMesoZoo) 
   call self%get_parameter(self%Half_Saturation_MesoZoo, 'Half_Saturation_MesoZoo', default=Half_Saturation_MesoZoo) 
   call self%get_parameter(self%HalfSatMort_MesoZoo, 'HalfSatMort_MesoZoo', default=HalfSatMort_MesoZoo) 
   call self%get_parameter(self%labilefraction, 'labilefraction', default=labilefraction) 
   call self%get_parameter(self%MaxgrazingrateMesoZoo, 'MaxgrazingrateMesoZoo', default=MaxgrazingrateMesoZoo) 
   call self%get_parameter(self%Messy_feeding_MesoZoo, 'Messy_feeding_MesoZoo', default=Messy_feeding_MesoZoo) 
   call self%get_parameter(self%Mortanoxic, 'Mortanoxic', default=Mortanoxic) 
   call self%get_parameter(self%NCrMesoZoo, 'NCrMesoZoo', default=NCrMesoZoo) 
   call self%get_parameter(self%NLin_Mort_MesoZoo, 'NLin_Mort_MesoZoo', default=NLin_Mort_MesoZoo) 
   call self%get_parameter(self%OCr, 'OCr', default=OCr) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', default=PNRedfield) 
   call self%get_parameter(self%Q10Zoo, 'Q10Zoo', default=Q10Zoo) 
   call self%get_parameter(self%SiNrDiatoms, 'SiNrDiatoms', default=SiNrDiatoms) 

   ! Register state variables 

   call self%register_state_variable(self%id_mes, 'MES'  & 
         , 'mmol C m-3', 'Mesozooplakton biomass' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_dependency(self%id_agg, 'Aggregates', 'm-3') 
   call self%register_state_dependency(self%id_bac, 'Bacterial biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cdi, 'Diatom biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cem, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cfl, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcl, 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
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
   call self%register_diagnostic_variable(self%id_bac_to_ZOO, 'bac_to_ZOO', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_bac_to_ZOOIntegrated, 'bac_to_ZOOIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_phy_to_ZOO, 'phy_to_ZOO', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_phy_to_ZOOIntegrated, 'phy_to_ZOOIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_POC_to_ZOO, 'POC_to_ZOO', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_POC_to_ZOOIntegrated, 'POC_to_ZOOIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Totalrespiration_ZOO, 'Totalrespiration_ZOO', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Totalrespiration_zooIntegrated, 'Totalrespiration_zooIntegrated', '-', & 
      '-', output=output_instantaneous) 

   return 

99 call self%fatal_error('Mesozoo', 'Error reading namelist ulg_Mesozoo') 

   end subroutine initialize 


   ! Right hand sides of Mesozoo model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  AGG,BAC,CDI,CEM,CFL,DCL,DCS,DIC,DNL,DNS,DOX,MIC,NDI,NEM,NFL,NHS,PHO,POC,PON,SID
      real(rk) ::  par,temp
      real(rk) ::  MES
      real(rk) ::   bac_to_ZOO,phy_to_ZOO,POC_to_ZOO,Totalrespiration_ZOO
      real(rk) ::   bac_to_ZOOIntegrated,phy_to_ZOOIntegrated,POC_to_ZOOIntegrated,Totalrespiration_zooIntegrated
      real(rk) ::   C_ZOOEgest	 + ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   C_ZOOIntake	 + ! mmol C m-3, Zooplankton carbon intake
      real(rk) ::   C_ZOOMessyfeeding	 + ! mmol C m-3, Zooplankton messy feeding to the DOM in carbon
      real(rk) ::   C_ZOOMort	 + ! mmol C m-3, Zooplankton mortality flux in carbon
      real(rk) ::   C_ZOOResp	 + ! flux, Zooplankton respiration
      real(rk) ::   FluxPrey_carbon	 + ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   grazing_carbonMesoZoo	 + ! mmol C m-3, Grazing in carbon by mesozooplankaton
      real(rk) ::   grazing_nitrogenZoo	 + ! mmol N m-3, Grazing in nitrogen all zooplankaton
      real(rk) ::   N_ZOOEgest	 + ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   N_ZOOExcr	 + ! mmol N m-3, Zooplankton excretion of ammonium
      real(rk) ::   N_ZOOIntake	 + ! mmol N m-3, Zooplnkton nitrogen intake
      real(rk) ::   N_ZOOMessyfeeding	 + ! mmol N m-3, Zooplankton messy feeding to the DOM in nitrogen
      real(rk) ::   N_ZOOMort	 + ! mmol N m-3, Zooplankton mortality flux in nitrogen
      real(rk) ::   NCrfoodMesoZoo	 + ! mmol N mmol C-1, N/C ratio in food of mesozooplankton
      real(rk) ::   tf	 + ! -, Temperature factor
      real(rk) ::   ZOOGrowth	 + ! mmol C m-3, Zooplankton growth flux
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_agg,AGG)       ! Aggregates
   _GET_(self%id_bac,BAC)       ! Bacterial biomass
   _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
   _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
   _GET_(self%id_cfl,CFL)       ! Small flagellate biomass in carbon
   _GET_(self%id_dcl,DCL)       ! Labile detritus concentration in carbon
   _GET_(self%id_dcs,DCS)       ! Semi-labile detritus concentration in carbon
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dnl,DNL)       ! Labile detritus concentration in nitrogen
   _GET_(self%id_dns,DNS)       ! Semi-labile detritus concentration in nitrogen
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
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
   ! TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor(temp,Q10Zoo)
   ! ZOOPLANKTON 
   ! Grazing rate of zooplankton (grazing_carbonZoo(I),NCrfoodZoo(I),mmolC/day) 
   CALL GRAZING_RATE(tf,MaxgrazingrateMesoZoo,Half_Saturation_MesoZoo,Capt_eff_MesoZoo_Flagellates, Capt_eff_MesoZoo_Emiliana,Capt_eff_MesoZoo_Diatoms,Capt_eff_MesoZoo_MicroZoo,Capt_eff_MesoZoo_MesoZoo,Capt_eff_MesoZoo_POM,Capt_eff_MesoZoo_BAC,MES,grazing_carbonMesoZoo,NCrfoodMesoZoo,FluxPrey_carbon,i,j,k) ! REMOVE (POSSIBLY)
   ! Grazing rate of zooplankton on nitrogen(grazing_nitrogenZoo,NCrfoodZoo(I),mmolN/day) 
    grazing_nitrogenZoo = grazing_carbonMesoZoo*NCrfoodMesoZoo
   ! Ingestion rate of zooplankton C_ZOOIntake (mmolC/m3/day),N_ZOOIntake mmolN/m3/day) 
   ! Intake of N and C, N_ZOOIntake (mmolN/m3/day), C_ZOOIntake (mmolC/m3/day) are the sum 
   ! of grazing on phytoplanktonm zooplankton, detritus and bacteria, less messy feeding losses 
   ! altrhought in fact some of these losses occur after passage through the gut. 
    C_ZOOIntake = grazing_carbonMesoZoo*(1. - self%Messy_feeding_MesoZoo)
    N_ZOOIntake = grazing_nitrogenZoo*(1. - self%Messy_feeding_MesoZoo)
   !Zooplankton messy feeding C_ZOOMessyfeeding,N_ZOOMessyfeeding 
    C_ZOOMessyfeeding = grazing_carbonMesoZoo*self%Messy_feeding_MesoZoo
    N_ZOOMessyfeeding = grazing_nitrogenZoo*self%Messy_feeding_MesoZoo
   ! Growth rate of zooplankton computed according to Anderson and Hensen 
   CALL ZOO_GROWTH_RATE(Ass_Eff_OnNitrogen,Ass_Eff_OnCarbon,NCrMesoZoo,efficiency_growth_MesoZoo,NCrfoodMesoZoo,C_ZOOIntake,N_ZOOIntake,N_ZOOExcr,ZOOGrowth) ! REMOVE (POSSIBLY)
   ! Zooplankton respiration(C_ZOOResp, mmolC/m3/day) 
    C_ZOOResp=self%Ass_Eff_OnCarbon*C_ZOOIntake-ZOOGrowth
   ! Egestion rate of zooplankton(C_ZOOEgest,N_ZOOEgest,mmol/day) 
    C_ZOOEgest= (1-self%Ass_Eff_OnCarbon)*C_ZOOIntake
    N_ZOOEgest = (1-self%Ass_Eff_OnNitrogen)*N_ZOOIntake
   ! Zooplankton mortality rate (C_ZOOMort,N_ZOOMort, /day) 
   ! Attention do not forget to add OXYGEN 
   CALL MORTALITY_RATE(HalfSatMort_MesoZoo,NLin_Mort_MesoZoo,expmortMesoZoo,DOXsatmort,Mortanoxic,tf,MES,DOX,C_ZOOMort) ! REMOVE (POSSIBLY)
   ! Mortality in nitrogen units 
    N_ZOOMort  = C_ZOOMort * self%NCrMesoZoo
   ! ADJUSTING THE RATE OF CHANGE 
   ! zoooplankton C increases by intake of preys, 
             ! it decreases by egestion,respiration,mortality,predation (computed below) with  ZOOGrowth =Intake -respiration - egestion ! REMOVE (POSSIBLY)
             !and intake = grazing - messyfeeding ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_mes,1.0*( ZOOGrowth)) 
   _ADD_SOURCE_(self%id_mes,-1.0*( C_ZOOMort)) 
   ! Detritus is formed by the non-assimilated zooplankton grazing, 
   ! when phytoplankton dies (see the phytoplankton subroutine), when zooplankton dies 
   _ADD_SOURCE_(self%id_poc,1.0*( C_ZOOEgest+C_ZOOMort)) 
   _ADD_SOURCE_(self%id_pon,1.0*( N_ZOOEgest+N_ZOOMort)) 
   ! The number of aggregates increases with PON 
   ! if diatoms aggregate 
             SELECT CASE (SinkingVelocityType) ! REMOVE (POSSIBLY)
             CASE ('aggregation') ! REMOVE (POSSIBLY)
               !dPAGG(i,j,k) = dPAGG(i,j,k) + (N_ZOOEgest+N_ZOOMort)*AGG(i,j,k)/(PON+PNS(i,j,K) ) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,1.0*( (N_ZOOEgest+N_ZOOMort)*AGG/(PON))) 
#ifdef nanquest 
               if (isnan(dPAGG(I,J,K))) then ! REMOVE (POSSIBLY)
                 write(*,*) '** NAN QUEST ** in Calcmesozoo' ! REMOVE (POSSIBLY)
                 write(*,*) 'i,j,k,PONI(I,J,K),AGGI(I,J,K)',i,j,k,TRB(I,J,K,PON),trb(I,J,K,AGG) ! REMOVE (POSSIBLY)
                 call flush(6) ! REMOVE (POSSIBLY)
                 stop ! REMOVE (POSSIBLY)
            endif 
#endif 
             END SELECT ! REMOVE (POSSIBLY)
   ! Dissolved orgaqnic matter is formed by messy feeding, a part (labilefraction) increases the labile pool 
   ! while the remaining (1 -labilefrac) increases the semi-labile fraction. 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%labilefraction*C_ZOOMessyfeeding)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( self%labilefraction*N_ZOOMessyfeeding)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( (1.0 - self%labilefraction)*C_ZOOMessyfeeding)) 
   _ADD_SOURCE_(self%id_dns,1.0*( (1.0 - self%labilefraction)*N_ZOOMessyfeeding)) 
   ! Ammonium is excreyed by zooplankton 
   _ADD_SOURCE_(self%id_nhs,1.0*( N_ZOOExcr)) 
   _ADD_SOURCE_(self%id_pho,1.0*( N_ZOOExcr*self%PNRedfield)) 
   ! Grazing on phytoplankton 
   _ADD_SOURCE_(self%id_cfl,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Flagellates*CFL)) 
   _ADD_SOURCE_(self%id_cem,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Emiliana*CEM)) 
   _ADD_SOURCE_(self%id_cdi,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms*CDI)) 
   _ADD_SOURCE_(self%id_nfl,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Flagellates*NFL)) 
   _ADD_SOURCE_(self%id_nem,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Emiliana*NEM)) 
   _ADD_SOURCE_(self%id_ndi,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms*NDI)) 
   !When eating diatoms, mesozoo ejects silicate as silicious_Detritus 
   _ADD_SOURCE_(self%id_sid,1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms*NDI*self%SiNrDiatoms)) 
   ! Grazing on zoooplankton 
   _ADD_SOURCE_(self%id_mic,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_MicroZoo*MIC)) 
   !  Grazing on detritus 
   _ADD_SOURCE_(self%id_poc,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_POM*POC)) 
   _ADD_SOURCE_(self%id_pon,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_POM*PON)) 
   ! Grazing on Bacteria 
   _ADD_SOURCE_(self%id_bac,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_BAC*BAC)) 
   ! POMNOS decreases by the grazing 
   ! if diatoms form aggregate 
             SELECT CASE (SinkingVelocityType) ! REMOVE (POSSIBLY)
             CASE ('aggregation') ! REMOVE (POSSIBLY)
               !dDAGG(i,j,k) = dDAGG(i,j,k) + grazing_carbonMesoZoo/FluxPrey_carbon*Capt_eff_MesoZoo_POM*PONI(i,j,k)*AGG(i,j,k)/(PON(i,j,k)+PNS(i,j,K)) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,-1.0*( grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_POM*AGG)) 
#ifdef nanquest 
               if (isnan(dDAGG(I,J,K))) then ! REMOVE (POSSIBLY)
                 write (*,*) '** NAN QUEST ** in CalcMESO' ! REMOVE (POSSIBLY)
                 write (*,*) 'i,j,k,PONI(I,J,K),AGGI(I,J,K)',i,j,k,TRB(I,J,K,PON),trb(I,J,K,AGG) ! REMOVE (POSSIBLY)
                 call flush(6) ! REMOVE (POSSIBLY)
                 stop ! REMOVE (POSSIBLY)
             endif 
#endif 
             END SELECT ! REMOVE (POSSIBLY)
   ! DOX decreases due to respiration of zooplankton 
   _ADD_SOURCE_(self%id_dox,-1.0*( C_ZOOResp*self%OCr)) 
   _ADD_SOURCE_(self%id_dic,1.0*( C_ZOOResp)) 
   ! CO2 production and consumption 
             ! dDIC=dDIC+C_ZOOResp ! REMOVE (POSSIBLY)
   ! CALL UpdateCO2(I, - GrowthPHY) 
#ifdef biodiag1 
   !  Diagnostics Store the fluxes 
          Totalrespiration_ZOO = Totalrespiration_ZOO + C_ZOOResp
#endif 
#ifdef biodiagtrophic 
   !  Diagnostics Store the fluxes 
phy_to_ZOO = phy_to_ZOO + grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Flagellates*CFL+ grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Emiliana *CEM+ grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms *CDI
          bac_to_ZOO = bac_to_ZOO + grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_BAC * BAC
          POC_to_ZOO = POC_to_ZOO + grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_POM * POC
#endif 
           end if ! REMOVE (POSSIBLY)
         END DO ! REMOVE (POSSIBLY)
       END DO ! REMOVE (POSSIBLY)
     END DO ! REMOVE (POSSIBLY)
   ! OUTPUT VARIABLES 
#ifdef biodiag1 
   _SET_DIAGNOSTIC_(self%id_Totalrespiration_ZOO, Totalrespiration_ZOO)
#endif 
#ifdef biodiagtrophic 
   _SET_DIAGNOSTIC_(self%id_phy_to_ZOO, phy_to_ZOO)
   _SET_DIAGNOSTIC_(self%id_bac_to_ZOO, bac_to_ZOO)
   _SET_DIAGNOSTIC_(self%id_POC_to_ZOO, POC_to_ZOO)
#endif 
   _LOOP_END_

   end subroutine do

