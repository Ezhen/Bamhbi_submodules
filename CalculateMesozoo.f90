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
      type (type_state_variable_id)         :: id_bac,id_cdi,id_cem,id_cfl,id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_mic,id_ndi,id_nem,id_nfl,id_nhs,id_pho,id_poc,id_pon,id_sid
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_bac_to_ZOO,id_phy_to_ZOO,id_POC_to_ZOO,id_Totalrespiration_ZOO
      type (type_diagnostic_variable_id)    :: id_bac_to_ZOOIntegrated,id_phy_to_ZOOIntegrated,id_POC_to_ZOOIntegrated,id_Totalrespiration_zooIntegrated

!     Model parameters 
      real(rk)     :: Ass_Eff_OnCarbon, Ass_Eff_OnNitrogen, Capt_eff_MesoZoo_BAC
      real(rk)     :: Capt_eff_MesoZoo_Diatoms, Capt_eff_MesoZoo_Emiliana
      real(rk)     :: Capt_eff_MesoZoo_Flagellates, Capt_eff_MesoZoo_MesoZoo
      real(rk)     :: Capt_eff_MesoZoo_MicroZoo, Capt_eff_MesoZoo_POM
      real(rk)     :: DOXsatmort, efficiency_growth_MesoZoo, expmortMesoZoo
      real(rk)     :: Half_Saturation_MesoZoo, HalfSatMort_MesoZoo
      real(rk)     :: labilefraction, MaxgrazingrateMesoZoo, Messy_feeding_MesoZoo
      real(rk)     :: Mortanoxic, NCrBac, NCrMesoZoo, NLin_Mort_MesoZoo
      real(rk)     :: OCr, PNRedfield, Q10Zoo, SiNrDiatoms

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
                      NCrBac, NCrMesoZoo, NLin_Mort_MesoZoo, 	 & 
                      OCr, PNRedfield, Q10Zoo, SiNrDiatoms

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%Ass_Eff_OnCarbon, 'Ass_Eff_OnCarbon', '-', 'ZOO assimilation efficiency on C', default=0.64_rk) 
   call self%get_parameter(self%Ass_Eff_OnNitrogen, 'Ass_Eff_OnNitrogen', '-', 'ZOO assimilation efficiencies on N', default=0.77_rk) 
   call self%get_parameter(self%Capt_eff_MesoZoo_BAC, 'Capt_eff_MesoZoo_BAC', '-', 'Capture efficiency of MES on BAC', default=0.0_rk) 
   call self%get_parameter(self%Capt_eff_MesoZoo_Diatoms, 'Capt_eff_MesoZoo_Diatoms', '-', 'Capture efficiency of MES on DI', default=1.0_rk) 
   call self%get_parameter(self%Capt_eff_MesoZoo_Emiliana, 'Capt_eff_MesoZoo_Emiliana', '-', 'Capture efficiency of MES on EM', default=0.4_rk) 
   call self%get_parameter(self%Capt_eff_MesoZoo_Flagellates, 'Capt_eff_MesoZoo_Flagellates', '-', 'Capture efficiency of MES on FL', default=0.4_rk) 
   call self%get_parameter(self%Capt_eff_MesoZoo_MesoZoo, 'Capt_eff_MesoZoo_MesoZoo', '-', 'Capture efficiency of MES on MES', default=0.0_rk) 
   call self%get_parameter(self%Capt_eff_MesoZoo_MicroZoo, 'Capt_eff_MesoZoo_MicroZoo', '-', 'Capture efficiency of MES on MIC', default=1.0_rk) 
   call self%get_parameter(self%Capt_eff_MesoZoo_POM, 'Capt_eff_MesoZoo_POM', '-', 'Capture efficiency of MES on POM', default=0.8_rk) 
   call self%get_parameter(self%DOXsatmort, 'DOXsatmort', 'mmolO2 m-3', 'Percentage of saturation where metabolic respiration is half the one under oxygen satyrated conditions', default=7.8125_rk) 
   call self%get_parameter(self%efficiency_growth_MesoZoo, 'efficiency_growth_MesoZoo', '-', 'MES net growth efficiency on C', default=0.8_rk) 
   call self%get_parameter(self%expmortMesoZoo, 'expmortMesoZoo', '(?°', '? mortality rate of MES', default=2.0_rk) 
   call self%get_parameter(self%Half_Saturation_MesoZoo, 'Half_Saturation_MesoZoo', 'mmolC m-3', 'Half-saturation constant for MEC grazing', default=5.0_rk) 
   call self%get_parameter(self%HalfSatMort_MesoZoo, 'HalfSatMort_MesoZoo', 'mmolC m-3', '', default=1.0_rk) 
   call self%get_parameter(self%labilefraction, 'labilefraction', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%MaxgrazingrateMesoZoo, 'MaxgrazingrateMesoZoo', 'd-1', 'Maximum grazing rate of MES', default=1.2_rk) 
   call self%get_parameter(self%Messy_feeding_MesoZoo, 'Messy_feeding_MesoZoo', '-', 'Messy feeding fraction of MES grazing ', default=0.23_rk) 
   call self%get_parameter(self%Mortanoxic, 'Mortanoxic', 'd-1', 'Mortality rate in anoxia', default=0.25_rk) 
   call self%get_parameter(self%NCrBac, 'NCrBac', 'molN molC-1', 'N:C', default=0.196_rk) 
   call self%get_parameter(self%NCrMesoZoo, 'NCrMesoZoo', 'molN molC-1', 'N:C molar ratio in MES', default=0.21_rk) 
   call self%get_parameter(self%NLin_Mort_MesoZoo, 'NLin_Mort_MesoZoo', 'd-1', 'Maximum mortality rate of MES', default=0.3_rk) 
   call self%get_parameter(self%OCr, 'OCr', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%Q10Zoo, 'Q10Zoo', '-', 'Temperature factor Soetart et al., 2001', default=2.0_rk) 
   call self%get_parameter(self%SiNrDiatoms, 'SiNrDiatoms', 'molSi molN-1', 'Si:N ratio in DI', default=0.83_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_mes, 'MES'  & 
         , 'mmol C m-3', 'Mesozooplakton biomass' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
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
   class (type_ulg_Mesozoo), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  BAC,CDI,CEM,CFL,DCL,DCS,DIC,DNL,DNS,DOX,MIC,NDI,NEM,NFL,NHS,PHO,POC,PON,SID
      real(rk) ::  par,temp
      real(rk) ::  MES
      real(rk) ::   bac_to_ZOO,phy_to_ZOO,POC_to_ZOO,Totalrespiration_ZOO
      real(rk) ::   bac_to_ZOOIntegrated,phy_to_ZOOIntegrated,POC_to_ZOOIntegrated,Totalrespiration_zooIntegrated
      real(rk) ::   C_ZOOEgest	  ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   C_ZOOIntake	  ! mmol C m-3, Zooplankton carbon intake
      real(rk) ::   C_ZOOMessyfeeding	  ! mmol C m-3, Zooplankton messy feeding to the DOM in carbon
      real(rk) ::   C_ZOOMort	  ! mmol C m-3, Zooplankton mortality flux in carbon
      real(rk) ::   C_ZOOResp	  ! flux, Zooplankton respiration
      real(rk) ::   FluxPrey_carbon	  ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   FluxPrey_nitrogen	  ! mmol N m-3, Flux of consummed preys in nitrogen
      real(rk) ::   grazing_carbonMesoZoo	  ! mmol C m-3, Grazing in carbon by mesozooplankaton
      real(rk) ::   grazing_nitrogenZoo	  ! mmol N m-3, Grazing in nitrogen all zooplankaton
      real(rk) ::   NCrfoodZooref	  ! mol N mol C-1, Food threshold elemental ratio
      real(rk) ::   N_ZOOEgest	  ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   N_ZOOExcr	  ! mmol N m-3, Zooplankton excretion of ammonium
      real(rk) ::   N_ZOOIntake	  ! mmol N m-3, Zooplnkton nitrogen intake
      real(rk) ::   N_ZOOMessyfeeding	  ! mmol N m-3, Zooplankton messy feeding to the DOM in nitrogen
      real(rk) ::   N_ZOOMort	  ! mmol N m-3, Zooplankton mortality flux in nitrogen
      real(rk) ::   NCrfoodMesoZoo	  ! mmol N mmol C-1, N/C ratio in food of mesozooplankton
      real(rk) ::   tf	  ! -, Temperature factor
      real(rk) ::   ZOOGrowth	  ! mmol C m-3, Zooplankton growth flux
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_mes,MES)       ! Mesozooplakton biomass
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
    
    tf = Q10Factor(temp,Q10Zoo)
    
   ! Flux of consummed preys in carbon 
    FluxPrey_carbon = self%Capt_eff_MesoZoo_Flagellates*CFL + self%Capt_eff_MesoZoo_Emiliana*CEM + self%Capt_eff_MesoZoo_Diatoms*CDI + self%Capt_eff_MesoZoo_MicroZoo*MIC+ self%Capt_eff_MesoZoo_MesoZoo*MES + self%Capt_eff_MesoZoo_POM*POC + self%Capt_eff_MesoZoo_BAC*BAC
    
   ! Flux of consummed preys in nitrogen 
    FluxPrey_nitrogen = self%Capt_eff_MesoZoo_Flagellates*NFL + self%Capt_eff_MesoZoo_Emiliana*NEM + self%Capt_eff_MesoZoo_Diatoms*NDI + self%Capt_eff_MesoZoo_MicroZoo*MIC*self%NCrMicroZoo + self%Capt_eff_MesoZoo_MesoZoo*MES*self%NCrMesoZoo+self%Capt_eff_MesoZoo_POM*PON + self%Capt_eff_MesoZoo_BAC*BAC*self%NCrBac
    
   ! N:C molar ratio of the consumed food 
    NCrfoodMesoZoo=FluxPrey_nitrogen/FluxPrey_carbon
    
   ! Grazing rate in carbon 
    grazing_carbonMesoZoo = tf*self%MaxgrazingrateMesoZoo*Michaelis(FluxPrey_carbon,Half_Saturation_MesoZoo)*MES
    
   ! Grazing rate of zooplankton on nitrogen 
    grazing_nitrogenZoo = grazing_carbonMesoZoo*NCrfoodMesoZoo
    
   ! Ingestion rate of zooplankton  
    C_ZOOIntake = grazing_carbonMesoZoo*(1. - self%Messy_feeding_MesoZoo)
    N_ZOOIntake = grazing_nitrogenZoo*(1. - self%Messy_feeding_MesoZoo)
    
   !Zooplankton messy feeding  
    C_ZOOMessyfeeding = grazing_carbonMesoZoo*self%Messy_feeding_MesoZoo
    N_ZOOMessyfeeding = grazing_nitrogenZoo*self%Messy_feeding_MesoZoo
    
   ! Food threshold elemental ratio 
    NCrfoodZooref = self%Ass_Eff_OnCarbon*self%efficiency_growth_MesoZoo/self%Ass_Eff_OnNitrogen*self%NCrMesoZoo
    
   ! Growth and extrection as a function of N:C ratio of food 
    if (NCrfoodMesoZoo < NCrfoodZooref) then 
      ZOOGrowth = self%Ass_Eff_OnNitrogen*N_ZOOIntake/self%NCrMesoZoo
      N_ZOOExcr = 0.0
    else 
      ZOOGrowth = self%Ass_Eff_OnCarbon*self%efficiency_growth_MesoZoo*C_ZOOIntake
      N_ZOOExcr = C_ZOOIntake*(self%Ass_Eff_OnNitrogen*NCrfoodMesoZoo-self%Ass_Eff_OnCarbon*self%efficiency_growth_MesoZoo*self%NCrMesoZoo)
    endif    
    
   ! Zooplankton respiration 
    C_ZOOResp = self%Ass_Eff_OnCarbon*C_ZOOIntake-ZOOGrowth
    
   ! Zooplankton egestion 
    C_ZOOEgest = ZOOEgestion(Ass_Eff_OnCarbon,C_ZOOIntake)
    N_ZOOEgest = ZOOEgestion(Ass_Eff_OnNitrogen,N_ZOOIntake)
    
   phy_to_ZOO = phy_to_ZOO + grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Flagellates*CFL+ grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Emiliana *CEM+ grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms *CDI
   bac_to_ZOO = grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_BAC * BAC
   POC_to_ZOO = grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_POM * POC

   ! Zooplankton mortality rate  
    C_ZOOMort = Mortality_consument(HalfSatMort_MesoZoo,NLin_Mort_MesoZoo,expmortMesoZoo,DOXsatmort,Mortanoxic,tf,MES,DOX)
    N_ZOOMort  = C_ZOOMort * self%NCrMesoZoo
    
   ! Zoooplankton C increases by intake of preys (growth minus messy feeding) and decreases by egestion (where is it?), respiration, mortality 

   _ADD_SOURCE_(self%id_mes, ZOOGrowth - C_ZOOMort) 
    
   ! Particulate detritus pool if formed from non-assimilated grazed food and dead zooplatkton 
   ! Dissolved orgaqnic matter is formed by messy feeding taking into account labile/semi-labile partitioning; Ammonium is excreyed by zooplankton
   ! Grazing on phytoplankton, detritus and bacteria; When eating diatoms, zooplankton, ejects silicate as silicious_detritus

   _ADD_SOURCE_(self%id_poc, C_ZOOEgest + C_ZOOMort - grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_POM*POC) 
   _ADD_SOURCE_(self%id_pon, N_ZOOEgest + N_ZOOMort - grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_POM*PON) 
   _ADD_SOURCE_(self%id_dcl, self%labilefraction*C_ZOOMessyfeeding) 
   _ADD_SOURCE_(self%id_dnl, self%labilefraction*N_ZOOMessyfeeding)
   _ADD_SOURCE_(self%id_dcs, (1.0 - self%labilefraction)*C_ZOOMessyfeeding) 
   _ADD_SOURCE_(self%id_dns, (1.0 - self%labilefraction)*N_ZOOMessyfeeding)
   _ADD_SOURCE_(self%id_nhs, N_ZOOExcr) 
   _ADD_SOURCE_(self%id_pho, N_ZOOExcr*self%PNRedfield)
   _ADD_SOURCE_(self%id_cfl, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Flagellates*CFL)
   _ADD_SOURCE_(self%id_cem, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Emiliana*CEM)
   _ADD_SOURCE_(self%id_cdi, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms*CDI)
   _ADD_SOURCE_(self%id_nfl, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Flagellates*NFL) 
   _ADD_SOURCE_(self%id_nem, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Emiliana*NEM) 
   _ADD_SOURCE_(self%id_ndi, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms*NDI)
   _ADD_SOURCE_(self%id_sid, grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_Diatoms*NDI*self%SiNrDiatoms)
   _ADD_SOURCE_(self%id_mic, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_MicroZoo*MIC)
   _ADD_SOURCE_(self%id_bac, -grazing_carbonMesoZoo/FluxPrey_carbon*self%Capt_eff_MesoZoo_BAC*BAC)  
   _ADD_SOURCE_(self%id_dox, -self%OCr*C_ZOOResp) 
   _ADD_SOURCE_(self%id_dic, C_ZOOResp) 


   _SET_DIAGNOSTIC_(self%id_Totalrespiration_ZOO, C_ZOOResp)
   _SET_DIAGNOSTIC_(self%id_phy_to_ZOO, phy_to_ZOO)
   _SET_DIAGNOSTIC_(self%id_bac_to_ZOO, bac_to_ZOO)
   _SET_DIAGNOSTIC_(self%id_POC_to_ZOO, POC_to_ZOO)

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_Mesozoo 
