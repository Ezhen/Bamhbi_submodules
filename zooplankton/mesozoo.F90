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

   module fabm_ulg_mesozoo 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_mesozoo 
      type (type_state_variable_id)         :: id_mes
      type (type_state_variable_id)         :: id_bac,id_cdi,id_cem,id_cfl,id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_mic,id_ndi,id_nem,id_nfl,id_nhs,id_pho,id_poc,id_pon,id_sid
      type (type_dependency_id)             :: id_temp 
      type (type_diagnostic_variable_id)    :: id_zoo_on_bac,id_zoo_on_phy,id_zoo_on_poc,id_respiration_zoo

!     Model parameters 
      real(rk)     :: doxsatmort, eff_ass_zoo_c, eff_ass_zoo_n
      real(rk)     :: eff_gr_mes_c, eff_mes_bac, eff_mes_dia, eff_mes_emi
      real(rk)     :: eff_mes_fla, eff_mes_mes, eff_mes_mic, eff_mes_pom
      real(rk)     :: f_dl_dom, gmax_mes, ks_mort_mes, ks_prey_mec
      real(rk)     :: mess_prey_mes, mo_anox_pred, moexp_mes, momax_mes
      real(rk)     :: q10_zoo, r_n_c_bac, r_n_c_mes, r_n_c_mic, r_o2_c_resp
      real(rk)     :: r_p_n_redfield, r_si_n_dia


      contains 

      procedure :: initialize 
      procedure :: do
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

   contains
   ! Initialise the Mesozoo model

   subroutine initialize(self,configunit)
   class (type_ulg_mesozoo), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

!     Model parameters 
      real(rk)     :: doxsatmort, eff_ass_zoo_c, eff_ass_zoo_n
      real(rk)     :: eff_gr_mes_c, eff_mes_bac, eff_mes_dia, eff_mes_emi
      real(rk)     :: eff_mes_fla, eff_mes_mes, eff_mes_mic, eff_mes_pom
      real(rk)     :: f_dl_dom, gmax_mes, ks_mort_mes, ks_prey_mec
      real(rk)     :: mess_prey_mes, mo_anox_pred, moexp_mes, momax_mes
      real(rk)     :: q10_zoo, r_n_c_bac, r_n_c_mes, r_n_c_mic, r_o2_c_resp
      real(rk)     :: r_p_n_redfield, r_si_n_dia


   namelist /ulg_mesozoo/ doxsatmort, 	 & 
                      eff_ass_zoo_c, eff_ass_zoo_n, 	 & 
                      eff_gr_mes_c, eff_mes_bac, eff_mes_dia, 	 & 
                      eff_mes_emi, eff_mes_fla, eff_mes_mes, 	 & 
                      eff_mes_mic, eff_mes_pom, f_dl_dom, 	 & 
                      gmax_mes, ks_mort_mes, ks_prey_mec, 	 & 
                      mess_prey_mes, mo_anox_pred, moexp_mes, 	 & 
                      momax_mes, q10_zoo, r_n_c_bac, r_n_c_mes,  & 
                      r_n_c_mic, r_o2_c_resp, r_p_n_redfield, 	 & 
                      r_si_n_dia

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%doxsatmort, 'doxsatmort', 'mmolO2 m-3 (?)', 'Perc. of sat. where metabolic respiration is 1/2 the one under O2 sat.', default=7.8125_rk) 
   call self%get_parameter(self%eff_ass_zoo_c, 'eff_ass_zoo_c', '-', 'ZOO assimilation efficiency on C', default=0.64_rk) 
   call self%get_parameter(self%eff_ass_zoo_n, 'eff_ass_zoo_n', '-', 'ZOO assimilation efficiencies on N', default=0.77_rk) 
   call self%get_parameter(self%eff_gr_mes_c, 'eff_gr_mes_c', '-', 'MES net growth efficiency on C', default=0.8_rk) 
   call self%get_parameter(self%eff_mes_bac, 'eff_mes_bac', '-', 'Capture efficiency of MES on BAC', default=0.0_rk) 
   call self%get_parameter(self%eff_mes_dia, 'eff_mes_dia', '-', 'Capture efficiency of MES on DI', default=1.0_rk) 
   call self%get_parameter(self%eff_mes_emi, 'eff_mes_emi', '-', 'Capture efficiency of MES on EM', default=0.4_rk) 
   call self%get_parameter(self%eff_mes_fla, 'eff_mes_fla', '-', 'Capture efficiency of MES on FL', default=0.4_rk) 
   call self%get_parameter(self%eff_mes_mes, 'eff_mes_mes', '-', 'Capture efficiency of MES on MES', default=0.0_rk) 
   call self%get_parameter(self%eff_mes_mic, 'eff_mes_mic', '-', 'Capture efficiency of MES on MIC', default=1.0_rk) 
   call self%get_parameter(self%eff_mes_pom, 'eff_mes_pom', '-', 'Capture efficiency of MES on POM', default=0.8_rk) 
   call self%get_parameter(self%f_dl_dom, 'f_dl_dom', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%gmax_mes, 'gmax_mes', 'd-1', 'Maximum grazing rate of MES', default=1.2_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%ks_mort_mes, 'ks_mort_mes', 'mmolC m-3', '', default=1.0_rk) 
   call self%get_parameter(self%ks_prey_mec, 'ks_prey_mec', 'mmolC m-3', 'Half-saturation constant for MEC grazing', default=5.0_rk) 
   call self%get_parameter(self%mess_prey_mes, 'mess_prey_mes', '-', 'Messy feeding fraction of MES grazing ', default=0.23_rk) 
   call self%get_parameter(self%mo_anox_pred, 'mo_anox_pred', 'd-1', 'Mortality rate in anoxia', default=0.25_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%moexp_mes, 'moexp_mes', '-', 'Order of the non-linearity of mortality rate for MES', default=2.0_rk) 
   call self%get_parameter(self%momax_mes, 'momax_mes', 'd-1', 'Maximum mortality rate of MES', default=0.3_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%q10_zoo, 'q10_zoo', '-', 'Temperature factor Soetart et al., 2001', default=2.0_rk) 
   call self%get_parameter(self%r_n_c_bac, 'r_n_c_bac', 'molN molC-1', 'N:C', default=0.196_rk) 
   call self%get_parameter(self%r_n_c_mes, 'r_n_c_mes', 'molN molC-1', 'N:C molar ratio in MES', default=0.21_rk) 
   call self%get_parameter(self%r_n_c_mic, 'r_n_c_mic', 'molN molC-1', 'N:C molar ratio in MIC', default=0.18_rk) 
   call self%get_parameter(self%r_o2_c_resp, 'r_o2_c_resp', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%r_si_n_dia, 'r_si_n_dia', 'molSi molN-1', 'Si:N ratio in DI', default=0.83_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_mes, 'MES', 'mmol C m-3', 'Mesozooplakton biomass', minimum=0.0e-7_rk) 

   call self%register_state_dependency(self%id_bac, 'BAC', 'Bacterial biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cdi, 'CDI', 'Diatom biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cem, 'CEM', 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_cfl, 'CFL', 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcl, 'DCL', 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'DCS', 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'DIC', 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'DNL', 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'DNS', 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'DOX', 'Dissolved oxygen concentration', 'mmol O2 m-3') 
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
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_mes)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_mes, scale_factor=self%r_n_c_mes)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_mes, scale_factor=self%r_n_c_mes*self%r_p_n_redfield)

    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_zoo_on_bac, 'zoo_on_bac', 'mmol C m-3 d-1', & 
      'Zooplankton grazing on bacteria', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_zoo_on_phy, 'zoo_on_phy', 'mmol C m-3 d-1', & 
      'Zooplankton grazing on zooplankton', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_zoo_on_poc, 'zoo_on_poc', 'mmol C m-3 d-1', & 
      'Zooplankton grazing on POC', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_respiration_zoo, 'respiration_zoo', 'mmol C m-3 d-1', & 
      'Total Respiration of all zooplankton', output=output_instantaneous) 
   return 

99 call self%fatal_error('Mesozoo', 'Error reading namelist ulg_mesozoo') 

   end subroutine initialize 


   ! Right hand sides of Mesozoo model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_mesozoo), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  BAC,CDI,CEM,CFL,DCL,DCS,DIC,DNL,DNS,DOX,MIC,NDI,NEM,NFL,NHS,PHO,POC,PON,SID
      real(rk) ::  temp
      real(rk) ::  MES
      real(rk) ::   Egestion_C	  ! mmol C m-3, Zooplankton POM egestion in carbon
      real(rk) ::   Egestion_N	  ! mmol N m-3, Zooplankton POM egestion in nitrogen
      real(rk) ::   Excretion	  ! mmol N m-3, Zooplankton excretion of ammonium
      real(rk) ::   Grazing_C	  ! mmol C m-3, Grazing in carbon by microzooplankaton
      real(rk) ::   Grazing_N	  ! mmol N m-3, Grazing in nitrogen all zooplankaton
      real(rk) ::   Growth	  ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   Intake_C	  ! mmol C m-3, Zooplankton carbon intake
      real(rk) ::   Intake_N	  ! mmol N m-3, Zooplnkton nitrogen intake
      real(rk) ::   Messy_feeding_C	  ! mmol C m-3, Zooplankton messy feeding to the DOM in carbon
      real(rk) ::   Messy_feeding_N	  ! mmol N m-3, Zooplankton messy feeding to the DOM in nitrogen
      real(rk) ::   Mortality_C	  ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Mortality_N	  ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   Prey_C	  ! mmol C m-3, Flux of ingested preys in carbon 
      real(rk) ::   Prey_N	  ! mmol N m-3, Flux of consummed preys in nitrogen
      real(rk) ::   Ratio_N_C	  ! mol N mol C-1, N/C ratio in small flagellates
      real(rk) ::   Ratio_N_C_thr	  ! mol N mol C-1, Food threshold elemental ratio
      real(rk) ::   Respiration_C	  ! mmol C m-3, Zooplankton respiration flux
      real(rk) ::   tf	  ! -, Temperature factor

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_mes,MES)       ! Mesozooplakton biomass
   _GET_(self%id_bac,BAC)       ! Bacterial biomass
   _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
   _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
   _GET_(self%id_cfl,CFL)       ! Large flagellate biomass in carbon
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
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor(temp,self%q10_zoo)
    
   ! Flux of consummed preys in carbon 
    Prey_C = self%eff_mes_fla*CFL + self%eff_mes_emi*CEM + self%eff_mes_dia*CDI + self%eff_mes_mic*MIC+ self%eff_mes_mes*MES + self%eff_mes_pom*POC + self%eff_mes_bac*BAC
    
   ! Flux of consummed preys in nitrogen 
    Prey_N = self%eff_mes_fla*NFL + self%eff_mes_emi*NEM + self%eff_mes_dia*NDI + self%eff_mes_mic*MIC*self%r_n_c_mic + self%eff_mes_mes*MES*self%r_n_c_mes+self%eff_mes_pom*PON + self%eff_mes_bac*BAC*self%r_n_c_bac
    
   ! N:C molar ratio of the consumed food 
    Ratio_N_C=Prey_N/Prey_C
    
   ! Grazing rate in carbon 
    Grazing_C = tf*self%gmax_mes*Michaelis(Prey_C,self%ks_prey_mec)*MES
    
   ! Grazing rate of zooplankton on nitrogen 
    Grazing_N = Grazing_C*Ratio_N_C
    
   ! Ingestion rate of zooplankton  
    Intake_C = Grazing_C*(1.0 - self%mess_prey_mes)
    Intake_N = Grazing_N*(1.0 - self%mess_prey_mes)
    
   !Zooplankton messy feeding  
    Messy_feeding_C = Grazing_C*self%mess_prey_mes
    Messy_feeding_N = Grazing_N*self%mess_prey_mes
    
   ! Food threshold elemental ratio 
    Ratio_N_C_thr = self%eff_ass_zoo_c*self%eff_gr_mes_c/self%eff_ass_zoo_n*self%r_n_c_mes
    
   ! Growth and extrection as a function of N:C ratio of food 
    if (Ratio_N_C < Ratio_N_C_thr) then 
      Growth = self%eff_ass_zoo_n*Intake_N/self%r_n_c_mes
      Excretion = 0.0
    else 
      Growth = self%eff_ass_zoo_c*self%eff_gr_mes_c*Intake_C
      Excretion = Intake_C*(self%eff_ass_zoo_n*Ratio_N_C-self%eff_ass_zoo_c*self%eff_gr_mes_c*self%r_n_c_mes)
    endif    
    
   ! Zooplankton respiration 
    Respiration_C = self%eff_ass_zoo_c*Intake_C-Growth
    
   ! Zooplankton egestion 
    Egestion_C = zoo_egestion(self%eff_ass_zoo_c,Intake_C)
    Egestion_N = zoo_egestion(self%eff_ass_zoo_n,Intake_N)

   ! Zooplankton mortality rate  
    Mortality_C = mortality_consument(self%ks_mort_mes,self%momax_mes,self%moexp_mes,self%doxsatmort,self%mo_anox_pred,tf,MES,DOX)
    Mortality_N  = Mortality_C * self%r_n_c_mes
    
   ! Zoooplankton C increases by intake of preys (growth minus messy feeding) and decreases by egestion (where is it?), respiration, mortality 

   _ADD_SOURCE_(self%id_mes, Growth - Mortality_C) 
    
   ! Particulate detritus pool if formed from non-assimilated grazed food and dead zooplatkton 
   ! Dissolved orgaqnic matter is formed by messy feeding taking into account labile/semi-labile partitioning; Ammonium is excreyed by zooplankton
   ! Grazing on phytoplankton, detritus and bacteria; When eating diatoms, zooplankton, ejects silicate as silicious_detritus

   _ADD_SOURCE_(self%id_poc, Egestion_C + Mortality_C - Grazing_C/Prey_C*self%eff_mes_pom*POC) 
   _ADD_SOURCE_(self%id_pon, Egestion_N + Mortality_N - Grazing_C/Prey_C*self%eff_mes_pom*PON) 
   _ADD_SOURCE_(self%id_dcl, self%f_dl_dom*Messy_feeding_C) 
   _ADD_SOURCE_(self%id_dnl, self%f_dl_dom*Messy_feeding_N)
   _ADD_SOURCE_(self%id_dcs, (1.0 - self%f_dl_dom)*Messy_feeding_C) 
   _ADD_SOURCE_(self%id_dns, (1.0 - self%f_dl_dom)*Messy_feeding_N)
   _ADD_SOURCE_(self%id_nhs, Excretion) 
   _ADD_SOURCE_(self%id_pho, Excretion*self%r_p_n_redfield)
   _ADD_SOURCE_(self%id_cfl, -Grazing_C/Prey_C*self%eff_mes_fla*CFL)
   _ADD_SOURCE_(self%id_cem, -Grazing_C/Prey_C*self%eff_mes_emi*CEM)
   _ADD_SOURCE_(self%id_cdi, -Grazing_C/Prey_C*self%eff_mes_dia*CDI)
   _ADD_SOURCE_(self%id_nfl, -Grazing_C/Prey_C*self%eff_mes_fla*NFL) 
   _ADD_SOURCE_(self%id_nem, -Grazing_C/Prey_C*self%eff_mes_emi*NEM) 
   _ADD_SOURCE_(self%id_ndi, -Grazing_C/Prey_C*self%eff_mes_dia*NDI)
   _ADD_SOURCE_(self%id_sid, Grazing_C/Prey_C*self%eff_mes_dia*NDI*self%r_si_n_dia)
   _ADD_SOURCE_(self%id_mic, -Grazing_C/Prey_C*self%eff_mes_mic*MIC)
   _ADD_SOURCE_(self%id_bac, -Grazing_C/Prey_C*self%eff_mes_bac*BAC)  
   _ADD_SOURCE_(self%id_dox, -self%r_o2_c_resp*Respiration_C) 
   _ADD_SOURCE_(self%id_dic, Respiration_C) 


   _SET_DIAGNOSTIC_(self%id_respiration_zoo, Respiration_C)
   _SET_DIAGNOSTIC_(self%id_zoo_on_phy, Grazing_C/Prey_C*self%eff_mes_fla*CFL+ Grazing_C/Prey_C*self%eff_mes_emi *CEM+ Grazing_C/Prey_C*self%eff_mes_dia *CDI)
   _SET_DIAGNOSTIC_(self%id_zoo_on_bac, Grazing_C/Prey_C*self%eff_mes_bac * BAC)
   _SET_DIAGNOSTIC_(self%id_zoo_on_poc, Grazing_C/Prey_C*self%eff_mes_pom * POC)

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_mesozoo 
