#include "fabm_driver.h" 
 
!#########################################################################################
!                              3Dchemicals.F90
!
! 1-D ecosystem model - Chemical Processus
!
! Chemical equations of the suboxic and anoxic layers. All the reduced substances of the anoxic layer (Fe2+,Mn2+
! HS-) are gathered in a state variable called ODU (Oxygen demand units) as in Soetaert et al., 1996.
! Oxydants are considered to be not limiting. It is assumed that everywhere in the water columm,we will have an oxidant of particulate prganic matter
! (DOX,NO3,iron oxide, manganese oxyde,sulfate). This reaction produces ODU and ODU diffuses upwards to be oxydized by either DOX or NOs
! In reality, ODU can be also oxidyzed by MNO2 or FEOOH but these reactions are very fast and produce FE2+ and Mn2+ which are then oxydized by
! nitrate. The rate of the first eraction 9oxydation by mno2 or feooh is very high (about one order of magnitude higher than thge oxidation by nitrate and oxygen
! and we can consider that it is almost instantaneous and short circuiting this reaction. Alos, the oxidation of ODU will be considered as a direct
! oxidation by NOs and oxygen producing oxydant that are modelled. Heterotrophic bacteria are assumed to respire either oxygen,nitrate or ODU.
! according to environmental chemical conditions. The respiration of oxygen is computed by multiplying the
! respiration rate  of carbon by an oxygen li;itation function (Monod), the respiration of nitrate is computed
! by multiplying the respiration of carbon bu a nitrate limit funct (Monod) and a inhibition funvtion of
! oxygen concentration (inverse monod), the respiration of other oxydant is not computed since these oxydants are not
! modelled but the production of ODU is computed.
! Ammonium can be oxidized by oxygen producing nitrate : this is the nitrification rate. Nitrification depends
! on the oxygen concentration according to a Monod kinetics. The ammonium diffusing upwrds from the deep anoxic layer
! is oxidized by nitrate producing nitrogen gas that escapes to the atmosphere.
! Implementation: Marilaure Gregoire,                NIOO-CEME
! Translation into FABM : Evgeny Ivanov
!######################################################################

   module fabm_ulg_Chemical 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Chemical 
      type (type_state_variable_id)         :: id_dox,id_nhs,id_nos,id_odu
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_ANAMMOX,id_Nitrification,id_Oxidation_by_nitrate,id_Oxidation_by_oxygen
      type (type_diagnostic_variable_id)    :: id_ANAMMOXIntegrated,id_NitrificationIntegrated,id_Oxidation_by_nitrateIntegrated,id_Oxidation_by_oxygenIntegrated

!     Model parameters 
      real(rk)     :: kinoxnhsdox, kinoxnhsodu, kinoxodudox, ksoxnhsdox
      real(rk)     :: ksoxodudox, ksoxodunos, NODUr, NOsNHsr, ONoxnhsr
      real(rk)     :: OODUr, Q10chem, Roxnhs, Roxnhsnos, Roxodu
      real(rk)     :: Roxodunos

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
   ! Initialise the Chemical model

   subroutine initialize(self,configunit)
   class (type_ulg_Chemical), intent(inout), target :: self
   integer,                        intent(in)          :: configunit


   namelist /ulg_Chemical/ kinoxnhsdox, 	 & 
                      kinoxnhsodu, kinoxodudox, ksoxnhsdox, 	 & 
                      ksoxodudox, ksoxodunos, NODUr, NOsNHsr, 	 & 
                      ONoxnhsr, OODUr, Q10chem, Roxnhs, 	 & 
                      Roxnhsnos, Roxodu, Roxodunos

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%kinoxnhsdox, 'kinoxnhsdox', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in NHS oxidation by NOS', default=8.0_rk) 
   call self%get_parameter(self%kinoxnhsodu, 'kinoxnhsodu', '(?)', '', default=0.5_rk) 
   call self%get_parameter(self%kinoxodudox, 'kinoxodudox', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in ODU oxidation by NOS', default=5.0_rk) 
   call self%get_parameter(self%ksoxnhsdox, 'ksoxnhsdox', 'mmolO2 m-3', 'Half-sat. constant for O2 lim. in NHS oxidation by O2', default=3.0_rk) 
   call self%get_parameter(self%ksoxodudox, 'ksoxodudox', 'mmolO2 m-3', 'Half-sat. constant for O2 lim. in ODU oxidation', default=1.0_rk) 
   call self%get_parameter(self%ksoxodunos, 'ksoxodunos', 'mmolN m-3', 'Half-sat. constant for NOS lim. in ODU oxidation by NOS', default=2.0_rk) 
   call self%get_parameter(self%NODUr, 'NODUr', 'molNOS molODU-1', 'NOS:ODU ratio in ODU oxidation', default=0.8_rk) 
   call self%get_parameter(self%NOsNHsr, 'NOsNHsr', 'molNOS molNHS-1', 'NOS:NHS ratio in NHS oxidation', default=0.6_rk) 
   call self%get_parameter(self%ONoxnhsr, 'ONoxnhsr', 'molO2 molNS-1', 'O2:NHS ratio in NHS oxidation in nitrification', default=2.0_rk) 
   call self%get_parameter(self%OODUr, 'OODUr', 'molO2 molODU-1', 'O2:ODU ratio in ODU oxidation', default=1.0_rk) 
   call self%get_parameter(self%Q10chem, 'Q10chem', '-', 'Temperature factor for chemical processes', default=2.0_rk) 
   call self%get_parameter(self%Roxnhs, 'Roxnhs', 'd-1', 'Maximum NHS oxidation rate of NHS by O2', default=0.03_rk) 
   call self%get_parameter(self%Roxnhsnos, 'Roxnhsnos', 'd-1', 'Maximum NHS oxidation rate by NOS', default=0.05_rk) 
   call self%get_parameter(self%Roxodu, 'Roxodu', 'd-1', 'Maximum ODU oxidation rate by O2', default=0.1_rk) 
   call self%get_parameter(self%Roxodunos, 'Roxodunos', 'd-1', 'Maximum ODU oxidation rate by NOS', default=0.05_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_dox, 'DOX'  & 
         , 'mmol O2 m-3', 'Dissolved oxygen concentration' & 
         minimum=0.0e-7_rk)
   call self%register_state_variable(self%id_nhs, 'NHS'  & 
         , 'mmol N m-3', 'Ammonium concentration' & 
         minimum=0.0e-7_rk)
   call self%register_state_variable(self%id_nos, 'NOS'  & 
         , 'mmol N m-3', 'Nitrate concentration' & 
         minimum=0.0e-7_rk)
   call self%register_state_variable(self%id_odu, 'ODU'  & 
         , 'mmol ODU m-3', 'Oxygen demand unit concentration' & 
         minimum=0.0e-7_rk)

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_ANAMMOX, 'ANAMMOX', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_ANAMMOXIntegrated, 'ANAMMOXIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Nitrification, 'Nitrification', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_NitrificationIntegrated, 'NitrificationIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Oxidation_by_nitrate, 'Oxidation_by_nitrate', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Oxidation_by_nitrateIntegrated, 'Oxidation_by_nitrateIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Oxidation_by_oxygen, 'Oxidation_by_oxygen', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Oxidation_by_oxygenIntegrated, 'Oxidation_by_oxygenIntegrated', '-', & 
      '-', output=output_instantaneous) 

   return 

99 call self%fatal_error('Chemical', 'Error reading namelist ulg_Chemical') 

   end subroutine initialize 


   ! Right hand sides of Chemical model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_Chemical), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  par,temp
      real(rk) ::  DOX,NHS,NOS,ODU
      real(rk) ::   ANAMMOX,Nitrification,Oxidation_by_nitrate,Oxidation_by_oxygen
      real(rk) ::   ANAMMOXIntegrated,NitrificationIntegrated,Oxidation_by_nitrateIntegrated,Oxidation_by_oxygenIntegrated
      real(rk) ::   Ammonium_Oxidation_rate_by_Nitrate	  ! mmol N m-3 d-1, Ammonium oxidation rate by nitrate
      real(rk) ::   inhib	  ! -, Inhibiting function
      real(rk) ::   lim	  ! -, Limiting Michaelis-Menten function
      real(rk) ::   Nitrification_Rate	  ! mmol N m-3 d-1, Ammonium oxidation rate by oxygen
      real(rk) ::   ODU_Oxidation_Rate_by_nitrate	  ! ODU m-3 d-1, ODU oxidation rate by nitrate
      real(rk) ::   ODU_Oxidation_Rate_by_oxygen	  ! O2 m-3 d-1, ODU oxidation rate by oxygen
      real(rk) ::   tf	  ! -, Temperature factor
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_nos,NOS)       ! Nitrate concentration
   _GET_(self%id_odu,ODU)       ! Oxygen demand unit concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor(temp,Q10chem)
    
   ! Compute ammonium oxidation by oxygen 
    Nitrification_Rate = tf * self%Roxnhs * Michaelis(DOX,ksoxnhsdox)
    
   ! Compute Ammonium oxidation by nitrate 
    Ammonium_Oxidation_rate_by_Nitrate = tf*self%Roxnhsnos*inhibition(ODU,kinoxnhsodu) * Inhibition(DOX,kinoxnhsdox) * Michaelis(NOS,0.3_wp)
    
   ! Compute ODU oxidation by oxygen 
    ODU_Oxidation_Rate_by_oxygen = tf * self%Roxodu * Michaelis(DOX,ksoxodudox)
    
   ! Compute ODU oxidation by nitrate 
    ODU_Oxidation_Rate_by_nitrate = tf * self%Roxodunos * Inhibition(DOX,kinoxodudox) * Michaelis(NOS,ksoxodunos)

   Oxidation_by_oxygen =  ODU_Oxidation_Rate_by_oxygen* ODU*self%OODUr + Nitrification_Rate*NHS*self%ONoxnhsr
   Oxidation_by_nitrate = ODU_Oxidation_Rate_by_nitrate*ODU*self%NODUr + Ammonium_Oxidation_rate_by_Nitrate*NHS*self%NOsNHsr
   Nitrification = Nitrification_Rate*NHS
   ANAMMOX = Ammonium_Oxidation_rate_by_Nitrate*NHS

   ! Oxidation of NHS by O2 produces NO3, oxidation of NHS by NO3 produces N2, oxidation of ODU by NO3 produces oxydant
   _ADD_SOURCE_(self%id_nos, Nitrification_Rate*NHS - Ammonium_Oxidation_rate_by_Nitrate*NHS*self%NOsNHsr - ODU_Oxidation_Rate_by_nitrate*ODU*self%NODUr) 
   _ADD_SOURCE_(self%id_nhs, -NHS * (Nitrification_Rate + Ammonium_Oxidation_rate_by_Nitrate)) 
   _ADD_SOURCE_(self%id_dox, -Nitrification_Rate*NHS*self%ONoxnhsr - ODU_Oxidation_Rate_by_oxygen*ODU*self%OODUr) 
   _ADD_SOURCE_(self%id_odu, -ODU * (ODU_Oxidation_Rate_by_oxygen*ODU + ODU_Oxidation_Rate_by_nitrate)) 

   _SET_DIAGNOSTIC_(self%id_Oxidation_by_nitrate, Oxidation_by_nitrate)
   _SET_DIAGNOSTIC_(self%id_Oxidation_by_oxygen, Oxidation_by_oxygen)
   _SET_DIAGNOSTIC_(self%id_Nitrification, Nitrification)
   _SET_DIAGNOSTIC_(self%id_ANAMMOX, ANAMMOX)

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_Chemical 
