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
      real(rk) :: kinoxnhsdox
      real(rk) :: kinoxnhsodu
      real(rk) :: kinoxodudox
      real(rk) :: ksoxnhsdox
      real(rk) :: ksoxodudox
      real(rk) :: ksoxodunos
      real(rk) :: NODUr
      real(rk) :: NOsNHsr
      real(rk) :: ONoxnhsr
      real(rk) :: OODUr
      real(rk) :: Q10chem
      real(rk) :: Roxnhs
      real(rk) :: Roxnhsnos
      real(rk) :: Roxodu
      real(rk) :: Roxodunos

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

   real(rk)     :: kinoxnhsdox=8.0
   real(rk)     :: kinoxnhsodu=0.5
   real(rk)     :: kinoxodudox=5.0
   real(rk)     :: ksoxnhsdox=3.0
   real(rk)     :: ksoxodudox=1.0
   real(rk)     :: ksoxodunos=2.0
   real(rk)     :: NODUr=0.8
   real(rk)     :: NOsNHsr=0.6
   real(rk)     :: ONoxnhsr=2.0
   real(rk)     :: OODUr=1.0
   real(rk)     :: Q10chem=2.0
   real(rk)     :: Roxnhs=0.03/daytosecond
   real(rk)     :: Roxnhsnos=0.05/daytosecond
   real(rk)     :: Roxodu=0.1/daytosecond
   real(rk)     :: Roxodunos=0.05/daytosecond

   namelist /ulg_Chemical/ kinoxnhsdox, 	 & 
                      kinoxnhsodu, kinoxodudox, ksoxnhsdox, 	 & 
                      ksoxodudox, ksoxodunos, NODUr, NOsNHsr, 	 & 
                      ONoxnhsr, OODUr, Q10chem, Roxnhs, 	 & 
                      Roxnhsnos, Roxodu, Roxodunos

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%kinoxnhsdox, 'kinoxnhsdox', default=kinoxnhsdox) 
   call self%get_parameter(self%kinoxnhsodu, 'kinoxnhsodu', default=kinoxnhsodu) 
   call self%get_parameter(self%kinoxodudox, 'kinoxodudox', default=kinoxodudox) 
   call self%get_parameter(self%ksoxnhsdox, 'ksoxnhsdox', default=ksoxnhsdox) 
   call self%get_parameter(self%ksoxodudox, 'ksoxodudox', default=ksoxodudox) 
   call self%get_parameter(self%ksoxodunos, 'ksoxodunos', default=ksoxodunos) 
   call self%get_parameter(self%NODUr, 'NODUr', default=NODUr) 
   call self%get_parameter(self%NOsNHsr, 'NOsNHsr', default=NOsNHsr) 
   call self%get_parameter(self%ONoxnhsr, 'ONoxnhsr', default=ONoxnhsr) 
   call self%get_parameter(self%OODUr, 'OODUr', default=OODUr) 
   call self%get_parameter(self%Q10chem, 'Q10chem', default=Q10chem) 
   call self%get_parameter(self%Roxnhs, 'Roxnhs', default=Roxnhs) 
   call self%get_parameter(self%Roxnhsnos, 'Roxnhsnos', default=Roxnhsnos) 
   call self%get_parameter(self%Roxodu, 'Roxodu', default=Roxodu) 
   call self%get_parameter(self%Roxodunos, 'Roxodunos', default=Roxodunos) 

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
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  par,temp
      real(rk) ::  DOX,NHS,NOS,ODU
      real(rk) ::   ANAMMOX,Nitrification,Oxidation_by_nitrate,Oxidation_by_oxygen
      real(rk) ::   ANAMMOXIntegrated,NitrificationIntegrated,Oxidation_by_nitrateIntegrated,Oxidation_by_oxygenIntegrated
      real(rk) ::   Ammonium_Oxidation_rate_by_Nitrate	 + ! mmol N m-3 d-1, Ammonium oxidation rate by nitrate
      real(rk) ::   inhib	 + ! ?, Inhibiting function
      real(rk) ::   lim	 + ! mmol m-3 s-1, Limiting Michaelis-Menten function
      real(rk) ::   Nitrification_Rate	 + ! mmol N m-3 d-1, Ammonium oxidation rate by oxygen
      real(rk) ::   ODU_Oxidation_Rate_by_nitrate	 + ! ODU m-3 d-1, ODU oxidation rate by nitrate
      real(rk) ::   ODU_Oxidation_Rate_by_oxygen	 + ! O2 m-3 d-1, ODU oxidation rate by oxygen
      real(rk) ::   tf	 + ! -, Temperature factor
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
   ! TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor (temp,Q10chem)
   ! Compute ammonium oxidation by oxygen 
    lim = michaelis(DOX,ksoxnhsdox)
    Nitrification_Rate = tf*self%Roxnhs*lim
   ! Compute Ammonium oxidation by nitrate 
             !lim   = inhibition(ODU,kinoxnhsodu) ! REMOVE (POSSIBLY)
    lim   = michaelis(NOS,0.3_wp)
    inhib = inhibition(DOX,kinoxnhsdox)
    Ammonium_Oxidation_rate_by_Nitrate = tf*self%Roxnhsnos*inhib*lim*inhibition(ODU,kinoxnhsodu)  ! BUG TEST 2019-01-21
   ! Compute ODU oxidation by oxygen 
    lim = michaelis(DOX,ksoxodudox)
    ODU_Oxidation_Rate_by_oxygen = tf*self%Roxodu*lim
   ! Compute ODU oxidation by nitrate 
    lim = michaelis(NOS,ksoxodunos)
    inhib = inhibition(DOX,kinoxodudox)
    ODU_Oxidation_Rate_by_nitrate = tf*self%Roxodunos*inhib*lim
   ! ADJUSTING THE RATE OF CHANGE 
   ! oxidation of Ammonium by oxygen consumes oxygen and ammonium and produces nitrate 
   _ADD_SOURCE_(self%id_nos,1.0*( Nitrification_Rate*NHS)) 
   _ADD_SOURCE_(self%id_nhs,-1.0*( Nitrification_Rate*NHS)) 
   _ADD_SOURCE_(self%id_dox,-1.0*( Nitrification_Rate*NHS*self%ONoxnhsr)) 
   ! Oxidation of ammonium by nitrate consumes nitrate and ammonium and produces N2 which is definitely lost for the system 
#ifdef testcons 
    Ammonium_Oxidation_rate_by_Nitrate=0.
#endif 
   _ADD_SOURCE_(self%id_nos,-1.0*( Ammonium_Oxidation_rate_by_Nitrate*NHS*self%NOsNHsr)) 
   _ADD_SOURCE_(self%id_nhs,-1.0*( Ammonium_Oxidation_rate_by_Nitrate*NHS)) 
   !oxidation of ODU by oxygen consumes oxygen and ODU and produces oxydant not modelled 
   _ADD_SOURCE_(self%id_dox,-1.0*( ODU_Oxidation_Rate_by_oxygen*ODU*self%OODUr)) 
   _ADD_SOURCE_(self%id_odu,-1.0*( ODU_Oxidation_Rate_by_oxygen*ODU)) 
   ! oxidation of ODU by nitrate consumes nitrate and ODU and produces oxydant not modelled 
#ifdef testcons 
    ODU_Oxidation_Rate_by_nitrate=0.
#endif 
   _ADD_SOURCE_(self%id_odu,-1.0*( ODU_Oxidation_Rate_by_nitrate*ODU)) 
   _ADD_SOURCE_(self%id_nos,-1.0*( ODU_Oxidation_Rate_by_nitrate*ODU*self%NODUr)) 
   ! Diagnostics Compute auxillary variables 
#ifdef biodiag2 
          Oxidation_by_oxygen =  ODU_Oxidation_Rate_by_oxygen* ODU*self%OODUr + Nitrification_Rate*NHS*self%ONoxnhsr
          Oxidation_by_nitrate = ODU_Oxidation_Rate_by_nitrate*ODU*self%NODUr + Ammonium_Oxidation_rate_by_Nitrate*NHS*self%NOsNHsr
#endif 
#ifdef biodiag1 
          Nitrification=Nitrification_Rate*NHS
          ANAMMOX=Ammonium_Oxidation_rate_by_Nitrate*NHS
#endif 
           end if ! REMOVE (POSSIBLY)
         end do ! REMOVE (POSSIBLY)
       end do ! REMOVE (POSSIBLY)
     end do ! REMOVE (POSSIBLY)
   ! OUTPUT VARIABLES 
#ifdef biodiag2 
   _SET_DIAGNOSTIC_(self%id_Oxidation_by_nitrate, Oxidation_by_nitrate)
   _SET_DIAGNOSTIC_(self%id_Oxidation_by_oxygen, Oxidation_by_oxygen)
#endif 
#ifdef biodiag1 
   _SET_DIAGNOSTIC_(self%id_Nitrification, Nitrification)
   _SET_DIAGNOSTIC_(self%id_ANAMMOX, ANAMMOX)
#endif 
   ! Diagnostics Averaged over entire water column 
   _LOOP_END_

   end subroutine do

