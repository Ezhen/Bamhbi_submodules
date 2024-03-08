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

   module fabm_ulg_chemical 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_chemical 
      type (type_state_variable_id)         :: id_dox,id_nhs,id_nos,id_odu
      type (type_dependency_id)             :: id_temp,id_salt,id_rho
      type (type_diagnostic_variable_id)    :: id_anammox,id_nitrification,id_oxidation_nitrate,id_oxidation_oxygen

!     Model parameters 
      real(rk)     :: ki_nhs_o2, ki_nhs_odu, ki_odu_o2, ks_nhs_o2
      real(rk)     :: ks_odu_nos, ks_odu_o2, q10_che, r_nos_nhs_oxid
      real(rk)     :: r_nos_odu_oxid, r_o2_nhs_nitr, r_o2_odu_oxid
      real(rk)     :: r_p_n_redfield, rox_nhs_nos, rox_nhs_o2
      real(rk)     :: rox_odu_nos, rox_odu_o2

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
   ! Initialise the Chemical model

   subroutine initialize(self,configunit)
   class (type_ulg_chemical), intent(inout), target :: self
   integer,                        intent(in)          :: configunit



   namelist /ulg_chemical/ ki_nhs_o2, 	 & 
                      ki_nhs_odu, ki_odu_o2, ks_nhs_o2, 	 & 
                      ks_odu_nos, ks_odu_o2, q10_che, 	 & 
                      r_nos_nhs_oxid, r_nos_odu_oxid, 	 & 
                      r_o2_nhs_nitr, r_o2_odu_oxid, 	 & 
                      r_p_n_redfield, rox_nhs_nos, rox_nhs_o2, 	 & 
                      rox_odu_nos, rox_odu_o2

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%ki_nhs_o2, 'ki_nhs_o2', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in NHS oxidation by NOS', default=8.0_rk) 
   call self%get_parameter(self%ki_nhs_odu, 'ki_nhs_odu', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in NHS oxidation by ODU', default=0.5_rk) 
   call self%get_parameter(self%ki_odu_o2, 'ki_odu_o2', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in ODU oxidation by NOS', default=5.0_rk) 
   call self%get_parameter(self%ks_nhs_o2, 'ks_nhs_o2', 'mmolO2 m-3', 'Half-sat. constant for O2 lim. in NHS oxidation by O2', default=3.0_rk) 
   call self%get_parameter(self%ks_odu_nos, 'ks_odu_nos', 'mmolN m-3', 'Half-sat. constant for NOS lim. in ODU oxidation by NOS', default=2.0_rk) 
   call self%get_parameter(self%ks_odu_o2, 'ks_odu_o2', 'mmolO2 m-3', 'Half-sat. constant for O2 lim. in ODU oxidation', default=1.0_rk) 
   call self%get_parameter(self%q10_che, 'q10_che', '-', 'Temperature factor for chemical processes', default=2.0_rk) 
   call self%get_parameter(self%r_nos_nhs_oxid, 'r_nos_nhs_oxid', 'molNOS molNHS-1', 'NOS:NHS ratio in NHS oxidation', default=0.6_rk) 
   call self%get_parameter(self%r_nos_odu_oxid, 'r_nos_odu_oxid', 'molNOS molODU-1', 'NOS:ODU ratio in ODU oxidation', default=0.8_rk) 
   call self%get_parameter(self%r_o2_nhs_nitr, 'r_o2_nhs_nitr', 'molO2 molNS-1', 'O2:NHS ratio in NHS oxidation in nitrification', default=2.0_rk) 
   call self%get_parameter(self%r_o2_odu_oxid, 'r_o2_odu_oxid', 'molO2 molODU-1', 'O2:ODU ratio in ODU oxidation', default=1.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%rox_nhs_nos, 'rox_nhs_nos', 'd-1', 'Maximum NHS oxidation rate by NOS', default=0.05_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rox_nhs_o2, 'rox_nhs_o2', 'd-1', 'Maximum NHS oxidation rate of NHS by O2', default=0.03_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rox_odu_nos, 'rox_odu_nos', 'd-1', 'Maximum ODU oxidation rate by NOS', default=0.05_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rox_odu_o2, 'rox_odu_o2', 'd-1', 'Maximum ODU oxidation rate by O2', default=0.1_rk, scale_factor=one_pr_day)

   ! Register state variables 

   call self%register_state_variable(self%id_dox, 'DOX', 'mmol O2 m-3', 'Dissolved oxygen concentration', minimum=0.0e-7_rk)
   call self%register_state_variable(self%id_nhs, 'NHS', 'mmol N m-3', 'Ammonium concentration',minimum=0.0e-7_rk)
   call self%register_state_variable(self%id_nos, 'NOS', 'mmol N m-3', 'Nitrate concentration', minimum=0.0e-7_rk)
   call self%register_state_variable(self%id_odu, 'ODU', 'mmol ODU m-3', 'Oxygen demand unit concentration', minimum=0.0e-7_rk)

    ! Register environmental dependencies 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 
   call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
   call self%register_dependency(self%id_rho,standard_variables%density)

    ! Add to aggregate variables 
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_dic)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nos)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nhs)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_nos, scale_factor=self%r_p_n_redfield)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_nhs, scale_factor=self%r_p_n_redfield)

    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_anammox, 'anammox', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_nitrification, 'nitrification', 'mmol N m-3 d-1', & 
      'Nitrification', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_oxidation_nitrate, 'oxidation_nitrate', 'mmol N m-3 d-1', & 
      'Oxudation by nitrate', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_oxidation_oxygen, 'oxidation_oxygen', 'mmol O m-3 d-1', & 
      'Oxydation by oxygen', output=output_instantaneous) 
   return 

99 call self%fatal_error('Chemical', 'Error reading namelist ulg_chemical') 

   end subroutine initialize 


   ! Right hand sides of Chemical model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_chemical), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  temp,salt,dens 
      real(rk) ::  DOX,NHS,NOS,ODU
      real(rk) ::   Rate_nitrification	  ! mmol N m-3 d-1, Ammonium oxidation rate by oxygen
      real(rk) ::   Rate_oxidation_NHS_NOS	  ! mmol N m-3 d-1, Ammonium oxidation rate by nitrate
      real(rk) ::   Rate_oxidation_ODU_NO3	  ! ODU m-3 d-1, ODU oxidation rate by nitrate
      real(rk) ::   Rate_oxidation_ODU_O2	  ! O2 m-3 d-1, ODU oxidation rate by oxygen
      real(rk) ::   tf	  ! -, Temperature factor

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_nos,NOS)       ! Nitrate concentration
   _GET_(self%id_odu,ODU)       ! Oxygen demand unit concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor(temp,self%q10_che)
    
   ! Compute ammonium oxidation by oxygen 
    Rate_nitrification = tf * self%rox_nhs_o2 * Michaelis(DOX,self%ks_nhs_o2)
    
   ! Compute Ammonium oxidation by nitrate 
    Rate_oxidation_NHS_NOS = tf * self%rox_nhs_nos * Inhibition(ODU,self%ki_nhs_odu) * Inhibition(DOX,self%ki_nhs_o2) * Michaelis(NOS,0.3)
    
   ! Compute ODU oxidation by oxygen 
    Rate_oxidation_ODU_O2 = tf * self%rox_odu_o2 * Michaelis(DOX,self%ks_odu_o2)
    
   ! Compute ODU oxidation by nitrate 
    Rate_oxidation_ODU_NO3 = tf * self%rox_odu_nos * Inhibition(DOX,self%ki_odu_o2) * Michaelis(NOS,self%ks_odu_nos)

   ! Oxidation of NHS by O2 produces NO3, oxidation of NHS by NO3 produces N2, oxidation of ODU by NO3 produces oxydant
   _ADD_SOURCE_(self%id_nos, Rate_nitrification * NHS - Rate_oxidation_NHS_NOS * NHS * self%r_nos_nhs_oxid - Rate_oxidation_ODU_NO3 * ODU * self%r_nos_odu_oxid) 
   _ADD_SOURCE_(self%id_nhs, -NHS * (Rate_nitrification + Rate_oxidation_NHS_NOS)) 
   _ADD_SOURCE_(self%id_dox, -Rate_nitrification * NHS * self%r_o2_nhs_nitr - Rate_oxidation_ODU_O2 * ODU * self%r_o2_odu_oxid) 
   _ADD_SOURCE_(self%id_odu, -ODU * (Rate_oxidation_ODU_O2 * ODU + Rate_oxidation_ODU_NO3)) 

   _SET_DIAGNOSTIC_(self%id_oxidation_nitrate, Rate_oxidation_ODU_NO3 * ODU * self%r_nos_odu_oxid + Rate_oxidation_NHS_NOS * NHS * self%r_nos_nhs_oxid)
   _SET_DIAGNOSTIC_(self%id_oxidation_oxygen, Rate_oxidation_ODU_O2 * ODU * self%r_o2_odu_oxid + Rate_nitrification * NHS * self%r_o2_nhs_nitr)
   _SET_DIAGNOSTIC_(self%id_nitrification, Rate_nitrification * NHS)
   _SET_DIAGNOSTIC_(self%id_anammox, Rate_oxidation_NHS_NOS * NHS)

   _LOOP_END_

   end subroutine do


   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_ulg_chemical),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

  real(rk)                 :: temp,wnd,salt,o2
  real(rk),parameter       :: secs_pr_day=86400.0_rk
  real(rk)                 :: Dox_sat

   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,temp)            ! local temperature
   _GET_(self%id_salt,salt)            ! local temperature
   _GET_(self%id_rho,rho)		! local water density
   _GET_HORIZONTAL_(self%id_wind,wind)

   _GET_(self%id_dox,dox)

   Dox_sat = o2_saturation(temp,salt,rho)
   Piston  = piston_wanninkhof(temp,wind)

   _ADD_SURFACE_FLUX_(self%id_o2, (Dox_sat - dox) * Piston )

   _HORIZONTAL_LOOP_END_
   end subroutine do_surface


   elemental real(rk) function piston_wanninkhof(t, w)
   ! Piston-Wanninkhof formula
  real(rk), intent :: t, w ! temp, wind
  real(wp) :: schmidt      ! The Schmidt number

  schmidt = 1953.4 - 128.00*t + 3.9918*t*t - 0.050091*t*t*t
  piston_wanninkhof = 0.266*w*w/SQRT(schmidt/600.) / 360000.

   end function piston_wanninkhof


   elemental real(rk) function o2_saturation(t, s, rho)
   ! Weiss formula
  real(rk), intent :: t, s, rho	! temp, salt, dens
  real(rk)         :: tk, sat_ox
  real(wp), parameter :: &
       a1 = -173.9894,      &
       a2 = 255.5907,       &
       a3 = 146.4813,       &
       a4 = -22.204,        &
       b1 = -0.037362,      &
       b2 = 0.016504,       &
       b3 = -0.0020564

  tk = (t + 273.15)*0.01

  ! Saturated oxygen concentration, in umol/kg & convert to in umol/dm3
  sat_ox = EXP(a1 + a2/tk + a3*LOG(tk) + a4*tk + s*(b1+b2*tk+b3*tk*tk))
  o2_saturation = sat_ox * rho * 0.001

   end function o2_saturation


   end module fabm_ulg_chemical 
