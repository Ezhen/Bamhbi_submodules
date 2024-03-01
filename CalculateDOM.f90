#include "fabm_driver.h" 
 
!#########################################################################################
!                              3Ddom.F90
!
! 1-D ecosystem model - Biological model of DOM
!
! Anderson, 1992,  Modelling the influence of food C:N ratio, and respiration on growth
! and nitrogen excretion in marine zooplankton and bacteria,
! Journal of Plankton Research, vol. 14, n 12, pp. 1645-1671, 1992
!
! Anderson and Williams,Modelling the Seasonal
! Cycle of Dissolved Organic Carbonj at Station E1 in the English Channel, Estuarine, Coastal
! and Shelf Science 1998, 46, 93-109.
!
! Anderson amd Pondhaven (2003),Non-redfield carbon and nitrogen cycling in the Sarasso Sea :
! pelagic imbalances and export flux, in press in DSR
!
! Implementation: Marilaure Gregoire,           NIOO-CEME
! Translation into FABM : Evgeny Ivanov,     ULg / MAST
!
!--------------------------------------------------------------------
! Contains the pelagic submodel, for DOM described in Anderson and Williams(1998),ECS and also
! in Anderson and Ponhaven (2003), DSR I. The DOM pool is divided between labile (DOCL, DONL)
! and semi-labile (DOCSL, DONSL) pools. each with dynamic C/N. The labile pool (turnover of hours to days)
! represents molecules which can be rapidly taken up by bacteria, whereas semi-labile pool (turnover of months)
!is characterized by complex bounding structures and must first be hydrolysed by exoenzymes before uptake by bacteria can occur
! (Anderson and Williams, 1998, 1999). DOM isproduced by messy feeding, phytoplankton and bacteria mortality (these terms are
! computed in the zopoplankton, phytoplankton and bacteria subroutines)
! and detrital turnover>
!--------------------------------------------------------------------

   module fabm_ulg_DOM 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_DOM 
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dnl,id_dns,id_poc,id_pon
      type (type_state_variable_id)         :: id_bac,id_dox
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: 
      type (type_diagnostic_variable_id)    :: 

!     Model parameters 
      real(rk)     :: csatdocsl, hydPOCmax, hydPONmax, ksatOxygenHydrolysis
      real(rk)     :: labilefraction, maxhydrDOCSL, Q10chem

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
   ! Initialise the DOM model

   subroutine initialize(self,configunit)
   class (type_ulg_DOM), intent(inout), target :: self
   integer,                        intent(in)          :: configunit


   namelist /ulg_DOM/ csatdocsl, 	 & 
                      hydPOCmax, hydPONmax, 	 & 
                      ksatOxygenHydrolysis, labilefraction, 	 & 
                      maxhydrDOCSL, Q10chem

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%csatdocsl, 'csatdocsl', 'mmolC m-3', 'Half-sat. constant for DSC uptake by BAC', default=417.0_rk) 
   call self%get_parameter(self%hydPOCmax, 'hydPOCmax', 'd-1', 'POC hydrolysis rate', default=0.04_rk) 
   call self%get_parameter(self%hydPONmax, 'hydPONmax', 'd-1', 'PON hydrolysis rate', default=0.055_rk) 
   call self%get_parameter(self%ksatOxygenHydrolysis, 'ksatOxygenHydrolysis', 'mmolO2 m-3', 'Half-saturation constant for oxic hydrolysis rate', default=2.7_rk) 
   call self%get_parameter(self%labilefraction, 'labilefraction', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%maxhydrDOCSL, 'maxhydrDOCSL', 'd-1', 'Maximum DSL hydrolysis', default=4.0_rk) 
   call self%get_parameter(self%Q10chem, 'Q10chem', '-', 'Temperature factor for chemical processes', default=2.0_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_dcl, 'DCL'  & 
         , 'mmol C m-3', 'Labile detritus concentration in carbon' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_dcs, 'DCS'  & 
         , 'mmol C m-3', 'Semi-labile detritus concentration in carbon' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_dnl, 'DNL'  & 
         , 'mmol N m-3', 'Labile detritus concentration in nitrogen' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_dns, 'DNS'  & 
         , 'mmol C m-3', 'Semi-labile detritus concentration in nitrogen' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_poc, 'POC'  & 
         , 'mmol C m-3', 'Particulate organic carbon concentration' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_pon, 'PON'  & 
         , 'mmol N m-3', 'Particulate organic nitrogen concentration' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_dependency(self%id_bac, 'Bacterial biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 

   return 

99 call self%fatal_error('DOM', 'Error reading namelist ulg_DOM') 

   end subroutine initialize 


   ! Right hand sides of DOM model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_DOM), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  BAC,DOX
      real(rk) ::  par,temp
      real(rk) ::  DCL,DCS,DNL,DNS,POC,PON
      real(rk) ::   
      real(rk) ::   
      real(rk) ::   DOCSLHYDR	  ! mmol C m-3, Hydrolysis flux of DSL to DL in carbon
      real(rk) ::   DONSLHYDR	  ! mmol N m-3, Hydrolysis flux of DSL to DL in nitrogen
      real(rk) ::   hydroPOMLim	  ! -, Limiting function on POM hydrolysis
      real(rk) ::   POCHYDR	  ! mmol C m-3, Hydrolysis flux of POC into DSL and DL 
      real(rk) ::   PONHYDR	  ! mmol N m-3, Hydrolysis flux of PON into DSL and DL
      real(rk) ::   tf	  ! -, Temperature factor
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_dcl,DCL)       ! Labile detritus concentration in carbon
   _GET_(self%id_dcs,DCS)       ! Semi-labile detritus concentration in carbon
   _GET_(self%id_dnl,DNL)       ! Labile detritus concentration in nitrogen
   _GET_(self%id_dns,DNS)       ! Semi-labile detritus concentration in nitrogen
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration
   _GET_(self%id_bac,BAC)       ! Bacterial biomass
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor(temp,Q10chem)
    
   ! Compute hydrolysis rates of particulate detritus 
    hydroPOMLim=(DOX+1.0)/(DOX+self%ksatOxygenHydrolysis +1.0)
    POCHYDR = tf*hydroPOMLim*self%hydPOCmax*POC
    PONHYDR = tf*hydroPOMLim*self%hydPONmax*PON
    
   ! Compute hydrolysis rates of semi-labile detritus 
    DOCSLHYDR = tf*self%maxhydrDOCSL*BAC*Michaelis((DCS,csatdocsl)
    DONSLHYDR = tf*self%maxhydrDOCSL*BAC*DNS/DCS*Michaelis(DCS,csatdocsl)
    
   ! Assign hydrolysis fluxes
   _ADD_SOURCE_(self%id_dcl, self%labilefraction*POCHYDR + DOCSLHYDR) 
   _ADD_SOURCE_(self%id_dnl, self%labilefraction*PONHYDR + DONSLHYDR) 
   _ADD_SOURCE_(self%id_dcs, (1.0 - self%labilefraction)*POCHYDR - DOCSLHYDR) 
   _ADD_SOURCE_(self%id_dns, (1.0 - self%labilefraction)*PONHYDR - DONSLHYDR) 
   _ADD_SOURCE_(self%id_poc,-POCHYDR) 
   _ADD_SOURCE_(self%id_pon,-PONHYDR) 
   
   _ADD_SOURCE_(self%id_dns,-1.0*( )) 

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_DOM 
