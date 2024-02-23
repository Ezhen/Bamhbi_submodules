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
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dnl,id_dns
      type (type_state_variable_id)         :: id_agg,id_bac,id_dox,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: 
      type (type_diagnostic_variable_id)    :: 

!     Model parameters 
      real(rk) :: csatdocsl
      real(rk) :: hydPOCmax
      real(rk) :: hydPONmax
      real(rk) :: ksatOxygenHydrolysis
      real(rk) :: labilefraction
      real(rk) :: maxhydrDOCSL
      real(rk) :: Q10chem

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

   real(rk)     :: csatdocsl=417.0
   real(rk)     :: hydPOCmax=0.04/daytosecond
   real(rk)     :: hydPONmax=0.055/daytosecond
   real(rk)     :: ksatOxygenHydrolysis=2.7
   real(rk)     :: labilefraction=0.7
   real(rk)     :: maxhydrDOCSL=4./daytosecond
   real(rk)     :: Q10chem=2.0

   namelist /ulg_DOM/ csatdocsl, 	 & 
                      hydPOCmax, hydPONmax, 	 & 
                      ksatOxygenHydrolysis, labilefraction, 	 & 
                      maxhydrDOCSL, Q10chem

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%csatdocsl, 'csatdocsl', default=csatdocsl) 
   call self%get_parameter(self%hydPOCmax, 'hydPOCmax', default=hydPOCmax) 
   call self%get_parameter(self%hydPONmax, 'hydPONmax', default=hydPONmax) 
   call self%get_parameter(self%ksatOxygenHydrolysis, 'ksatOxygenHydrolysis', default=ksatOxygenHydrolysis) 
   call self%get_parameter(self%labilefraction, 'labilefraction', default=labilefraction) 
   call self%get_parameter(self%maxhydrDOCSL, 'maxhydrDOCSL', default=maxhydrDOCSL) 
   call self%get_parameter(self%Q10chem, 'Q10chem', default=Q10chem) 

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
   call self%register_state_dependency(self%id_agg, 'Aggregates', 'm-3') 
   call self%register_state_dependency(self%id_bac, 'Bacterial biomass', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_poc, 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'Particulate organic nitrogen concentration', 'mmol N m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 

   return 

99 call self%fatal_error('DOM', 'Error reading namelist ulg_DOM') 

   end subroutine initialize 


   ! Right hand sides of DOM model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  AGG,BAC,DOX,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  DCL,DCS,DNL,DNS
      real(rk) ::   
      real(rk) ::   
      real(rk) ::   DOCSLHYDR	 + ! mmol C m-3, Hydrolysis flux of DSL to DL in carbon
      real(rk) ::   DONSLHYDR	 + ! mmol N m-3, Hydrolysis flux of DSL to DL in nitrogen
      real(rk) ::   hydroPOC	 + ! mmol C m-3, Hydrolysis rate of POC
      real(rk) ::   hydroPOMLim	 + ! ?, ?
      real(rk) ::   hydroPON	 + ! mmol N m-3, Hydrolysis rate of PON
      real(rk) ::   POCHYDR	 + ! mmol C m-3, Hydrolysis flux of POC into DSL and DL 
      real(rk) ::   PONHYDR	 + ! mmol N m-3, Hydrolysis flux of PON into DSL and DL
      real(rk) ::   tf	 + ! -, Temperature factor
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_agg,AGG)       ! Aggregates
   _GET_(self%id_bac,BAC)       ! Bacterial biomass
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
   ! TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor(temp,Q10chem)
   ! Dissolved Organic matter divided in a labile part DOC and DON and a semi-labile part DOCSL and DONSL 
   ! the C:N ratio of the DOM is not constant 
   ! The DOM increases as a result of hydrolysis of the POM (POCHYDR,PONHYDR). This hydrolysis produced labile DOM with a coefficient 
   ! labilefrac and a semi-labile part with a coefficient (1-labilefrac) 
   ! The rate of hydrolysis is assumed to be proportional to the POC and PON concentration and the rate of hydrolysis is constant 
   ! but different for the POC and PON. It is assumed as in Anderson and Williams that the PON is faster decomposed 
   ! We have to consider that the decomposition of POM is not so active in anoxic conditions (Devol et al, 2001, LO,! Van Mooy et al., 2002 Geochimica et Cosmochimica) 
    hydroPOMLim=(DOX+1.0)/(DOX+self%ksatOxygenHydrolysis +1.0)
    hydroPOC=hydroPOMLim*self%hydPOCmax
    hydroPON=hydroPOMLim*self%hydPONmax
    POCHYDR = tf*hydroPOC*POC
    PONHYDR = tf*hydroPON*PON
             !hydroPCS=hydroPOMLim*hydPCSmax ! REMOVE (POSSIBLY)
             !hydroPNS=hydroPOMLim*hydPNSmax ! REMOVE (POSSIBLY)
             !PCSHYDR = tf*hydroPCS*PCSI(i,j,k) ! REMOVE (POSSIBLY)
             !PNSHYDR = tf*hydroPNS*PNSI(i,j,k) ! REMOVE (POSSIBLY)
   ! The semi-labile DOM is decomposed into labile dom according to the model proposed in Anderson and Williams (1998) 
   ! and Anderson and Pondhaven (2003). 
    DOCSLHYDR = tf*self%maxhydrDOCSL*DCS/(DCS + self%csatdocsl)*BAC
    DONSLHYDR = tf*self%maxhydrDOCSL*DCS/(DCS + self%csatdocsl)*BAC*DNS/DCS
   ! ADJUSTING THE RATE OF CHANGE 
   ! The hydrolysis of POM increases the pool of DOm and decreases the pool of POM 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%labilefraction*POCHYDR ! + self%labilefraction_fromS*PCSHYDR)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( self%labilefraction*PONHYDR ! + self%labilefraction_fromS*PNSHYDR)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( (1.0 - self%labilefraction)*POCHYDR !  +  (1.0 - self%labilefraction_fromS)*PCSHYDR)) 
   _ADD_SOURCE_(self%id_dns,1.0*( (1.0 - self%labilefraction)*PONHYDR !  +  (1.0 - self%labilefraction_fromS)*PNSHYDR)) 
   _ADD_SOURCE_(self%id_poc,-1.0*( POCHYDR)) 
   _ADD_SOURCE_(self%id_pon,-1.0*( PONHYDR)) 
            !dDPCS(i,j,k) = dDPCS(i,j,k) + PCSHYDR ! REMOVE (POSSIBLY)
            !dDPNS(i,j,k) = dDPNS(i,j,k) + PNSHYDR ! REMOVE (POSSIBLY)
#ifdef PCScheck 
             if  ( ((I.eq.14).and.(J.eq.4)).and.((K.eq.2).OR.(K.EQ.3))) Then ! REMOVE (POSSIBLY)
               write (*,*) "*********PELAGIC*********" ! REMOVE (POSSIBLY)
               write (*,*) " i , j ,k ", i , j ,k ! REMOVE (POSSIBLY)
               write (*,*) " DOX(i,j,k)",DOX," NHS(i,j,k)",NHS,"NOS(i,j,k)",NOS ! REMOVE (POSSIBLY)
               write (*,*) "*   *   *   *   *   *   *" ! REMOVE (POSSIBLY)
          endif 
#endif 
   ! The hydrolysis of the POM decreases the number of aggregates POMNOS 
   ! If diatoms can form aggregates : 
             SELECT CASE (SinkingVelocityType) ! REMOVE (POSSIBLY)
             CASE ('aggregation') ! REMOVE (POSSIBLY)
              !dDAGG(i,j,k) = dDAGG(i,j,k) + (PONHYDR+PNSHYDR)*AGG(i,j,k)/(PON+PNS(i,j,k)) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,-1.0*( (PONHYDR)*AGG/(PON))) 
#ifdef nanquest 
               if (isnan(dDAGG(I,J,K))) then ! REMOVE (POSSIBLY)
                 write (*,*) '** NAN QUEST ** in CalcDOM' ! REMOVE (POSSIBLY)
                 write (*,*) 'i,j,k,PON(I,J,K),AGG(I,J,K)',i,j,k,TRB(I,J,K,PON),trb(I,J,K,AGG) ! REMOVE (POSSIBLY)
                 write (*,*) 'PONHYDR',PONHYDR ! REMOVE (POSSIBLY)
                 call flush(6) ! REMOVE (POSSIBLY)
                 stop ! REMOVE (POSSIBLY)
            endif 
#endif 
             END SELECT ! REMOVE (POSSIBLY)
   ! The hydrolysis of the semi labile DOM increases the pool of labile DOM and decreases the pool of semi-labile DOM 
   _ADD_SOURCE_(self%id_dcl,1.0*( DOCSLHYDR)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( DONSLHYDR)) 
   _ADD_SOURCE_(self%id_dcs,-1.0*( DOCSLHYDR)) 
   _ADD_SOURCE_(self%id_dns,-1.0*( DONSLHYDR)) 
           end if ! REMOVE (POSSIBLY)
         end do ! REMOVE (POSSIBLY)
       end do ! REMOVE (POSSIBLY)
     end do ! REMOVE (POSSIBLY)
   ! OUTPUT VARIABLES 
   _LOOP_END_

   end subroutine do

