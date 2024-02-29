#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DBACTERIA.F90
!
! 1-D ecosystem model - Biological model of bacteria
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
! Implementation: Marilaure Gregoire,       NIOO-CEME
! Translation into FABM : Evgeny Ivanov                                                                 
! Contains the pelagic submodel, for bacteria described in Anderson (1992),JPR and also
! in Anderson and Ponhaven (2003), DSR I. This model is a nitrogen-carbon balanced model.
! It is assumed that the organic matter is composed of nitrogenous (proteins, amino acids)
! and non nitrogenous compounds (carbohydrate, lipids). Fixed C/N ratio is aasigned to bacteria
! The cycling of C and N by bacteria is described by elemental stoichiometry (Anderson,1992, Anderson and Williams (1998).
! Bacteria act as either remineralizers or consumers of ammonium, depending on the relative imbalance
! in the C/N ratios of their biomass relative to the DOM they consume, mediated by the C gross growth efficiency of utilization
! The hypothesis is that bacteria preferentially use nitrogeneous compounds for growth and carbon compounds for
! respiration (nore energy in non nitrogenous substrats). Bacteria growth, excretion and respiration are calculated
! from elemental stoichiometry. This method assumes that labile DOC amd DON are the primary growth substrates with ammonium supplementing
! DOM when the C/N ratio of DOM is high.
!######################################################################

   module fabm_ulg_Bacteria 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Bacteria 
      type (type_state_variable_id)         :: id_bac
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nhs,id_nos,id_odu,id_pho
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_bacteria_anoxrem,id_bacteria_oxygenconsumption,id_Bacteria_Respiration,id_denitrification,id_Uptake_DOCL
      type (type_diagnostic_variable_id)    :: id_bacteria_anoxremIntegrated,id_bacteria_oxygenconsumptionIntegrated,id_Bacteria_RespirationIntegrated,id_DenitrificationIntegrated,id_Uptake_DOCLIntegrated

!     Model parameters 
      real(rk) :: bactgrowtheff
      real(rk) :: csatamm
      real(rk) :: csatdocl
      real(rk) :: csatpo4
      real(rk) :: Halfsaturation_Iron
      real(rk) :: IronCsurf
      real(rk) :: kinanoxremdox
      real(rk) :: kinanoxremnos
      real(rk) :: kindenidox
      real(rk) :: ksdeninos
      real(rk) :: ksremindox
      real(rk) :: labilefraction
      real(rk) :: maxgrowthbac
      real(rk) :: mortbac
      real(rk) :: NCr
      real(rk) :: NCrBac
      real(rk) :: OCr
      real(rk) :: ODU_solid
      real(rk) :: ODUCr
      real(rk) :: Param1IronCurve
      real(rk) :: Param2IronCurve
      real(rk) :: PNRedfield
      real(rk) :: Q10bac

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
   ! Initialise the Bacteria model

   subroutine initialize(self,configunit)
   class (type_ulg_Bacteria), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   real(rk)     :: bactgrowtheff=0.17
   real(rk)     :: csatamm=0.5
   real(rk)     :: csatdocl=25.0
   real(rk)     :: csatpo4=0.5/16.0
   real(rk)     :: Halfsaturation_Iron=100.0
   real(rk)     :: IronCsurf=10.0
   real(rk)     :: kinanoxremdox=0.0005
   real(rk)     :: kinanoxremnos=0.0005
   real(rk)     :: kindenidox=0.5
   real(rk)     :: ksdeninos=0.3
   real(rk)     :: ksremindox=3.0
   real(rk)     :: labilefraction=0.7
   real(rk)     :: maxgrowthbac=13.3/daytosecond
   real(rk)     :: mortbac=0.05/daytosecond
   real(rk)     :: NCr=0.8
   real(rk)     :: NCrBac=1./5.1
   real(rk)     :: OCr=1.0
   real(rk)     :: ODU_solid=0.2
   real(rk)     :: ODUCr=1.0
   real(rk)     :: Param1IronCurve=25000.0
   real(rk)     :: Param2IronCurve=50.0
   real(rk)     :: PNRedfield=1.0/16.0
   real(rk)     :: Q10bac=2.0

   namelist /ulg_Bacteria/ bactgrowtheff, 	 & 
                      csatamm, csatdocl, csatpo4, 	 & 
                      Halfsaturation_Iron, IronCsurf, 	 & 
                      kinanoxremdox, kinanoxremnos, 	 & 
                      kindenidox, ksdeninos, ksremindox, 	 & 
                      labilefraction, maxgrowthbac, mortbac, 	 & 
                      NCr, NCrBac, OCr, ODU_solid, ODUCr, 	 & 
                      Param1IronCurve, Param2IronCurve, 	 & 
                      PNRedfield, Q10bac

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%bactgrowtheff, 'bactgrowtheff', default=bactgrowtheff) 
   call self%get_parameter(self%csatamm, 'csatamm', default=csatamm) 
   call self%get_parameter(self%csatdocl, 'csatdocl', default=csatdocl) 
   call self%get_parameter(self%csatpo4, 'csatpo4', default=csatpo4) 
   call self%get_parameter(self%Halfsaturation_Iron, 'Halfsaturation_Iron', default=Halfsaturation_Iron) 
   call self%get_parameter(self%IronCsurf, 'IronCsurf', default=IronCsurf) 
   call self%get_parameter(self%kinanoxremdox, 'kinanoxremdox', default=kinanoxremdox) 
   call self%get_parameter(self%kinanoxremnos, 'kinanoxremnos', default=kinanoxremnos) 
   call self%get_parameter(self%kindenidox, 'kindenidox', default=kindenidox) 
   call self%get_parameter(self%ksdeninos, 'ksdeninos', default=ksdeninos) 
   call self%get_parameter(self%ksremindox, 'ksremindox', default=ksremindox) 
   call self%get_parameter(self%labilefraction, 'labilefraction', default=labilefraction) 
   call self%get_parameter(self%maxgrowthbac, 'maxgrowthbac', default=maxgrowthbac) 
   call self%get_parameter(self%mortbac, 'mortbac', default=mortbac) 
   call self%get_parameter(self%NCr, 'NCr', default=NCr) 
   call self%get_parameter(self%NCrBac, 'NCrBac', default=NCrBac) 
   call self%get_parameter(self%OCr, 'OCr', default=OCr) 
   call self%get_parameter(self%ODU_solid, 'ODU_solid', default=ODU_solid) 
   call self%get_parameter(self%ODUCr, 'ODUCr', default=ODUCr) 
   call self%get_parameter(self%Param1IronCurve, 'Param1IronCurve', default=Param1IronCurve) 
   call self%get_parameter(self%Param2IronCurve, 'Param2IronCurve', default=Param2IronCurve) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', default=PNRedfield) 
   call self%get_parameter(self%Q10bac, 'Q10bac', default=Q10bac) 

   ! Register state variables 

   call self%register_state_variable(self%id_bac, 'BAC'  & 
         , 'mmol C m-3', 'Bacterial biomass' & 
         minimum=0.0e-7_rk)
   call self%register_state_dependency(self%id_dcl, 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_nhs, 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nos, 'Nitrate concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_odu, 'Oxygen demand unit concentration', 'mmol ODU m-3') 
   call self%register_state_dependency(self%id_pho, 'Phosphorus', 'mmol P m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_bacteria_anoxrem, 'bacteria_anoxrem', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_bacteria_anoxremIntegrated, 'bacteria_anoxremIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_bacteria_oxygenconsumption, 'bacteria_oxygenconsumption', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_bacteria_oxygenconsumptionIntegrated, 'bacteria_oxygenconsumptionIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Bacteria_Respiration, 'Bacteria_Respiration', 'mmol C m-3 d-1', & 
      'Respiration of bacteria', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Bacteria_RespirationIntegrated, 'Bacteria_RespirationIntegrated', 'mmol C m-3 d-1', & 
      'Respiration of bacteria (vertically-integrated)', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_denitrification, 'denitrification', 'mmol N m-3 d-1', & 
      'Denitrification rate', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_DenitrificationIntegrated, 'DenitrificationIntegrated', 'mmol N m-3 d-1', & 
      'Denitrification rate (vertically-integrated)', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Uptake_DOCL, 'Uptake_DOCL', 'mmol C m-3 d-1', & 
      'Uptake of labile DOC', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Uptake_DOCLIntegrated, 'Uptake_DOCLIntegrated', 'mmol C m-3 d-1', & 
      'Uptake of labile DOC (vertically-integrated)', output=output_instantaneous) 

   return 

99 call self%fatal_error('Bacteria', 'Error reading namelist ulg_Bacteria') 

   end subroutine initialize 


   ! Right hand sides of Bacteria model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,ODU,PHO
      real(rk) ::  par,temp
      real(rk) ::  BAC
      real(rk) ::   bacteria_anoxrem,bacteria_oxygenconsumption,Bacteria_Respiration,denitrification,Uptake_DOCL
      real(rk) ::   bacteria_anoxremIntegrated,bacteria_oxygenconsumptionIntegrated,Bacteria_RespirationIntegrated,DenitrificationIntegrated,Uptake_DOCLIntegrated
      real(rk) ::   BACExcr	 + ! mmol N m-3, Bacteria excretion flux of ammonium
      real(rk) ::   BACGrowth	 + ! mmol C m-3, Bacterial growth
      real(rk) ::   BACResp	 + ! mmol C m-3, Bacteria respiration flux
      real(rk) ::   bacteria_anoxrem_local	 + ! mmol ODU m-3, Bacterial anoxic remineralisation
      real(rk) ::   bacteria_oxygenconsumption_local	 + ! mmol O2 m-3, Bacterial respiration
      real(rk) ::   C_BACMort	 + ! mmol C m-3, Bacteria mortality flux in carbon
      real(rk) ::   denitrif	 + ! mmol N m-3, Denitrification flux
      real(rk) ::   Iron	 + ! ?, ?
      real(rk) ::   Limitation_By_Iron	 + ! ?, ?
      real(rk) ::   NutLim	 + ! ?, ?
      real(rk) ::   N_BACMort	 + ! mmol N m-3, Bacteria mortality flux in nitrogen
      real(rk) ::   tf	 + ! -, Temperature factor
      real(rk) ::   testratio	 + ! ?, ?
      real(rk) ::   Uptake_DOCL_local	 + ! mmol C m-3, Bacteria uptake of DOC
      real(rk) ::   Uptake_DONL_local	 + ! mmol N m-3, Bacteria uptake of DON
      real(rk) ::   Uptake_NHS_local	 + ! mmol N m-3, Bacteria uptake of ammonium
      real(rk) ::   Uptake_Potential_NHS	 + ! ?, ?
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_dcl,DCL)       ! Labile detritus concentration in carbon
   _GET_(self%id_dcs,DCS)       ! Semi-labile detritus concentration in carbon
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dnl,DNL)       ! Labile detritus concentration in nitrogen
   _GET_(self%id_dns,DNS)       ! Semi-labile detritus concentration in nitrogen
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_nos,NOS)       ! Nitrate concentration
   _GET_(self%id_odu,ODU)       ! Oxygen demand unit concentration
   _GET_(self%id_pho,PHO)       ! Phosphorus

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
   !TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor(temp,Q10bac)
   ! FORMER SUBROUTINE BAC GROWTH RATE  
   !Compute Iron Limitation 
    Iron=self%IronCsurf+self%Param1IronCurve/(self%self%Param2IronCurve*sqrt(2.*3.14159265359))*exp(-(depth-275.)**2/(2.*self%self%Param2IronCurve**2))
    Limitation_By_Iron = Iron/(Iron+Halfsaturation_Iron)															
    Uptake_DOCL_local = tf*self%maxgrowthbac*DCL/(DCL+self%csatdocl)*BAC
    Uptake_DONL_local = Uptake_DOCL_local*(DNL/DCL)
    NutLim = min(NHS/(self%csatamm+NHS)),PHO/(self%csatpo4+PHO)) 	
    Uptake_Potential_NHS = self%maxgrowthbac*NutLim*BAC*NCrBac					
   ! Test if bacterial growth is limited by carbon (DOClabile) or nitrogen (ammonium + DONLabile) 
   ! Define a threshold value and compare the Uptake_Potential_NHS, to this value 
    testratio = Uptake_DOCL_local*(Uptake_DONL_local/Uptake_DOCL_local - self%bactgrowtheff*NCrBac)		
   ! Growth rate of bacteria (BACGROWTH, in mmolC/m3/day) 
             if (Uptake_Potential_NHS > (-testratio)) then ! REMOVE (POSSIBLY)
   ! In this case we are in a situation of carbon limitation 
    BACGrowth = self%bactgrowtheff*Uptake_DOCL_local
   ! Growth rate computed taking into account the iron limitation 
    BACResp = Uptake_DOCL_local*(1.0 - self%bactgrowtheff)
   ! We have to test now, if NHS uptake is necessary to use all the DOCl 
               if (testratio > 0) then ! REMOVE (POSSIBLY)
   ! We are in case of remineralisation of ammonium through bacteria excretion and no net uptake of ammonium is necessary to use all the DOC 
    Uptake_NHS_local = 0
    BACExcr = testratio
               else ! REMOVE (POSSIBLY)
    Uptake_NHS_local = -testratio
    BACExcr = 0
            endif 
             else ! REMOVE (POSSIBLY)
   ! if we are in case of nitrogen limitation,it means that all the DON and the potential uptake of NHS is not sufficient to consume all the DOC 
    Uptake_NHS_local = Uptake_Potential_NHS
    BACGrowth = (Uptake_Potential_NHS + Uptake_DONL_local)/self%NCrBac
   ! Growth rate computed taking into account the iron limitation 
    BACExcr = 0
    BACResp = BACGrowth*(1.0/bactgrowthefficiency  - 1.0)
          endif 
   ! END OF FORMER SUBROUTINE BAC GROWTH RATE  
   ! Bacteria mortality rate (C_BACMort,N_BACMort, /day) 
    C_BACMort = self%mortbac*tf*BAC
   ! Mortality in nitrogen units 
    N_BACMort  = self%mortbac*tf*BAC*self%NCrBac
   ! Adjusting the rate of change: Bacteria C increases by intake of DOCl and decreases by mortality & predation 
   _ADD_SOURCE_(self%id_bac,1.0*( BACGrowth)) 
   _ADD_SOURCE_(self%id_bac,-1.0*( C_BACMort)) 
   ! Ammonium is excreyed by bacteria and can be taken up by bacteria 
   _ADD_SOURCE_(self%id_nhs,1.0*( BACExcr)) 
   _ADD_SOURCE_(self%id_nhs,-1.0*( Uptake_NHS_local)) 
   ! phosphore 
   _ADD_SOURCE_(self%id_pho,1.0*( BACExcr*self%PNRedfield)) 
   _ADD_SOURCE_(self%id_pho,-1.0*( Uptake_NHS_local*self%PNRedfield)) 
   ! As a result of respiration. Increases as a result of bactreial mortality which is fractionned bewteen the labile and semi-labile 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%labilefraction*C_BACMort)) 
   _ADD_SOURCE_(self%id_dcl,-1.0*( BACGrowth + BACResp)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( self%labilefraction*N_BACMort)) 
   _ADD_SOURCE_(self%id_dnl,-1.0*( Uptake_DONL_local)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( (1.0 - self%labilefraction)*C_BACMort)) 
   _ADD_SOURCE_(self%id_dns,1.0*( (1.0 - self%labilefraction)*N_BACMort)) 
    denitrif = BACResp*NOS/(NOS+self%ksdeninos)*(self%self%kindenidox/(DOX+self%self%kindenidox))*self%NCr
#ifdef testcons 
    denitrif=0
#endif 
    bacteria_oxygenconsumption_local=BACResp*DOX/(DOX+self%ksremindox) * self%OCr
    bacteria_anoxrem_local= BACResp*(self%self%kinanoxremnos/(NOS+self%self%kinanoxremnos))* (self%self%kinanoxremdox/(DOX+self%self%kinanoxremdox))*self%ODUCr
   ! Oxygen decreases due to bacterial respiration 
   _ADD_SOURCE_(self%id_dox,-1.0*( bacteria_oxygenconsumption_local)) 
   ! NOs decreases due to bacterial respiration 
   _ADD_SOURCE_(self%id_nos,-1.0*( denitrif)) 
   _ADD_SOURCE_(self%id_dic,1.0*( BACResp)) 
   ! ODU increases due to bacterial anoxic remineralisation 
   _ADD_SOURCE_(self%id_odu,1.0*( bacteria_anoxrem_local*(1. - Limitation_By_Iron*self%ODU_solid))) 
#ifdef biodiag1 
          Uptake_DOCL=Uptake_DOCL_local
          Bacteria_Respiration=BACResp
          denitrification=denitrif
#endif 
#ifdef biodiag2 
          bacteria_anoxrem=bacteria_anoxrem_local
          bacteria_oxygenconsumption= bacteria_oxygenconsumption_local
#endif 
        endif 
         end do ! REMOVE (POSSIBLY)
       end do ! REMOVE (POSSIBLY)
     end do ! REMOVE (POSSIBLY)
#ifdef biodiag1 
   _SET_DIAGNOSTIC_(self%id_Uptake_DOCL, Uptake_DOCL)
   _SET_DIAGNOSTIC_(self%id_Bacteria_Respiration, Bacteria_Respiration)
   _SET_DIAGNOSTIC_(self%id_denitrification, denitrification)
#endif 
#ifdef biodiag2 
   _SET_DIAGNOSTIC_(self%id_bacteria_oxygenconsumption, bacteria_oxygenconsumption)
   _SET_DIAGNOSTIC_(self%id_bacteria_anoxrem, bacteria_anoxrem)
#endif 
   _LOOP_END_

   end subroutine do

