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
      real(rk)     :: bactgrowtheff, csatamm, csatdocl, csatpo4
      real(rk)     :: Halfsaturation_Iron, IronCsurf, kinanoxremdox
      real(rk)     :: kinanoxremnos, kindenidox, ksdeninos, ksremindox
      real(rk)     :: labilefraction, maxgrowthbac, mortbac, NCr
      real(rk)     :: NCrBac, OCr, ODU_solid, ODUCr, Param1IronCurve
      real(rk)     :: Param2IronCurve, PNRedfield, Q10bac

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
   call self%get_parameter(self%bactgrowtheff, 'bactgrowtheff', '-', 'BAC gross growth efficiency on C', default=0.17_rk) 
   call self%get_parameter(self%csatamm, 'csatamm', 'mmolN m-3', 'Half-sat. constant for NHS uptake by BAC', default=0.5_rk) 
   call self%get_parameter(self%csatdocl, 'csatdocl', 'mmolC m-3', 'Half-sat. constant for DLC uptake by BAC', default=25.0_rk) 
   call self%get_parameter(self%csatpo4, 'csatpo4', 'mmolP m-3', 'Half-saturation constant for PO4 uptake by BAC', default=0.031_rk) 
   call self%get_parameter(self%Halfsaturation_Iron, 'Halfsaturation_Iron', 'mmolFe m-3', 'Half-sat. constant for iron lim. in solid ODU formation', default=100.0_rk) 
   call self%get_parameter(self%IronCsurf, 'IronCsurf', 'mmolFe m-3', 'Concentration of iron in surface water', default=10.0_rk) 
   call self%get_parameter(self%kinanoxremdox, 'kinanoxremdox', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in anoxic remineralization', default=0.0005_rk) 
   call self%get_parameter(self%kinanoxremnos, 'kinanoxremnos', 'mmolN m-3', 'Half-sat. constant for NOS inhibition in anoxic remineralization', default=0.0005_rk) 
   call self%get_parameter(self%kindenidox, 'kindenidox', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in denitrification', default=0.5_rk) 
   call self%get_parameter(self%ksdeninos, 'ksdeninos', 'mmolN m-3', 'Half-sat. constant for NOS lim. in denitrif.', default=0.3_rk) 
   call self%get_parameter(self%ksremindox, 'ksremindox', 'mmolO2 m-3', 'Half-sat. constant for O2 lim. in oxic min.', default=3.0_rk) 
   call self%get_parameter(self%labilefraction, 'labilefraction', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%maxgrowthbac, 'maxgrowthbac', 'd-1', 'Maximum labile DOC or NHS uptake', default=0.000154 (?)_rk) 
   call self%get_parameter(self%mortbac, 'mortbac', 'd-1', 'Bacteria natural mortality', default=0.05_rk) 
   call self%get_parameter(self%NCr, 'NCr', 'molN molC-1', 'N:C ratio of denitrification', default=0.8_rk) 
   call self%get_parameter(self%NCrBac, 'NCrBac', 'molN molC-1', 'N:C', default=0.196_rk) 
   call self%get_parameter(self%OCr, 'OCr', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%ODU_solid, 'ODU_solid', '-', 'Percentage of solid ODU formation', default=0.2_rk) 
   call self%get_parameter(self%ODUCr, 'ODUCr', 'molODU molC-1', 'ODU:C ratio in anoxic remineralization', default=1.0_rk) 
   call self%get_parameter(self%Param1IronCurve, 'Param1IronCurve', '-', 'Parameter of the curve simulating the iron concentration', default=25000.0_rk) 
   call self%get_parameter(self%Param2IronCurve, 'Param2IronCurve', '-', 'Parameter of the curve simulating the iron c', default=50.0_rk) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%Q10bac, 'Q10bac', '-', 'Temperature factor for BAC', default=2.0_rk) 

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
   class (type_ulg_Bacteria), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,ODU,PHO
      real(rk) ::  par,temp
      real(rk) ::  BAC
      real(rk) ::   bacteria_anoxrem,bacteria_oxygenconsumption,Bacteria_Respiration,denitrification,Uptake_DOCL
      real(rk) ::   bacteria_anoxremIntegrated,bacteria_oxygenconsumptionIntegrated,Bacteria_RespirationIntegrated,DenitrificationIntegrated,Uptake_DOCLIntegrated
      real(rk) ::   BACExcr	  ! mmol N m-3, Bacteria excretion flux of ammonium
      real(rk) ::   BACGrowth	  ! mmol C m-3, Bacterial growth
      real(rk) ::   BACResp	  ! mmol C m-3, Bacteria respiration flux
      real(rk) ::   bacteria_anoxrem_local	  ! mmol ODU m-3, Bacterial anoxic remineralisation
      real(rk) ::   bacteria_oxygenconsumption_local	  ! mmol O2 m-3, Bacterial respiration
      real(rk) ::   C_BACMort	  ! mmol C m-3, Bacteria mortality flux in carbon
      real(rk) ::   denitrif	  ! mmol N m-3, Denitrification flux
      real(rk) ::   Iron	  ! mmol Fe m-3, Iron concentration
      real(rk) ::   Limitation_By_Iron	  ! -, Limitation by Iron
      real(rk) ::   NutLim	  ! -, Nutrient limitation
      real(rk) ::   N_BACMort	  ! mmol N m-3, Bacteria mortality flux in nitrogen
      real(rk) ::   tf	  ! -, Temperature factor
      real(rk) ::   testratio	  ! mmol N m-3, Value showing if BAC growth limited by C or N
      real(rk) ::   Uptake_DOCL_local	  ! mmol C m-3, Bacteria uptake of DOC
      real(rk) ::   Uptake_DONL_local	  ! mmol N m-3, Bacteria uptake of DON
      real(rk) ::   Uptake_NHS_local	  ! mmol N m-3, Bacteria uptake of ammonium
      real(rk) ::   Uptake_Potential_NHS	  ! mmol N m-3, Bacteria potential uptake of ammonium
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_bac,BAC)       ! Bacterial biomass
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
    
    tf = Q10Factor(temp,Q10bac) 
    
   !Compute Iron Limitation 
    Iron = self%IronCsurf + self%Param1IronCurve/(self%self%Param2IronCurve*sqrt(2.*3.1416))*exp(-(depth-275.0)**2/(2.*self%self%Param2IronCurve**2))
    Limitation_By_Iron = Michaelis(Iron,Halfsaturation_Iron)
   													 
    Uptake_DOCL_local = tf*self%maxgrowthbac*BAC*Michaelis(DCL,csatdocl)
    Uptake_DONL_local = Uptake_DOCL_local*(DNL/DCL)
    Nutlim = min(Michaelis(NHS,csatamm),Michaelis(PHO,csatpo4))
    Uptake_Potential_NHS = self%maxgrowthbac*NutLim*BAC*NCrBac		
   	 
   ! Test if bacterial growth is limited by carbon (DOClabile) or nitrogen (ammonium + DONLabile) 
    testratio = Uptake_DOCL_local*(Uptake_DONL_local/Uptake_DOCL_local - self%bactgrowtheff*NCrBac)	
    
   ! Growth rate of bacteria depends on a threshold value compared with Uptake_Potential_NHS 
    if (Uptake_Potential_NHS > (-testratio)) then 
      ! In this case we are in a situation of carbon limitation 
      BACGrowth = self%bactgrowtheff*Uptake_DOCL_local
      ! Growth rate computed taking into account the iron limitation 
      BACResp = Uptake_DOCL_local*(1.0 - self%bactgrowtheff)
      ! Test if NHS uptake is necessary to use all the DOCl 
      if (testratio > 0) then 
        ! We are in case of remineralisation of ammonium through bacteria excretion and no net uptake of ammonium is necessary to use all the DOC 
        Uptake_NHS_local = 0
        BACExcr = testratio
      else 
        Uptake_NHS_local = -testratio
        BACExcr = 0
      endif 
    else 
      ! if we are in case of nitrogen limitation,it means that all the DON and the potential uptake of NHS is not sufficient to consume all the DOC 
      Uptake_NHS_local = Uptake_Potential_NHS
      BACGrowth = (Uptake_Potential_NHS + Uptake_DONL_local)/self%NCrBac
      ! Growth rate computed taking into account the iron limitation 
      BACExcr = 0
      BACResp = BACGrowth*(1.0/bactgrowthefficiency  - 1.0)
    endif 
    
   ! Bacteria mortality rate (C_BACMort,N_BACMort, /day) 
    C_BACMort = self%mortbac*tf*BAC
    N_BACMort = self%mortbac*tf*BAC*self%NCrBac
    
    denitrif = BACResp*NOS/(NOS+self%ksdeninos)*(self%self%kindenidox/(DOX+self%self%kindenidox))*self%NCr
    bacteria_oxygenconsumption_local = BACResp*DOX/(DOX+self%ksremindox) * self%OCr
    bacteria_anoxrem_local = BACResp*(self%self%kinanoxremnos/(NOS+self%self%kinanoxremnos))* (self%self%kinanoxremdox/(DOX+self%self%kinanoxremdox))*self%ODUCr
    
   ! Carbon content increases by intake of nutrients and decreases by mortality (and predation) 
   _ADD_SOURCE_(self%id_bac, BACGrowth - C_BACMort)
   ! Ammonium is excreted or can be taken up 
   _ADD_SOURCE_(self%id_nhs, BACExcr - Uptake_NHS_local) 
   _ADD_SOURCE_(self%id_pho, (BACExcr - Uptake_NHS_local)*self%PNRedfield) 
   _ADD_SOURCE_(self%id_dcl, self%labilefraction*C_BACMort - BACGrowth - BACResp) 
   _ADD_SOURCE_(self%id_dnl, self%labilefraction*N_BACMort - Uptake_DONL_local) 
   _ADD_SOURCE_(self%id_dcs, (1.0 - self%labilefraction)*C_BACMort) 
   _ADD_SOURCE_(self%id_dns, (1.0 - self%labilefraction)*N_BACMort) 
   _ADD_SOURCE_(self%id_dox, - bacteria_oxygenconsumption_local) 
   ! NOS decreases due to bacterial respiration 
   _ADD_SOURCE_(self%id_nos, -denitrif) 
   _ADD_SOURCE_(self%id_dic, BACResp) 
   ! ODU increases due to bacterial anoxic remineralisation 
   _ADD_SOURCE_(self%id_odu, (1.0 - Limitation_By_Iron*self%ODU_solid)*bacteria_anoxrem_local) 

   _SET_DIAGNOSTIC_(self%id_Uptake_DOCL, Uptake_DOCL_local)
   _SET_DIAGNOSTIC_(self%id_Bacteria_Respiration, BACResp)
   _SET_DIAGNOSTIC_(self%id_denitrification, denitrif)
   _SET_DIAGNOSTIC_(self%id_bacteria_oxygenconsumption, bacteria_anoxrem_local)
   _SET_DIAGNOSTIC_(self%id_bacteria_anoxrem, bacteria_oxygenconsumption_local)

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_Bacteria 
