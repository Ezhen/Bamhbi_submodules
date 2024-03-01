#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DemilianaF90
!
! 1-D ecosystem model - Biological model of Tett
!
! Tett, P., 1998. Parameterising a microplankton model.
! Department of Biological Sciences, Napier University,
! Report ISBN 0 902703 60 9, 60 pp.
!
! Sharples Tett (1994). Modeling the effect of physical! variability on the midwater chlorophyll maximum.
! Journal of marine research 52: 219-238
!
! Implementation: Marilaure Gregoire,             NIOO-CEME
! Translation into FABM: E. Ivanov, ULg / MAST
!
! Contains the pelagic submodel, as used in Soetaert et al., 2001.
!
!--------------------------------------------------------------------*

   module fabm_ulg_Emiliana 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Emiliana 
      type (type_state_variable_id)         :: id_cem,id_nem
      type (type_state_variable_id)         :: id_cfl,id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nfl,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeEmiliana,id_Nitrogen_Uptake_Emiliana,id_NPP,id_PhytoNitrateReduction,id_TotalRespirationEmiliana
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeEmilianaIntegrated,id_Nitrogen_Uptake_EmilianaIntegrated,id_NPPIntegrated,id_PhytoNitrateReductionIntegrated,id_TotalRespirationEmilianaIntegrated

!     Model parameters 
      real(rk)     :: alphaPIEmiliana, extradocphyexcr, GrowthRespEmiliana
      real(rk)     :: kinNHsPhy, ksNHsEmiliana, ksNOsEmiliana
      real(rk)     :: ksPO4Emiliana, labileextradocphyexcr, labilefraction
      real(rk)     :: leakagephy, MaxChlNrEmiliana, MaxNCrEmiliana
      real(rk)     :: MinChlNrEmiliana, MinNCrEmiliana, MortalityEmiliana
      real(rk)     :: mortphydom, MuMaxEmiliana, NHsMaxUptakeEmiliana
      real(rk)     :: NosMaxUptakeEmiliana, OCr, ONoxnhsr, PNRedfield
      real(rk)     :: PO4MaxUptakeEmiliana, Q10Phy, QuantumYieldEmiliana
      real(rk)     :: RespirationEmiliana

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
   ! Initialise the Emiliana model

   subroutine initialize(self,configunit)
   class (type_ulg_Emiliana), intent(inout), target :: self
   integer,                        intent(in)          :: configunit


   namelist /ulg_Emiliana/ alphaPIEmiliana, 	 & 
                      extradocphyexcr, GrowthRespEmiliana, 	 & 
                      kinNHsPhy, ksNHsEmiliana, ksNOsEmiliana, 	 & 
                      ksPO4Emiliana, labileextradocphyexcr, 	 & 
                      labilefraction, leakagephy, 	 & 
                      MaxChlNrEmiliana, MaxNCrEmiliana, 	 & 
                      MinChlNrEmiliana, MinNCrEmiliana, 	 & 
                      MortalityEmiliana, mortphydom, 	 & 
                      MuMaxEmiliana, NHsMaxUptakeEmiliana, 	 & 
                      NosMaxUptakeEmiliana, OCr, ONoxnhsr, 	 & 
                      PNRedfield, PO4MaxUptakeEmiliana, 	 & 
                      Q10Phy, QuantumYieldEmiliana, 	 & 
                      RespirationEmiliana

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%alphaPIEmiliana, 'alphaPIEmiliana', 'm2 W-1 d-1', 'Initial slope of photosynthesis-light curve for EM', default=0.3_rk) 
   call self%get_parameter(self%extradocphyexcr, 'extradocphyexcr', '-', 'Extra-photosynthetic DOC excretion', default=0.05_rk) 
   call self%get_parameter(self%GrowthRespEmiliana, 'GrowthRespEmiliana', '-', 'Part of primary production used for respiration by EM ', default=0.1_rk) 
   call self%get_parameter(self%kinNHsPhy, 'kinNHsPhy', 'mmolN m-3', 'Inhib. constant of NHS for NOS uptake by PHY', default=0.5_rk) 
   call self%get_parameter(self%ksNHsEmiliana, 'ksNHsEmiliana', 'mmolN m-3', 'Half-saturation constant for NHS uptake by EM', default=0.05_rk) 
   call self%get_parameter(self%ksNOsEmiliana, 'ksNOsEmiliana', 'mmolN m-3', 'Half-saturation constant for NOS uptake by EM', default=0.05_rk) 
   call self%get_parameter(self%ksPO4Emiliana, 'ksPO4Emiliana', 'mmolP m-3', 'Half-saturation constant for PO4 uptake by EM', default=0.02_rk) 
   call self%get_parameter(self%labileextradocphyexcr, 'labileextradocphyexcr', '-', 'Labile fraction phytoxcreted DOC', default=0.65_rk) 
   call self%get_parameter(self%labilefraction, 'labilefraction', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%leakagephy, 'leakagephy', '-', 'Phytoplankton leakage fraction', default=0.02_rk) 
   call self%get_parameter(self%MaxChlNrEmiliana, 'MaxChlNrEmiliana', 'g Chla molN-1', 'Maximum Chl:N ratio in EM', default=2.0_rk) 
   call self%get_parameter(self%MaxNCrEmiliana, 'MaxNCrEmiliana', 'molN molC-1', 'Maximum N:C ratio in EM', default=0.2_rk) 
   call self%get_parameter(self%MinChlNrEmiliana, 'MinChlNrEmiliana', 'g Chla molN-1', 'Minimum Chl:N ratio in EM', default=1.0_rk) 
   call self%get_parameter(self%MinNCrEmiliana, 'MinNCrEmiliana', 'molN molC-1', 'Minimum N:C ratio in EM', default=0.05_rk) 
   call self%get_parameter(self%MortalityEmiliana, 'MortalityEmiliana', 'd-1', 'Mortality rate of EM', default=0.03_rk) 
   call self%get_parameter(self%mortphydom, 'mortphydom', '-', 'DOM fraction of phytoplankton mortality', default=0.34_rk) 
   call self%get_parameter(self%MuMaxEmiliana, 'MuMaxEmiliana', 'd-1', 'Maximum specific growth rate of EM', default=2.5_rk) 
   call self%get_parameter(self%NHsMaxUptakeEmiliana, 'NHsMaxUptakeEmiliana', 'molN molC-1 d-1', 'Maximal NHS uptake rate by EM', default=1.5_rk) 
   call self%get_parameter(self%NosMaxUptakeEmiliana, 'NosMaxUptakeEmiliana', 'molN molC-1 d-1', 'Maximal NOS uptake rate by EM', default=1.5_rk) 
   call self%get_parameter(self%OCr, 'OCr', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%ONoxnhsr, 'ONoxnhsr', 'molO2 molNS-1', 'O2:NHS ratio in NHS oxidation in nitrification', default=2.0_rk) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%PO4MaxUptakeEmiliana, 'PO4MaxUptakeEmiliana', 'molP molC-1 d-1', 'Maximal PO4 uptake rate by EM', default=0.09375_rk) 
   call self%get_parameter(self%Q10Phy, 'Q10Phy', '-', 'Temperature factor', default=2.0_rk) 
   call self%get_parameter(self%QuantumYieldEmiliana, 'QuantumYieldEmiliana', 'mmolC (mg Chl dW m-2)-1', 'Maximum quantum yield of EM', default=0.6_rk) 
   call self%get_parameter(self%RespirationEmiliana, 'RespirationEmiliana', 'd-1', 'Basal respiration rate of EM', default=0.009_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_cem, 'CEM'  & 
         , 'mmol C m-3', 'Small flagellate biomass in carbon' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_nem, 'NEM'  & 
         , 'mmol N m-3', 'Small flagellate biomass in nitrogen' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_dependency(self%id_cfl, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcl, 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_nfl, 'Large flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nhs, 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nos, 'Nitrate concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'Particulate organic nitrogen concentration', 'mmol N m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_Carbon_UptakeEmiliana, 'Carbon_UptakeEmiliana', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Carbon_UptakeEmilianaIntegrated, 'Carbon_UptakeEmilianaIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Nitrogen_Uptake_Emiliana, 'Nitrogen_Uptake_Emiliana', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Nitrogen_Uptake_EmilianaIntegrated, 'Nitrogen_Uptake_EmilianaIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_NPP, 'NPP', 'mmol N m-3 d-1', & 
      ' Primary production', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_NPPIntegrated, 'NPPIntegrated', 'mmol N m-3 d-1', & 
      ' Primary production (vertically-integrated)', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_PhytoNitrateReduction, 'PhytoNitrateReduction', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_PhytoNitrateReductionIntegrated, 'PhytoNitrateReductionIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_TotalRespirationEmiliana, 'TotalRespirationEmiliana', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_TotalRespirationEmilianaIntegrated, 'TotalRespirationEmilianaIntegrated', '-', & 
      '-', output=output_instantaneous) 

   return 

99 call self%fatal_error('Emiliana', 'Error reading namelist ulg_Emiliana') 

   end subroutine initialize 


   ! Right hand sides of Emiliana model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_Emiliana), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  CFL,DCL,DCS,DIC,DNL,DNS,DOX,NFL,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CEM,NEM
      real(rk) ::   Carbon_UptakeEmiliana,Nitrogen_Uptake_Emiliana,NPP,PhytoNitrateReduction,TotalRespirationEmiliana
      real(rk) ::   Carbon_UptakeEmilianaIntegrated,Nitrogen_Uptake_EmilianaIntegrated,NPPIntegrated,PhytoNitrateReductionIntegrated,TotalRespirationEmilianaIntegrated
      real(rk) ::   Ammonium_UpPHY	  ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   C_PHYMort	  ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Carbon_UptakePHY	  ! mmol C m-3, C assimilation of phytoplankton
      real(rk) ::   ChlCrEmiliana	  ! g Chla mol C-1, Chl/C ratio in small flagellates
      real(rk) ::   DOC_extra_excr	  ! mmol C d-1, Phytoplankton extra excretion
      real(rk) ::   DOC_leakage	  ! mmol C d-1, Phytoplankton passive leakage rate for carbon
      real(rk) ::   DON_leakage	  ! mmol N d-1, Phytoplankton passive leakage rate for nitrogen
      real(rk) ::   GrowthPHY	  ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   LightLimitationEmiliana	  ! -, Light limitation for small flagellates
      real(rk) ::   Mu_Nitrogen	  ! ?, ?
      real(rk) ::   Mu_Silicate	  ! ?, ?
      real(rk) ::   N_PHYMort	  ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   NCrat	  ! ?, ?
      real(rk) ::   NCrEmiliana	  ! mol N mol C-1, N/C ratio in small flagellates
      real(rk) ::   Nitrate_UpPHY	  ! mmol N m-3, Nitrate uptake of phytoplankton
      real(rk) ::   Nitrogen_UpPHY	  ! mmol N m-3, Nitrogen uptake of phytoplankton
      real(rk) ::   Nutrient_UpPHY	  ! mmol m-3, Nutrient uptake of phytoplankton
      real(rk) ::   NutrientLimitationEmiliana	  ! -, Nutrient limitation for small flagellates
      real(rk) ::   Phosphate_upEmiliana	  ! mmol P m-3 , Small flagellates phosphate uptake
      real(rk) ::   PHYMort	  ! mmol m-3, Phytoplankton mortality rate
      real(rk) ::   SiCrat	  ! ?, ?
      real(rk) ::   tf	  ! -, Temperature factor
      real(rk) ::   TotalRespirationPHY	  ! mmol C m-3, Total phytoplankton respiration (basal & activity)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
   _GET_(self%id_nem,NEM)       ! Small flagellate biomass in nitrogen
   _GET_(self%id_cfl,CFL)       ! Small flagellate biomass in carbon
   _GET_(self%id_dcl,DCL)       ! Labile detritus concentration in carbon
   _GET_(self%id_dcs,DCS)       ! Semi-labile detritus concentration in carbon
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dnl,DNL)       ! Labile detritus concentration in nitrogen
   _GET_(self%id_dns,DNS)       ! Semi-labile detritus concentration in nitrogen
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_nfl,NFL)       ! Large flagellate biomass in nitrogen
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_nos,NOS)       ! Nitrate concentration
   _GET_(self%id_pho,PHO)       ! Phosphorus
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor(temp,Q10Phy)
   ! Calculate ratios in phytoplankton 
    NCrEmiliana  = NEM/CEM
    ChlCrEmiliana = ChlCrPHY(NCratio,MaxNCr,MinNCr,MinChlNr,MaxChlNr)
    
             ! Chlorophyll(i,j,k,3) = ChlCrEmiliana * CEM	! REMOVE (POSSIBLY) 
    
   ! Nitrate uptake rate 
    Nitrate_UpPHY = NOuptake(NCrEmiliana,tf,MaxNCrEmiliana,NosMaxUptakeEmiliana) * Michaelis(NOS,ksNOsEmiliana) * Inhibition(NHS,kinNHsPhy) * CEM
    
   ! Ammonium uptake rate 
    Ammonium_UpPHY = NUT_UPTAKE_RATE(NCrEmiliana,(NHS-0.0),tf,MaxNCrEmiliana,NHsMaxUptakeEmiliana,ksNHsEmiliana) * CEM
    
   ! Phosphate uptake rate 
    Phosphate_upEmiliana = NUT_UPTAKE_RATE(NCrEmiliana,(PHO-0.0),tf,MaxNCrEmiliana,PO4MaxUptakeEmiliana,ksPO4Emiliana) * CEM
    
    
   ! Potential nitrogen uptake 
    Nitrogen_UpPHY = Ammonium_UpPHY + Nitrate_UpPHY
    
   ! Nutrient uptake 
    Nutrient_UpPHY = min(Nitrogen_UpPHY,Phosphate_upEmiliana/self%PNRedfield)
    
   ! Compute actual N:C and Si:C ratios 
    NCrat = Ratio_PHYT(NCrEmiliana,MaxNCrEmiliana,MinNCrEmiliana)
    SiCrat = Ratio_PHYT(SiCrEmiliana,MaxSiCrEmiliana,MinSiCrEmiliana)
    
    
   ! Compute nutrient and light limitation 
    NutrientLimitationEmiliana = Nutr_LIM(MinNCrEmiliana,MinSiCrEmiliana,NCrat,SiCrat)
          LightLimitationEmiliana = 1.0-exp(-self%alphaPIEmiliana*PAR/self%MuMaxEmiliana)
    
   ! Compute carbon uptake 
    Carbon_UptakePHY=self%MuMaxEmiliana*LightLimitationEmiliana*NutrientLimitationEmiliana*CEM*tf
    
   ! Compute respiration 
    TotalRespirationPHY=Carbon_UptakePHY*self%GrowthRespEmiliana + self%RespirationEmiliana*CEM*tf
    
   ! Compute growth 
    GrowthPHY=Carbon_UptakePHY-TotalRespirationPHY
    
   ! Compute the extra DOC excretion from the difference of growth in nutrient limited (actual NC ratio) and nutrient saturated (max NC ratio) 
    DOC_extra_excr = CEM * tf * self%extradocphyexcr * self%MuMaxEmiliana * ExtraEXCR_term(LightLimitationEmiliana,MinNCrEmiliana,MaxNCrEmiliana,NCrat)        
    
   ! Compute the leakage 
    DOC_leakage = self%leakagephy * Carbon_UptakePHY
    DON_leakage = self%leakagephy * abs(Nutrient_UpPHY)
    
   ! Compute mortality 
    C_PHYMort  = PHYMORT(MortalityEmiliana,tf) * CEM
    N_PHYMort  = PHYMORT(MortalityEmiliana,tf) * NEM
    
   ! Carbon in Emiliana increases by growth and decreases by leakage and mortality 
   _ADD_SOURCE_(self%id_cem,1.0*( GrowthPHY)) 
   _ADD_SOURCE_(self%id_cem,-1.0*( C_PHYMort + DOC_leakage)) 
   _ADD_SOURCE_(self%id_nem,1.0*( Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_nem,-1.0*( N_PHYMort + DON_leakage)) 
    
   ! IF CN ratio of phytoplankton higher than CNmin, than nitrogen is taken up, unless it gets excreted 
             IF (Nutrient_UpPHY.gt.0) THEN 
   _ADD_SOURCE_(self%id_nos,-1.0*( Nutrient_UpPHY*Nitrate_UpPHY/Nitrogen_UpPHY)) 
   _ADD_SOURCE_(self%id_dox,1.0*( Nutrient_UpPHY*Nitrate_UpPHY/Nitrogen_UpPHY*self%ONoxnhsr)) 
   _ADD_SOURCE_(self%id_nhs,-1.0*( Nutrient_UpPHY*Ammonium_UpPHY/Nitrogen_UpPHY)) 
   _ADD_SOURCE_(self%id_pho,-1.0*( Nutrient_UpPHY*self%PNRedfield)) 
             ELSE 
   _ADD_SOURCE_(self%id_nhs,1.0*( - Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_pho,1.0*( - Nutrient_UpPHY*self%PNRedfield)) 
             END IF 
    
   ! Mortality increases the pool of POM and DOM with proper partitioning and leakage adds to the labile pool 
   _ADD_SOURCE_(self%id_poc,1.0*( (1.0 - self%mortphydom)*C_PHYMort)) 
   _ADD_SOURCE_(self%id_pon,1.0*( (1.0 - self%mortphydom)*N_PHYMort)) 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%mortphydom*C_PHYMort*self%labilefraction + DOC_leakage)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( self%mortphydom*N_PHYMort*self%labilefraction + DON_leakage)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( self%mortphydom*C_PHYMort*(1.0 - self%labilefraction))) 
   _ADD_SOURCE_(self%id_dns,1.0*( self%mortphydom*N_PHYMort*(1.0 - self%labilefraction))) 
    
   ! Extra-excretion and leakage add to the labile and semi-labile detritus 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%self%leakagephy*DOC_extra_excr + self%labileextradocphyexcr*(1.0 - self%self%leakagephy)*DOC_extra_excr)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( (1.0 - self%labileextradocphyexcr)*(1.0 - self%leakagephy)*DOC_extra_excr)) 
    
   ! The oxygen increases due to photosynthesis 
   _ADD_SOURCE_(self%id_dox,1.0*( (GrowthPHY + DOC_extra_excr)*self%OCr)) 
    
   ! CO2 production and consumption 
   _ADD_SOURCE_(self%id_dic,-1.0*( GrowthPHY + DOC_extra_excr)) 
    
#ifdef primaryprod 
          NPP = NPP + Carbon_UptakePHY-TotalRespirationPHY
#endif 
#ifdef biodiag2 
          Nitrogen_Uptake_Emiliana = Ammonium_UpPHY + Nitrate_UpPHY
          Carbon_UptakeEmiliana    = Carbon_UptakePHY
          TotalRespirationEmiliana = TotalRespirationPHY
#endif 
#ifdef biodiag1 
          PhytoNitrateReduction = PhytoNitrateReduction+Nutrient_UpPHY*ratio(Nitrate_UpPHY,Nitrogen_UpPHY)*self%ONoxnhsr
#endif 

   ! Averaged over entire water column Diagnostics 
#ifdef biodiag2 
   _SET_DIAGNOSTIC_(self%id_Nitrogen_Uptake_Emiliana, Nitrogen_Uptake_Emiliana)
   _SET_DIAGNOSTIC_(self%id_Carbon_UptakeEmiliana, Carbon_UptakeEmiliana)
   _SET_DIAGNOSTIC_(self%id_TotalRespirationEmiliana, TotalRespirationEmiliana)
#endif 
#ifdef biodiag1 
   _SET_DIAGNOSTIC_(self%id_PhytoNitrateReduction, PhytoNitrateReduction)
#endif 
#ifdef primaryprod 
   _SET_DIAGNOSTIC_(self%id_NPP, NPP)
#endif 
   _LOOP_END_

   end subroutine do


   end module fabm_ulg_Emiliana 
