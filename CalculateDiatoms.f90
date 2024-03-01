#include "fabm_driver.h" 
 
!#########################################################################################
!                              3Ddiatoms.F90
!
! 1-D ecosystem model - Biological model of Tett

! References to the microplankton model:
! Tett, P., 1998. Parameterising a microplankton model.
! Department of Biological Sciences, Napier University,
! Report ISBN 0 902703 60 9, 60 pp.
!
! Sharples Tett (1994). Modeling the effect of physical! variability on the midwater chlorophyll maximum.
! Journal of marine research 52: 219-238
!
! Implementation: Marilaure Gregoire,                 NIOO-CEME
! Translation in FABM: Evgeny Ivanov,   Universite de Liege, MAST
!
! Contains the pelagic submodel, as used in Soetaert et al., 2001.
!######################################################################

   module fabm_ulg_Diatoms 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Diatoms 
      type (type_state_variable_id)         :: id_cdi,id_ndi,id_sid,id_sio
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeDiatoms,id_Nitrogen_Uptake_Diatoms,id_NPP,id_PhytoNitrateReduction,id_Silicate_upDiatoms,id_TotalRespirationDiatoms
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeDiatomsIntegrated,id_Nitrogen_Uptake_DiatomsIntegrated,id_Silicate_upDiatomsIntegrated,id_TotalRespirationDiatomsIntegrated

!     Model parameters 
      real(rk)     :: alphaPIDiatoms, extradocphyexcr, GrowthRespDiatoms
      real(rk)     :: kdis_Silicious_Detritus, kinNHsPhy, ksNHsDiatoms
      real(rk)     :: ksNOsDiatoms, ksPO4Diatoms, ksSiDiatoms
      real(rk)     :: labileextradocphyexcr, labilefraction, leakagephy
      real(rk)     :: MaxChlNrDiatoms, MaxNCrDiatoms, MinChlNrDiatoms
      real(rk)     :: MinNCrDiatoms, MortalityDiatoms, mortphydom
      real(rk)     :: MuMaxDiatoms, NHsMaxUptakeDiatoms, NosMaxUptakeDiatoms
      real(rk)     :: OCr, ONoxnhsr, PNRedfield, PO4MaxUptakeDiatoms
      real(rk)     :: Q10PhyDiatoms, Q10SilicateDiss, QuantumYieldDiatoms
      real(rk)     :: RespirationDiatoms, SiMaxUptakeDiatoms, SiNrDiatoms

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
   ! Initialise the Diatoms model

   subroutine initialize(self,configunit)
   class (type_ulg_Diatoms), intent(inout), target :: self
   integer,                        intent(in)          :: configunit


   namelist /ulg_Diatoms/ alphaPIDiatoms, 	 & 
                      extradocphyexcr, GrowthRespDiatoms, 	 & 
                      kdis_Silicious_Detritus, kinNHsPhy, 	 & 
                      ksNHsDiatoms, ksNOsDiatoms, 	 & 
                      ksPO4Diatoms, ksSiDiatoms, 	 & 
                      labileextradocphyexcr, labilefraction, 	 & 
                      leakagephy, MaxChlNrDiatoms, 	 & 
                      MaxNCrDiatoms, MinChlNrDiatoms, 	 & 
                      MinNCrDiatoms, MortalityDiatoms, 	 & 
                      mortphydom, MuMaxDiatoms, 	 & 
                      NHsMaxUptakeDiatoms, 	 & 
                      NosMaxUptakeDiatoms, OCr, ONoxnhsr, 	 & 
                      PNRedfield, PO4MaxUptakeDiatoms, 	 & 
                      Q10PhyDiatoms, Q10SilicateDiss, 	 & 
                      QuantumYieldDiatoms, RespirationDiatoms, 	 & 
                      SiMaxUptakeDiatoms, SiNrDiatoms

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%alphaPIDiatoms, 'alphaPIDiatoms', 'm2 W-1 d-1', 'Initial slope of photosynthesis-light curve for DI', default=0.3312_rk) 
   call self%get_parameter(self%extradocphyexcr, 'extradocphyexcr', '-', 'Extra-photosynthetic DOC excretion', default=0.05_rk) 
   call self%get_parameter(self%GrowthRespDiatoms, 'GrowthRespDiatoms', '-', 'Part of primary production used for respiration by DI', default=0.1_rk) 
   call self%get_parameter(self%kdis_Silicious_Detritus, 'kdis_Silicious_Detritus', 'd-1', 'Rate of dissolution of silicious detritus', default=0.08_rk) 
   call self%get_parameter(self%kinNHsPhy, 'kinNHsPhy', 'mmolN m-3', 'Inhib. constant of NHS for NOS uptake by PHY', default=0.5_rk) 
   call self%get_parameter(self%ksNHsDiatoms, 'ksNHsDiatoms', 'mmolN m-3', 'Half-saturation constant for NHS uptake by DI', default=1.0_rk) 
   call self%get_parameter(self%ksNOsDiatoms, 'ksNOsDiatoms', 'mmolN m-3', 'Half-saturation constant for NOS uptake by DI', default=1.0_rk) 
   call self%get_parameter(self%ksPO4Diatoms, 'ksPO4Diatoms', 'mmolP m-3', 'Half-saturation constant for PO4 uptake by DI', default=0.1_rk) 
   call self%get_parameter(self%ksSiDiatoms, 'ksSiDiatoms', 'mmolSi m-3', 'Half-saturation constant for SiOs uptake by DI', default=3.5_rk) 
   call self%get_parameter(self%labileextradocphyexcr, 'labileextradocphyexcr', '-', 'Labile fraction phytoxcreted DOC', default=0.65_rk) 
   call self%get_parameter(self%labilefraction, 'labilefraction', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%leakagephy, 'leakagephy', '-', 'Phytoplankton leakage fraction', default=0.02_rk) 
   call self%get_parameter(self%MaxChlNrDiatoms, 'MaxChlNrDiatoms', 'g Chla molN-1', 'Maximum Chl:N ratio in DI', default=2.0_rk) 
   call self%get_parameter(self%MaxNCrDiatoms, 'MaxNCrDiatoms', 'molN molC-1', 'Maximum N:C ratio in DI', default=0.2_rk) 
   call self%get_parameter(self%MinChlNrDiatoms, 'MinChlNrDiatoms', 'g Chla molN-1', 'Minimum Chl:N ratio in DI', default=1.0_rk) 
   call self%get_parameter(self%MinNCrDiatoms, 'MinNCrDiatoms', 'molN molC-1', 'Minimum N:C ratio in DI', default=0.05_rk) 
   call self%get_parameter(self%MortalityDiatoms, 'MortalityDiatoms', 'd-1', 'Mortality rate of DI', default=0.03_rk) 
   call self%get_parameter(self%mortphydom, 'mortphydom', '-', 'DOM fraction of phytoplankton mortality', default=0.34_rk) 
   call self%get_parameter(self%MuMaxDiatoms, 'MuMaxDiatoms', 'd-1', 'Maximum specific growth rate of DI', default=3.5_rk) 
   call self%get_parameter(self%NHsMaxUptakeDiatoms, 'NHsMaxUptakeDiatoms', 'molN molC-1 d-1', 'Maximal NHS uptake rate by DI', default=1.0_rk) 
   call self%get_parameter(self%NosMaxUptakeDiatoms, 'NosMaxUptakeDiatoms', 'molN molC-1 d-1', 'Maximal NOS uptake rate by DI', default=1.0_rk) 
   call self%get_parameter(self%OCr, 'OCr', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%ONoxnhsr, 'ONoxnhsr', 'molO2 molNS-1', 'O2:NHS ratio in NHS oxidation in nitrification', default=2.0_rk) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 
   call self%get_parameter(self%PO4MaxUptakeDiatoms, 'PO4MaxUptakeDiatoms', 'molP molC-1 d-1', 'Maximal PO4 uptake rate by DI', default=0.0625_rk) 
   call self%get_parameter(self%Q10PhyDiatoms, 'Q10PhyDiatoms', '-', 'Temperature factor for DI', default=1.8_rk) 
   call self%get_parameter(self%Q10SilicateDiss, 'Q10SilicateDiss', '-', 'Temperature factor for chemical processes', default=3.3_rk) 
   call self%get_parameter(self%QuantumYieldDiatoms, 'QuantumYieldDiatoms', 'mmolC (mg Chl dW m-2)-1', 'Maximum quantum yield of DI', default=0.8_rk) 
   call self%get_parameter(self%RespirationDiatoms, 'RespirationDiatoms', 'd-1', 'Basal respiration rate of DI', default=0.009_rk) 
   call self%get_parameter(self%SiMaxUptakeDiatoms, 'SiMaxUptakeDiatoms', 'molSi molC-1 d-1', 'Maximal SiOs uptake rate by DI', default=0.5_rk) 
   call self%get_parameter(self%SiNrDiatoms, 'SiNrDiatoms', 'molSi molN-1', 'Si:N ratio in DI', default=0.83_rk) 

   ! Register state variables 

   call self%register_state_variable(self%id_cdi, 'CDI'  & 
         , 'mmol C m-3', 'Diatom biomass in carbon' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_ndi, 'NDI'  & 
         , 'mmol N m-3', 'Diatom biomass in nitrogen' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_sid, 'SID'  & 
         , 'mmol Si m-3', 'Detrital silicate concentration' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_sio, 'SIO'  & 
         , 'mmol Si m-3', 'Silicilic acid concentration' & 
         minimum=0.0e-7_rk)
   call self%register_state_dependency(self%id_dcl, 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_nhs, 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nos, 'Nitrate concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'Particulate organic nitrogen concentration', 'mmol N m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_Carbon_UptakeDiatoms, 'Carbon_UptakeDiatoms', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Carbon_UptakeDiatomsIntegrated, 'Carbon_UptakeDiatomsIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Nitrogen_Uptake_Diatoms, 'Nitrogen_Uptake_Diatoms', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Nitrogen_Uptake_DiatomsIntegrated, 'Nitrogen_Uptake_DiatomsIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_NPP, 'NPP', 'mmol N m-3 d-1', & 
      ' Primary production', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_PhytoNitrateReduction, 'PhytoNitrateReduction', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Silicate_upDiatoms, 'Silicate_upDiatoms', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Silicate_upDiatomsIntegrated, 'Silicate_upDiatomsIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_TotalRespirationDiatoms, 'TotalRespirationDiatoms', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_TotalRespirationDiatomsIntegrated, 'TotalRespirationDiatomsIntegrated', '-', & 
      '-', output=output_instantaneous) 

   return 

99 call self%fatal_error('Diatoms', 'Error reading namelist ulg_Diatoms') 

   end subroutine initialize 


   ! Right hand sides of Diatoms model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_ulg_Diatoms), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CDI,NDI,SID,SIO
      real(rk) ::   Carbon_UptakeDiatoms,Nitrogen_Uptake_Diatoms,NPP,PhytoNitrateReduction,Silicate_upDiatoms,TotalRespirationDiatoms
      real(rk) ::   Carbon_UptakeDiatomsIntegrated,Nitrogen_Uptake_DiatomsIntegrated,Silicate_upDiatomsIntegrated,TotalRespirationDiatomsIntegrated
      real(rk) ::   Ammonium_UpPHY	  ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   C_PHYMort	  ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Carbon_UptakePHY	  ! mmol C m-3, C assimilation of phytoplankton
      real(rk) ::   ChlCrDiatoms	  ! g Chla mol C-1, Chl/C ratio in large flagellates
      real(rk) ::   DOC_extra_excr	  ! mmol C d-1, Phytoplankton extra excretion
      real(rk) ::   DOC_leakage	  ! mmol C d-1, Phytoplankton passive leakage rate for carbon
      real(rk) ::   DON_leakage	  ! mmol N d-1, Phytoplankton passive leakage rate for nitrogen
      real(rk) ::   GrowthPHY	  ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   LightLimitationDiatoms	  ! -, Light limitation diatoms
      real(rk) ::   MaxSiCrDiatoms	  ! mol Si mol C-1, Maximum Si/C ratio in diatoms
      real(rk) ::   MinSiCrDiatoms	  ! mol Si mol C-1, Minimum Si/C ratio in diatoms
      real(rk) ::   Mu_Nitrogen	  ! ?, ?
      real(rk) ::   Mu_Silicate	  ! ?, ?
      real(rk) ::   NCrat	  ! ?, ?
      real(rk) ::   N_PHYMort	  ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   NCrDiatoms	  ! mol N mol C-1, N/C ratio in diatoms
      real(rk) ::   Nitrate_UpPHY	  ! mmol N m-3, Nitrate uptake of phytoplankton
      real(rk) ::   Nitrogen_UpPHY	  ! mmol N m-3, Nitrogen uptake of phytoplankton
      real(rk) ::   Nutrient_UpPHY	  ! mmol m-3, Nutrient uptake of phytoplankton
      real(rk) ::   NutrientLimitationDiatoms	  ! -, Nutrient limitation diatoms
      real(rk) ::   Phosphate_upDiatoms	  ! mmol P m-3, Diatoms phosphate uptake
      real(rk) ::   PHYMort	  ! mmol m-3, Phytoplankton mortality rate
      real(rk) ::   SiCrat	  ! ?, ?
      real(rk) ::   SiCrDiatoms	  ! mol Si mol C-1, Si/C ratio in diatoms
      real(rk) ::   Silicate_upDia	  ! mmol Si m-3, Diatoms silicate uptake
      real(rk) ::   tf	  ! -, Temperature factor
      real(rk) ::   tfsilicate	  ! -, Silicate dissolution adjusted for temperature
      real(rk) ::   TotalRespirationPHY	  ! mmol C m-3, Total phytoplankton respiration (basal & activity)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
   _GET_(self%id_ndi,NDI)       ! Diatom biomass in nitrogen
   _GET_(self%id_sid,SID)       ! Detrital silicate concentration
   _GET_(self%id_sio,SIO)       ! Silicilic acid concentration
   _GET_(self%id_dcl,DCL)       ! Labile detritus concentration in carbon
   _GET_(self%id_dcs,DCS)       ! Semi-labile detritus concentration in carbon
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dnl,DNL)       ! Labile detritus concentration in nitrogen
   _GET_(self%id_dns,DNS)       ! Semi-labile detritus concentration in nitrogen
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_nos,NOS)       ! Nitrate concentration
   _GET_(self%id_pho,PHO)       ! Phosphorus
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
    
    tf = Q10Factor (temp,Q10PhyDiatoms)
    tfsilicate = Q10Factor(temp,Q10SilicateDiss)
    
   ! Calculate ratios in phytoplankton 
    NCrDiatoms  = NDI/CDI
    ChlCrDiatoms = ChlCrPHY(NCratio,MaxNCr,MinNCr,MinChlNr,MaxChlNr)
    
             ! chlorophyll(i,j,k,2) = ChlCrDiatoms * CDI		! REMOVE (POSSIBLY) 
    
   ! Nitrate uptake rate 
    Nitrate_UpPHY = NOuptake(NCrDiatoms,tf,MaxNCrDiatoms,NosMaxUptakeDiatoms) * Michaelis(NOS,ksNOsDiatoms) * Inhibition(NHS,kinNHsPhy) * CDI
    
   ! Ammonium uptake rate 
    Ammonium_UpPHY = NUT_UPTAKE_RATE(NCrDiatoms,(NHS-0.0),tf,MaxNCrDiatoms,NHsMaxUptakeDiatoms,ksNHsDiatoms) * CDI
    
   ! Phosphate uptake rate 
    Phosphate_upDiatoms = NUT_UPTAKE_RATE(NCrDiatoms,(PHO-0.0),tf,MaxNCrDiatoms,PO4MaxUptakeDiatoms,ksPO4Diatoms) * CDI
    
   ! Compute silicate uptake 
    MaxSiCrDiatoms = self%MaxNCrDiatoms*self%SiNrDiatoms
    MinSiCrDiatoms = self%MinNCrDiatoms*self%SiNrDiatoms
    SiCrDiatoms = NCrDiatoms*self%SiNrDiatoms
    Silicate_upDia = NUT_UPTAKE_RATE(SiCrDiatoms,(SIO-0.0),tf,MaxSiCrDiatoms,SiMaxUptakeDiatoms,ksSiDiatoms) * CDI
    
   ! Potential nitrogen uptake 
    Nitrogen_UpPHY = Ammonium_UpPHY + Nitrate_UpPHY
    
   ! Nutrient uptake 
    Nutrient_UpPHY = min(Nitrogen_UpPHY,Silicate_upDia/SiNrDiatoms,Phosphate_upDiatoms/self%PNRedfield)
    
   ! Compute actual N:C and Si:C ratios 
    NCrat = Ratio_PHYT(NCrDiatoms,MaxNCrDiatoms,MinNCrDiatoms)
    SiCrat = Ratio_PHYT(SiCrDiatoms,MaxSiCrDiatoms,MinSiCrDiatoms)
    
   ! Compute nutrient and light limitation 
    NutrientLimitationDiatoms = Nutr_LIM(MinNCrDiatoms,MinSiCrDiatoms,NCrat,SiCrat)
          LightLimitationDiatoms = 1.-exp(-self%alphaPIDiatoms*PAR/self%MuMaxDiatoms)
    
   ! Compute carbon uptake 
    Carbon_UptakePHY = self%MuMaxDiatoms*LightLimitationDiatoms*NutrientLimitationDiatoms*CDI*tf
    
   ! Compute respiration 
    TotalRespirationPHY = Carbon_UptakePHY*self%GrowthRespDiatoms + self%RespirationDiatoms*CDI*tf
    
   ! Compute growth 
    GrowthPHY = Carbon_UptakePHY-TotalRespirationPHY
    
   ! Compute the extra DOC excretion from the difference of growth in nutrient limited (actual NC ratio) and nutrient saturated (max NC ratio)        
    DOC_extra_excr = CDI * tf * self%extradocphyexcr * self%MuMaxDiatoms * ExtraEXCR_term(LightLimitationDiatoms,MinNCrDiatoms,MaxNCrDiatoms,NCrat)    
    
   ! Compute the leakage 
    DOC_leakage = self%leakagephy * Carbon_UptakePHY
    DON_leakage = self%leakagephy * abs(Nutrient_UpPHY)
    
   ! Compute mortality 
    C_PHYMort  = PHYMORT(MortalityDiatoms,tf) * CDI
    N_PHYMort  = PHYMORT(MortalityDiatoms,tf) * NDI
    
   !Carbon in Diatoms increases by growth and decreases by leakage and mortality 
   _ADD_SOURCE_(self%id_cdi,1.0*( GrowthPHY)) 
   _ADD_SOURCE_(self%id_cdi,-1.0*( C_PHYMort + DOC_leakage)) 
   _ADD_SOURCE_(self%id_ndi,1.0*( Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_ndi,-1.0*( N_PHYMort + DON_leakage)) 
    
   ! IF CN ratio of phytoplankton higher than CNmin, than nitrogen is taken up, unless it gets excreted 
             IF (Nutrient_UpPHY.gt.0) THEN 
   _ADD_SOURCE_(self%id_nos,-1.0*( Nutrient_UpPHY*Nitrate_UpPHY/Nitrogen_UpPHY)) 
   _ADD_SOURCE_(self%id_dox,1.0*( Nutrient_UpPHY*Nitrate_UpPHY/Nitrogen_UpPHY*self%ONoxnhsr)) 
   _ADD_SOURCE_(self%id_nhs,-1.0*( Nutrient_UpPHY*Ammonium_UpPHY/Nitrogen_UpPHY)) 
   _ADD_SOURCE_(self%id_pho,-1.0*( Nutrient_UpPHY*self%PNRedfield)) 
             ELSE 
   _ADD_SOURCE_(self%id_nhs,1.0*( - Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_pho,1.0*( - Nutrient_UpPHY*self%PNRedfield)) 
             ENDIF 
    
   ! Update Silicate concentration 
    Silicate_upDia = (Nutrient_UpPHY - DON_leakage)*self%SiNrDiatoms
   _ADD_SOURCE_(self%id_sio,-1.0*( Silicate_upDia)) 
   _ADD_SOURCE_(self%id_sio,1.0*( tfsilicate*self%kdis_Silicious_Detritus*SID)) 
   _ADD_SOURCE_(self%id_sid,-1.0*( tfsilicate*self%kdis_Silicious_Detritus*SID)) 
   _ADD_SOURCE_(self%id_sid,1.0*( N_PHYMort*self%SiNrDiatoms)) 
    
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
   _ADD_SOURCE_(self%id_dic,-1.0*(GrowthPHY+DOC_extra_excr)) 
    
#ifdef primaryprod 
          NPP = NPP+Carbon_UptakePHY-TotalRespirationPHY
#endif 
#ifdef biodiag2 
          Nitrogen_Uptake_Diatoms = Nutrient_UpPHY
          Carbon_UptakeDiatoms    = Carbon_UptakePHY
          TotalRespirationDiatoms = TotalRespirationPHY
          Silicate_upDiatoms      = Silicate_upDia
#endif 
#ifdef biodiag1 
          PhytoNitrateReduction   = PhytoNitrateReduction+Nutrient_UpPHY*ratio(Nitrate_upPHY,Nitrogen_UpPHY)*self%ONoxnhsr
#endif 

   ! Averaged over entire water column Diagnostics 
#ifdef biodiag2 
   _SET_DIAGNOSTIC_(self%id_Nitrogen_Uptake_Diatoms, Nitrogen_Uptake_Diatoms)
   _SET_DIAGNOSTIC_(self%id_Silicate_upDiatoms, Silicate_upDiatoms)
   _SET_DIAGNOSTIC_(self%id_Carbon_UptakeDiatoms, Carbon_UptakeDiatoms)
   _SET_DIAGNOSTIC_(self%id_TotalRespirationDiatoms, TotalRespirationDiatoms)
#endif 
   _LOOP_END_

   end subroutine do


   end module fabm_ulg_Diatoms 
