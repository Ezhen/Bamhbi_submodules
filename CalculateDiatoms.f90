#include "fabm_driver.h" 
 
!#########################################################################################
!                              3Ddiatoms.F90
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
! Implementation: Marilaure Gregoire,                 NIOO-CEME
! Translation in FABM: Evgeny Ivanov,   Universite de Liege, MAST
!
!######################################################################

   module fabm_ulg_Diatoms 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Diatoms 
      type (type_state_variable_id)         :: id_cdi,id_ndi,id_sid,id_sio
      type (type_state_variable_id)         :: id_agg,id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeDiatoms,id_Nitrogen_Uptake_Diatoms,id_NPP,id_PhytoNitrateReduction,id_Silicate_upDiatoms,id_TotalRespirationDiatoms
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeDiatomsIntegrated,id_Nitrogen_Uptake_DiatomsIntegrated,id_Silicate_upDiatomsIntegrated,id_TotalRespirationDiatomsIntegrated

!     Model parameters 
      real(rk) :: alphaPIDiatoms
      real(rk) :: extradocphyexcr
      real(rk) :: GrowthRespDiatoms
      real(rk) :: kdis_Silicious_Detritus
      real(rk) :: kinNHsPhy
      real(rk) :: ksNHsDiatoms
      real(rk) :: ksNOsDiatoms
      real(rk) :: ksPO4Diatoms
      real(rk) :: ksSiDiatoms
      real(rk) :: labileextradocphyexcr
      real(rk) :: labilefraction
      real(rk) :: leakagephy
      real(rk) :: MaxChlNrDiatoms
      real(rk) :: MaxNCrDiatoms
      real(rk) :: MinChlNrDiatoms
      real(rk) :: MinNCrDiatoms
      real(rk) :: MortalityDiatoms
      real(rk) :: mortphydom
      real(rk) :: MuMaxDiatoms
      real(rk) :: NHsMaxUptakeDiatoms
      real(rk) :: NosMaxUptakeDiatoms
      real(rk) :: OCr
      real(rk) :: ONoxnhsr
      real(rk) :: PNRedfield
      real(rk) :: PO4MaxUptakeDiatoms
      real(rk) :: Q10PhyDiatoms
      real(rk) :: Q10SilicateDiss
      real(rk) :: QuantumYieldDiatoms
      real(rk) :: RespirationDiatoms
      real(rk) :: SiMaxUptakeDiatoms
      real(rk) :: SiNrDiatoms

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

   real(rk)     :: alphaPIDiatoms=0.3312/daytosecond
   real(rk)     :: extradocphyexcr=0.05
   real(rk)     :: GrowthRespDiatoms=0.1
   real(rk)     :: kdis_Silicious_Detritus=0.08/daytosecond
   real(rk)     :: kinNHsPhy=0.5
   real(rk)     :: ksNHsDiatoms=1.0
   real(rk)     :: ksNOsDiatoms=1.0
   real(rk)     :: ksPO4Diatoms=0.1
   real(rk)     :: ksSiDiatoms=3.5
   real(rk)     :: labileextradocphyexcr=0.65
   real(rk)     :: labilefraction=0.7
   real(rk)     :: leakagephy=0.02
   real(rk)     :: MaxChlNrDiatoms=2.0
   real(rk)     :: MaxNCrDiatoms=0.2
   real(rk)     :: MinChlNrDiatoms=1.0
   real(rk)     :: MinNCrDiatoms=0.05
   real(rk)     :: MortalityDiatoms=0.03/daytosecond
   real(rk)     :: mortphydom=0.34
   real(rk)     :: MuMaxDiatoms=3.50/daytosecond
   real(rk)     :: NHsMaxUptakeDiatoms=1.0/daytosecond
   real(rk)     :: NosMaxUptakeDiatoms=1.0/daytosecond
   real(rk)     :: OCr=1.0
   real(rk)     :: ONoxnhsr=2.0
   real(rk)     :: PNRedfield=1.0/16.0
   real(rk)     :: PO4MaxUptakeDiatoms=1.0/16.0/daytosecond
   real(rk)     :: Q10PhyDiatoms=1.8
   real(rk)     :: Q10SilicateDiss=3.3
   real(rk)     :: QuantumYieldDiatoms=0.8
   real(rk)     :: RespirationDiatoms=0.009/daytosecond
   real(rk)     :: SiMaxUptakeDiatoms=0.50/daytosecond
   real(rk)     :: SiNrDiatoms=5./6.

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
   call self%get_parameter(self%alphaPIDiatoms, 'alphaPIDiatoms', default=alphaPIDiatoms) 
   call self%get_parameter(self%extradocphyexcr, 'extradocphyexcr', default=extradocphyexcr) 
   call self%get_parameter(self%GrowthRespDiatoms, 'GrowthRespDiatoms', default=GrowthRespDiatoms) 
   call self%get_parameter(self%kdis_Silicious_Detritus, 'kdis_Silicious_Detritus', default=kdis_Silicious_Detritus) 
   call self%get_parameter(self%kinNHsPhy, 'kinNHsPhy', default=kinNHsPhy) 
   call self%get_parameter(self%ksNHsDiatoms, 'ksNHsDiatoms', default=ksNHsDiatoms) 
   call self%get_parameter(self%ksNOsDiatoms, 'ksNOsDiatoms', default=ksNOsDiatoms) 
   call self%get_parameter(self%ksPO4Diatoms, 'ksPO4Diatoms', default=ksPO4Diatoms) 
   call self%get_parameter(self%ksSiDiatoms, 'ksSiDiatoms', default=ksSiDiatoms) 
   call self%get_parameter(self%labileextradocphyexcr, 'labileextradocphyexcr', default=labileextradocphyexcr) 
   call self%get_parameter(self%labilefraction, 'labilefraction', default=labilefraction) 
   call self%get_parameter(self%leakagephy, 'leakagephy', default=leakagephy) 
   call self%get_parameter(self%MaxChlNrDiatoms, 'MaxChlNrDiatoms', default=MaxChlNrDiatoms) 
   call self%get_parameter(self%MaxNCrDiatoms, 'MaxNCrDiatoms', default=MaxNCrDiatoms) 
   call self%get_parameter(self%MinChlNrDiatoms, 'MinChlNrDiatoms', default=MinChlNrDiatoms) 
   call self%get_parameter(self%MinNCrDiatoms, 'MinNCrDiatoms', default=MinNCrDiatoms) 
   call self%get_parameter(self%MortalityDiatoms, 'MortalityDiatoms', default=MortalityDiatoms) 
   call self%get_parameter(self%mortphydom, 'mortphydom', default=mortphydom) 
   call self%get_parameter(self%MuMaxDiatoms, 'MuMaxDiatoms', default=MuMaxDiatoms) 
   call self%get_parameter(self%NHsMaxUptakeDiatoms, 'NHsMaxUptakeDiatoms', default=NHsMaxUptakeDiatoms) 
   call self%get_parameter(self%NosMaxUptakeDiatoms, 'NosMaxUptakeDiatoms', default=NosMaxUptakeDiatoms) 
   call self%get_parameter(self%OCr, 'OCr', default=OCr) 
   call self%get_parameter(self%ONoxnhsr, 'ONoxnhsr', default=ONoxnhsr) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', default=PNRedfield) 
   call self%get_parameter(self%PO4MaxUptakeDiatoms, 'PO4MaxUptakeDiatoms', default=PO4MaxUptakeDiatoms) 
   call self%get_parameter(self%Q10PhyDiatoms, 'Q10PhyDiatoms', default=Q10PhyDiatoms) 
   call self%get_parameter(self%Q10SilicateDiss, 'Q10SilicateDiss', default=Q10SilicateDiss) 
   call self%get_parameter(self%QuantumYieldDiatoms, 'QuantumYieldDiatoms', default=QuantumYieldDiatoms) 
   call self%get_parameter(self%RespirationDiatoms, 'RespirationDiatoms', default=RespirationDiatoms) 
   call self%get_parameter(self%SiMaxUptakeDiatoms, 'SiMaxUptakeDiatoms', default=SiMaxUptakeDiatoms) 
   call self%get_parameter(self%SiNrDiatoms, 'SiNrDiatoms', default=SiNrDiatoms) 

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
   call self%register_state_dependency(self%id_agg, 'Aggregates', 'm-3') 
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
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  AGG,DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CDI,NDI,SID,SIO
      real(rk) ::   Carbon_UptakeDiatoms,Nitrogen_Uptake_Diatoms,NPP,PhytoNitrateReduction,Silicate_upDiatoms,TotalRespirationDiatoms
      real(rk) ::   Carbon_UptakeDiatomsIntegrated,Nitrogen_Uptake_DiatomsIntegrated,Silicate_upDiatomsIntegrated,TotalRespirationDiatomsIntegrated
      real(rk) ::   Ammonium_UpPHY	 + ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   C_PHYMort	 + ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Carbon_UptakePHY	 + ! mmol C m-3, C assimilation of phytoplankton
      real(rk) ::   ChlCrDiatoms	 + ! g Chla mol C-1, Chl/C ratio in large flagellates
      real(rk) ::   DOC_extra_excr	 + ! mmol C d-1, Phytoplankton extra excretion
      real(rk) ::   DOC_leakage	 + ! mmol C d-1, Phytoplankton passive leakage rate for carbon
      real(rk) ::   DON_leakage	 + ! mmol N d-1, Phytoplankton passive leakage rate for nitrogen
      real(rk) ::   GrowthPHY	 + ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   LightLimitationDiatoms	 + ! ?, Light limitation diatoms
      real(rk) ::   MaxSiCrDiatoms	 + ! mol Si mol C-1, Maximum Si/C ratio in diatoms
      real(rk) ::   MinSiCrDiatoms	 + ! mol Si mol C-1, Minimum Si/C ratio in diatoms
      real(rk) ::   Mu_Nitrogen	 + ! ?, ?
      real(rk) ::   Mu_Silicate	 + ! ?, ?
      real(rk) ::   NCrat	 + ! ?, ?
      real(rk) ::   N_PHYMort	 + ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   NCrDiatoms	 + ! mol N mol C-1, N/C ratio in diatoms
      real(rk) ::   Nitrate_UpPHY	 + ! mmol N m-3, Nitrate uptake of phytoplankton
      real(rk) ::   Nitrogen_UpPHY	 + ! mmol N m-3, Nitrogen uptake of phytoplankton
      real(rk) ::   Nutrient_UpPHY	 + ! mmol m-3, Nutrient uptake of phytoplankton
      real(rk) ::   NutrientLimitationDiatoms	 + ! ?, Nutrient limitation diatoms
      real(rk) ::   Phosphate_upDiatoms	 + ! mmol P m-3, Diatoms phosphate uptake
      real(rk) ::   PHYMort	 + ! mmol m-3, Phytoplankton mortality rate
      real(rk) ::   SiCrat	 + ! ?, ?
      real(rk) ::   SiCrDiatoms	 + ! mol Si mol C-1, Si/C ratio in diatoms
      real(rk) ::   Silicate_upDia	 + ! mmol Si m-3, Diatoms silicate uptake
      real(rk) ::   tf	 + ! -, Temperature factor
      real(rk) ::   tfsilicate	 + ! ?, Silicate dissolution adjusted for temperature
      real(rk) ::   TotalRespirationPHY	 + ! mmol C m-3, Total phytoplankton respiration (basal & activity)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_agg,AGG)       ! Aggregates
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
   ! TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor (temp,Q10PhyDiatoms)
    tfsilicate = Q10Factor (temp,Q10SilicateDiss)
   ! PHYTOPLANKTON 
   ! N/C ratio 
    NCrDiatoms  = Ratio (NDI,CDI)
             !CNrDiatoms  = Ratio (CDI,NDI) ! not used ! REMOVE (POSSIBLY)
   ! Chlorophyll to carbon ratio (ChlCrPHY, mg Chl/mol C) 
             CALL CHL_C_RATIO(NCrDiatoms,MaxNCrDiatoms,MinNCrDiatoms,MinChlNrDiatoms,MaxChlNrDiatoms,ChlCrDiatoms) ! REMOVE (POSSIBLY)
    chlorophyll(i,j,k,2) = ChlCrDiatoms * CDI
   ! Nitrate uptake rates of phytoplankton (NO3_upPHY, mmol N/m3/day) 
   CALL NO_UPTAKE_RATE(NCrDiatoms,NOS,NHS,CDI,tf,MaxNCrDiatoms,NosMaxUptakeDiatoms,ksNOsDiatoms,kinNHsPhy, Nitrate_UpPHY) ! REMOVE (POSSIBLY)
   ! Ammonium uptake rates of phytoplankton (NH3_upPHY, mmolN /m3/day) 
   ! (and excretion if NC ratio too high) 
   CALL NUT_UPTAKE_RATE(NCrDiatoms,NHS,CDI,tf,MaxNCrDiatoms,NHsMaxUptakeDiatoms,ksNHsDiatoms,0.0_wp, Ammonium_UpPHY) ! REMOVE (POSSIBLY)
             if (IncludeSIlicate) THEN ! REMOVE (POSSIBLY)
    MaxSiCrDiatoms = self%MaxNCrDiatoms*self%SiNrDiatoms
    MinSiCrDiatoms = self%MinNCrDiatoms*self%SiNrDiatoms
    SiCrDiatoms = NCrDiatoms*self%SiNrDiatoms
   !Compute Silicate uptake 
   CALL NUT_UPTAKE_RATE(SiCrDiatoms,SIO,CDI,tf,MaxSiCrDiatoms,SiMaxUptakeDiatoms,ksSiDiatoms,0.0_wp, Silicate_upDia) ! REMOVE (POSSIBLY)
             ELSE ! REMOVE (POSSIBLY)
    MaxSiCrDiatoms = 1.
    MinSiCrDiatoms = 0.
    SiCrDiatoms = 1.
             ENDIF ! REMOVE (POSSIBLY)
   CALL NUT_UPTAKE_RATE(NCrDiatoms,PHO,CDI,tf,MaxNCrDiatoms,PO4MaxUptakeDiatoms,ksPO4Diatoms,0.0_wp, Phosphate_upDiatoms) ! REMOVE (POSSIBLY)
    Nitrogen_UpPHY=Ammonium_UpPHY + Nitrate_UpPHY
    Nutrient_UpPHY=min(Nitrogen_UpPHY,Silicate_upDia/SiNrDiatoms,Phosphate_upDiatoms/self%PNRedfield)
   ! Compute the diatoms growth taking into account the potential limitimg effect of silicate 
   ! ################################# FORMER SUBROUTINE PHYT GROWTH RATE ########################################### 
   ! Potential nutrient controlled growth rate; due to droop kinetics, the actual NC ratio can slightly surpass the maximal or minimal ratio 
             IF(NCrDiatoms > MaxNCrDiatoms) THEN ! REMOVE (POSSIBLY)
    NCrat = self%MaxNCrDiatoms
             ELSEIF (NCrDiatoms < MinNCrDiatoms) THEN ! REMOVE (POSSIBLY)
    NCrat = self%MinNCrDiatoms
             ELSE ! REMOVE (POSSIBLY)
    NCrat = NCrDiatoms
             ENDIF ! REMOVE (POSSIBLY)
   ! the actual SiC ratio can slightly surpass the maximal or minimal ratio 
             IF(SiCrDiatoms > MaxSiCrDiatoms) THEN ! REMOVE (POSSIBLY)
    SiCrat = MaxSiCrDiatoms
             ELSEIF (SiCrDiatoms < MinSiCrDiatoms) THEN ! REMOVE (POSSIBLY)
    SiCrat = MinSiCrDiatoms
             ELSE ! REMOVE (POSSIBLY)
    SiCrat = SiCrDiatoms
             ENDIF ! REMOVE (POSSIBLY)
   ! Tett model 
    Mu_Nitrogen = 1. - self%MinNCrDiatoms  / NCrat
    Mu_Silicate = 1. - MinSiCrDiatoms / SiCrat
    NutrientLimitationDiatoms = min(Mu_Nitrogen,Mu_Silicate)
          LightLimitationDiatoms = 1.-exp(-self%alphaPIDiatoms*PAR/self%MuMaxDiatoms)
    Carbon_UptakePHY=self%MuMaxDiatoms*LightLimitationDiatoms*NutrientLimitationDiatoms*CDI*tf
    TotalRespirationPHY=Carbon_UptakePHY*self%GrowthRespDiatoms + self%RespirationDiatoms*CDI*tf
    GrowthPHY=Carbon_UptakePHY-TotalRespirationPHY
   ! Compute the extra DOC excretion as a fraction (extradocphyexcr) of the difference of growth in nutrient limited (actual NC ratio) and nutrient saturated (max NC ratio) as in VDM et al 2004 L&O 
    DOC_extra_excr=abs(extradocphyexcr*self%MuMaxDiatoms*CDI*tf*(min(LightLimitationDiatoms,(1. - self%self%MinNCrDiatoms  / self%MaxNCrDiatoms))- min(LightLimitationDiatoms,(1. - self%self%MinNCrDiatoms  / NCrat))))
   ! ################################# END OF FORMER SUBROUTINE PHYT GROWTH RATE ########################################### 
   !Compute the leakage and extra DOC excretion, DOC_extra_excr 
    DOC_leakage = self%leakagephy*Carbon_UptakePHY
    DON_leakage  =self%leakagephy*abs(Nutrient_UpPHY)
   ! Phytoplankton mortality rate (PHYMort, /day) 
             CALL PHYMORT_RATE(MortalityDiatoms,tf, PHYmort) ! REMOVE (POSSIBLY)
   ! Phytoplankton mortality flux C_PHYMort,N_PHYMort (in mmol C/m3/day or mmol N/m3/day) 
    C_PHYMort  = PHYMort * CDI
    N_PHYMort  = PHYMort * NDI
   ! ADJUSTING THE RATE OF CHANGE 
   ! phytoplankton C increases by growth, 
   ! it decreases by zooplankton grazing (see the zooplankton subroutine) and phytoplankton mortality 
   _ADD_SOURCE_(self%id_cdi,1.0*(  GrowthPHY)) 
   _ADD_SOURCE_(self%id_cdi,-1.0*( C_PHYMort+DOC_leakage)) 
   ! phytoplankton N increases by N uptake, 
   ! it decreases by zooplankton grazing and phytoplankton mortality and leakage 
             ! NetNGrowthDiatoms=  Nutrient_UpPHY - N_PHYMort-DON_leakage ! luc commented this, it doesn't seem to be used ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_ndi,1.0*(  Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_ndi,-1.0*(  N_PHYMort+DON_leakage)) 
   ! IF CN ratio of phytoplankton not lower than CNmin, than nitrogen is taken up 
   ! Nitrate is taken up by phytoplankton 
             IF (Nutrient_UpPHY.gt.0) THEN ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_nos,-1.0*( Nutrient_UpPHY*Nitrate_UpPHY/Nitrogen_UpPHY)) 
   _ADD_SOURCE_(self%id_dox,1.0*( Nutrient_UpPHY*Nitrate_UpPHY/Nitrogen_UpPHY*self%ONoxnhsr)) 
   ! Ammonium is taken up by phytoplankton 
   _ADD_SOURCE_(self%id_nhs,-1.0*( Nutrient_UpPHY*Ammonium_UpPHY/Nitrogen_UpPHY)) 
   _ADD_SOURCE_(self%id_pho,-1.0*( Nutrient_UpPHY*self%PNRedfield)) 
             ELSE ! REMOVE (POSSIBLY)
   ! IF CN ratio of phytoplankton is lower than CNmin, than ammonium is excreted 
   _ADD_SOURCE_(self%id_nhs,1.0*( - Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_pho,1.0*( - Nutrient_UpPHY*self%PNRedfield)) 
             ENDIF ! REMOVE (POSSIBLY)
   ! UPdate Silicate concentration 
             if (IncludeSIlicate) then ! REMOVE (POSSIBLY)
    Silicate_upDia= (Nutrient_UpPHY - DON_leakage)*self%SiNrDiatoms
   _ADD_SOURCE_(self%id_sio,-1.0*( Silicate_upDia)) 
   _ADD_SOURCE_(self%id_sio,1.0*( tfsilicate*self%kdis_Silicious_Detritus*SID)) 
   _ADD_SOURCE_(self%id_sid,-1.0*(tfsilicate*self%kdis_Silicious_Detritus*SID)) 
   _ADD_SOURCE_(self%id_sid,1.0*( N_PHYMort*self%SiNrDiatoms)) 
             ENDIF ! REMOVE (POSSIBLY)
   ! As in Anderson and Pondhaven (2003), the phytoplanton mortality increases the pool of POM and DOM with a coefficient of mortdom 
   ! for the DOM pool a part is considered as labile (labilefraction) and another part as semi labile (1-labilefraction) 
   _ADD_SOURCE_(self%id_poc,1.0*( (1.0 - self%mortphydom)*C_PHYMort)) 
   _ADD_SOURCE_(self%id_pon,1.0*( (1.0 - self%mortphydom)*N_PHYMort)) 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%mortphydom*C_PHYMort*self%labilefraction)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( self%mortphydom*N_PHYMort*self%labilefraction)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( self%mortphydom*C_PHYMort*(1.0 - self%labilefraction))) 
   _ADD_SOURCE_(self%id_dns,1.0*( self%mortphydom*N_PHYMort*(1.0 - self%labilefraction))) 
   ! As in Anderson and Pondhaven (2003), the DCLI and DNLI concentration increases also due to phytoplankton leakage which is considered 
   ! to peoduced only labile DOC and DON 
   _ADD_SOURCE_(self%id_dcl,1.0*( DOC_leakage)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( DON_leakage)) 
   ! As in Anderson and Pondhaven (2003), the phytoplankton extra DOC excretion, DOC_extra_excr, increases the pool of DOC (only DOC and not DON) 
   ! this extra-excretion is composed of a part (leakage) which is assimilated to leakage and is thus only labile 
   ! the remaining part (1-leakage) is composed of a labile (labileextradocexcr) and semi labile part (1-labileextradocexcr) 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%self%leakagephy*DOC_extra_excr + self%labileextradocphyexcr*(1.0 - self%self%leakagephy)*DOC_extra_excr)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( (1.0 - self%labileextradocphyexcr)*(1.0 - self%leakagephy)*DOC_extra_excr)) 
             SELECT CASE (SinkingVelocityType) ! REMOVE (POSSIBLY)
             CASE ('aggregation') ! REMOVE (POSSIBLY)
   ! The number of aggregates, POMNOS including PON increases and decreases with PON 
               !     dPAGGI(i,j,k) = dPAGGI(i,j,k) + (1.0 - mortphydom)*N_PHYMort*AGGI(i,j,k)/(PONI(i,j,k) +PNSI(i,j,k)) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,1.0*( (1.0 - self%mortphydom)*N_PHYMort*AGG/(PON))) 
#ifdef nanquest 
               if (isnan(dPAGG(I,J,K))) then ! REMOVE (POSSIBLY)
                 write (*,*) '** NAN QUEST ** in CalcDIATOMS' ! REMOVE (POSSIBLY)
                 write (*,*) 'i,j,k,PON(I,J,K),AGG(I,J,K)',i,j,k,TRB(I,J,K,PON),trb(I,J,K,AGG) ! REMOVE (POSSIBLY)
                 call flush(6) ! REMOVE (POSSIBLY)
                 stop ! REMOVE (POSSIBLY)
            endif 
#endif 
             END SELECT ! REMOVE (POSSIBLY)
   ! The oxygen concentration increases due to photosynthesis 
   _ADD_SOURCE_(self%id_dox,1.0*( (GrowthPHY + DOC_extra_excr)*self%OCr)) 
   ! CO2 production and consumption 
   _ADD_SOURCE_(self%id_dic,-1.0*(GrowthPHY+DOC_extra_excr)) 
   !Store Diagnostics the fluxes 
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
           end if ! REMOVE (POSSIBLY)
         end do ! REMOVE (POSSIBLY)
       end do ! REMOVE (POSSIBLY)
     end do ! REMOVE (POSSIBLY)
   ! OUTPUT VARIABLES 
   ! Averaged over entire water column Diagnostics 
#ifdef biodiag2 
   _SET_DIAGNOSTIC_(self%id_Nitrogen_Uptake_Diatoms, Nitrogen_Uptake_Diatoms)
   _SET_DIAGNOSTIC_(self%id_Silicate_upDiatoms, Silicate_upDiatoms)
   _SET_DIAGNOSTIC_(self%id_Carbon_UptakeDiatoms, Carbon_UptakeDiatoms)
   _SET_DIAGNOSTIC_(self%id_TotalRespirationDiatoms, TotalRespirationDiatoms)
#endif 
   _LOOP_END_

   end subroutine do

