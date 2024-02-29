#include "fabm_driver.h" 
 
!#########################################################################################
!                              3DflagellatesF90
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
! Translation into FABM : Evgeny Ivanov, ULg / MAST
!
!--------------------------------------------------------------------
! Contains the pelagic submodel, as used in Soetaert et al., 2001.
! References to the microplankton model:
! Tett, P., 1998. Parameterising a microplankton model.
! Department of Biological Sciences, Napier University,
! Report ISBN 0 902703 60 9, 60 pp.                                                                                                                                                                                                                                           !
! Sharples Tett (1994). Modeling the effect of physical! variability on the midwater chlorophyll maximum.
! Journal of marine research 52: 219-238
!
!--------------------------------------------------------------------*

   module fabm_ulg_Flagellates 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Flagellates 
      type (type_state_variable_id)         :: id_cfl,id_nfl
      type (type_state_variable_id)         :: id_agg,id_cem,id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nem,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeFlagellates,id_Nitrogen_Uptake_Flagellates,id_NPP,id_PhytoNitrateReduction,id_TotalRespirationFlagellates
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeFlagellatesIntegrated,id_Nitrogen_Uptake_FlagellatesIntegrated,id_TotalRespirationFlagellatesIntegrated

!     Model parameters 
      real(rk) :: alphaPIFlagellates
      real(rk) :: extradocphyexcr
      real(rk) :: GrowthRespFlagellates
      real(rk) :: kinNHsPhy
      real(rk) :: ksNHsFlagellates
      real(rk) :: ksNOsFlagellates
      real(rk) :: ksPO4Flagellates
      real(rk) :: labileextradocphyexcr
      real(rk) :: labilefraction
      real(rk) :: leakagephy
      real(rk) :: MaxChlNrFlagellates
      real(rk) :: MaxNCrFlagellates
      real(rk) :: MinChlNrFlagellates
      real(rk) :: MinNCrFlagellates
      real(rk) :: MortalityFlagellates
      real(rk) :: mortphydom
      real(rk) :: MuMaxFlagellates
      real(rk) :: NHsMaxUptakeFlagellates
      real(rk) :: NosMaxUptakeFlagellates
      real(rk) :: OCr
      real(rk) :: ONoxnhsr
      real(rk) :: PNRedfield
      real(rk) :: PO4MaxUptakeFlagellates
      real(rk) :: Q10Phy
      real(rk) :: QuantumYieldFlagellates
      real(rk) :: RespirationFlagellates

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
   ! Initialise the Flagellates model

   subroutine initialize(self,configunit)
   class (type_ulg_Flagellates), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   real(rk)     :: alphaPIFlagellates=0.2153/daytosecond
   real(rk)     :: extradocphyexcr=0.05
   real(rk)     :: GrowthRespFlagellates=0.1
   real(rk)     :: kinNHsPhy=0.5
   real(rk)     :: ksNHsFlagellates=3.0
   real(rk)     :: ksNOsFlagellates=3.0
   real(rk)     :: ksPO4Flagellates=0.2
   real(rk)     :: labileextradocphyexcr=0.65
   real(rk)     :: labilefraction=0.7
   real(rk)     :: leakagephy=0.02
   real(rk)     :: MaxChlNrFlagellates=2.0
   real(rk)     :: MaxNCrFlagellates=0.2
   real(rk)     :: MinChlNrFlagellates=1.0
   real(rk)     :: MinNCrFlagellates=0.05
   real(rk)     :: MortalityFlagellates=0.03/daytosecond
   real(rk)     :: mortphydom=0.34
   real(rk)     :: MuMaxFlagellates=1./daytosecond
   real(rk)     :: NHsMaxUptakeFlagellates=0.50/daytosecond
   real(rk)     :: NosMaxUptakeFlagellates=0.50/daytosecond
   real(rk)     :: OCr=1.0
   real(rk)     :: ONoxnhsr=2.0
   real(rk)     :: PNRedfield=1.0/16.0
   real(rk)     :: PO4MaxUptakeFlagellates=0.5/16.0/daytosecond
   real(rk)     :: Q10Phy=2.0
   real(rk)     :: QuantumYieldFlagellates=0.6
   real(rk)     :: RespirationFlagellates=0.009/daytosecond

   namelist /ulg_Flagellates/ alphaPIFlagellates, 	 & 
                      extradocphyexcr, GrowthRespFlagellates, 	 & 
                      kinNHsPhy, ksNHsFlagellates, 	 & 
                      ksNOsFlagellates, ksPO4Flagellates, 	 & 
                      labileextradocphyexcr, labilefraction, 	 & 
                      leakagephy, MaxChlNrFlagellates, 	 & 
                      MaxNCrFlagellates, MinChlNrFlagellates, 	 & 
                      MinNCrFlagellates, MortalityFlagellates, 	 & 
                      mortphydom, MuMaxFlagellates, 	 & 
                      NHsMaxUptakeFlagellates, 	 & 
                      NosMaxUptakeFlagellates, OCr, ONoxnhsr, 	 & 
                      PNRedfield, PO4MaxUptakeFlagellates, 	 & 
                      Q10Phy, QuantumYieldFlagellates, 	 & 
                      RespirationFlagellates

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%alphaPIFlagellates, 'alphaPIFlagellates', default=alphaPIFlagellates) 
   call self%get_parameter(self%extradocphyexcr, 'extradocphyexcr', default=extradocphyexcr) 
   call self%get_parameter(self%GrowthRespFlagellates, 'GrowthRespFlagellates', default=GrowthRespFlagellates) 
   call self%get_parameter(self%kinNHsPhy, 'kinNHsPhy', default=kinNHsPhy) 
   call self%get_parameter(self%ksNHsFlagellates, 'ksNHsFlagellates', default=ksNHsFlagellates) 
   call self%get_parameter(self%ksNOsFlagellates, 'ksNOsFlagellates', default=ksNOsFlagellates) 
   call self%get_parameter(self%ksPO4Flagellates, 'ksPO4Flagellates', default=ksPO4Flagellates) 
   call self%get_parameter(self%labileextradocphyexcr, 'labileextradocphyexcr', default=labileextradocphyexcr) 
   call self%get_parameter(self%labilefraction, 'labilefraction', default=labilefraction) 
   call self%get_parameter(self%leakagephy, 'leakagephy', default=leakagephy) 
   call self%get_parameter(self%MaxChlNrFlagellates, 'MaxChlNrFlagellates', default=MaxChlNrFlagellates) 
   call self%get_parameter(self%MaxNCrFlagellates, 'MaxNCrFlagellates', default=MaxNCrFlagellates) 
   call self%get_parameter(self%MinChlNrFlagellates, 'MinChlNrFlagellates', default=MinChlNrFlagellates) 
   call self%get_parameter(self%MinNCrFlagellates, 'MinNCrFlagellates', default=MinNCrFlagellates) 
   call self%get_parameter(self%MortalityFlagellates, 'MortalityFlagellates', default=MortalityFlagellates) 
   call self%get_parameter(self%mortphydom, 'mortphydom', default=mortphydom) 
   call self%get_parameter(self%MuMaxFlagellates, 'MuMaxFlagellates', default=MuMaxFlagellates) 
   call self%get_parameter(self%NHsMaxUptakeFlagellates, 'NHsMaxUptakeFlagellates', default=NHsMaxUptakeFlagellates) 
   call self%get_parameter(self%NosMaxUptakeFlagellates, 'NosMaxUptakeFlagellates', default=NosMaxUptakeFlagellates) 
   call self%get_parameter(self%OCr, 'OCr', default=OCr) 
   call self%get_parameter(self%ONoxnhsr, 'ONoxnhsr', default=ONoxnhsr) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', default=PNRedfield) 
   call self%get_parameter(self%PO4MaxUptakeFlagellates, 'PO4MaxUptakeFlagellates', default=PO4MaxUptakeFlagellates) 
   call self%get_parameter(self%Q10Phy, 'Q10Phy', default=Q10Phy) 
   call self%get_parameter(self%QuantumYieldFlagellates, 'QuantumYieldFlagellates', default=QuantumYieldFlagellates) 
   call self%get_parameter(self%RespirationFlagellates, 'RespirationFlagellates', default=RespirationFlagellates) 

   ! Register state variables 

   call self%register_state_variable(self%id_cfl, 'CFL'  & 
         , 'mmol C m-3', 'Small flagellate biomass in carbon' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_nfl, 'NFL'  & 
         , 'mmol N m-3', 'Large flagellate biomass in nitrogen' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_dependency(self%id_agg, 'Aggregates', 'm-3') 
   call self%register_state_dependency(self%id_cem, 'Small flagellate biomass in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcl, 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_nem, 'Small flagellate biomass in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nhs, 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nos, 'Nitrate concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_pho, 'Phosphorus', 'mmol P m-3') 
   call self%register_state_dependency(self%id_poc, 'Particulate organic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_pon, 'Particulate organic nitrogen concentration', 'mmol N m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux) 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 


    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_Carbon_UptakeFlagellates, 'Carbon_UptakeFlagellates', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Carbon_UptakeFlagellatesIntegrated, 'Carbon_UptakeFlagellatesIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Nitrogen_Uptake_Flagellates, 'Nitrogen_Uptake_Flagellates', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_Nitrogen_Uptake_FlagellatesIntegrated, 'Nitrogen_Uptake_FlagellatesIntegrated', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_NPP, 'NPP', 'mmol N m-3 d-1', & 
      ' Primary production', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_PhytoNitrateReduction, 'PhytoNitrateReduction', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_TotalRespirationFlagellates, 'TotalRespirationFlagellates', '-', & 
      '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_TotalRespirationFlagellatesIntegrated, 'TotalRespirationFlagellatesIntegrated', '-', & 
      '-', output=output_instantaneous) 

   return 

99 call self%fatal_error('Flagellates', 'Error reading namelist ulg_Flagellates') 

   end subroutine initialize 


   ! Right hand sides of Flagellates model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  AGG,CEM,DCL,DCS,DIC,DNL,DNS,DOX,NEM,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CFL,NFL
      real(rk) ::   Carbon_UptakeFlagellates,Nitrogen_Uptake_Flagellates,NPP,PhytoNitrateReduction,TotalRespirationFlagellates
      real(rk) ::   Carbon_UptakeFlagellatesIntegrated,Nitrogen_Uptake_FlagellatesIntegrated,TotalRespirationFlagellatesIntegrated
      real(rk) ::   Ammonium_UpPHY	 + ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   C_PHYMort	 + ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Carbon_UptakePHY	 + ! mmol C m-3, C assimilation of phytoplankton
      real(rk) ::   ChlCrFlagellates	 + ! g Chla mol C-1, Chl/C ratio in large flagellates
      real(rk) ::   DOC_extra_excr	 + ! mmol C d-1, Phytoplankton extra excretion
      real(rk) ::   DOC_leakage	 + ! mmol C d-1, Phytoplankton passive leakage rate for carbon
      real(rk) ::   DON_leakage	 + ! mmol N d-1, Phytoplankton passive leakage rate for nitrogen
      real(rk) ::   GrowthPHY	 + ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   LightLimitationFlagellates	 + ! ?, Light limitation for flagellates
      real(rk) ::   Mu_Nitrogen	 + ! ?, ?
      real(rk) ::   Mu_Silicate	 + ! ?, ?
      real(rk) ::   N_PHYMort	 + ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   NCrat	 + ! ?, ?
      real(rk) ::   NCrFlagellates	 + ! mol N mol C-1, N/C ratio in large flagellates
      real(rk) ::   Nitrate_UpPHY	 + ! mmol N m-3, Nitrate uptake of phytoplankton
      real(rk) ::   Nitrogen_UpPHY	 + ! mmol N m-3, Nitrogen uptake of phytoplankton
      real(rk) ::   Nutrient_UpPHY	 + ! mmol m-3, Nutrient uptake of phytoplankton
      real(rk) ::   NutrientLimitationFlagellates	 + ! ?, Nutrient limitation for flagellates
      real(rk) ::   Phosphate_upFlagellates	 + ! mmol P m-3, Phosphate uptake by large flagellates
      real(rk) ::   PHYMort	 + ! mmol m-3, Phytoplankton mortality rate
      real(rk) ::   SiCrat	 + ! ?, ?
      real(rk) ::   tf	 + ! -, Temperature factor
      real(rk) ::   TotalRespirationPHY	 + ! mmol C m-3, Total phytoplankton respiration (basal & activity)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_agg,AGG)       ! Aggregates
   _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
   _GET_(self%id_dcl,DCL)       ! Labile detritus concentration in carbon
   _GET_(self%id_dcs,DCS)       ! Semi-labile detritus concentration in carbon
   _GET_(self%id_dic,DIC)       ! Dissolved inorganic carbon concentration
   _GET_(self%id_dnl,DNL)       ! Labile detritus concentration in nitrogen
   _GET_(self%id_dns,DNS)       ! Semi-labile detritus concentration in nitrogen
   _GET_(self%id_dox,DOX)       ! Dissolved oxygen concentration
   _GET_(self%id_nem,NEM)       ! Small flagellate biomass in nitrogen
   _GET_(self%id_nhs,NHS)       ! Ammonium concentration
   _GET_(self%id_nos,NOS)       ! Nitrate concentration
   _GET_(self%id_pho,PHO)       ! Phosphorus
   _GET_(self%id_poc,POC)       ! Particulate organic carbon concentration
   _GET_(self%id_pon,PON)       ! Particulate organic nitrogen concentration

   ! Retrieve current environmental conditions.
    _GET_(self%id_par,par)              ! local photosynthetically active radiation
    _GET_(self%id_temp,temp)            ! local temperature
   ! TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor(temp,Q10Phy)
   ! PHYTOPLANKTON 
   ! N/C ratio 
    NCrFlagellates  = Ratio(NFL,CFL)
   ! Chlorophyll to carbon ratio (ChlCrPHY, mg Chl/mol C) 
   CALL CHL_C_RATIO(NCrFlagellates,MaxNCrFlagellates,MinNCrFlagellates,MinChlNrFlagellates,MaxChlNrFlagellates,ChlCrFlagellates) ! REMOVE (POSSIBLY)
    chlorophyll(i,j,k,1) = ChlCrFlagellates * CFL
   ! Nitrate uptake rates of phytoplankton (NO3_upPHY, mmol N/m3/day) 
   CALL NO_UPTAKE_RATE(NCrFlagellates,NOS,NHS,CFL,tf,MaxNCrFlagellates,NosMaxUptakeFlagellates,ksNOsFlagellates,kinNHsPhy, Nitrate_UpPHY) ! REMOVE (POSSIBLY)
   ! Ammonium uptake rates of phytoplankton (NH3_upPHY, mmolN /m3/day) 
   ! (and excretion if NC ratio too high) 
   CALL NUT_UPTAKE_RATE(NCrFlagellates,NHS,CFL,tf,MaxNCrFlagellates,NHsMaxUptakeFlagellates,ksNHsFlagellates,0.0_wp, Ammonium_UpPHY) ! REMOVE (POSSIBLY)
   CALL NUT_UPTAKE_RATE(NCrFlagellates,PHO,CFL,tf,MaxNCrFlagellates,PO4MaxUptakeFlagellates,ksPO4Flagellates,0.0_wp, Phosphate_upFlagellates) ! REMOVE (POSSIBLY)
   ! Potential Nitrogen uptake 
    Nitrogen_UpPHY=Ammonium_UpPHY + Nitrate_UpPHY
   ! Nutrient uptake 
    Nutrient_UpPHY=min(Nitrogen_UpPHY,Phosphate_upFlagellates/self%PNRedfield)
             ! Nutrient_UpPHY=max(Nutrient_UpPHY,0)                                                                                                                                                                                                                         ! REMOVE (POSSIBLY)
   ! ################################# FORMER SUBROUTINE PHYT GROWTH RATE ########################################### 
   ! Potential nutrient controlled growth rate; due to droop kinetics, the actual NC ratio can slightly surpass the maximal or minimal ratio 
             IF(NCrFlagellates > MaxNCrFlagellates) THEN ! REMOVE (POSSIBLY)
    NCrat = self%MaxNCrFlagellates
             ELSEIF (NCrFlagellates < MinNCrFlagellates) THEN ! REMOVE (POSSIBLY)
    NCrat = self%MinNCrFlagellates
             ELSE ! REMOVE (POSSIBLY)
    NCrat = NCrFlagellates
             ENDIF ! REMOVE (POSSIBLY)
   ! the actual SiC ratio can slightly surpass the maximal or minimal ratio 
             IF(SiCrFlagellates > MaxSiCrFlagellates) THEN ! REMOVE (POSSIBLY)
    SiCrat = MaxSiCrFlagellates
             ELSEIF (SiCrFlagellates < MinSiCrFlagellates) THEN ! REMOVE (POSSIBLY)
    SiCrat = MinSiCrFlagellates
             ELSE ! REMOVE (POSSIBLY)
    SiCrat = SiCrFlagellates
             ENDIF ! REMOVE (POSSIBLY)
   ! Tett model 
    Mu_Nitrogen = 1. - self%MinNCrFlagellates  / NCrat
    Mu_Silicate = 1. - MinSiCrFlagellates / SiCrat
    NutrientLimitationFlagellates = min(Mu_Nitrogen,Mu_Silicate)
          LightLimitationFlagellates = 1.-exp(-self%alphaPIFlagellates*PAR/self%MuMaxFlagellates)
    Carbon_UptakePHY=self%MuMaxFlagellates*LightLimitationFlagellates*NutrientLimitationFlagellates*CFL*tf
    TotalRespirationPHY=Carbon_UptakePHY*self%GrowthRespFlagellates + self%RespirationFlagellates*CFL*tf
    GrowthPHY=Carbon_UptakePHY-TotalRespirationPHY
   ! Compute the extra DOC excretion as a fraction (extradocphyexcr) of the difference of growth in nutrient limited (actual NC ratio) and nutrient saturated (max NC ratio) as in VDM et al 2004 L&O 
    DOC_extra_excr=abs(extradocphyexcr*self%MuMaxFlagellates*CFL*tf*(min(LightLimitationFlagellates,(1. - self%self%MinNCrFlagellates  / self%MaxNCrFlagellates))- min(LightLimitationFlagellates,(1. - self%self%MinNCrFlagellates  / NCrat))))
   ! ################################# END OF FORMER SUBROUTINE PHYT GROWTH RATE ########################################### 
   !Compute the leakage    and extra DOC excretion, DOC_extra_excr 
    DOC_leakage = self%leakagephy*Carbon_UptakePHY
    DON_leakage  =self%leakagephy*abs(Nutrient_upPHY)
   ! Phytoplankton mortality rate (PHYMort, /day) 
             CALL PHYMORT_RATE(MortalityFlagellates,tf, PHYmort) ! REMOVE (POSSIBLY)
   ! Phytoplankton mortality flux C_PHYmort,N_PHYMort (in mmol C/m3/day or mmol N/m3/day) 
    C_PHYMort  = PHYMort * CFL
    N_PHYMort  = PHYMort * NFL
   ! ADJUSTING THE RATE OF CHANGE 
   ! phytoplankton C increases by growth, 
   ! it decreases by zooplankton grazing (see the zooplankton subroutine) and phytoplankton mortality 
   _ADD_SOURCE_(self%id_cfl,1.0*( GrowthPHY)) 
   _ADD_SOURCE_(self%id_cfl,-1.0*( C_PHYMort + DOC_leakage)) 
   ! phytoplankton N increases by N uptake, 
   ! it decreases by zooplankton grazing and phytoplankton mortality 
   _ADD_SOURCE_(self%id_nfl,1.0*( Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_nfl,-1.0*( N_PHYMort +  DON_leakage)) 
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
   ! As in Anderson and Pondhaven (2003), the phytoplanton mortality increases the pool of POM and DOM with a coefficient of mortdom 
   ! for the DOM pool a part is considered as labile (labilefraction) and another part as semi labile (1-labilefraction) 
   _ADD_SOURCE_(self%id_poc,1.0*( (1.0 - self%mortphydom)*C_PHYMort)) 
   _ADD_SOURCE_(self%id_pon,1.0*( (1.0 - self%mortphydom)*N_PHYMort)) 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%mortphydom*C_PHYMort*self%labilefraction)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( self%mortphydom*N_PHYMort*self%labilefraction)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( self%mortphydom*C_PHYMort*(1.0 - self%labilefraction))) 
   _ADD_SOURCE_(self%id_dns,1.0*( self%mortphydom*N_PHYMort*(1.0 - self%labilefraction))) 
   ! As in Anderson and Pondhaven (2003), the DOCL and DONL concentration increases also due to phytoplankton leakage which is considered 
   ! to produce only labile DOC and DON 
   _ADD_SOURCE_(self%id_dcl,1.0*( DOC_leakage)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( DON_leakage)) 
   ! As in Anderson and Pondhaven (2003), the phytoplankton extra DOC excretion, DOC_extra_excr, increases the pool of DOC (only DOC and not DON) 
   ! this extra-excretion is composed of a part (leakage) which is assimilated to leakage and is thus only labile 
   ! the remaining part (1-leakage) is composed of a labile (labileextradocexcr) and semi labile part (1-labileextradocexcr) 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%self%leakagephy*DOC_extra_excr + self%labileextradocphyexcr*(1.0 - self%self%leakagephy)*DOC_extra_excr)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( (1.0 - self%labileextradocphyexcr)*(1.0 - self%leakagephy)*DOC_extra_excr)) 
             SELECT CASE (SinkingVelocityType) ! REMOVE (POSSIBLY)
             CASE ('aggregation') ! REMOVE (POSSIBLY)
   ! The number of aggregates, POMNOS including  PON increases and decreases with PON 
               !      dPAGGI(i,j,k) = dPAGGI(i,j,k) + (1.0 - mortphydom)*N_PHYMort*AGGI(i,j,k)/(PONI(I,J,K)+PNSI(I,J,K)) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,1.0*( (1.0 - self%mortphydom)*N_PHYMort*AGG/(PON))) 
#ifdef nanquest 
               if (isnan(dPAGG(I,J,K))) then ! REMOVE (POSSIBLY)
                 write (*,*) '** NAN QUEST ** in Calcflagellates' ! REMOVE (POSSIBLY)
                 write (*,*) 'i,j,k,PON(I,J,K),AGG(I,J,K)',i,j,k,PON,AGG ! REMOVE (POSSIBLY)
                 stop ! REMOVE (POSSIBLY)
            endif 
#endif 
             END SELECT ! REMOVE (POSSIBLY)
   ! The oxygen concentration increases due to photosynthesis 
   _ADD_SOURCE_(self%id_dox,1.0*( (GrowthPHY + DOC_extra_excr)*self%OCr)) 
   ! CO2 production and consumption 
   _ADD_SOURCE_(self%id_dic,-1.0*( GrowthPHY + DOC_extra_excr)) 
   !Diagnostics Store the fluxes 
#ifdef primaryprod 
   ! initialized here because FL is called first ... 
          NPP = Carbon_UptakePHY-TotalRespirationPHY
#endif 
#ifdef biodiag2 
          Nitrogen_Uptake_Flagellates = Ammonium_UpPHY + Nitrate_UpPHY
          Carbon_UptakeFlagellates=Carbon_UptakePHY
          TotalRespirationFlagellates=TotalRespirationPHY
#endif 
#ifdef biodiag1 
          PhytoNitrateReduction=Nutrient_UpPHY*ratio(Nitrate_UpPHY,Nitrogen_UpPHY)*self%ONoxnhsr
#endif 
        endif 
         END DO ! REMOVE (POSSIBLY)
       END DO ! REMOVE (POSSIBLY)
     END DO ! REMOVE (POSSIBLY)
   ! OUTPUT VARIABLES 
   ! Diagnostics Averaged over entire water column 
#ifdef biodiag2 
   _SET_DIAGNOSTIC_(self%id_Nitrogen_Uptake_Flagellates, Nitrogen_Uptake_Flagellates)
   _SET_DIAGNOSTIC_(self%id_Carbon_UptakeFlagellates, Carbon_UptakeFlagellates)
   _SET_DIAGNOSTIC_(self%id_TotalRespirationFlagellates, TotalRespirationFlagellates)
#endif 
   _LOOP_END_

   end subroutine do

