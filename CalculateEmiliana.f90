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
! Translation into FABM: L. Ivanov, ULg / MAST
!
! Contains the pelagic submodel, as used in Soetaert et al., 2001.
! References to the microplankton model:
! Tett, P., 1998. Parameterising a microplankton model.
! Department of Biological Sciences, Napier University,
! Report ISBN 0 902703 60 9, 60 pp.
!
! Sharples Tett (1994). Modeling the effect of physical! variability on the midwater chlorophyll maximum.
! Journal of marine research 52: 219-238
!--------------------------------------------------------------------*

   module fabm_ulg_Emiliana 
 
   use fabm_types 
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_Emiliana 
      type (type_state_variable_id)         :: id_cem,id_nem
      type (type_state_variable_id)         :: id_agg,id_cfl,id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nfl,id_nhs,id_nos,id_pho,id_poc,id_pon
      type (type_dependency_id)             :: id_par,id_temp 
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeEmiliana,id_Nitrogen_Uptake_Emiliana,id_NPP,id_PhytoNitrateReduction,id_TotalRespirationEmiliana
      type (type_diagnostic_variable_id)    :: id_Carbon_UptakeEmilianaIntegrated,id_Nitrogen_Uptake_EmilianaIntegrated,id_NPPIntegrated,id_PhytoNitrateReductionIntegrated,id_TotalRespirationEmilianaIntegrated

!     Model parameters 
      real(rk) :: alphaPIEmiliana
      real(rk) :: extradocphyexcr
      real(rk) :: GrowthRespEmiliana
      real(rk) :: kinNHsPhy
      real(rk) :: ksNHsEmiliana
      real(rk) :: ksNOsEmiliana
      real(rk) :: ksPO4Emiliana
      real(rk) :: labileextradocphyexcr
      real(rk) :: labilefraction
      real(rk) :: leakagephy
      real(rk) :: MaxChlNrEmiliana
      real(rk) :: MaxNCrEmiliana
      real(rk) :: MinChlNrEmiliana
      real(rk) :: MinNCrEmiliana
      real(rk) :: MortalityEmiliana
      real(rk) :: mortphydom
      real(rk) :: MuMaxEmiliana
      real(rk) :: NHsMaxUptakeEmiliana
      real(rk) :: NosMaxUptakeEmiliana
      real(rk) :: OCr
      real(rk) :: ONoxnhsr
      real(rk) :: PNRedfield
      real(rk) :: PO4MaxUptakeEmiliana
      real(rk) :: Q10Phy
      real(rk) :: QuantumYieldEmiliana
      real(rk) :: RespirationEmiliana

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

   real(rk)     :: alphaPIEmiliana=0.3/daytosecond
   real(rk)     :: extradocphyexcr=0.05
   real(rk)     :: GrowthRespEmiliana=0.1
   real(rk)     :: kinNHsPhy=0.5
   real(rk)     :: ksNHsEmiliana=0.05
   real(rk)     :: ksNOsEmiliana=0.05
   real(rk)     :: ksPO4Emiliana=0.02
   real(rk)     :: labileextradocphyexcr=0.65
   real(rk)     :: labilefraction=0.7
   real(rk)     :: leakagephy=0.02
   real(rk)     :: MaxChlNrEmiliana=2.0
   real(rk)     :: MaxNCrEmiliana=0.2
   real(rk)     :: MinChlNrEmiliana=1.0
   real(rk)     :: MinNCrEmiliana=0.05
   real(rk)     :: MortalityEmiliana=0.03/daytosecond
   real(rk)     :: mortphydom=0.34
   real(rk)     :: MuMaxEmiliana=2.5/daytosecond
   real(rk)     :: NHsMaxUptakeEmiliana=1.5/daytosecond
   real(rk)     :: NosMaxUptakeEmiliana=1.5/daytosecond
   real(rk)     :: OCr=1.0
   real(rk)     :: ONoxnhsr=2.0
   real(rk)     :: PNRedfield=1.0/16.0
   real(rk)     :: PO4MaxUptakeEmiliana=1.5/16.0/daytosecond
   real(rk)     :: Q10Phy=2.0
   real(rk)     :: QuantumYieldEmiliana=0.6
   real(rk)     :: RespirationEmiliana=0.009/daytosecond

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
   call self%get_parameter(self%alphaPIEmiliana, 'alphaPIEmiliana', default=alphaPIEmiliana) 
   call self%get_parameter(self%extradocphyexcr, 'extradocphyexcr', default=extradocphyexcr) 
   call self%get_parameter(self%GrowthRespEmiliana, 'GrowthRespEmiliana', default=GrowthRespEmiliana) 
   call self%get_parameter(self%kinNHsPhy, 'kinNHsPhy', default=kinNHsPhy) 
   call self%get_parameter(self%ksNHsEmiliana, 'ksNHsEmiliana', default=ksNHsEmiliana) 
   call self%get_parameter(self%ksNOsEmiliana, 'ksNOsEmiliana', default=ksNOsEmiliana) 
   call self%get_parameter(self%ksPO4Emiliana, 'ksPO4Emiliana', default=ksPO4Emiliana) 
   call self%get_parameter(self%labileextradocphyexcr, 'labileextradocphyexcr', default=labileextradocphyexcr) 
   call self%get_parameter(self%labilefraction, 'labilefraction', default=labilefraction) 
   call self%get_parameter(self%leakagephy, 'leakagephy', default=leakagephy) 
   call self%get_parameter(self%MaxChlNrEmiliana, 'MaxChlNrEmiliana', default=MaxChlNrEmiliana) 
   call self%get_parameter(self%MaxNCrEmiliana, 'MaxNCrEmiliana', default=MaxNCrEmiliana) 
   call self%get_parameter(self%MinChlNrEmiliana, 'MinChlNrEmiliana', default=MinChlNrEmiliana) 
   call self%get_parameter(self%MinNCrEmiliana, 'MinNCrEmiliana', default=MinNCrEmiliana) 
   call self%get_parameter(self%MortalityEmiliana, 'MortalityEmiliana', default=MortalityEmiliana) 
   call self%get_parameter(self%mortphydom, 'mortphydom', default=mortphydom) 
   call self%get_parameter(self%MuMaxEmiliana, 'MuMaxEmiliana', default=MuMaxEmiliana) 
   call self%get_parameter(self%NHsMaxUptakeEmiliana, 'NHsMaxUptakeEmiliana', default=NHsMaxUptakeEmiliana) 
   call self%get_parameter(self%NosMaxUptakeEmiliana, 'NosMaxUptakeEmiliana', default=NosMaxUptakeEmiliana) 
   call self%get_parameter(self%OCr, 'OCr', default=OCr) 
   call self%get_parameter(self%ONoxnhsr, 'ONoxnhsr', default=ONoxnhsr) 
   call self%get_parameter(self%PNRedfield, 'PNRedfield', default=PNRedfield) 
   call self%get_parameter(self%PO4MaxUptakeEmiliana, 'PO4MaxUptakeEmiliana', default=PO4MaxUptakeEmiliana) 
   call self%get_parameter(self%Q10Phy, 'Q10Phy', default=Q10Phy) 
   call self%get_parameter(self%QuantumYieldEmiliana, 'QuantumYieldEmiliana', default=QuantumYieldEmiliana) 
   call self%get_parameter(self%RespirationEmiliana, 'RespirationEmiliana', default=RespirationEmiliana) 

   ! Register state variables 

   call self%register_state_variable(self%id_cem, 'CEM'  & 
         , 'mmol C m-3', 'Small flagellate biomass in carbon' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_variable(self%id_nem, 'NEM'  & 
         , 'mmol N m-3', 'Small flagellate biomass in nitrogen' & 
         minimum=0.0e-7_rk, vertical_movement=self%2.0_rk) 
   call self%register_state_dependency(self%id_agg, 'Aggregates', 'm-3') 
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
   class (type_uhh_dinoflag), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  AGG,CFL,DCL,DCS,DIC,DNL,DNS,DOX,NFL,NHS,NOS,PHO,POC,PON
      real(rk) ::  par,temp
      real(rk) ::  CEM,NEM
      real(rk) ::   Carbon_UptakeEmiliana,Nitrogen_Uptake_Emiliana,NPP,PhytoNitrateReduction,TotalRespirationEmiliana
      real(rk) ::   Carbon_UptakeEmilianaIntegrated,Nitrogen_Uptake_EmilianaIntegrated,NPPIntegrated,PhytoNitrateReductionIntegrated,TotalRespirationEmilianaIntegrated
      real(rk) ::   Ammonium_UpPHY	 + ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   C_PHYMort	 + ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Carbon_UptakePHY	 + ! mmol C m-3, C assimilation of phytoplankton
      real(rk) ::   ChlCrEmiliana	 + ! g Chla mol C-1, Chl/C ratio in small flagellates
      real(rk) ::   DOC_extra_excr	 + ! mmol C d-1, Phytoplankton extra excretion
      real(rk) ::   DOC_leakage	 + ! mmol C d-1, Phytoplankton passive leakage rate for carbon
      real(rk) ::   DON_leakage	 + ! mmol N d-1, Phytoplankton passive leakage rate for nitrogen
      real(rk) ::   GrowthPHY	 + ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   LightLimitationEmiliana	 + ! ?, Light limitation for small flagellates
      real(rk) ::   N_PHYMort	 + ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   NCrEmiliana	 + ! mol N mol C-1, N/C ratio in small flagellates
      real(rk) ::   Nitrate_UpPHY	 + ! mmol N m-3, Nitrate uptake of phytoplankton
      real(rk) ::   Nitrogen_UpPHY	 + ! mmol N m-3, Nitrogen uptake of phytoplankton
      real(rk) ::   Nutrient_UpPHY	 + ! mmol m-3, Nutrient uptake of phytoplankton
      real(rk) ::   NutrientLimitationEmiliana	 + ! ?, Nutrient limitation for small flagellates
      real(rk) ::   Phosphate_upEmiliana	 + ! mmol P m-3 , Small flagellates phosphate uptake
      real(rk) ::   PHYMort	 + ! mmol m-3, Phytoplankton mortality rate
      real(rk) ::   tf	 + ! -, Temperature factor
      real(rk) ::   TotalRespirationPHY	 + ! mmol C m-3, Total phytoplankton respiration (basal & activity)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_agg,AGG)       ! Aggregates
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
   ! TEMPERATURE EFFECT on rates (all rates are defined at 20 dg C) 
    tf = Q10Factor(temp,Q10Phy)
   ! PHYTOPLANKTON 
   ! N/C ratio 
    NCrEmiliana  = Ratio (NEM,CEM)
   ! Chlorophyll to carbon ratio (ChlCrPHY, mg Chl/mol C) 
   CALL CHL_C_RATIO(NCrEmiliana,MaxNCrEmiliana,MinNCrEmiliana,MinChlNrEmiliana,MaxChlNrEmiliana,ChlCrEmiliana) ! REMOVE (POSSIBLY)
    Chlorophyll(i,j,k,3) = ChlCrEmiliana * CEM
   ! Nitrate uptake rates of phytoplankton (NO3_upPHY, mmol N/m3/day) 
   CALL NO_UPTAKE_RATE(NCrEmiliana,NOS,NHS,CEM,tf,MaxNCrEmiliana,NosMaxUptakeEmiliana,ksNOsEmiliana,kinNHsPhy, Nitrate_UpPHY) ! REMOVE (POSSIBLY)
   ! Ammonium uptake rates of phytoplankton (NH3_upPHY, mmolN /m3/day) 
   ! (and excretion if NC ratio too high) 
   CALL NUT_UPTAKE_RATE(NCrEmiliana,NHS,CEM,tf,MaxNCrEmiliana,NHsMaxUptakeEmiliana,ksNHsEmiliana,0.0_wp, Ammonium_UpPHY) ! REMOVE (POSSIBLY)
   CALL NUT_UPTAKE_RATE(NCrEmiliana,PHO,CEM,tf,MaxNCrEmiliana,PO4MaxUptakeEmiliana,ksPO4Emiliana,0.0_wp, Phosphate_upEmiliana) ! REMOVE (POSSIBLY)
   ! Potential Nitrogen uptake 
    Nitrogen_UpPHY=Ammonium_UpPHY + Nitrate_UpPHY
   ! Nutrient uptake 
    Nutrient_UpPHY=min(Nitrogen_UpPHY,Phosphate_upEmiliana/self%PNRedfield)
             !Nutrient_UpPHY=max(Nutrient_UpPHY,0) ! REMOVE (POSSIBLY)
   CALL GROWTH_RATE(ChlCrEmiliana,QuantumYieldEmiliana,alphaPIEmiliana,par(i,j,k),LightLimitationEmiliana,NutrientLimitationEmiliana,MaxNCrEmiliana,MinNCrEmiliana,NCrEmiliana,1.0_wp,0.0_wp,1.0_wp,MuMaxEmiliana,RespirationEmiliana,GrowthRespEmiliana,extradocphyexcr,tf,CEM,GrowthPHY,Carbon_UptakePHY,TotalRespirationPHY,DOC_extra_excr) ! REMOVE (POSSIBLY)
   !Compute the leakage    and extra DOC excretion, DOC_extra_excr 
    DOC_leakage = self%leakagephy*Carbon_UptakePHY
    DON_leakage  =self%leakagephy*abs(Nutrient_UpPHY)
   ! Phytoplankton mortality rate (PHYMort, /day) 
             CALL PHYMORT_RATE(MortalityEmiliana,tf, PHYmort) ! REMOVE (POSSIBLY)
   ! Phytoplankton mortality flux C_PHYMort,N_PHYMort (in mmol C/m3/day or mmol N/m3/day) 
    C_PHYMort  = PHYMort * CEM
    N_PHYMort  = PHYMort * NEM
   ! ADJUSTING THE RATE OF CHANGE 
   ! phytoplankton C increases by growth, 
   ! it decreases by zooplankton grazing (see the zooplankton subroutine) and phytoplankton mortality 
   _ADD_SOURCE_(self%id_cem,1.0*( GrowthPHY)) 
   _ADD_SOURCE_(self%id_cem,-1.0*( C_PHYMort+DOC_leakage)) 
   ! phytoplankton N increases by N uptake, 
   ! it decreases by zooplankton grazing and phytoplankton mortality 
   _ADD_SOURCE_(self%id_nem,1.0*( Nutrient_UpPHY)) 
   _ADD_SOURCE_(self%id_nem,-1.0*( N_PHYMort+DON_leakage)) 
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
             END IF ! REMOVE (POSSIBLY)
   ! As in Anderson and Pondhaven (2003), the phytoplanton mortality increases the pool of POM and DOM with a coefficient of mortdom 
   ! for the DOM pool a part is considered as labile (labilefraction) and another part as semi labile (1-labilefraction) 
   _ADD_SOURCE_(self%id_poc,1.0*( (1.0 - self%mortphydom)*C_PHYMort)) 
   _ADD_SOURCE_(self%id_pon,1.0*( (1.0 - self%mortphydom)*N_PHYMort)) 
   _ADD_SOURCE_(self%id_dcl,1.0*( self%mortphydom*C_PHYMort*self%labilefraction)) 
   _ADD_SOURCE_(self%id_dnl,1.0*( self%mortphydom*N_PHYMort*self%labilefraction)) 
   _ADD_SOURCE_(self%id_dcs,1.0*( self%mortphydom*C_PHYMort*(1.0 - self%labilefraction))) 
   _ADD_SOURCE_(self%id_dns,1.0*( self%mortphydom*N_PHYMort*(1.0 - self%labilefraction))) 
   ! As in Anderson and Pondhaven (2003), the DOCL and DONL concentration increases also due to phytoplankton leakage which is considered 
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
   ! The number of aggregates, POMNOS including  PON increases and decreases with PON 
               !      dPAGGI(i,j,k) = dPAGGI(i,j,k) + (1.0 - mortphydom)*N_PHYMort*AGGI(i,j,k)/(PONI(i,j,k)+PNSI(i,j,k)) ! REMOVE (POSSIBLY)
   _ADD_SOURCE_(self%id_agg,1.0*( (1.0 - self%mortphydom)*N_PHYMort*AGG/(PON))) 
#ifdef nanquest 
               if (isnan(dPAGG(I,J,K))) then ! REMOVE (POSSIBLY)
                 write (*,*) '** NAN QUEST ** in CalcEMILIANNA' ! REMOVE (POSSIBLY)
                 write (*,*) 'i,j,k,PON(I,J,K),trb(I,J,K,AGG)' ! REMOVE (POSSIBLY)
                 write (*,*) i,j,k,trb(I,J,K,PON),trb(I,J,K,AGG) ! REMOVE (POSSIBLY)
                 call flush(6) ! REMOVE (POSSIBLY)
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
           end if ! REMOVE (POSSIBLY)
         end do  ! REMOVE (POSSIBLY)
       end do ! REMOVE (POSSIBLY)
     end do ! REMOVE (POSSIBLY)
   ! OUTPUT VARIABLES 
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

