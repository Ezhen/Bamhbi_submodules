#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module fabm_ulg_bamhbi_split_utilities
!
! !DESCRIPTION:
! This module holds common functions and data types for the bamhbi model
!
! !USES:
   use fabm_types

   implicit none
   
   contains

   elemental real(rk) function Q10Factor(Temperature, Q10)
   ! Temperature dependence according to a Q10 formulation
      real(rk), intent :: Temperature 	 		! Temperature (C)
      real(rk), intent :: Q10	 	 		! Temperature factor
      IF (Q10 <= 0.) THEN
        Q10Factor = 0.
      ELSE
        Q10Factor = EXP(LOG(Q10)*(Temperature-20.0)/10.0)
      ENDIF
   end function Q10Factor


   elemental real(rk) function Michaelis(limitingelement,halfsaturationconstant)
   ! Limiting Michaelis-Menten function
      real(rk), intent :: limitingelement 	 	! Limiting element (quantity/concentration)
      real(rk), intent :: halfsaturationconstant	! Half-saturation constant for this element
      Michaelis = limitingelement/(limitingelement + halfsaturationconstant)
   end function Michaelis


   elemental real(rk) function Ratio(t1,t2)
   ! Ratio
      real(rk), intent :: t1 	 			! First tracer
      real(rk), intent :: t2	 			! Second tracer
      Ratio = t1/t2
   end function Ratio


   elemental real(rk) function Inhibition(inhibitingelement,inhibitionconstant)
   ! Inibiting function computed as 1- michaelis menten function
      real(rk), intent :: inhibitingelement 	 	! Inhibiting element (quantity/concentration)
      real(rk), intent :: inhibitionconstant	 	! Half-saturation constant for inhibition
      Inhibition = inhibitionconstant/(inhibitingelement + inhibitionconstant)
   end function Inhibition


   elemental real(rk) function zoo_egestion(Asszoo,grazing)
   ! Calculate egestion by zooplankton or gelatinous
      real(rk), intent :: Asszoo 	 	! Assimilation efficiency
      real(rk), intent :: grazing	 	! mmol m-3, Intake of C or N
      zoo_egestion = (1-Asszoo)*grazing
   end function zoo_egestion


   elemental real(rk) function mortality_consument(HalfSatMort_Zoo,mortnonlinzoo,expmortzoo,DOXsatmortality,mortalityanoxic,tf,ZOO,DOX)
   ! Mortality rate : sum of a linear mortality, a quadratic mortality and an oxygen dependent term
      real(rk), intent :: HalfSatMort_Zoo 	! mmol C m-3, Half-saturation constant of mortalirty
      real(rk), intent :: mortnonlinzoo 	! mmolC m-3 d-1, Non-linear mortality rate
      real(rk), intent :: expmortzoo 	 	! -, Order of the non-linearity of mortality rate
      real(rk), intent :: DOXsatmortality 	! mmol O m-3, Percentage of saturation where metabolic respiration is half the one under oxygen satyrated conditions
      real(rk), intent :: mortalityanoxic 	! d-1, Anoxic mortality of zooplankton 
      real(rk), intent :: tf 	 		! -, Temperature function
      real(rk), intent :: ZOO 	 		! mmol C m-3, Zooplankton biomass in carbon
      real(rk), intent :: DOX 	 		! mmol O m-3, Oxygen concentration
      real(rk) :: foxygen 			! Oxygen limitation function
      foxygen = Michaelis(DOX,DOXsatmortality)
      mortality_consument = tf*(mortnonlinzoo*ZOO**expmortzoo/(ZOO+HalfSatMort_Zoo)+(1-foxygen)*mortalityanoxic*ZOO)
   end function mortality_consument


   elemental real(rk) function respiration_jelly(tf,Asszoo,Effzoo,bzoo,ZOO,grazing_carbonZoo)
  ! Calculates the respiration rate of zooplankton as a sum of a basal respiration and an activity respiration
      real(rk), intent :: tf 			! -, Temperature function
      real(rk), intent :: Asszoo 	 	! -, Assimilation efficiency
      real(rk), intent :: Effzoo 	 	! mmol m-3, Inhibition constant
      real(rk), intent :: bzoo 	 		! d-1, Basal respiration rate
      real(rk), intent :: ZOO 	 		! mmol C m-3, Zooplankton biomass in carbon
      real(rk), intent :: grazing_carbonZoo 	! mmol C m-3, Intake of carbon
      real(rk) :: basalrespiration 		! mmol C m-3, Basal respiration
      real(rk) :: activityrespiration 	 	! mmol C m-3, Activity respiration
  ! Basal respiration
      basalrespiration = bzoo*tf*Zoo
  ! Activity respiration
      activityrespiration = (1.0-Effzoo)*Asszoo*grazing_carbonZoo
  ! Total respiration
      respiration_jelly = basalrespiration + activityrespiration
   end function respiration_jelly


   elemental real(rk) function ratio_chl_c_phyt(NCratio,MaxNCr,MinNCr,MinChlNr,MaxChlNr)
  ! Calculates the Chlorophyll to carbon ratio of the phytoplankton (mg Chl/mol C), based on the N / C content
      real(rk), intent :: NCratio 		! mol N mol C-1, Actual N/C ratio
      real(rk), intent :: MaxNCr 		! mol N mol C-1, Maximum N/C ratio
      real(rk), intent :: MinNCr 		! mol N mol C-1, Minimum N/C ratio
      real(rk), intent :: MinChlNr 		! g Chla mol N-1, Minimum Chl/C ratio
      real(rk), intent :: MaxChlNr 		! g Chla mol N-1, Minimum Chl/C ratio
      real(rk), intent :: RangeChl_N 	 	! Range between maximum and minimum N/C ratios
   IF (NCratio >= MaxNCr) THEN
      ratio_chl_c_phyt = MaxNCr  * MaxChlNr
   ELSEIF (NCratio <= MinNCr) THEN
      ratio_chl_c_phyt= MinNCr  * MinChlNr
   ELSE
      RangeChl_N = MaxChlNr - MinChlNr
      ratio_chl_c_phyt = NCratio* (MinChlNr +rangeChl_n * (NCratio- MinNCr)/(MaxNCr-MinNCr))
   ENDIF
   end function ratio_chl_c_phyt


   elemental real(rk) function uptake_nitrate_phyt(NCratio,nitrate,ammonium,tf,MaxNCr,umax,ks,half_inhib_amm)
  ! Calculates the uptake rate of nitrate by the phytoplankton
      real(rk), intent :: NCratio 	 	! mol N mol C-1, Actual N/C ratio
      real(rk), intent :: tf 			! -, Temperature function
      real(rk), intent :: MaxNCr 	 	! mol N mol C-1, Maximum N/C ratio
      real(rk), intent :: umax 	 		! mol N d-1 (?), Maximal NO3 uptake rate
  IF (NCratio >= MaxNCr) THEN
      uptake_nitrate_phyt=0.0
  ELSE
      uptake_nitrate_phyt= umax*tf* (1.0-(NCratio/MaxNCr))
  ENDIF
   end function uptake_nitrate_phyt


   elemental real(rk) function uptake_nutrient_phyt(NutCratio,Nut0,tf,MaxNutCr,umax,ks)
  ! Calculates the uptake rate of a nutrient by phytoplankton (if negative then excretion)
      real(rk), intent :: NutCratio 		! -, Actual Nut/C ratio
      real(rk), intent :: tf 			! -, Temperature function
      real(rk), intent :: MaxNutCr 	 	! -, Maximal Nut/C ratio
      real(rk), intent :: umax 	 		! mol N d-1 (?), Maximal nutrient uptake rate by a phytoplankton
      real(rk), intent :: ks 	 		! mmol m-3, Half-saturation constant for nutrient uptake
      real(rk), intent :: Nut0 	 		! mmol m-3, Nutrient concentration minus inhibition constant
  IF(NutCratio > MaxNutCr) THEN
      uptake_nutrient_phyt =  umax*tf * (1.0-(NutCratio/ MaxNutCr))
  ELSEIF Nut0 > 0.0 THEN
      uptake_nutrient_phyt =  umax*tf * (1.0-(NutCratio/MaxNutCr)) * Nut0/(ks+Nut0)
  ELSE
      uptake_nutrient_phyt = 0.0
  ENDIF
   end function uptake_nutrient_phyt


   elemental real(rk) function ratio_adjustment(ActualRatio,MaxRatio,MinRatio)
  ! Correction of the ratio of elements in phytoplankton
      real(rk), intent :: ActualRatio 		! -, Actual E1:E2 ratio
      real(rk), intent :: MaxRatio		! -, Maximum E1:E2 ratio
      real(rk), intent :: MinRatio		! -, Minimum E1:E2 ratio
  IF (ActualRatio > MaxRatio) THEN
        ratio_adjustment = MaxRatio
  ELSE IF (ActualRatio < MinRatio) THEN
        ratio_adjustment = MinRatio
  ELSE
        ratio_adjustment = ActualRatio
  ENDIF
   end function ratio_adjustment


   elemental real(rk) function limitation_nutrient(MinNCrPHYT,MinSiCrPHYT,NCratio,SiCratio)
  ! Compute limitation by nutrients
      real(rk), intent :: MinNCrPHYT 		! -, Minimal N:C ratio
      real(rk), intent :: MinSiCrPHYT		! -, Minimal Si:C ratio
      real(rk), intent :: NCratio		! -, Actual N:C ratio
      real(rk), intent :: SiCratio		! -, Actual Si:C ratio
      real(rk) :: Mu_Nitrogen 			! -, Limitation by nitrogen
      real(rk) :: Mu_Silicate 	 		! -, Limitation by silicate
      Mu_Nitrogen = 1.0 - MinNCrPHYT  / NCratio
      Mu_Silicate = 1.0 - MinSiCrPHYT / SiCratio
      limitation_nutrient = min(Mu_Nitrogen,Mu_Silicate)
   end function limitation_nutrient



   elemental real(rk) function extra_excretion(LightLim,MinNCrPHYT,MaxNCrPHYT,NCratio)
  ! Compute a term for DOC extra-excretion
      real(rk), intent :: LightLim 		! -, Minimal N:C ratio
      real(rk), intent :: MinNCrPHYT		! -, Minimal Si:C ratio
      real(rk), intent :: MaxNCrPHYT		! -, Minimal Si:C ratio
      real(rk), intent :: NCratio		! -, Actual N:C ratio
      real(rk) :: Lim1				! -, Limitation N1
      real(rk) :: Lim2				! -, Limitation N2
      Lim1 = min(LightLim,(1.0 - MinNCrPHYT  / MaxNCrPHYT))
      Lim2 = min(LightLim,(1.0 - MinNCrPHYT  / NCratio))
      extra_excretion = abs(Lim1 - Lim2)
   end function extra_excretion


   elemental real(rk) function mortality_phyt(Mortality,tf)
      real(rk), intent :: Mortality 	 	! d-1, First-order mortality rate (P)
      real(rk)  :: tf 	 			! -, Temperature function
      mortality_phyt = Mortality * tf
   end function mortality_phyt

   end module fabm_ulg_bamhbi_split_utilities
