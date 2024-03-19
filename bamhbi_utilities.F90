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

  pure real(rk) function yy(a,x)
  IMPLICIT NONE
  real(rk),intent(in)        :: a,x
  yy=x**2/(a**2+x**2)
  RETURN
  END function yy

   pure real(rk) function Q10Factor(Temperature, Q10)
   implicit none
   ! Temperature dependence according to a Q10 formulation
      real(rk),intent(in) :: Temperature 	 		! Temperature (C)
      real(rk),intent(in) :: Q10	 	 		! Temperature factor
      IF (Q10 <= 0.) THEN
        Q10Factor = 0.
      ELSE
        Q10Factor = EXP(LOG(Q10)*(Temperature-20.0)/10.0)
      ENDIF
   return
   end function Q10Factor


   pure real(rk) function Michaelis(limitingelement,halfsaturationconstant)
   implicit none
   ! Limiting Michaelis-Menten function
      real(rk),intent(in) :: limitingelement 	 	! Limiting element (quantity/concentration)
      real(rk),intent(in) :: halfsaturationconstant	! Half-saturation constant for this element
      Michaelis = limitingelement/(limitingelement + halfsaturationconstant)
   return
   end function Michaelis


   pure real(rk) function Ratio(t1,t2)
   implicit none
   ! Ratio
      real(rk),intent(in) :: t1 	 			! First tracer
      real(rk),intent(in) :: t2	 			! Second tracer
      Ratio = t1/t2
   return
   end function Ratio


   pure real(rk) function Inhibition(inhibitingelement,inhibitionconstant)
   implicit none
   ! Inibiting function computed as 1- michaelis menten function
      real(rk),intent(in) :: inhibitingelement 	 	! Inhibiting element (quantity/concentration)
      real(rk),intent(in) :: inhibitionconstant	 	! Half-saturation constant for inhibition
      Inhibition = inhibitionconstant/(inhibitingelement + inhibitionconstant)
   return
   end function Inhibition


   pure real(rk) function zoo_egestion(Asszoo,grazing)
   ! Calculate egestion by zooplankton or gelatinous
      real(rk),intent(in) :: Asszoo 	 	! Assimilation efficiency
      real(rk),intent(in) :: grazing	 	! mmol m-3, Intake of C or N
      zoo_egestion = (1-Asszoo)*grazing
   end function zoo_egestion


   pure real(rk) function mortality_consument(HalfSatMort_Zoo,mortnonlinzoo,expmortzoo,DOXsatmortality,mortalityanoxic,tf,ZOO,DOX)
   implicit none
   ! Mortality rate : sum of a linear mortality, a quadratic mortality and an oxygen dependent term
      real(rk),intent(in) :: HalfSatMort_Zoo 	! mmol C m-3, Half-saturation constant of mortalirty
      real(rk),intent(in) :: mortnonlinzoo 	! mmolC m-3 d-1, Non-linear mortality rate
      real(rk),intent(in) :: expmortzoo 	 	! -, Order of the non-linearity of mortality rate
      real(rk),intent(in) :: DOXsatmortality 	! mmol O m-3, Percentage of saturation where metabolic respiration is half the one under oxygen satyrated conditions
      real(rk),intent(in) :: mortalityanoxic 	! d-1, Anoxic mortality of zooplankton 
      real(rk),intent(in) :: tf 	 		! -, Temperature function
      real(rk),intent(in) :: ZOO 	 		! mmol C m-3, Zooplankton biomass in carbon
      real(rk),intent(in) :: DOX 	 		! mmol O m-3, Oxygen concentration
      real(rk) :: foxygen 			! Oxygen limitation function
      foxygen = Michaelis(DOX,DOXsatmortality)
      mortality_consument = tf*(mortnonlinzoo*ZOO**expmortzoo/(ZOO+HalfSatMort_Zoo)+(1-foxygen)*mortalityanoxic*ZOO)
   return
   end function mortality_consument


   pure real(rk) function respiration_jelly(tf,Asszoo,Effzoo,bzoo,ZOO,grazing_carbonZoo)
   implicit none
  ! Calculates the respiration rate of zooplankton as a sum of a basal respiration and an activity respiration
      real(rk),intent(in) :: tf 			! -, Temperature function
      real(rk),intent(in) :: Asszoo 	 	! -, Assimilation efficiency
      real(rk),intent(in) :: Effzoo 	 	! mmol m-3, Inhibition constant
      real(rk),intent(in) :: bzoo 	 		! d-1, Basal respiration rate
      real(rk),intent(in) :: ZOO 	 		! mmol C m-3, Zooplankton biomass in carbon
      real(rk),intent(in) :: grazing_carbonZoo 	! mmol C m-3, Intake of carbon
      real(rk) :: basalrespiration 		! mmol C m-3, Basal respiration
      real(rk) :: activityrespiration 	 	! mmol C m-3, Activity respiration
  ! Basal respiration
      basalrespiration = bzoo*tf*Zoo
  ! Activity respiration
      activityrespiration = (1.0-Effzoo)*Asszoo*grazing_carbonZoo
  ! Total respiration
      respiration_jelly = basalrespiration + activityrespiration
   return
   end function respiration_jelly


   pure real(rk) function ratio_chl_c_phyt(NCratio,MaxNCr,MinNCr,MinChlNr,MaxChlNr)
   implicit none
  ! Calculates the Chlorophyll to carbon ratio of the phytoplankton (mg Chl/mol C), based on the N / C content
      real(rk),intent(in) :: NCratio 		! mol N mol C-1, Actual N/C ratio
      real(rk),intent(in) :: MaxNCr 		! mol N mol C-1, Maximum N/C ratio
      real(rk),intent(in) :: MinNCr 		! mol N mol C-1, Minimum N/C ratio
      real(rk),intent(in) :: MinChlNr 		! g Chla mol N-1, Minimum Chl/C ratio
      real(rk),intent(in) :: MaxChlNr 		! g Chla mol N-1, Minimum Chl/C ratio
      real(rk) :: RangeChl_N 	 	! Range between maximum and minimum N/C ratios
   IF (NCratio >= MaxNCr) THEN
      ratio_chl_c_phyt = MaxNCr  * MaxChlNr
   ELSE IF (NCratio <= MinNCr) THEN
      ratio_chl_c_phyt= MinNCr  * MinChlNr
   ELSE
      RangeChl_N = MaxChlNr - MinChlNr
      ratio_chl_c_phyt = NCratio* (MinChlNr +rangeChl_n * (NCratio- MinNCr)/(MaxNCr-MinNCr))
   ENDIF
   return
   end function ratio_chl_c_phyt


   pure real(rk) function uptake_nitrate_phyt(NCratio,tf,MaxNCr,umax)
   implicit none
  ! Calculates the uptake rate of nitrate by the phytoplankton
      real(rk),intent(in) :: NCratio 	 	! mol N mol C-1, Actual N/C ratio
      real(rk),intent(in) :: tf 			! -, Temperature function
      real(rk),intent(in) :: MaxNCr 	 	! mol N mol C-1, Maximum N/C ratio
      real(rk),intent(in) :: umax 	 		! mol N d-1 (?), Maximal NO3 uptake rate
  IF (NCratio >= MaxNCr) THEN
      uptake_nitrate_phyt=0.0
  ELSE
      uptake_nitrate_phyt= umax*tf* (1.0-(NCratio/MaxNCr))
  ENDIF
   return
   end function uptake_nitrate_phyt


   pure real(rk) function uptake_nutrient_phyt(NutCratio,Nut0,tf,MaxNutCr,umax,ks)
   implicit none
  ! Calculates the uptake rate of a nutrient by phytoplankton (if negative then excretion)
      real(rk),intent(in) :: NutCratio 		! -, Actual Nut/C ratio
      real(rk),intent(in) :: tf 			! -, Temperature function
      real(rk),intent(in) :: MaxNutCr 	 	! -, Maximal Nut/C ratio
      real(rk),intent(in) :: umax 	 		! mol N d-1 (?), Maximal nutrient uptake rate by a phytoplankton
      real(rk),intent(in) :: ks 	 		! mmol m-3, Half-saturation constant for nutrient uptake
      real(rk),intent(in) :: Nut0 	 		! mmol m-3, Nutrient concentration minus inhibition constant
  IF (NutCratio > MaxNutCr) THEN
      uptake_nutrient_phyt =  umax*tf * (1.0-(NutCratio/ MaxNutCr))
  ELSE IF (Nut0 > 0.0) THEN
      uptake_nutrient_phyt =  umax*tf * (1.0-(NutCratio/MaxNutCr)) * Nut0/(ks+Nut0)
  ELSE
      uptake_nutrient_phyt = 0.0
  ENDIF
   return
   end function uptake_nutrient_phyt


   pure real(rk) function ratio_adjustment(ActualRatio,MaxRatio,MinRatio)
   implicit none
  ! Correction of the ratio of elements in phytoplankton
      real(rk),intent(in) :: ActualRatio 		! -, Actual E1:E2 ratio
      real(rk),intent(in) :: MaxRatio		! -, Maximum E1:E2 ratio
      real(rk),intent(in) :: MinRatio		! -, Minimum E1:E2 ratio
  IF (ActualRatio > MaxRatio) THEN
        ratio_adjustment = MaxRatio
  ELSE IF (ActualRatio < MinRatio) THEN
        ratio_adjustment = MinRatio
  ELSE
        ratio_adjustment = ActualRatio
  ENDIF
   return
   end function ratio_adjustment


   pure real(rk) function limitation_by_nutrient(MinNCrPHYT,MinSiCrPHYT,NCratio,SiCratio)
  ! Compute limitation by nutrients
      real(rk),intent(in) :: MinNCrPHYT 		! -, Minimal N:C ratio
      real(rk),intent(in) :: MinSiCrPHYT		! -, Minimal Si:C ratio
      real(rk),intent(in) :: NCratio		! -, Actual N:C ratio
      real(rk),intent(in) :: SiCratio		! -, Actual Si:C ratio
      real(rk) :: Mu_Nitrogen 			! -, Limitation by nitrogen
      real(rk) :: Mu_Silicate 	 		! -, Limitation by silicate
      Mu_Nitrogen = 1.0 - MinNCrPHYT  / NCratio
      Mu_Silicate = 1.0 - MinSiCrPHYT / SiCratio
      limitation_by_nutrient = min(Mu_Nitrogen,Mu_Silicate)
   return
   end function limitation_by_nutrient



   pure real(rk) function extra_excretion(LightLim,MinNCrPHYT,MaxNCrPHYT,NCratio)
   implicit none
  ! Compute a term for DOC extra-excretion
      real(rk),intent(in) :: LightLim 		! -, Minimal N:C ratio
      real(rk),intent(in) :: MinNCrPHYT		! -, Minimal Si:C ratio
      real(rk),intent(in) :: MaxNCrPHYT		! -, Minimal Si:C ratio
      real(rk),intent(in) :: NCratio		! -, Actual N:C ratio
      real(rk) :: Lim1				! -, Limitation N1
      real(rk) :: Lim2				! -, Limitation N2
      Lim1 = min(LightLim,(1.0 - MinNCrPHYT  / MaxNCrPHYT))
      Lim2 = min(LightLim,(1.0 - MinNCrPHYT  / NCratio))
      extra_excretion = abs(Lim1 - Lim2)
   return
   end function extra_excretion


   pure real(rk) function mortality_phyt(Mortality,tf)
   implicit none
      real(rk),intent(in)  :: Mortality 	 	! d-1, First-order mortality rate (P)
      real(rk),intent(in)  :: tf 	 			! -, Temperature function
      mortality_phyt = Mortality * tf
   return
   end function mortality_phyt

   end module fabm_ulg_bamhbi_split_utilities
