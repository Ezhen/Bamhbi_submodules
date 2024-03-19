#include "fabm_driver.h" 
 
!#########################################################################################
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

   module fabm_ulg_bacteria 
 
   use fabm_types 
   use fabm_ulg_bamhbi_split_utilities
 
   implicit none 
 
!  default: all is private. 
   private 
 
! PUBLIC DERIVED TYPES: 
   type,extends(type_base_model),public :: type_ulg_bacteria 
      type (type_state_variable_id)         :: id_bac
      type (type_state_variable_id)         :: id_dcl,id_dcs,id_dic,id_dnl,id_dns,id_dox,id_nhs,id_nos,id_odu,id_pho
      type (type_dependency_id)             :: id_temp,id_depth 
      type (type_diagnostic_variable_id)    :: id_remineralization_anoxic_bac,id_oxygen_consumption_bac,id_respiration_bac,id_denitrification,id_uptake_dcl_bac

!     Model parameters 
      real(rk)     :: eff_gr_bac_c, f_dl_dom, f_solid_odu, i1_curve
      real(rk)     :: i2_curve, iron, ki_anox_nos, ki_anox_o2
      real(rk)     :: ki_denit_o2, ks_denitr_nos, ks_dls_bac, ks_nhs_bac
      real(rk)     :: ks_odu_iron, ks_oxic_o2, ks_po4_bac, mo_bac
      real(rk)     :: mumax_bac, q10_bac, r_n_c_bac, r_n_c_denit
      real(rk)     :: r_o2_c_resp, r_odu_c_anox, r_p_n_redfield

   contains 

      procedure :: initialize 
      procedure :: do 

   end type

   contains

   ! Initialise the Bacteria model

   subroutine initialize(self,configunit)

   class (type_ulg_bacteria), intent(inout), target :: self
   integer,                   intent(in)            :: configunit

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk

!     Model parameters 
      real(rk)     :: eff_gr_bac_c, f_dl_dom, f_solid_odu, i1_curve
      real(rk)     :: i2_curve, iron, ki_anox_nos, ki_anox_o2
      real(rk)     :: ki_denit_o2, ks_denitr_nos, ks_dls_bac, ks_nhs_bac
      real(rk)     :: ks_odu_iron, ks_oxic_o2, ks_po4_bac, mo_bac
      real(rk)     :: mumax_bac, q10_bac, r_n_c_bac, r_n_c_denit
      real(rk)     :: r_o2_c_resp, r_odu_c_anox, r_p_n_redfield

   namelist /ulg_bacteria/ eff_gr_bac_c, 	 & 
                      f_dl_dom, f_solid_odu, i1_curve, 	 & 
                      i2_curve, iron, ki_anox_nos, ki_anox_o2, 	 & 
                      ki_denit_o2, ks_denitr_nos, ks_dls_bac, 	 & 
                      ks_nhs_bac, ks_odu_iron, ks_oxic_o2, 	 & 
                      ks_po4_bac, mo_bac, mumax_bac, q10_bac, 	 & 
                      r_n_c_bac, r_n_c_denit, r_o2_c_resp, 	 & 
                      r_odu_c_anox, r_p_n_redfield

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=ulg_bacteria,err=99)

   ! Store parameter values in our own derived type 
   ! NB: all rates must be provided in values per day, 
   ! and are converted here to values per second. 
   call self%get_parameter(self%eff_gr_bac_c, 'eff_gr_bac_c', '-', 'BAC gross growth efficiency on C', default=0.17_rk) 
   call self%get_parameter(self%f_dl_dom, 'f_dl_dom', '-', 'Labile fraction of PHY- and nonPHY-produced DOM', default=0.7_rk) 
   call self%get_parameter(self%f_solid_odu, 'f_solid_odu', '-', 'Percentage of solid ODU formation', default=0.2_rk) 
   call self%get_parameter(self%i1_curve, 'i1_curve', '-', 'Parameter of the curve simulating the iron concentration', default=25000.0_rk) 
   call self%get_parameter(self%i2_curve, 'i2_curve', '-', 'Parameter of the curve simulating the iron c', default=50.0_rk) 
   call self%get_parameter(self%iron, 'iron', 'mmolFe m-3', 'Concentration of iron in surface water', default=10.0_rk) 
   call self%get_parameter(self%ki_anox_nos, 'ki_anox_nos', 'mmolN m-3', 'Half-sat. constant for NOS inhibition in anoxic remineralization', default=0.0005_rk) 
   call self%get_parameter(self%ki_anox_o2, 'ki_anox_o2', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in anoxic remineralization', default=0.0005_rk) 
   call self%get_parameter(self%ki_denit_o2, 'ki_denit_o2', 'mmolO2 m-3', 'Half-sat. constant for O2 inhibition in denitrification', default=0.5_rk) 
   call self%get_parameter(self%ks_denitr_nos, 'ks_denitr_nos', 'mmolN m-3', 'Half-sat. constant for NOS lim. in denitrif.', default=0.3_rk) 
   call self%get_parameter(self%ks_dls_bac, 'ks_dls_bac', 'mmolC m-3', 'Half-sat. constant for DLC uptake by BAC', default=25.0_rk) 
   call self%get_parameter(self%ks_nhs_bac, 'ks_nhs_bac', 'mmolN m-3', 'Half-sat. constant for NHS uptake by BAC', default=0.5_rk) 
   call self%get_parameter(self%ks_odu_iron, 'ks_odu_iron', 'mmolFe m-3', 'Half-sat. constant for iron lim. in solid ODU formation', default=100.0_rk) 
   call self%get_parameter(self%ks_oxic_o2, 'ks_oxic_o2', 'mmolO2 m-3', 'Half-sat. constant for O2 lim. in oxic min.', default=3.0_rk) 
   call self%get_parameter(self%ks_po4_bac, 'ks_po4_bac', 'mmolP m-3', 'Half-saturation constant for PO4 uptake by BAC', default=0.031_rk) 
   call self%get_parameter(self%mo_bac, 'mo_bac', 'd-1', 'Bacteria natural mortality', default=0.05_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%mumax_bac, 'mumax_bac', 'd-1', 'Maximum labile DOC or NHS uptake by BAC', default=0.000154_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%q10_bac, 'q10_bac', '-', 'Temperature factor for BAC', default=2.0_rk) 
   call self%get_parameter(self%r_n_c_bac, 'r_n_c_bac', 'molN molC-1', 'N:C', default=0.196_rk) 
   call self%get_parameter(self%r_n_c_denit, 'r_n_c_denit', 'molN molC-1', 'N:C ratio of denitrification', default=0.8_rk) 
   call self%get_parameter(self%r_o2_c_resp, 'r_o2_c_resp', 'molO2 molC-1', 'O2:C ratio of respiration process', default=1.0_rk) 
   call self%get_parameter(self%r_odu_c_anox, 'r_odu_c_anox', 'molODU molC-1', 'ODU:C ratio in anoxic remineralization', default=1.0_rk) 
   call self%get_parameter(self%r_p_n_redfield, 'r_p_n_redfield', 'molP molN-1', 'N:P Redfield ratio in PHY', default=0.0625_rk) 


   ! Register state variables 

   call self%register_state_variable(self%id_bac, 'BAC', 'mmol C m-3', 'Bacterial biomass',minimum=0.0e-7_rk)

   call self%register_state_dependency(self%id_dcl, 'DCL', 'Labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dcs, 'DCS', 'Semi-labile detritus concentration in carbon', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dic, 'DIC', 'Dissolved inorganic carbon concentration', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dnl, 'DNL', 'Labile detritus concentration in nitrogen', 'mmol N m-3') 
   call self%register_state_dependency(self%id_dns, 'DNS', 'Semi-labile detritus concentration in nitrogen', 'mmol C m-3') 
   call self%register_state_dependency(self%id_dox, 'DOX', 'Dissolved oxygen concentration', 'mmol O2 m-3') 
   call self%register_state_dependency(self%id_nhs, 'NHS', 'Ammonium concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_nos, 'NOS', 'Nitrate concentration', 'mmol N m-3') 
   call self%register_state_dependency(self%id_odu, 'ODU', 'Oxygen demand unit concentration', 'mmol ODU m-3') 
   call self%register_state_dependency(self%id_pho, 'PHO', 'Phosphorus', 'mmol P m-3') 

    ! Register environmental dependencies 
   call self%register_dependency(self%id_temp, standard_variables%temperature) 
   call self%register_dependency(self%id_depth,standard_variables%pressure)  

    ! Add to aggregate variables 
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_bac)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_bac, scale_factor=self%r_n_c_bac)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_bac, scale_factor=self%r_n_c_bac*self%r_p_n_redfield)

    ! Register diagnostic variables 
   call self%register_diagnostic_variable(self%id_remineralization_anoxic_bac, 'remineralization_anoxic_bac', '-', '-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_oxygen_consumption_bac, 'oxygen_consumption_bac', '-','-', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_respiration_bac, 'respiration_bac', 'mmol C m-3 d-1', 'Respiration of bacteria', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_denitrification, 'denitrification', 'mmol N m-3 d-1','Denitrification rate', output=output_instantaneous) 
   call self%register_diagnostic_variable(self%id_uptake_dcl_bac, 'uptake_dcl_bac', 'mmol C m-3 d-1','Uptake of labile DOC by bacteria', output=output_instantaneous) 

   return 

99 call self%fatal_error('Bacteria', 'Error reading namelist ulg_bacteria') 

   end subroutine initialize 


   ! Right hand sides of Bacteria model
   subroutine do(self,_ARGUMENTS_DO_)

   class (type_ulg_bacteria), intent(in) :: self

   _DECLARE_ARGUMENTS_DO_

      real(rk) ::  DCL,DCS,DIC,DNL,DNS,DOX,NHS,NOS,ODU,PHO
      real(rk) ::  temp, depth
      real(rk) ::  BAC
      real(rk) ::   Denitrification	  ! mmol N m-3, Denitrification flux
      real(rk) ::   Excretion	  ! mmol N m-3, Zooplankton excretion of ammonium
      real(rk) ::   Growth	  ! mmol C m-3 d-1, Phytoplankton growth
      real(rk) ::   Iron	  ! mmol Fe m-3, Iron concentration
      real(rk) ::   Limitation_iron	  ! -, Limitation by Iron
      real(rk) ::   Limitation_nutrient	  ! -, Nutrient limitation for small flagellates
      real(rk) ::   Mortality_C	  ! mmol C m-3, Phytoplankton mortality flux
      real(rk) ::   Mortality_N	  ! mmol N m-3, Phytoplankton mortality flux
      real(rk) ::   Remineralization_anoxic_loc	  ! mmol ODU m-3, Bacterial anoxic remineralisation
      real(rk) ::   Respiration	  ! mmol C m-3, Bacteria respiration flux
      real(rk) ::   Respiration_loc	  ! mmol O2 m-3, Bacterial respiration
      real(rk) ::   Uptake_DCL	  ! mmol C m-3, Bacteria uptake of DOC
      real(rk) ::   Uptake_DNL	  ! mmol N m-3, Bacteria uptake of DON
      real(rk) ::   Uptake_NHS	  ! mmol N m-3, Ammonium uptake of phytoplankton
      real(rk) ::   Uptake_NHS_pot	  ! mmol N m-3, Bacteria potential uptake of ammonium
      real(rk) ::   testratio	  ! mmol N m-3, Value showing if BAC growth limited by C or N
      real(rk) ::   tf	  ! -, Temperature factor

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
    _GET_(self%id_temp,temp)            ! local temperature
    _GET_(self%id_depth,depth)            ! depth
    
    tf = Q10Factor(temp,self%q10_bac) 
    
   !Compute Iron Limitation 
    Iron = self%iron + self%i1_curve / (self%i2_curve * sqrt(2.*3.1416)) * exp(-(depth - 275.0)**2 / (2. * self%i2_curve**2))
    Limitation_iron = Michaelis(Iron,self%ks_odu_iron)

    Uptake_DCL = tf * self%mumax_bac * BAC * Michaelis(DCL,self%ks_dls_bac)
    Uptake_DNL = Uptake_DCL * (DNL/DCL)
    Limitation_nutrient = min(Michaelis(NHS,self%ks_nhs_bac),Michaelis(PHO,self%ks_po4_bac))
    Uptake_NHS_pot = self%mumax_bac * Limitation_nutrient * BAC * self%r_n_c_bac

   ! Test if bacterial growth is limited by carbon (DOClabile) or nitrogen (ammonium + DONLabile) 
    testratio = Uptake_DCL * (Uptake_DNL/Uptake_DCL - self%eff_gr_bac_c*self%r_n_c_bac)
    
   ! Growth rate of bacteria depends on a threshold value compared with Uptake_Potential_NHS 
    if (Uptake_NHS_pot > (-testratio)) then 
      ! In this case we are in a situation of carbon limitation 
      Growth = self%eff_gr_bac_c*Uptake_DCL
      ! Growth rate computed taking into account the iron limitation 
      Respiration = Uptake_DCL*(1.0 - self%eff_gr_bac_c)
      ! Test if NHS uptake is necessary to use all the DOCl 
      if (testratio > 0) then 
        ! We are in case of remineralisation of ammonium through bacteria excretion and no net uptake of ammonium is necessary to use all the DOC 
        Uptake_NHS = 0
        Excretion = testratio
      else 
        Uptake_NHS = -testratio
        Excretion = 0
      endif 
    else 
      ! if we are in case of nitrogen limitation,it means that all the DON and the potential uptake of NHS is not sufficient to consume all the DOC 
      Uptake_NHS = Uptake_NHS_pot
      Growth = (Uptake_NHS_pot + Uptake_DNL)/self%r_n_c_bac
      ! Growth rate computed taking into account the iron limitation 
      Excretion = 0
      Respiration = Growth * (1.0/self%eff_gr_bac_c  - 1.0)
    endif 
    
   ! Bacteria mortality rate (C_BACMort,N_BACMort, /day) 
    Mortality_C = self%mo_bac * tf * BAC
    Mortality_N = Mortality_C * self%r_n_c_bac
    
    Denitrification = Respiration * Michaelis(NOS,self%ks_denitr_nos) * Inhibition(DOX,self%ki_denit_o2) * self%r_n_c_denit
    Respiration_loc = Respiration * Michaelis(DOX,self%ks_oxic_o2) * self%r_o2_c_resp
    Remineralization_anoxic_loc = Respiration * Inhibition(NOS,self%ki_anox_nos) * Inhibition(DOX,self%ki_anox_o2) * self%r_odu_c_anox
    
   ! Carbon content increases by intake of nutrients and decreases by mortality (and predation) 
   _ADD_SOURCE_(self%id_bac, Growth - Mortality_C)
   ! Ammonium is excreted or can be taken up 
   _ADD_SOURCE_(self%id_nhs, Excretion - Uptake_NHS) 
   _ADD_SOURCE_(self%id_pho, (Excretion - Uptake_NHS) *self%r_p_n_redfield) 
   _ADD_SOURCE_(self%id_dcl, self%f_dl_dom * Mortality_C - Growth - Respiration) 
   _ADD_SOURCE_(self%id_dnl, self%f_dl_dom * Mortality_N - Uptake_DNL) 
   _ADD_SOURCE_(self%id_dcs, (1.0 - self%f_dl_dom) * Mortality_C) 
   _ADD_SOURCE_(self%id_dns, (1.0 - self%f_dl_dom) * Mortality_N) 
   _ADD_SOURCE_(self%id_dox, -Respiration_loc) 
   ! NOS decreases due to bacterial respiration 
   _ADD_SOURCE_(self%id_nos, -Denitrification) 
   _ADD_SOURCE_(self%id_dic, Respiration) 
   ! ODU increases due to bacterial anoxic remineralisation 
   _ADD_SOURCE_(self%id_odu, (1.0 - Limitation_iron * self%f_solid_odu) * Remineralization_anoxic_loc) 

   _SET_DIAGNOSTIC_(self%id_uptake_dcl_bac, Uptake_DCL)
   _SET_DIAGNOSTIC_(self%id_respiration_bac, Respiration)
   _SET_DIAGNOSTIC_(self%id_denitrification, Denitrification)
   _SET_DIAGNOSTIC_(self%id_oxygen_consumption_bac, Remineralization_anoxic_loc)
   _SET_DIAGNOSTIC_(self%id_remineralization_anoxic_bac, Respiration_loc)

   _LOOP_END_

   end subroutine do


   end module fabm_ulg_bacteria 
