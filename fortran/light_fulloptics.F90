#include "fabm_driver.h"

! This is a FABM implementation of the two-band 
! spectral model used in the BAMHBI
!

module fabm_ulg_light

   use fabm_types
   use fabm_ulg_bamhbi_split_utilities

   implicit none

   private

   type, extends(type_base_model), public :: type_ulg_light
      type (type_state_variable_id)         :: id_poc
      type (type_state_variable_id)         :: id_cdi,id_ndi
      type (type_state_variable_id)         :: id_cem,id_nem
      type (type_state_variable_id)         :: id_cfl,id_nfl

      type (type_surface_dependency_id) :: id_swr0 ! Surface shortwave radiation
      type (type_dependency_id)         :: id_dz   ! Cell thickness
      type (type_dependency_id)         :: id_ext  ! Attentuation coefficient for PAR
      type (type_dependency_id)         :: id_z_r
      type (type_dependency_id)         :: id_salt
      !type (type_dependency_id)         :: id_TChl
      type (type_dependency_id)         :: id_ab_s, id_ab_l
      type (type_dependency_id)         :: id_bs_s, id_bs_l

      type (type_horizontal_dependency_id) :: id_lat, id_lon
      type (type_global_dependency_id)  :: id_yearday

      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)         :: id_par  ! Photosynthetically active radiation
      type (type_surface_diagnostic_variable_id) :: id_par0 ! Surface photosynthetically active radiation

      ! Parameters
      real(rk) :: ab_cdom_itc_s, ab_cdom_slp_s, ab_water_s
      real(rk) :: bs_water_s, bs_water_l
      real(rk) :: ab_cdom_itc_l, ab_cdom_slp_l, ab_water_l
      real(rk) :: bs_dia_s, bs_dia_l, bs_emi_s, bs_emi_l, bs_fla_s, bs_fla_l
      real(rk) :: ab_chl_s1, ab_chl_s2, ab_chl_l1, ab_chl_l2 
      real(rk) :: ab_poc_s, ab_poc_l, bs_poc_s, bs_poc_l
      real(rk) :: light_a, light_b
      real(rk) :: rmax_chl_n_dia, rmax_n_c_dia, rmin_chl_n_dia, rmin_n_c_dia
      real(rk) :: rmax_chl_n_emi, rmax_n_c_emi, rmin_chl_n_emi, rmin_n_c_emi
      real(rk) :: rmax_chl_n_fla, rmax_n_c_fla, rmin_chl_n_fla, rmin_n_c_fla

   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column

   end type type_ulg_light

   !type (type_bulk_standard_variable), parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg chla m-3',aggregate_variable=.true.)
   !type (type_bulk_standard_variable), parameter :: ab_short = type_bulk_standard_variable(name='absorption_of_swr',units='m-1',aggregate_variable=.true.)
   !type (type_bulk_standard_variable), parameter :: ab_long = type_bulk_standard_variable(name='absorption_of_lwr',units='m-1',aggregate_variable=.true.)
   !type (type_bulk_standard_variable), parameter :: bs_short = type_bulk_standard_variable(name='backscattering_of_swr',units='m-1',aggregate_variable=.true.)
   !type (type_bulk_standard_variable), parameter :: bs_long = type_bulk_standard_variable(name='backscattering_of_lwr',units='m-1',aggregate_variable=.true.)

contains

   subroutine initialize(self, configunit)
      class (type_ulg_light), intent(inout), target :: self
      integer,                 intent(in)            :: configunit
      call self%get_parameter(self%ab_chl_s1,  'ab_chl_s1',  'm-1','A*CHL^B', default=0.029_rk)						! a_chl_A_short
      call self%get_parameter(self%ab_chl_s2,  'ab_chl_s2',  'm-1','A*CHL^B', default=0.6_rk)						! a_chl_B_short
      call self%get_parameter(self%ab_chl_l1,  'ab_chl_l1',  'm-1', 'A*CHL^B', default=0.0066_rk)						! a_chl_A_long
      call self%get_parameter(self%ab_chl_l2,  'ab_chl_l2',  'm-1', 'A*CHL^B', default=0.8_rk)						! a_chl_B_long
      call self%get_parameter(self%ab_poc_s,  'ab_poc_s',  'm2 g-1', 'Absorption by organic matter', default=0.05_rk)			! a_poc_short
      call self%get_parameter(self%ab_poc_l,  'ab_poc_l',  'm2 g-1', 'Absorption by organic matter', default=0.0195_rk)			! a_poc_long

      call self%get_parameter(self%bs_dia_s,  'bs_dia_s',  'm2 (mg POC_DIA)-1','Backscattering by DI, shortwave', default=5.8e-06_rk)
      call self%get_parameter(self%bs_dia_l,  'bs_dia_l',  'm2 (mg POC_DIA)-1', 'Backscattering by DI, longwave', default=4.41e-06_rk)
      call self%get_parameter(self%bs_emi_s,  'bs_emi_s',  'm2 (mg POC_EMI)-1','Backscattering by EM, shortwave', default=1.008e-05_rk)
      call self%get_parameter(self%bs_emi_l,  'bs_emi_l',  'm2 (mg POC_EMI)-1', 'Backscattering by EM, longwave', default=7.56e-06_rk)
      call self%get_parameter(self%bs_fla_s,  'bs_fla_s',  'm2 (mg POC_FLA)-1','Backscattering by FL, shortwave', default=8.17e-06_rk)
      call self%get_parameter(self%bs_fla_l,  'bs_fla_l',  'm2 (mg POC_FLA)-1', 'Backscattering by FL, longwave', default=6.12e-06_rk)
      call self%get_parameter(self%bs_poc_s,  'bs_poc_s',  'm-1', 'Backscattering by POM, shortwave', default=0.0055_rk)
      call self%get_parameter(self%bs_poc_l,  'bs_poc_l',  'm-1', 'Backscattering by POM, longwave', default=0.005_rk)

      call self%get_parameter(self%ab_cdom_itc_s,  'ab_cdom_itc_s',  '?', '?', default=0.2522_rk)					! a_cdom_intercept_short
      call self%get_parameter(self%ab_cdom_slp_s,  'ab_cdom_slp_s',  '?', '?', default=-0.0122_rk)					! a_cdom_slope_short
      call self%get_parameter(self%ab_cdom_itc_l,  'ab_cdom_itc_l',  '?', '?', default=0.2155_rk)					! a_cdom_intercept_long
      call self%get_parameter(self%ab_cdom_slp_l,  'ab_cdom_slp_l',  '?', '?', default=-0.0113_rk)					! a_cdom_slope_long
      call self%get_parameter(self%ab_water_s,  'ab_water_s',  'm-1','Absorption by pure sea water', default=0.0196_rk)			! aa_w_short
      call self%get_parameter(self%bs_water_s,  'bs_water_s',  'm-1','Backscattering by pure sea water', default=0.0015_rk)		! bb_w_short
      call self%get_parameter(self%ab_water_l,  'ab_water_l',  'm-1', 'Absorption by pure sea water', default=0.24_rk)			! a_w_long
      call self%get_parameter(self%bs_water_l,  'bs_water_l',  'm-1', 'Backscattering by pure sea water', default=0.0007_rk)		! bb_w_long
      call self%get_parameter(self%light_a,  'light_a',  '-', 'Non-visible fraction of shortwave radiation', default=0.54_rk)		! LightAbsA
      call self%get_parameter(self%light_b,  'light_b',  '-', 'Fraction of longwave part of the light', default=0.63_rk)		! LightAbsB

      call self%get_parameter(self%rmax_chl_n_dia, 'rmax_chl_n_dia', 'g Chla molN-1', 'Maximum Chl:N ratio in DI', default=2.0_rk) 
      call self%get_parameter(self%rmax_n_c_dia, 'rmax_n_c_dia', 'molN molC-1', 'Maximum N:C ratio in DI', default=0.2_rk) 
      call self%get_parameter(self%rmin_chl_n_dia, 'rmin_chl_n_dia', 'g Chla molN-1', 'Minimum Chl:N ratio in DI', default=1.0_rk) 
      call self%get_parameter(self%rmin_n_c_dia, 'rmin_n_c_dia', 'molN molC-1', 'Minimum N:C ratio in DI', default=0.05_rk) 
      call self%get_parameter(self%rmax_chl_n_emi, 'rmax_chl_n_emi', 'g Chla molN-1', 'Maximum Chl:N ratio in EM', default=2.0_rk) 
      call self%get_parameter(self%rmax_n_c_emi, 'rmax_n_c_emi', 'molN molC-1', 'Maximum N:C ratio in EM', default=0.2_rk) 
      call self%get_parameter(self%rmin_chl_n_emi, 'rmin_chl_n_emi', 'g Chla molN-1', 'Minimum Chl:N ratio in EM', default=1.0_rk) 
      call self%get_parameter(self%rmin_n_c_emi, 'rmin_n_c_emi', 'molN molC-1', 'Minimum N:C ratio in EM', default=0.05_rk) 
      call self%get_parameter(self%rmax_chl_n_fla, 'rmax_chl_n_fla', 'g Chla molN-1', 'Maximum Chl:N ratio in FL', default=2.0_rk) 
      call self%get_parameter(self%rmax_n_c_fla, 'rmax_n_c_fla', 'mol N molC-1', 'Maximum N:C ratio in FL', default=0.2_rk) 
      call self%get_parameter(self%rmin_chl_n_fla, 'rmin_chl_n_fla', 'g Chla molN-1', 'Minimum Chl:N ratio in FL', default=1.0_rk) 
      call self%get_parameter(self%rmin_n_c_fla, 'rmin_n_c_fla', 'molN molC-1', 'Minimum N:C ratio in FL', default=0.05_rk) 

      ! Register diagnostic variables
      !call self%register_diagnostic_variable(self%id_swr, 'swr', 'W m-2', 'shortwave radiation', &
      !   standard_variable=standard_variables%downwelling_shortwave_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par, 'par', 'W m-2', 'photosynthetically active radiation', &
         standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par0, 'par0', 'W m-2', 'surface photosynthetically active radiation', &
         standard_variable=standard_variables%surface_downwelling_photosynthetic_radiative_flux, source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_swr0, standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_dz,   standard_variables%cell_thickness)
      call self%register_dependency(self%id_z_r, standard_variables%depth)
      call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
      ! for zendeg
      call self%register_dependency(self%id_lat,standard_variables%latitude)
      call self%register_dependency(self%id_lon,standard_variables%longitude)
      call self%register_dependency(self%id_yearday, standard_variables%number_of_days_since_start_of_the_year)

   end subroutine


   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_ulg_light), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: POC,CFL,NFL,CDI,NDI,NEM,CEM
      real(rk) :: par0, par
      real(rk) :: p1, p2
      real(rk) :: ext_l, ext_s
      real(rk) :: swr0
      real(rk) :: z_r, z_cdom, dz
      real(rk) :: ab_s, bs_s, ab_l, bs_l
      real(rk) :: cdom_est
      real(rk) :: salt , TChl
      real(rk) :: chla_fla, chla_dia, chla_emi
      !!
      real(rk) :: lat, lon
      real(rk) :: h_per_day, d_per_year, s_per_day
      real(rk) :: yearday, dt
      real(rk) :: nsec, days, th0, th02, th03
      real(rk) :: hour, sundec, thsun
      real(rk) :: nday_year, alat, coszen, zendeg

      real(rk), parameter :: hours_per_day = 24.0_rk
      real(rk), parameter :: days_per_year = 365.0_rk
      real(rk), parameter :: seconds_per_day = 86400.0_rk
      real(rk), parameter :: rad = 1.7453e-2

      _GET_SURFACE_(self%id_swr0,swr0)

      par0 = swr0 * (1.0_rk - self%light_a)

     ! Calculate light extinction swr and lwr

      _SET_SURFACE_DIAGNOSTIC_(self%id_par0,par0)

      p1=1 ; p2=1
      dt = 400  

      ! zendeg calculation
      _GET_GLOBAL_(self%id_yearday,yearday)

      _GET_HORIZONTAL_(self%id_lat,lat)
      _GET_HORIZONTAL_(self%id_lon,lon)

      nday_year = int(yearday + 1.0) !int(np.floor(i*dt/seconds_per_day)) + 1
      nsec = int(yearday*86400.0)

      days = int(nday_year - 1.0)

      th0 = 2. * 3.1415 * days / days_per_year	! in rad
      th02 = 2*th0
      th03 = 3*th0

      hour = ((nsec / seconds_per_day) - int(nsec / seconds_per_day)) * hours_per_day

      sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) - 0.006758*cos(th02) + 0.000907*sin(th02) - 0.002697*cos(th03) + 0.001480*sin(th03)

      thsun = ((hour-12.)*15. + lon) * rad
      alat = rad * lat
      coszen = sin(alat)*sin(sundec) + cos(alat)*cos(sundec)*cos(thsun)

      IF (coszen.lt.5.035e-4)  THEN
        coszen = 0.0
      END IF

      IF (coszen.gt.1.0) THEN
        zendeg = 0.0
      ELSE
        zendeg = acos(coszen)/rad
      END IF

      _DOWNWARD_LOOP_BEGIN_
          _GET_(self%id_dz,dz)     ! Layer height (m)

          ! Retrieve current (local) state variable values.
          _GET_(self%id_poc,POC)
          _GET_(self%id_cem,CEM)       ! Small flagellate biomass in carbon
          _GET_(self%id_nem,NEM)       ! Small flagellate biomass in nitrogen
          _GET_(self%id_cdi,CDI)       ! Diatom biomass in carbon
          _GET_(self%id_ndi,NDI)       ! Diatom biomass in nitrogen
          _GET_(self%id_cfl,CFL)       ! Large flagellate biomass in carbon
          _GET_(self%id_nfl,NFL)       ! Large flagellate biomass in nitrogen

          _GET_(self%id_z_r,z_r)	! depth of layer midpoints (m)
          _GET_(self%id_salt,salt)

          chla_fla = CFL * ratio_chl_c_phyt(Ratio(NFL,CFL),self%rmax_n_c_fla,self%rmin_n_c_fla,self%rmin_chl_n_fla,self%rmax_chl_n_fla) ! Chl flagellates
          chla_dia = CDI * ratio_chl_c_phyt(Ratio(NDI,CDI),self%rmax_n_c_dia,self%rmin_n_c_dia,self%rmin_chl_n_dia,self%rmax_chl_n_dia) ! Chl diatoms
          chla_emi = CEM * ratio_chl_c_phyt(Ratio(NEM,CEM),self%rmax_n_c_emi,self%rmin_n_c_emi,self%rmin_chl_n_emi,self%rmax_chl_n_emi) ! Chl emiliana
          TChl = chla_fla + chla_dia + chla_emi

          cdom_est = min(max(salt,10.0),18.5)

          z_cdom=1 ; if (z_r.gt.120.0) z_cdom=0

          ! shortwave adsorption
          ab_s =  self%ab_water_s &
               + (self%ab_chl_s1 * TChl**self%ab_chl_s2) &  
               + self%ab_poc_s * POC * 12.0 / 1000.0 &  
               + z_cdom * ( self%ab_cdom_itc_s + self%ab_cdom_slp_s * cdom_est) 

          ! (back)scattering shortwave
          bs_s =  self%bs_water_s &	
               + self%bs_dia_s * CDI * 12.0 + self%bs_fla_s * CFL * 12.0 + self%bs_emi_s * CEM * 12.0 &
               + self%bs_poc_s * POC * 12.0 / 1000.0 

          ! longwave adsorption
          ab_l =  self%ab_water_l & 
               + (self%ab_chl_l1 * TChl**self%ab_chl_l2) & 
               + self%ab_poc_l * POC * 12.0 / 1000.0 & 
               + z_cdom * ( self%ab_cdom_itc_l + self%ab_cdom_slp_l * cdom_est)

          ! (back)scattering longwave
          bs_l =  self%bs_water_l & 
               + self%bs_dia_l * CDI * 12.0 + self%bs_fla_l * CFL * 12.0 + self%bs_emi_l * CEM * 12.0 &
               +  self%bs_poc_l * POC * 12.0 / 1000.0 

          ext_s = (1.0 + 0.005*zendeg) * ab_s + 4.18 * (1.0 - 0.52 * exp(-10.8*ab_s)) * bs_s  ! extinction short			!0.005 * ZENDEG
          ext_l = (1.0 + 0.005*zendeg) * ab_l + 4.18 * (1.0 - 0.52 * exp(-10.8*ab_l)) * bs_l  ! extinction long			!0.005 * ZENDEG

          par = par0/dz * ( self%light_b*p1/ext_l * (1.0-exp(-ext_l*dz))  +  (1.0-self%light_b)*p2/ext_s * (1.0-exp(-ext_s*dz)) )
          p1 = p1*exp(-ext_l*dz)
          p2 = p2*exp(-ext_s*dz)

          _SET_DIAGNOSTIC_(self%id_par,par)		 ! Photosynthetically active radiation at layer center

      _DOWNWARD_LOOP_END_
   end subroutine do_column

end module fabm_ulg_light
