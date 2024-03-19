module ulg_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: ulg_model_factory

contains

   subroutine create(self,name,model)

      use fabm_ulg_chemical
      use fabm_ulg_diatoms
      use fabm_ulg_emiliana
      use fabm_ulg_flagellates
      use fabm_ulg_microzoo
      use fabm_ulg_mesozoo
      use fabm_ulg_gelatinous
      use fabm_ulg_noctiluca
      use fabm_ulg_bacteria
      use fabm_ulg_dom

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('chemical');     allocate(type_ulg_chemical::model)
         case ('diatoms');      allocate(type_ulg_diatoms::model)
         case ('emiliana');	allocate(type_ulg_emiliana::model)
         case ('flagellates');	allocate(type_ulg_flagellates::model)
         case ('microzoo'); 	allocate(type_ulg_microzoo::model)
         case ('mesozoo');	allocate(type_ulg_mesozoo::model)
         case ('gelatinous');	allocate(type_ulg_gelatinous::model)
         case ('noctiluca');	allocate(type_ulg_noctiluca::model)
         case ('bacteria');	allocate(type_ulg_bacteria::model)
         case ('dom');		allocate(type_ulg_dom::model)
      end select

   end subroutine

end module
