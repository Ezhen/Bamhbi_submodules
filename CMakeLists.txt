add_library(fabm_models_bamhbi OBJECT
            bamhbi_model_library.f90
            chem/chemical.f90
            phytoplankton/diatoms.f90
            phytoplankton/emiliana.f90
            phytoplankton/flagellates.f90
            zooplankton/microzoo.f90
            zooplankton/mesozoo.f90
            jellyfish/gelatinous.f90
            jellyfish/noctiluca.f90
            bacteria/bacteria.f90
            dom/dom.F90
           )

target_link_libraries(fabm_models_bamhbi PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
