add_library(fabm_models_ulg OBJECT
            ulg_model_library.F90
            bamhbi_utilities.F90
            chem/chemical.F90
            phytoplankton/diatoms.F90
            phytoplankton/emiliana.F90
            phytoplankton/flagellates.F90
            zooplankton/microzoo.F90
            zooplankton/mesozoo.F90
            jellyfish/gelatinous.F90
            jellyfish/noctiluca.F90
            bacteria/bacteria.F90
            dom/dom.F90
           )

target_link_libraries(fabm_models_ulg PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
