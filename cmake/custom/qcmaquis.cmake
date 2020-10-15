#.rst:
#
# Detects QCMaquis DMRG module.
#
# Variables used::
#
#   QCMaquis_ROOT         - User-set path to the QCMaquis DMRG module
#
# Variables defined::
#
#   QCMaquis_FOUND        - TRUE if the QCMaquis DMRG module was found
#

# User indicated to try using a pre-built QCMaquis
if(${QCMaquis_ROOT} STREQUAL "None")
else()
  message(STATUS "QCMaquis search activated in ${QCMaquis_ROOT}")
  find_package(MAQUIS_DMRG
               REQUIRED
               CONFIG
               PATHS "${QCMaquis_ROOT}/share"
              )
endif()
# Build QCMaquis as external package if pre-built version not found or not indicated
if(NOT MAQUIS_DMRG_FOUND)
  message(STATUS "QCMaquis DMRG not found. A pre-packaged version will be built.")
else()
  message(STATUS "Existing version of QCMaquis found at ${QCMaquis_ROOT} which will be used.")
  add_library(qcmaquis-suite ALIAS maquis_dmrg)
endif()

