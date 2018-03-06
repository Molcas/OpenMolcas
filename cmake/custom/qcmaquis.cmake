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
  message(STATUS "Existing version of QCMaquis found at ${QCMaquis_ROOT}\n   which will be used. Remember to set all environment variables (PATH, LD_LIBRARY_PATH, ...) correctly.\n   HINT: use `source ${QCMaquis_ROOT}/bin/qcmaquis.sh`")
endif()
