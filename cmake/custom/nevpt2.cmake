#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2017,2018, Leon Freitag                                *
#***********************************************************************
#                                                                      *
# Based on HDF5_QCM Makefile by Stefan Knecht                          *
#                                                                      *
#***********************************************************************
# CMakeLists.txt for NEVPT2                                            *
#***********************************************************************

# load External Project macro
include(ExternalProject)
# Set up compilation of QCMaquis components
set(CUSTOM_NEVPT2_LOCATION ${PROJECT_BINARY_DIR}/External/nevpt2_ext)
set(CUSTOM_QCMaquis_LOCATION ${PROJECT_BINARY_DIR}/External/qcmaquis)

# QCMaquis does not know profile
if(CMAKE_BUILD_TYPE MATCHES "profile")
  set(NEVPT2_BUILD_TYPE "release")
else()
  set(NEVPT2_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_ROOT})
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/custom)

set(OPENMOLCAS_TOOLS_DIR ${CMAKE_BINARY_DIR}/Tools/distributed-4rdm)
if(SINGLE_MOD_DIR)
  set(mod_dir ${MAIN_MOD_DIR}/_single)
else()
  set(mod_dir ${MAIN_MOD_DIR}/nevpt2)
endif()

# CMake does not support passing lists inside lists, so we need to
# replace the semicolons in the lists and pass them as normal strings
# and then replace the new separators with semicolons on the other side

string(REPLACE ";" "<->" LINALG_LIBRARIES_NEVPT "${LINALG_LIBRARIES}")
        list(APPEND NEVPT2CMakeArgs
  "-DCMAKE_BUILD_TYPE=${NEVPT2_BUILD_TYPE}"
  "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/External"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMake_C_FLAGS}"
  "-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>"
  "-DCMAKE_INSTALL_LIBDIR=lib"
  "-DENABLE_DEBUG_DMRG=OFF"
  "-DENABLE_DMRG=ON"
  "-DENABLE_OPENMP=${OPENMP}"
  "-DMOLCAS_MATHLIBS=${LINALG_LIBRARIES_NEVPT}"
  "-DENABLE_MOLCAS=ON"
  "-DMOLCAS_BUILD_DIR=${PROJECT_BINARY_DIR}"
  "-DCMAKE_Fortran_MODULE_DIRECTORY=${mod_dir}"
  "-DDMRG_INCLUDE=${HDF5_QCM_INCLUDE}"
  )

if (MPI AND GA)
# MPI and GA
list(APPEND GA_LIBRARIES "${MPI_Fortran_LIBRARIES}")
string(REPLACE ";" "<->" GA_LIBRARIES_NEVPT "${GA_LIBRARIES}")
  list(APPEND NEVPT2CMakeArgs
    "-DGA_LIBS=${GA_LIBRARIES_NEVPT}")
endif()

if (QCMaquis_ROOT)
list(APPEND NEVPT2CMakeArgs
  "-DMAQUIS_DMRG_DIR=${QCMaquis_ROOT}/share/cmake")
endif()


######################################
# git references for NEVPT2          #
######################################
set(reference_git_repo https://github.com/qcscine/nevpt2.git)
# set(reference_git_tag release-3.0) # uncomment before merging into master, since before merging we expect more patches into upstream

set(EP_PROJECT nevpt2_ext)

# Enabling source changes to keep ExternalProject happy
set (CMAKE_DISABLE_SOURCE_CHANGES OFF)

ExternalProject_Add(${EP_PROJECT}
                    PREFIX ${CUSTOM_NEVPT2_LOCATION}
                    GIT_REPOSITORY ${reference_git_repo}
                    GIT_TAG ${reference_git_tag}
                    CMAKE_ARGS "${NEVPT2CMakeArgs}"
                    INSTALL_DIR "${PROJECT_BINARY_DIR}/qcmaquis"
                   )

ExternalProject_Get_Property(${EP_PROJECT} install_dir)
set(TOOL_SUBDIR Tools/distributed-4rdm)
file(MAKE_DIRECTORY ${OPENMOLCAS_TOOLS_DIR})
ExternalProject_Add_Step(${EP_PROJECT} install_tools DEPENDEES install
                         COMMAND ${CMAKE_COMMAND} -E copy_if_different ${install_dir}/${TOOL_SUBDIR}/jobmanager.py
                                                                       ${install_dir}/${TOOL_SUBDIR}/prepare_rdm_template.sh
                                                                       ${install_dir}/${TOOL_SUBDIR}/submit-3rdm.sh
                                                                       ${install_dir}/${TOOL_SUBDIR}/submit-4rdm.sh
                                                                       ${OPENMOLCAS_TOOLS_DIR})

set (CMAKE_DISABLE_SOURCE_CHANGES ON)

# set variables for use in parent CMakeLists.txt
ExternalProject_Get_Property(${EP_PROJECT} BINARY_DIR)

set(NEVPT2_INCLUDE ${mod_dir} PARENT_SCOPE)
set(NEVPT2_LIBRARIES ${BINARY_DIR}/${CMAKE_FIND_LIBRARY_PREFIXES}qdnevpt2.a ${BINARY_DIR}/${CMAKE_FIND_LIBRARY_PREFIXES}AUXLIB_F.a PARENT_SCOPE)
