#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#***********************************************************************

#***********************************************************************
# CMakeLists.txt for building a bundled copy of libwignernj            *
#                                                                      *
# Used as a fallback when no system libwignernj is detected by         *
# find_package. libwignernj evaluates Wigner 3j/6j/9j symbols (and     *
# related coupling coefficients) exactly, via prime factorisation and  *
# multi-word integer arithmetic.                                       *
#***********************************************************************

# load External Project macro
include (ExternalProject)

set (CUSTOM_WIGNERNJ_LOCATION ${PROJECT_BINARY_DIR}/External/wignernj)

# libwignernj does not know the "profile" build type
if (CMAKE_BUILD_TYPE MATCHES "profile")
  set(WIGNERNJ_BUILD_TYPE "release")
else ()
  set(WIGNERNJ_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif ()

# libwignernj is built with OpenMolcas's baseline C and Fortran flags
# (CMAKE_C_FLAGS_DEFAULT / CMAKE_Fortran_FLAGS_DEFAULT) so that the
# Fortran .mod file it produces is ABI-compatible with the OpenMolcas
# objects that `use` it -- NAG, in particular, encodes the -kind= setting
# into every .mod and rejects a mismatch. The _DEFAULT baselines are used
# in preference to the full CMAKE_<lang>_FLAGS because the latter also
# carry flags that must not reach this sub-build: -Werror (a bundled
# library should not be held to OpenMolcas's zero-warning policy) and
# -fsyntax-only (set by the syntax-check CI, under which no object files
# are produced). This mirrors how the bundled Libxc is built.
set (WIGNERNJ_C_FLAGS ${CMAKE_C_FLAGS_DEFAULT})
if (CMAKE_C_COMPILER_ARG1)
  set (WIGNERNJ_C_FLAGS "${CMAKE_C_COMPILER_ARG1} ${WIGNERNJ_C_FLAGS}")
endif ()
set (WIGNERNJ_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEFAULT})
if (CMAKE_Fortran_COMPILER_ARG1)
  set (WIGNERNJ_Fortran_FLAGS "${CMAKE_Fortran_COMPILER_ARG1} ${WIGNERNJ_Fortran_FLAGS}")
endif ()
list (APPEND WIGNERNJCMakeArgs
      -DCMAKE_BUILD_TYPE=${WIGNERNJ_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
      -DCMAKE_INSTALL_LIBDIR=lib
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_C_FLAGS=${WIGNERNJ_C_FLAGS}
      -DCMAKE_Fortran_FLAGS=${WIGNERNJ_Fortran_FLAGS}
      -DCMAKE_POSITION_INDEPENDENT_CODE=ON
      -DBUILD_SHARED_LIBS=OFF
      -DBUILD_FORTRAN=ON
      -DBUILD_TESTS=OFF
      -DBUILD_CXX_TESTS=OFF
      -DBUILD_EXAMPLES=OFF
      -DBUILD_PYTHON=OFF
      -DBUILD_LTO=OFF
)

##############################################
# git references for the libwignernj library #
##############################################
set (reference_git_repo https://github.com/susilehtola/libwignernj.git)
set (reference_git_commit v0.6.0)
set (EP_PROJECT wignernj)

# Enabling source changes to keep ExternalProject happy
set (CMAKE_DISABLE_SOURCE_CHANGES OFF)

set (last_hash "None")
set (hash_file ${CUSTOM_WIGNERNJ_LOCATION}/${EP_PROJECT}.hash)
if (EXISTS ${hash_file})
  file (READ ${hash_file} last_hash)
  string (REGEX REPLACE "\n$" "" last_hash "${last_hash}")
endif ()
if (last_hash STREQUAL ${reference_git_commit})
  set (EP_SkipUpdate ON)
else ()
  set (EP_SkipUpdate OFF)
endif ()

# BUILD_BYPRODUCTS is required by the Ninja generator: the static
# libraries are linked into OpenMolcas executables by file path, so Ninja
# must know that this ExternalProject is the rule that produces them.
ExternalProject_Add (${EP_PROJECT}
                     PREFIX ${CUSTOM_WIGNERNJ_LOCATION}
                     CMAKE_ARGS "${WIGNERNJCMakeArgs}"
                     GIT_REPOSITORY ${reference_git_repo}
                     GIT_TAG ${reference_git_commit}
                     GIT_PROGRESS 1
                     UPDATE_DISCONNECTED ${EP_SkipUpdate}
                     INSTALL_DIR "${PROJECT_BINARY_DIR}"
                     BUILD_BYPRODUCTS
                       ${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}wignernj_f03${CMAKE_STATIC_LIBRARY_SUFFIX}
                       ${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}wignernj${CMAKE_STATIC_LIBRARY_SUFFIX}
)

ExternalProject_Add_Step (${EP_PROJECT} update_hash
                          COMMAND echo ${reference_git_commit} > ${hash_file}
                          DEPENDEES build
)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)

# set variables for use in the parent CMakeLists.txt
ExternalProject_Get_Property (${EP_PROJECT} install_dir)
set (WIGNERNJ_INCLUDE ${install_dir}/include/wignernj/fortran PARENT_SCOPE)
set (WIGNERNJ_LIBRARIES
     ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}wignernj_f03${CMAKE_STATIC_LIBRARY_SUFFIX}
     ${install_dir}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}wignernj${CMAKE_STATIC_LIBRARY_SUFFIX}
     m
     PARENT_SCOPE
)
