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
# Copyright (C) 2014-2016, Steven Vancoillie                           *
#               2014-2023, Ignacio Fdez. Galván                        *
#               2015, Stefan Knecht                                    *
#               2020, Gerardo Raggi                                    *
#***********************************************************************
# Global CMakeLists.txt file for building OpenMolcas with CMake.
#
# Steven Vancoillie - Spring 2014, on a Saturday afternoon, with coffee
# Ignacio Fdez. Galván made several improvements / additions
# Gerardo Raggi - AppleClang

cmake_minimum_required (VERSION 3.12)
cmake_policy (VERSION 3.12...3.21)

project (OpenMolcas
         LANGUAGES Fortran C
         VERSION ""
)

#set (CMAKE_VERBOSE_MAKEFILE ON)

include (${PROJECT_SOURCE_DIR}/cmake/cmake_helpers.cmake)

# Prevent in-source builds
# If an in-source build is attempted, you will still need to clean up a few files manually.
set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
get_filename_component (sourcedir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component (binarydir "${CMAKE_BINARY_DIR}" REALPATH)
if ("${sourcedir}" STREQUAL "${binarydir}")
  message (FATAL_ERROR "In-source builds in ${CMAKE_BINARY_DIR} are not "
          "allowed, please remove ./CMakeCache.txt and ./CMakeFiles/, create a "
          "separate build directory and run cmake from there.")
endif ()

# discover some necessary tools
find_program (GIT "git")
find_program (PERL "perl")
mark_as_advanced (FORCE GIT PERL)

FIND_PACKAGE (Python COMPONENTS Interpreter)

################################################################################
#                                                                              #
# Options                                                                      #
#                                                                              #
################################################################################

# The build type, specifies the different optimization types.
if (NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Type of build, options are: None (CFLAGS/FFLAGS can be used), Debug, Garble, RelWithDebInfo, Release, Fast."
       FORCE)
endif ()
set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS None Debug Garble RelWithDebInfo Release Fast)

# The place where OpenMolcas is installed with make install.
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "/opt/OpenMolcas" CACHE PATH "Location where OpenMolcas will be installed." FORCE)
endif ()
set (EXTRA_INSTALL_TARGETS "")
set (EXTRA_INSTALL_DIRS "")
set (EXTRA_INSTALL_IMPORTED_LIBS "")

option (INSTALL_TESTS "Install the test files and verify script." OFF)

set (DEBUG_DEFS "${DEBUG_DEFS}" CACHE STRING "Additional user-supplied definitions.")
separate_arguments (DEBUG_DEFS)

option (BUILD_SHARED_LIBS "Build dynamically-linked libraries." OFF)
option (BUILD_STATIC_LIBS "Build statically-linked libraries too." OFF)
if (BUILD_SHARED_LIBS)
  mark_as_advanced (CLEAR BUILD_STATIC_LIBS)
else ()
  mark_as_advanced (FORCE BUILD_STATIC_LIBS)
endif ()

option (MPI "Enable MPI parallellization." OFF)
option (GA "Use Global Arrays library." OFF)

option (OPENMP "Enable multi-threading." OFF)
option (IPO "Enable interprocedural optimization (if supported)." OFF)

option (BOUNDS "Enable bounds checking (only GNU Fortran)." OFF)

option (EXPERT "Show advanced compiler flags in CMake GUI)." OFF)

option (BIGOT "Do not allow any compiler warning (treat them as errors)." OFF)

set (LINALG_OPTIONS "Choose linear algebra library, options: Internal (default), Runtime, MKL, OpenBLAS, AOCL, Accelerate, Manual.")
if (DEFINED LINALG)
  set (LINALG "${LINALG}" CACHE STRING "${LINALG_OPTIONS}")
  if (NOT ";Internal;Runtime;MKL;OpenBLAS;AOCL;Accelerate;Manual;" MATCHES ";${LINALG};")
      message (FATAL_ERROR "Unsuported LINALG: ${LINALG}")
  endif ()
else ()
  set (LINALG "Internal" CACHE STRING "${LINALG_OPTIONS}")
endif ()
set_property (CACHE LINALG PROPERTY STRINGS Internal Runtime MKL OpenBLAS AOCL Accelerate Manual)

option (CUBLAS "Enable CUDA BLAS library." OFF)
option (NVBLAS "Enable NVidia BLAS library." OFF)

set (EXTERNAL_LIBXC "" CACHE STRING "Location for an external Libxc build")

option (FDE "Enable Frozen-density-embedding (FDE) interface." OFF)
option (GROMACS "Compile Gromacs interface." OFF)
option (BLOCK "Activate BLOCK-DMRG support." OFF)
option (CHEMPS2 "Activate CheMPS2-DMRG support." OFF)
option (DICE "Activate Dice-SHCI support." OFF)

option (GPERFTOOLS "Activate gperftools CPU profiling." OFF)
option (GCOV "Activate code coverage profiling (use with Debug/-O0)." OFF)
option (HDF5 "Activate HDF5 support for wavefunction format." ON)

option (TOOLS "Compile supported tools." OFF)

option (BUILD_TESTING "Build the unit_tests." ON)

if (";Debug;Garble;" MATCHES ";${CMAKE_BUILD_TYPE};")
  add_compile_definitions (_ADDITIONAL_RUNTIME_CHECK_)
endif ()

if (";Debug;Garble;RelWithDebInfo;" MATCHES ";${CMAKE_BUILD_TYPE};" OR
    CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" OR
    BIGOT)
  add_compile_definitions (_WARNING_WORKAROUND_)
endif ()

# External projects
option (DMRG     "Activate QCMaquis DMRG driver and library."                                                 OFF)
option (EFPLIB   "Enable EFPLib library for effective fragment potentials (requires External/efp submodule)." OFF)
option (GA_BUILD "Build Global Arrays."                                                                       OFF)
option (GEN1INT  "Enable Gen1Int library for 1-electron integrals."                                           OFF)
option (MolGUI   "Enable the Molcas Graphical User Interface."                                                OFF)
option (MSYM     "Activate MSYM support (requires External/libmsym submodule)."                               OFF)
option (NECI     "Activate NECI support (requires External/NECI submodule)."                                  OFF)
option (NEVPT2   "Activate (DMRG)-NEVPT2 support."                                                            OFF)
option (WFA      "Activate extended wavefunction analysis (requires External/libwfa submodule)."              OFF)

################################################################################
#                                                                              #
# Global settings                                                              #
#                                                                              #
################################################################################

# mark that we are running CMake (for subprocesses)
set (ENV{CMAKE_SESSION} "OpenMolcas")

# project name and supported languages
#=====================================
message ("Configuring compilers:")

if (WFA OR MolGUI OR DMRG)
  enable_language (CXX)
endif ()

#avoid using CXX compiler for linking
set (CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES "FALSE")

add_compile_definitions (_MOLCAS_)

#set cmake module search paths
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

# location of molcas-extra
if (NOT DEFINED EXTRA)
  if (DEFINED ENV{EXTRA})
    set (EXTRA $ENV{EXTRA})
  else ()
    set (EXTRA "")
  endif ()
endif ()
if (NOT "${EXTRA}" STREQUAL "")
  if (NOT IS_ABSOLUTE "${EXTRA}")
    set (EXTRA "${PROJECT_BINARY_DIR}/${EXTRA}")
  endif ()
  get_filename_component (EXTRA "${EXTRA}" ABSOLUTE)
  if (NOT EXISTS "${EXTRA}")
    message (WARNING "The molcas-extra directory ${EXTRA} does not exist.")
  else ()
    FILE (STRINGS "${EXTRA}/.molcashome" DIR_TYPE)
    if ("${DIR_TYPE}" STREQUAL "molcas-extra")
      set (EXTRA_DIR ${EXTRA})
      add_compile_definitions (_HAVE_EXTRA_)
    else ()
      message (FATAL_ERROR "The directory ${EXTRA} does not contain molcas-extra.")
    endif ()
  endif ()
endif ()
set (EXTRA "${EXTRA}" CACHE PATH "Location of the molcas-extra directory." FORCE)

# possible locations for source code
set (OPENMOLCAS_DIR ${CMAKE_CURRENT_LIST_DIR})
if (NOT ${PROJECT_SOURCE_DIR} STREQUAL ${OPENMOLCAS_DIR})
  set (basedirs ${PROJECT_SOURCE_DIR})
  set (EXTRA_DIR ${PROJECT_SOURCE_DIR})
  add_compile_definitions (_HAVE_EXTRA_)
endif ()
list (APPEND basedirs ${OPENMOLCAS_DIR})
if (NOT "${EXTRA}" STREQUAL "")
  if (NOT ";${basedirs};" MATCHES ";${EXTRA};")
    list (APPEND basedirs ${EXTRA})
  endif ()
endif ()
foreach (BASE_DIR ${basedirs})
  file (GLOB source_roots_tmp RELATIVE ${PROJECT_SOURCE_DIR} "${BASE_DIR}/src*")
  list (SORT source_roots_tmp)
  list (APPEND source_roots ${source_roots_tmp})
endforeach ()
list (REVERSE basedirs)

# get version from git if available, otherwise look in .molcasversion
#====================================================================

message ("Detecting Molcas version info:")

set (ENV{MOLCAS} ${PROJECT_BINARY_DIR})
set (ENV{OPENMOLCAS} ${OPENMOLCAS_DIR})
if (DEFINED EXTRA_DIR)
  set (ENV{MOLCAS_SOURCE} ${EXTRA_DIR})
else ()
  set (ENV{MOLCAS_SOURCE} ${OPENMOLCAS_DIR})
endif ()

execute_process (
  COMMAND ${GIT} "rev-parse" "--git-dir"
  WORKING_DIRECTORY ${OPENMOLCAS_DIR}
  OUTPUT_VARIABLE GIT_REVPARSE_OUTPUT
  ERROR_VARIABLE GIT_REVPARSE_ERROR
  RESULT_VARIABLE GIT_REVPARSE_RC
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if (GIT_REVPARSE_RC)
  set (OPENMOLCAS_GIT_REPO "${OPENMOLCAS_DIR}/.git")
else ()
  set (OPENMOLCAS_GIT_REPO "${GIT_REVPARSE_OUTPUT}")
  if (NOT IS_ABSOLUTE "${OPENMOLCAS_GIT_REPO}")
    set (OPENMOLCAS_GIT_REPO "${OPENMOLCAS_DIR}/${OPENMOLCAS_GIT_REPO}")
  endif ()
endif ()
set (MOLCAS_VERSION_FILE "${PROJECT_SOURCE_DIR}/.molcasversion")

if (DEFINED EXTRA_DIR)
  execute_process (
    COMMAND ${GIT} "rev-parse" "--git-dir"
    WORKING_DIRECTORY ${EXTRA_DIR}
    OUTPUT_VARIABLE GIT_REVPARSE_OUTPUT
    ERROR_VARIABLE GIT_REVPARSE_ERROR
    RESULT_VARIABLE GIT_REVPARSE_RC
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if (GIT_REVPARSE_RC)
    set (EXTRA_GIT_REPO "${EXTRA_DIR}/.git")
  else ()
    set (EXTRA_GIT_REPO "${GIT_REVPARSE_OUTPUT}")
    if (NOT IS_ABSOLUTE "${EXTRA_GIT_REPO}")
      set (EXTRA_GIT_REPO "${EXTRA_DIR}/${EXTRA_GIT_REPO}")
    endif ()
  endif ()
endif ()

if (EXISTS ${GIT})
  foreach (BASE_DIR ${basedirs})
    if ("${BASE_DIR}" STREQUAL "${OPENMOLCAS_DIR}")
      if (EXISTS ${OPENMOLCAS_GIT_REPO})
        set (MOLCAS_GIT_REPO ${OPENMOLCAS_GIT_REPO})
      else ()
        continue ()
      endif ()
      set (label "")
    elseif ("${BASE_DIR}" STREQUAL "${EXTRA_DIR}")
      if (EXISTS ${EXTRA_GIT_REPO})
        set (MOLCAS_GIT_REPO ${EXTRA_GIT_REPO})
      else ()
        continue ()
      endif ()
      set (label "-extra")
    endif ()
    # always set the HEAD as a source for configuring
    configure_file ("${MOLCAS_GIT_REPO}/HEAD" git-head${label} COPYONLY)
    # if it is a ref to a named branch, include the branch ref too
    file (READ "${MOLCAS_GIT_REPO}/HEAD" HEAD_STRING)
    string (STRIP "${HEAD_STRING}" HEAD_STRING)
    string (REPLACE "ref: " "" HEAD_REF "${HEAD_STRING}")
    if (EXISTS "${MOLCAS_GIT_REPO}/${HEAD_REF}")
      configure_file ("${MOLCAS_GIT_REPO}/${HEAD_REF}" git-head-ref${label} COPYONLY)
    endif ()
    # now, get the actual versions from the git describe tool
    execute_process (
      COMMAND ${GIT} "describe" "--always" "--dirty" "--match" "v*"
      WORKING_DIRECTORY ${BASE_DIR}
      OUTPUT_VARIABLE MOLCAS_VERSION
      ERROR_VARIABLE GIT_DESCRIBE_ERROR
      RESULT_VARIABLE GIT_DESCRIBE_RC
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if (GIT_DESCRIBE_RC)
      message (FATAL_ERROR "Failed to run git: ${GIT_DESCRIBE_ERROR}")
    endif ()
    if ("${BASE_DIR}" STREQUAL "${OPENMOLCAS_DIR}")
      set (OPENMOLCAS_VERSION ${MOLCAS_VERSION})
    elseif ("${BASE_DIR}" STREQUAL "${EXTRA_DIR}")
      set (EXTRA_VERSION ${MOLCAS_VERSION})
    endif ()
  endforeach (BASE_DIR)

  # used for external projects (git submodules)
  set (DEVELOPMENT_CODE)
  execute_process (
    COMMAND ${GIT} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    ERROR_QUIET
  )
  if (NOT ${GIT_BRANCH} STREQUAL "")
    string (STRIP ${GIT_BRANCH} GIT_BRANCH)
  endif ()
endif ()

if (DEFINED OPENMOLCAS_VERSION)
# Read the version from the .tags file, which should be populated when the code
# has been obtained from "git archive". Since there's no way to know if the
# code has been modified, add a question mark
elseif (EXISTS ${OPENMOLCAS_DIR}/.tags)
  file (STRINGS ${OPENMOLCAS_DIR}/.tags LINES)
  foreach (line ${LINES})
    string (FIND "${line}" "(" has_paren)
    string (FIND "${line}" "%" has_percent)
    string (REGEX MATCH "tag: (v[^,]*)," match ${line})
    # Parse the tag if present
    if (has_paren GREATER -1)
      if (NOT ${match} STREQUAL "")
        set (OPENMOLCAS_VERSION "${CMAKE_MATCH_1} ?")
      endif ()
    # Or just use the commit sha, making sure this is not the raw file
    elseif (has_percent EQUAL -1)
      string (SUBSTRING ${line} 0 14 OPENMOLCAS_VERSION)
      set (OPENMOLCAS_VERSION "o${OPENMOLCAS_VERSION} ?")
    endif ()
  endforeach ()
elseif (EXISTS ${MOLCAS_VERSION_FILE})
  file (STRINGS ${MOLCAS_VERSION_FILE} LINES)
  foreach (line ${LINES})
    string (FIND "${line}" ".o" is_open)
    string (FIND "${line}" ".x" is_extra)
    if (is_open GREATER -1)
      string (STRIP ${line} OPENMOLCAS_VERSION)
    elseif (is_extra GREATER -1)
      string (STRIP ${line} EXTRA_VERSION)
    endif ()
  endforeach ()
else ()
  set (OPENMOLCAS_VERSION "unknown")
  set (GIT_BRANCH "")
endif ()

message ("-- OPENMOLCAS_VERSION: ${OPENMOLCAS_VERSION}")
set (MOLCAS_VERSION ${OPENMOLCAS_VERSION})
if (DEFINED EXTRA_VERSION)
  message ("-- EXTRA_VERSION: ${EXTRA_VERSION}")
  set (MOLCAS_VERSION "${OPENMOLCAS_VERSION} & ${EXTRA_VERSION}")
endif ()

# set location of libraries and executables
#==========================================

set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# location for external projects (git submodules)
#================================================

set (extprojsdir "${OPENMOLCAS_DIR}/External")

# global include directory
#=========================

set (MAIN_INCLUDE_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include)
configure_file (${OPENMOLCAS_DIR}/src/Include/molcasversion.h.in ${MAIN_INCLUDE_DIR}/molcasversion.h)

foreach (BASE_DIR ${basedirs})
  if (EXISTS "${BASE_DIR}/src/Include")
    include_directories (${BASE_DIR}/src/Include)
  endif ()
endforeach ()
include_directories (${MAIN_INCLUDE_DIR})

# system information
#===================

message ("Detecting system info:")

# operating system
set (OS ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})

message ("-- OS: ${OS}")

# address mode
if (${CMAKE_SIZEOF_VOID_P} EQUAL 8)
  set (ADDRMODE 64)
  add_compile_definitions (_I8_)
else ()
  set (ADDRMODE 32)
endif ()

message ("-- ADDRMODE: ${ADDRMODE}")

# platform settings
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  add_compile_definitions (_LINUX_)
  if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" OR
      ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i686")
    set (PLATFORM "LINUX")
  elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set (PLATFORM "LINUX64")
  elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ia64")
    set (PLATFORM "LINUX64_IA")
  elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc_64" OR
          ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
    set (PLATFORM "PPC64")
  endif ()
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  add_compile_definitions (_LINUX_ _DARWIN_ _DARWIN_C_SOURCE)
  set (CMAKE_MACOSX_RPATH 0)
  set (PLATFORM "MacOS")
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
  add_compile_definitions (_LINUX_ _CYGWIN_)
  if (ADDRMODE EQUAL 32)
    set (PLATFORM "WIN")
  elseif (ADDRMODE EQUAL 64)
    set (PLATFORM "WIN64")
  endif ()
endif ()

if (NOT PLATFORM)
  if (DEFINED EXTRA_DIR)
    message (WARNING "Unsupported platform.")
  endif ()
else ()
  message ("-- PLATFORM: ${PLATFORM}")
endif ()

################################################################################
#                                                                              #
# Compiler-specific settings                                                   #
#                                                                              #
################################################################################

include (${PROJECT_SOURCE_DIR}/cmake/compilers.cmake)

# Generate the actual flags used in the build
#============================================

if (DMRG OR WFA)
  if (NOT CXXCID)
    message (FATAL_ERROR "unsupported CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
  endif ()
endif ()

if (NOT CCID)
  message (FATAL_ERROR "unsupported C compiler: ${CMAKE_C_COMPILER_ID}")
endif ()

if (NOT FCID)
  message (FATAL_ERROR "unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif ()

# flatten namespace on Mac OS X
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  set (CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace")
endif ()

# Exception trapping
if (CMAKE_BUILD_TYPE STREQUAL "Garble")
  add_compile_definitions (_GARBLE_ _FPE_TRAP_)
  # Disable FPE traps for compilers that do not support IEEE Fortran modules
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "5.0")
    remove_definitions (-D_FPE_TRAP_)
  elseif (NOT LINALG STREQUAL "Internal")
    message (WARNING
             "Garble is active, but an external LINALG library is used! "
             "The use of an external BLAS/LAPACK can result in floating point "
             "exceptions handled by the library being trapped instead. "
             "Please use -DLINALG=Internal for better chances for Garble to work correctly!")
    set (FFLAGS_${FCID}_GARBLE "${FFLAGS_${FCID}_GARBLE} ${${FCID}_TRAP}")
  else ()
    message (WARNING
             "Garble and floating point exception trapping are active! "
             "It is possible that intended exceptions within BLAS/LAPACK crash the program. "
             "Be warned!")
    set (FFLAGS_${FCID}_GARBLE "${FFLAGS_${FCID}_GARBLE} ${${FCID}_TRAP}")
  endif ()
endif ()

if (NOT CMAKE_GENERATOR STREQUAL "Ninja")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_PREPROCESS}")
endif ()

if (ADDRMODE EQUAL 64)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ILP64}")
endif ()

if (OPENMP)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_OPENMP}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_OPENMP}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP}")
else ()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_NO_OPENMP}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_NO_OPENMP}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_NO_OPENMP}")
endif ()

if (BOUNDS)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_BOUNDS}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_BOUNDS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_BOUNDS}")
  list (APPEND mma_util_defs "_MALLOC_INTERCEPT_")
endif ()

# Allow user-specified modifications
if (EXTRA_CXXFLAGS)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXXFLAGS}")
endif ()
if (EXTRA_CFLAGS)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_CFLAGS}")
endif ()
if (EXTRA_FFLAGS)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FFLAGS}")
endif ()

################################################################################
#                                                                              #
# Process options                                                              #
#                                                                              #
################################################################################

if ((BLOCK AND CHEMPS2) OR
    (BLOCK AND DMRG) OR
    (BLOCK AND DICE) OR
    (CHEMPS2 AND DMRG) OR
    (CHEMPS2 AND DICE) OR
    (DICE AND DMRG))
  message (FATAL_ERROR "BLOCK, CHEMPS2, DICE, and DMRG options are not compatible.")
endif ()

# Install directory
#==================
if (NOT IS_ABSOLUTE "${CMAKE_INSTALL_PREFIX}")
  set (CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}")
endif ()
get_filename_component (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)

# Parallel settings
#==================

if (MPI)
  message ("Configuring with MPI parallellization:")
  find_package (MPI REQUIRED)

  add_compile_definitions (_MOLCAS_MPP_)
  include_directories (${MPI_Fortran_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})

  message ("-- MPI_C_INCLUDE_PATH: ${MPI_C_INCLUDE_PATH}")
  message ("-- MPI_Fortran_INCLUDE_PATH: ${MPI_Fortran_INCLUDE_PATH}")
  message ("-- MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
  message ("-- MPI_Fortran_LIBRARIES: ${MPI_Fortran_LIBRARIES}")
  message ("-- MPIEXEC: ${MPIEXEC}")

  execute_process (
    COMMAND ${MPIEXEC} "--version"
    OUTPUT_VARIABLE MPIEXEC_VERSION
    ERROR_VARIABLE MPIEXEC_VERSION
  )
  if (MPIEXEC_VERSION MATCHES "OpenRTE")
    set (MPI_IMPLEMENTATION "openmpi")
  elseif (MPIEXEC_VERSION MATCHES "Intel\\(R\\) MPI")
    set (MPI_IMPLEMENTATION "impi")
  elseif (MPIEXEC_VERSION MATCHES "HYDRA")
    set (MPI_IMPLEMENTATION "mpich")
  else ()
    # check for MPICH 3, which has only "mpichversion" as an executable,
    # but "mpiexec --version" does not work
    execute_process (
    COMMAND "mpichversion"
    OUTPUT_VARIABLE MPIEXEC_VERSION
    ERROR_VARIABLE MPIEXEC_VERSION
    )
    if (MPIEXEC_VERSION MATCHES "MPICH")
      set (MPI_IMPLEMENTATION "mpich")
    else ()
      message (FATAL_ERROR "Unknown or unsupported MPI implementation:\n${MPIEXEC_VERSION}")
    endif ()
  endif ()
  message ("-- MPI_IMPLEMENTATION: ${MPI_IMPLEMENTATION}")

  # MPICH and GCC > 10.0 requires explicit MPI interface
  if (MPI_IMPLEMENTATION STREQUAL "mpich" AND
      CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER "10.0")
    add_compile_definitions (_NEED_EXPLICIT_MPI_INTERFACE_)
  endif ()

  list (APPEND EXTERNAL_LIBRARIES ${MPI_Fortran_LIBRARIES})
endif ()

# Support for HDF5 I/O library
#=============================

if (HDF5)
  if (ADDRMODE EQUAL 32)
    message ("-- HDF5 interface not compatible with 32bit installations")
  else ()
    message ("Configuring HDF5 support:")
    get_property (languages GLOBAL PROPERTY ENABLED_LANGUAGES)
    if (";${languages};" MATCHES ";CXX;")
      find_package (HDF5 COMPONENTS C CXX)
    else ()
      find_package (HDF5 COMPONENTS C)
    endif ()
    if (HDF5_FOUND)
      message ("-- HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
      message ("-- HDF5_C_LIBRARIES: ${HDF5_C_LIBRARIES}")
      if (HDF5_VERSION)
        if (DMRG AND HDF5_VERSION VERSION_EQUAL "1.10.0.1")
          message (FATAL_ERROR "This HDF5 version (${HDF5_VERSION}) is not supported by QCMaquis.")
        endif ()
      endif ()
      list (APPEND EXTERNAL_LIBRARIES ${HDF5_C_LIBRARIES})
      if (HDF5_CXX_LIBRARIES)
        message ("-- HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}")
        list (APPEND EXTERNAL_LIBRARIES ${HDF5_CXX_LIBRARIES})
      endif ()
    else ()
      message ("-- HDF5 not found, it will be deactivated")
    endif ()
  endif ()
endif ()

# BLAS/LAPACK settings
#=====================

message ("Configuring linear algebra libraries:")

# If one uses BLAS/LAPACK, you can use the CMake FindBLAS/FindLAPACK modules to
# set everything, which is determined by setting the BLA_VENDOR environment
# variable. Unfortunately, currently, the modules that comes with CMake doesn't
# seem to support any ilp64 version of MKL, as it only has e.g. Intel10_64lp
# For now, I'm using some non-portable thing.

set (MKLROOT "" CACHE PATH "MKL root directory.")
if (LINALG STREQUAL "MKL")
  message ("-- Using Intel Math Kernel Library (MKL)")

  # quick bypass for 32bit, don't bother with it (for now?)
  if (NOT ADDRMODE EQUAL 64)
    message (FATAL_ERROR "Use of MKL for 32bit installation not implemented.")
  endif ()

  mark_as_advanced (CLEAR MKLROOT)
  if (MKLROOT STREQUAL "")
    set (MKLROOT $ENV{MKLROOT} CACHE PATH "MKL root directory." FORCE)
    if (NOT MKLROOT)
      message (FATAL_ERROR
               "You must set environment variable MKLROOT, "
               "or specify -DMKLROOT=/path/to/mkl_root_dir "
               "when running cmake.")
    endif ()
  endif ()
  # at this point, MKLROOT should be defined and not empty
  message ("-- MKLROOT = ${MKLROOT}")

  # here we check if MKLROOT has changed after it was last used
  # succesfully, if so then we need to reprocess the checks here.
  if (NOT MKLROOT STREQUAL MKLROOT_LAST)
    # reset MKL paths
    set (MKL_INCLUDE_PATH "${MKLROOT}/include" CACHE PATH "location of MKL include files." FORCE)

    if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
      if (ADDRMODE EQUAL 64)
        set (libpath "${MKLROOT}/lib")
      elseif (ADDRMODE EQUAL 32)
        message (FATAL_ERROR
                 "PLATFORM=MacOS, but you use ADDRMODE=32. This is "
                 "either a VERY old, unsupported platform, or "
                 "there is a configuration error.")
      endif ()
    else ()
      if (ADDRMODE EQUAL 64)
        set (libpath "${MKLROOT}/lib/intel64")
      elseif (ADDRMODE EQUAL 32)
        set (libpath "${MKLROOT}/lib/ia32")
      endif ()
    endif ()
    set (MKL_LIBRARY_PATH ${libpath} CACHE PATH "location of MKL libraries." FORCE)
    # uncache variables
    unset (LIBMKL_BLAS CACHE)
    unset (LIBMKL_LAPACK CACHE)
    unset (LIBMKL_CORE CACHE)
    unset (LIBMKL_INTERFACE CACHE)
    unset (LIBMKL_SEQUENTIAL CACHE)
    unset (LIBMKL_THREADING CACHE)
    unset (LIBMKL_SCALAPACK CACHE)
    unset (LIBMKL_BLACS CACHE)
    # cache last used MKLROOT
    set (MKLROOT_LAST ${MKLROOT} CACHE INTERNAL "last value." FORCE)
  endif ()
  message ("-- MKL_INCLUDE_PATH = ${MKL_INCLUDE_PATH}")
  message ("-- MKL_LIBRARY_PATH = ${MKL_LIBRARY_PATH}")

  # core library
  if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    find_library (LIBMKL_CORE NAMES "libmkl_core.a"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  else ()
    find_library (LIBMKL_CORE NAMES "mkl_core"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()
  # compiler-specific library interface
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
      find_library (LIBMKL_BLAS NAMES "libmkl_blas95_ilp64.a"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
      find_library (LIBMKL_LAPACK NAMES "libmkl_lapack95_ilp64.a"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
      find_library (LIBMKL_INTERFACE NAMES "libmkl_intel_ilp64.a"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
    else ()
    find_library (LIBMKL_INTERFACE NAMES "mkl_gf_ilp64"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
    endif ()
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI") # apparently PGI uses this too
    find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()
  # sequential/compiler-specific threading interface
  if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    find_library (LIBMKL_SEQUENTIAL NAMES "libmkl_sequential.a"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  else ()
    find_library (LIBMKL_SEQUENTIAL NAMES "mkl_sequential"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set (LIBMKL_OMP_LINK_FLAGS "${FFLAGS_GNU_OPENMP}")
    find_library (LIBMKL_THREADING NAMES "mkl_gnu_thread"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set (LIBMKL_OMP_LINK_FLAGS "${FFLAGS_Intel_OPENMP}")
    find_library (LIBMKL_THREADING NAMES "mkl_intel_thread"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()
  # find scalapack/blacs for parallel lapack support
  find_library (LIBMKL_SCALAPACK NAMES "mkl_scalapack_ilp64"
                PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  if (MPI_IMPLEMENTATION STREQUAL "openmpi")
    find_library (LIBMKL_BLACS NAMES "mkl_blacs_openmpi_ilp64"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (MPI_IMPLEMENTATION STREQUAL "impi")
    find_library (LIBMKL_BLACS NAMES "mkl_blacs_intelmpi_ilp64"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (MPI_IMPLEMENTATION STREQUAL "mpich")
    find_library (LIBMKL_BLACS NAMES "mkl_blacs_intelmpi_ilp64"
                  PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()

  # generate actual library list with paths
  if (MPI AND LIBMKL_SCALAPACK)
    add_compile_definitions (_SCALAPACK_)
    list (APPEND MKL_LIBRARIES ${LIBMKL_SCALAPACK})
  endif ()
  if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    list (APPEND MKL_LIBRARIES ${LIBMKL_BLAS})
    list (APPEND MKL_LIBRARIES ${LIBMKL_LAPACK})
  endif ()
  list (APPEND MKL_LIBRARIES ${LIBMKL_INTERFACE})
  list (APPEND MKL_LIBRARIES ${LIBMKL_CORE})
  if (OPENMP)
    set (CMAKE_EXE_LINKER_FLAGS ${LIBMKL_OMP_LINK_FLAGS})
    list (APPEND MKL_LIBRARIES ${LIBMKL_THREADING})
  else ()
    list (APPEND MKL_LIBRARIES ${LIBMKL_SEQUENTIAL})
  endif ()
  if (MPI AND LIBMKL_BLACS)
    list (APPEND MKL_LIBRARIES ${LIBMKL_BLACS})
  endif ()

  list (APPEND linalg_util_defs "_MKL_")
  set (LINALG_LIBRARIES ${MKL_LIBRARIES})
else ()
  mark_as_advanced (FORCE MKLROOT)
endif ()

# AOCL settings, I don't think FindBLAS supports recent AOCL either
set (AOCLROOT "" CACHE PATH "AOCL path.")
if (LINALG STREQUAL "AOCL")
  message ("-- Using AMD Optimizing CPU Libraries (AOCL)")

  mark_as_advanced (CLEAR AOCLROOT)
  if (AOCLROOT STREQUAL "")
    set (AOCLROOT $ENV{AOCLROOT} CACHE PATH "AOCL root directory." FORCE)
    if (NOT AOCLROOT)
      message (FATAL_ERROR
               "You must set environment variable AOCLROOT, "
               "or specify -DAOCLROOT=/path/to/AOCL_dir "
               "when running cmake.")
    endif ()
  endif ()

  # Search libflame first
  find_library (LIBFLAME
    NAMES flame
    PATHS ${AOCLROOT}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
  )

  if (NOT LIBFLAME)
    message (FATAL_ERROR "LIBFLAME not found" )
  endif ()

  # Search libblis
  find_library (LIBBLIS
    NAMES blis
    PATHS ${AOCLROOT}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
  )

  if (NOT LIBBLIS)
    message (FATAL_ERROR "LIBBLIS not found" )
  endif ()

  unset (LIBAOCL CACHE)
  set (LIBAOCL ${LIBFLAME} ${LIBBLIS})

  if (NOT LIBAOCL)
    message (FATAL_ERROR
             "AOCL library not found, please check that "
             "the AOCLROOT variable is set and points to "
             "a valid AOCL installation directory.")
  endif ()

  list (APPEND linalg_util_defs "_AOCL_")
  set (LINALG_LIBRARIES ${LIBAOCL})
else ()
  mark_as_advanced (FORCE AOCLROOT)
endif ()

# OpenBLAS
set (OPENBLASROOT "" CACHE PATH "OpenBLAS root directory.")
if (LINALG STREQUAL "OpenBLAS")
  message ("-- Using OpenBLAS LAPACK+BLAS libraries")

  # find out what OPENBLASROOT is
  mark_as_advanced (CLEAR OPENBLASROOT)
  if (OPENBLASROOT STREQUAL "")
    set (OPENBLASROOT $ENV{OPENBLASROOT} CACHE PATH "OpenBLAS root directory." FORCE)
    if (NOT OPENBLASROOT)
      message (FATAL_ERROR
               "You must set environment variable OPENBLASROOT, "
               "or specify -DOPENBLASROOT=/path/to/openblas_dir "
               "when running cmake.")
    endif ()
  endif ()
  # at this point, OPENBLASROOT should be defined and not empty
  message ("-- OPENBLASROOT = ${OPENBLASROOT}")

  # here we check if OPENBLASROOT has changed after it was last used
  # succesfully, if so then we need to reprocess the checks here.
  if (NOT OPENBLASROOT STREQUAL OPENBLASROOT_LAST)
    # we need to redo the library search later, because OPENBLASROOT
    # changed since the last time we defined it.
    unset (LIBOPENBLAS CACHE)
    # check if the configuration header is available
    set (OPENBLAS_CONFIG "${OPENBLASROOT}/include/openblas_config.h")
    if (NOT EXISTS ${OPENBLAS_CONFIG})
      # for system-wide OpenBLAS installations the config path might be different
      # so try with an alternative path
      set (OPENBLAS_CONFIG "${OPENBLASROOT}/include/openblas/openblas_config.h")
    endif ()
    if (NOT EXISTS ${OPENBLAS_CONFIG})
      message (FATAL_ERROR
               "Could not find OpenBLAS config header in: ${OPENBLASROOT} "
               "(tried ${OPENBLASROOT}/include/openblas_config.h and "
               "please check if the OPENBLASROOT variable points to a "
               "valid OpenBLAS installation directory.")
    endif ()
    # check if the OpenBLAS installation was configured with 64bit integer support
    message ("-- Checking OpenBLAS for 64-bit integer interface...")
    include (CheckSymbolExists)
    unset (OPENBLAS_WITH_ILP64 CACHE)
    check_symbol_exists ("OPENBLAS_USE64BITINT" ${OPENBLAS_CONFIG} OPENBLAS_WITH_ILP64)
    if (ADDRMODE EQUAL 64 AND NOT OPENBLAS_WITH_ILP64)
      message (FATAL_ERROR
               "OpenBLAS was not configured for 64-bit integer interface, "
               "please build OpenBLAS with INTERFACE64=1 defined.")
    endif ()
    # save the last location to check if it changed between configurations
    set (OPENBLASROOT_LAST ${OPENBLASROOT} CACHE INTERNAL "last value." FORCE)
  endif ()

  # search for the OpenBLAS library
  find_library (LIBOPENBLAS
    NAMES openblas
    PATHS ${OPENBLASROOT}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
  )

  if (NOT LIBOPENBLAS)
    message (FATAL_ERROR
             "OpenBLAS library not found, please check that "
             "the OPENBLASROOT variable is set and points to "
             "a valid OpenBLAS installation directory.")
  endif ()

  # here we check if LIBOPENBLAS has changed after it was processed
  # succesfully, if not we do not need to rerun anything here.
  if (NOT LIBOPENBLAS STREQUAL LIBOPENBLAS_LAST)
    # check if the OpenBLAS library contains LAPACK functionality
    message ("-- Checking OpenBLAS for LAPACK functionality...")
    include (CheckFortranFunctionExists)
    set (CMAKE_REQUIRED_LIBRARIES ${LIBOPENBLAS})
    unset (OPENBLAS_WITH_LAPACK CACHE)
    check_fortran_function_exists ("dsyev" OPENBLAS_WITH_LAPACK)
    if (NOT OPENBLAS_WITH_LAPACK)
      unset (LIBOPENBLAS CACHE)
      message (FATAL_ERROR
               "LAPACK functionality missing from OpenBLAS library, "
               "please build OpenBLAS with NO_LAPACK=0 defined.")
    endif ()
    # save the last location to check if it changed between configurations
    set (LIBOPENBLAS_LAST ${LIBOPENBLAS} CACHE INTERNAL "last value." FORCE)
  endif ()

  set (LINALG_LIBRARIES ${LIBOPENBLAS})
else ()
  mark_as_advanced (FORCE OPENBLASROOT)
endif ()

if (LINALG STREQUAL "Accelerate")
  message ("-- Using XCode's veclib libraries")
  set (LINALG_LIBRARIES "-framework Accelerate")
endif ()

# Internal or Runtime: a git submodule is required
if (LINALG STREQUAL "Internal" OR LINALG STREQUAL "Runtime")
  set (proj ${extprojsdir}/lapack)
  if (EXISTS ${proj}/BLAS/SRC/ddot.f)
    set (internal_linalg ${proj}/SRC ${proj}/INSTALL ${proj}/BLAS/SRC)
  else ()
    # Note: initial space causes line not to break in output
    message (FATAL_ERROR " LAPACK+BLAS sources not available, run \"${GIT} submodule update --init ${proj}\"")
  endif ()
endif ()

# Internal LAPACK+BLAS libraries from source code
set (INTERNAL_LINALG "lapack" "blas")
find_source (LinAlg_internal ${source_roots})
if (LINALG STREQUAL "Internal")
  message ("-- Using internal LAPACK+BLAS libraries (SLOW!)")
  set (LINALG_LIBRARIES ${INTERNAL_LINALG})
  # create the BLAS+LAPACK targets
  set (util "LinAlg_internal")
  if (NOT DEFINED ${util}_src)
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
  set (linalgmod "${${util}_src}/la_constants_.F90" "${${util}_src}/la_xisnan_.F90")
  foreach (lib ${LINALG_LIBRARIES})
    set (source ${${util}_src}/${lib}.f ${${util}_src}/${lib}_f90.F90 ${linalgmod})
    add_library (${lib} ${source})
    list (APPEND ALL_SOURCES ${source})
    list (APPEND ${util}_incs ${internal_linalg})
    foreach (inc ${${util}_incs})
      target_include_directories (${lib} PRIVATE ${inc})
    endforeach ()
    if (BUILD_SHARED_LIBS)
      set_target_properties (${lib} PROPERTIES POSITION_INDEPENDENT_CODE "TRUE")
    endif ()
  endforeach ()
endif ()

# Runtime loading of LAPACK+BLAS libraries (experimental)
if (LINALG STREQUAL "Runtime")
  message ("-- Using Runtime LAPACK+BLAS libraries (Internal as fallback)")
  set (LINALG_LIBRARIES "dl")

  list (APPEND system_util_defs "_DELAYED_")
endif ()

if (LINALG STREQUAL "Manual")
  message ("-- Manual linear algebra libraries, using LINALG_LIBRARIES as specified")
  string (REPLACE ":" ";" LINALG_LIBRARIES "${LINALG_LIBRARIES}")
  set (LINALG_LIBRARIES "${LINALG_LIBRARIES}" CACHE STRING "Manual linear algebra libraries." FORCE)
else ()
  unset (LINALG_LIBRARIES CACHE)
endif ()

list (APPEND EXTERNAL_LIBRARIES ${LINALG_LIBRARIES})

message ("-- LINALG_LIBRARIES: ${LINALG_LIBRARIES}")

# GPU acceleration
#=================

set (CUDAROOT "" CACHE PATH "CUDA root directory.")
if (CUBLAS OR NVBLAS)
  message ("Configuring with CUDA BLAS library:")

  mark_as_advanced (CLEAR CUDAROOT)
  if (CUDAROOT STREQUAL "")
    set (CUDAROOT $ENV{CUDAROOT} CACHE PATH "CUDA root directory." FORCE)
    if (NOT CUDAROOT)
      message (FATAL_ERROR "CUDA requested, but environment variable CUDAROOT not set.")
    endif ()
  endif ()

  # the CUDA library itself is installed with the driver in its own
  # location, so we need to find it first
  find_library (LIBCUDA NAMES "cuda")
  list (APPEND CUDA_LIBRARIES ${LIBCUDA})

  if (CUBLAS)
    # compile Fortran wrappers as a library
    set (CUDA_INCLUDE_PATH "${CUDAROOT}/include")
    set (CUDA_FORTRAN_WRAPPERS "${CUDAROOT}/src/fortran_thunking.c")
    if (NOT EXISTS ${CUDA_FORTRAN_WRAPPERS})
      message (FATAL_ERROR "could not find CUDA fortran wrappers: ${CUDA_FORTRAN_WRAPPERS}")
    endif ()
    add_library (cublas_fortran OBJECT "${CUDA_FORTRAN_WRAPPERS}")
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_definitions (cublas_fortran PRIVATE "CUBLAS_GFORTRAN")
    endif ()
    target_include_directories (cublas_fortran PRIVATE ${CUDA_INCLUDE_PATH})
    list (APPEND utils_obj cublas_fortran)

    # CUDA needs to link to the C++ standard library, so enable CXX here to
    # be able to use it as a linker later on for the program executables
    enable_language (CXX)

    add_compile_definitions (_CUDA_BLAS_)

    find_library (LIBCUBLAS cublas
      PATHS ${CUDAROOT}
      PATH_SUFFIXES lib64
      NO_DEFAULT_PATH
    )

    list (APPEND CUDA_LIBRARIES ${LIBCUBLAS})
  else ()
    if (NOT LINALG)
      message (FATAL_ERROR "NVBLAS needs an external linear algebra library such as MKL or AOCL.")
    endif ()
    find_library (LIBNVBLAS nvblas
      PATHS ${CUDAROOT}
      PATH_SUFFIXES lib64
      NO_DEFAULT_PATH
    )
    find_library (LIBCUBLAS cublas
      PATHS ${CUDAROOT}
      PATH_SUFFIXES lib64
      NO_DEFAULT_PATH
    )

    list (APPEND CUDA_LIBRARIES ${LIBNVBLAS} ${LIBCUBLAS})
  endif ()

  message ("-- CUDA_LIBRARIES: ${CUDA_LIBRARIES}")

  list (APPEND EXTERNAL_LIBRARIES ${CUDA_LIBRARIES})
else ()
  mark_as_advanced (FORCE CUDAROOT)
endif ()

# Libxc settings
#===============
message ("Configuring Libxc:")

if (EXTERNAL_LIBXC)
  message ("-- External Libxc root specified at: ${EXTERNAL_LIBXC}")
  find_library (xc xc PATHS ${EXTERNAL_LIBXC} PATH_SUFFIXES lib lib32 lib64)
  find_library (xcf03 xcf03 PATHS ${EXTERNAL_LIBXC} PATH_SUFFIXES lib lib32 lib64)
  if (NOT xc OR NOT xcf03)
    message (FATAL_ERROR "Libxc could not be found in ${EXTERNAL_LIBXC}.")
  endif ()
  # find Libxc version
  file (READ ${EXTERNAL_LIBXC}/include/xc_version.h header)
  string (REGEX MATCH "#define XC_VERSION \"[0-9.]+\"" macrodef "${header}")
  string (REGEX MATCH "[0-9.]+" LIBXC_VERSION "${macrodef}")
  if (LIBXC_VERSION VERSION_LESS "5.2.0")
    message (FATAL_ERROR "Libxc version >= 5.2.0 is required (found ${LIBXC_VERSION})")
  endif ()
  set (Libxc_INCLUDE_DIRS ${EXTERNAL_LIBXC}/include)
  set (Libxc_LIBRARIES ${xcf03} ${xc})
else ()
  message ("-- Configuring built-in Libxc")
  include (ExternalProject)
  set (EP_PROJECT Libxc)
  set (EP_PATH ${PROJECT_BINARY_DIR}/External/${EP_PROJECT})
  set (EP_INSTALL ${PROJECT_BINARY_DIR}/External/${EP_PROJECT}_install)

  set (Libxc_REPO https://gitlab.com/libxc/libxc.git)
  set (Libxc_HASH 0a19bedb) # Release 6.2.2

  if (";Debug;Release;RelWithDebInfo;" MATCHES ";${CMAKE_BUILD_TYPE};")
    set (LIBXC_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  elseif (";Garble;" MATCHES ";${CMAKE_BUILD_TYPE};")
    set (LIBXC_BUILD_TYPE "Debug")
  else ()
    set (LIBXC_BUILD_TYPE "Release")
  endif ()

  # Note shared libraries won't work here because of installation issues.
  # TODO: Solving this will require creating a wrapper main project that includes
  #       both Libxc and OpenMolcas as ExternalProject. This could also simplify
  #       the handling of other codes like GlobalArrays or QCMaquis.
  list (APPEND LibxcCMakeArgs
        -D CMAKE_BUILD_TYPE=${LIBXC_BUILD_TYPE}
        -D CMAKE_INSTALL_PREFIX=${EP_INSTALL}
        -D BUILD_SHARED_LIBS=NO
        -D CMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -D CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -D CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS_DEFAULT}
        -D CMAKE_INSTALL_LIBDIR=lib
        -D ENABLE_FORTRAN=ON
        -D ENABLE_PYTHON=OFF
        -D DISABLE_FXC=ON
        -D DISABLE_KXC=ON
        -D DISABLE_LXC=ON
        -D BUILD_TESTING=NO
  )
  # "Hidden" option to pass extra arguments to Libxc CMake
  list (APPEND LibxcCMakeArgs ${Libxc_EXTRA_CMakeArgs})

  # Workaround for ExternalProject's shortcoming
  # (see https://gitlab.kitware.com/cmake/cmake/-/issues/16419)
  # UPDATE_DISCONNECTED is ON if the current commit has already
  # been built, and OFF when the commit changes, but it requires
  # to run cmake for detection
  set (last_hash "None")
  set (hash_file ${EP_PATH}/${EP_PROJECT}.hash)
  if (EXISTS ${hash_file})
    file (READ ${hash_file} last_hash)
    string (REGEX REPLACE "\n$" "" last_hash "${last_hash}")
  endif ()
  if (last_hash STREQUAL ${Libxc_HASH})
    set (LibxcSkipUpdate ON)
  else ()
    set (LibxcSkipUpdate OFF)
  endif ()

  ExternalProject_Add (${EP_PROJECT}
                       PREFIX ${EP_PATH}
                       GIT_REPOSITORY ${Libxc_REPO}
                       GIT_TAG ${Libxc_HASH}
                       CMAKE_ARGS ${LibxcCMakeArgs}
                       UPDATE_DISCONNECTED ${LibxcSkipUpdate}
  )

  ExternalProject_Add_Step (${EP_PROJECT} update_hash
                            COMMAND echo ${Libxc_HASH} > ${hash_file}
                            DEPENDEES build
  )

  link_directories(${EP_INSTALL}/lib)

  set (Libxc_INCLUDE_DIRS ${EP_INSTALL}/include)
  set (Libxc_LIBRARIES xcf03 xc)

endif ()
message ("-- Libxc_INCLUDE_DIRS: ${Libxc_INCLUDE_DIRS}")
message ("-- Libxc_LIBRARIES: ${Libxc_LIBRARIES}")
list (APPEND EXTERNAL_LIBRARIES ${Libxc_LIBRARIES})

# Some compilers seem to have issues with nested modules
set (_comp_list PGI NAG)
if (${FCID} IN_LIST _comp_list)
  include_directories (${Libxc_INCLUDE_DIRS})
else ()
  foreach (prog
    dft_util
    nq_util
  )
    list (APPEND ${prog}_incs ${Libxc_INCLUDE_DIRS})
  endforeach ()
endif ()

# Global Arrays settings
#=======================

if (GA AND MPI)
  # Build from source
  if (GA_BUILD)
    message ("Configuring built-in Global Arrays library:")

    if (LINALG STREQUAL "Runtime")
      message (FATAL_ERROR "GA_BUILD does not support LINALG=Runtime.")
    endif ()

    include (ExternalProject)
    set (EP_PROJECT GlobalArrays)
    set (EP_PATH ${PROJECT_BINARY_DIR}/External/${EP_PROJECT})
    set (EP_INSTALL ${PROJECT_BINARY_DIR}/External/${EP_PROJECT}_install)

    set (GA_REPO https://github.com/GlobalArrays/ga.git)
    set (GA_HASH da7a53ff)

    target_files (LINALG_LIBRARIES_FILES ${LINALG_LIBRARIES})
    string (REPLACE ";" "|" LINALG_LIBRARIES_FILES "${LINALG_LIBRARIES_FILES}")
    list (APPEND GACMakeArgs
          -D CMAKE_INSTALL_PREFIX=${EP_INSTALL}
          -D CMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE
          -D BUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
          -D ENABLE_FORTRAN=ON
          -D ENABLE_CXX=OFF
          -D ENABLE_TESTS=OFF
          -D ENABLE_BLAS=ON
          -D LINALG_PREFIX=${PROJECT_BINARY_DIR}
          -D BLAS_LIBRARIES=${LINALG_LIBRARIES_FILES}
    )
    if (ADDRMODE EQUAL 64)
      list (APPEND GACMakeArgs
            -D ENABLE_I8=ON
            -D LINALG_REQUIRED_COMPONENTS=ilp64
      )
    else ()
      list (APPEND GACMakeArgs
            -D ENABLE_I8=OFF
            -D LINALG_REQUIRED_COMPONENTS=lp64
      )
    endif ()

    # Same workaround as for Libxc (see above)
    set (last_hash "None")
    set (hash_file ${EP_PATH}/${EP_PROJECT}.hash)
    if (EXISTS ${hash_file})
      file (READ ${hash_file} last_hash)
      string (REGEX REPLACE "\n$" "" last_hash "${last_hash}")
    endif ()
    if (last_hash STREQUAL ${GA_HASH})
      set (GASkipUpdate ON)
    else ()
      set (GASkipUpdate OFF)
    endif ()

    ExternalProject_Add (${EP_PROJECT}
                         PREFIX ${EP_PATH}
                         GIT_REPOSITORY ${GA_REPO}
                         GIT_TAG ${GA_HASH}
                         CMAKE_ARGS ${GACMakeArgs}
                         UPDATE_DISCONNECTED ${GASkipUpdate}
                         LIST_SEPARATOR "|"
    )

    if (NOT GASkipUpdate)
      ExternalProject_add_patches (${EP_PROJECT} "cmake/custom/ga.patch")
    endif ()

    ExternalProject_Add_Step (${EP_PROJECT} update_hash
                              COMMAND echo ${GA_HASH} > ${hash_file}
                              DEPENDEES build
    )

    if (LINALG STREQUAL "Internal")
      add_dependencies (${EP_PROJECT} ${LINALG_LIBRARIES})
    endif ()

    set (GA_INCLUDE_PATH ${EP_INSTALL}/include/ga)
    set (GA_LIBRARIES ga armci)
    # This is not a very portable way of finding and installing
    # the library files, but I can't find anything better
    if (BUILD_SHARED_LIBS)
      set (lib_ext ".so")
    else ()
      set (lib_ext ".a")
    endif ()
    foreach (lib ${GA_LIBRARIES})
      set (${lib}_file ${EP_INSTALL}/lib/lib${lib}${lib_ext})
      add_library (${lib} UNKNOWN IMPORTED)
      set_target_properties (${lib} PROPERTIES IMPORTED_LOCATION ${${lib}_file})
    endforeach ()
    list (APPEND EXTRA_INSTALL_IMPORTED_LIBS ${GA_LIBRARIES})

  # External precompiled GA library
  else ()
    message ("Configuring external Global Arrays library:")

    find_path (GA_INCLUDE_PATH ga.h
      PATHS ENV GAROOT
      PATH_SUFFIXES include include/ga
      NO_DEFAULT_PATH
    )
    if (NOT GA_INCLUDE_PATH)
      message (FATAL_ERROR "GA_INCLUDE_PATH not found, make sure GAROOT environment variable is set.")
    endif ()
    find_library (LIBGA ga
      PATHS ENV GAROOT
      PATH_SUFFIXES lib
      NO_DEFAULT_PATH
    )
    if (NOT LIBGA)
      message (FATAL_ERROR "LIBGA not found, make sure GAROOT environment variable is set.")
    else ()
      add_library (ga UNKNOWN IMPORTED)
      set_target_properties (ga PROPERTIES IMPORTED_LOCATION ${LIBGA})
      list (APPEND GA_LIBRARIES ga)
    endif ()
    find_library (LIBARMCI armci
      PATHS ENV GAROOT
      PATH_SUFFIXES lib
      NO_DEFAULT_PATH
    )
    if (NOT LIBARMCI)
      message (FATAL_ERROR "LIBARMCI not found, make sure GAROOT environment variable is set.")
    else ()
      add_library (armci UNKNOWN IMPORTED)
      set_target_properties (armci PROPERTIES IMPORTED_LOCATION ${LIBARMCI})
      list (APPEND GA_LIBRARIES armci)
    endif ()
    # Hackish: try to find out if ARMCI needs IBVERBS by looking for ibv_open_device
    execute_process (
      COMMAND "nm" "${LIBARMCI}"
      COMMAND "grep" "-iq" "ibv_open_device"
      RESULT_VARIABLE FAILURE
    )
    if (FAILURE)
      message ("-> no libibverbs dependency found for GA.")
    else ()
      message ("-> detected ibverbs dependency for GA.")
      find_library (LIBIBVERBS ibverbs)
      if (NOT LIBARMCI)
        message (WARNING "LIBIBVERBS not found, but probably needed by GA.")
      else ()
        target_link_libraries (armci INTERFACE ${LIBIBVERBS})
      endif ()
    endif ()
    set (GA_FOUND "TRUE")
  endif ()

  # Final common settings

  foreach (lib ${GA_LIBRARIES})
    target_link_libraries (${lib} INTERFACE ${MPI_LIBRARIES} ${LINALG_LIBRARIES})
  endforeach ()

  add_compile_definitions (_GA_)

  message ("-- GA_INCLUDE_PATH: ${GA_INCLUDE_PATH}")
  message ("-- GA_LIBRARIES: ${GA_LIBRARIES}")

  list (APPEND EXTERNAL_LIBRARIES ${GA_LIBRARIES})
endif ()

# Gromacs
#========

set (GROMACS_DIR "" CACHE PATH "Gromacs installation path.")
if (GROMACS)
  mark_as_advanced (CLEAR GROMACS_DIR)

  message ("Configuring Gromacs:")

  # Try to auto-detect the location of Gromacs
  if (GROMACS_DIR STREQUAL "")
    if (DEFINED ENV{GMXLDLIB})
      get_filename_component (GMX_DIR $ENV{GMXLDLIB} PATH)
      get_filename_component (GMX_DIR ${GMX_DIR} PATH)
      set (GROMACS_DIR ${GMX_DIR} CACHE PATH "Gromacs installation path." FORCE)
      unset (GMX_DIR)
    else ()
      message (FATAL_ERROR "You must supply a path for GROMACS_DIR.")
    endif ()
  endif ()

  # The lib directory is set by Gromacs
  if (DEFINED ENV{GMXLDLIB} AND NOT GROMACS_DIR)
    set (GMX_LIBDIR $ENV{GMXLDLIB})
  else ()
    set (GMX_LIBDIR "${GROMACS_DIR}/lib")
  endif ()

  find_library (LIBGMX
    NAMES gromacs_d gromacs
    PATHS ${GMX_LIBDIR}
    NO_DEFAULT_PATH
  )
  if (NOT LIBGMX)
    message (FATAL_ERROR "LIBGMX not found, make sure the GROMACS_DIR variable is set correctly.")
  else ()
    list (APPEND GMX_LIBRARIES ${LIBGMX})
  endif ()
  unset (LIBGMX CACHE)

  # For the include directory we search the header file
  find_path (INCGMX
    NAMES gromacs/mmslave.h
    PATHS ${GROMACS_DIR}/include
    NO_DEFAULT_PATH
  )
  if (NOT INCGMX)
    message (FATAL_ERROR "INCGMX not found, make sure the GROMACS_DIR variable is set correctly.")
  else ()
    list (APPEND GMX_INCLUDE_PATH ${INCGMX})
  endif ()
  unset (INCGMX CACHE)

  foreach (prog
    gateway_util
    espf_util
  )
    list (APPEND ${prog}_incs ${GMX_INCLUDE_PATH})
    list (APPEND ${prog}_defs "_GROMACS_")
  endforeach ()

  message ("-- GMX_LIBRARIES: ${GMX_LIBRARIES}")
  message ("-- GMX_INCLUDE_PATH: ${GMX_INCLUDE_PATH}")

  list (APPEND EXTERNAL_LIBRARIES ${GMX_LIBRARIES})
else ()
  mark_as_advanced (FORCE GROMACS_DIR)
  message ("Gromacs interface DISABLED")
endif ()

# G. Chan group's DMRG (Block code) support
#==========================================

if (BLOCK)
  mark_as_advanced (CLEAR BLOCK_DIR)

  message ("Configuring Block-DMRG support:")

  # Since CMake build is not available in current Block code,
  # suppose Block has been installed at ${BLOCK_DIR} for the time.
  set (BLOCK_DIR ${BLOCK_DIR} CACHE PATH "Block installation path." FORCE)
  if ("${BLOCK_DIR}" STREQUAL "")
    message (FATAL_ERROR "You must supply a path for BLOCK_DIR.")
  endif ()

  if (EXISTS ${BLOCK_DIR}/libqcdmrg.so)
    foreach (prog
      wfn_util
      caspt2
      rasscf
    )
      list (APPEND ${prog}_defs "_ENABLE_BLOCK_DMRG_")
    endforeach ()

    link_directories (${BLOCK_DIR})
    add_library (block_dmrg_lib_name SHARED IMPORTED)
    set_target_properties (block_dmrg_lib_name PROPERTIES IMPORTED_LOCATION ${BLOCK_DIR}/libqcdmrg.so)

    message ("-- BLOCK_DIR: ${BLOCK_DIR}")

    list (APPEND EXTERNAL_LIBRARIES block_dmrg_lib_name)
  else ()
    set (BLOCK OFF)
    message ("-- BLOCK interface ignored (no BLOCK installation found at ${BLOCK_DIR})")
    message ("BLOCK interface DISABLED")
  endif ()
else ()
  mark_as_advanced (FORCE BLOCK_DIR)
  message ("BLOCK interface DISABLED")
endif ()

# S. Wouters' DMRG (CheMPS2 code) support
#========================================

if (CHEMPS2)
  mark_as_advanced (CLEAR CHEMPS2_DIR)

  message ("Configuring CheMPS2-DMRG support:")

  if (NOT HDF5)
    message (FATAL_ERROR "CheMPS2 interface requires enabling HDF5.")
  endif ()

  find_program (CHEMPS2_EXE chemps2)
  if (NOT CHEMPS2_EXE)
    message (WARNING "No CheMPS2 installation found, you will have to make sure it is available at run time.")
  endif ()

  foreach (prog
    wfn_util
    caspt2
    rasscf
  )
    list (APPEND ${prog}_defs "_ENABLE_CHEMPS2_DMRG_")
  endforeach ()
  message ("-- CHEMPS2_EXE: ${CHEMPS2_EXE}")
  unset (CHEMPS2_EXE CACHE)
else ()
  mark_as_advanced (FORCE CHEMPS2_DIR)
  message ("CHEMPS2 interface DISABLED")
endif ()

# Dice-SHCI support
#==================

if (DICE)
  mark_as_advanced (CLEAR DICE_DIR)

  message ("Configuring Dice-SHCI support:")

  find_program (DICE_EXE Dice)
  if (NOT DICE_EXE)
    message (WARNING "No Dice installation found, you will have to make sure it is available at run time.")
  endif ()

  foreach (prog
    wfn_util
    rasscf
  )
    list (APPEND ${prog}_defs "_ENABLE_DICE_SHCI_")
  endforeach ()
  message ("-- DICE_EXE: ${DICE_EXE}")
  unset (DICE_EXE CACHE)
else ()
  mark_as_advanced (FORCE DICE_DIR)
  message ("DICE interface DISABLED")
endif ()

# GPerfTools: CPU profiling
#==========================

if (GPERFTOOLS)
  message ("Configuring GPerfTools:")
  find_library (LIBPROFILER NAMES "profiler" "libprofiler.so.0")
  if (NOT LIBPROFILER)
    message (FATAL_ERROR "gperftools CPU profiler library not found.")
  endif ()
  message ("-- LIBPROFILER: ${LIBPROFILER}")

  list (APPEND EXTERNAL_LIBRARIES ${LIBPROFILER})
endif ()

# GCov: code coverage
#====================

if (GCOV)
  message ("Configuring for code coverage:")
  if (NOT ${CXXCID} STREQUAL "GNU" OR NOT ${CCID} STREQUAL "GNU" OR NOT ${FCID} STREQUAL "GNU")
    message (WARNING "gcov code coverage only supported with GCC.")
  else ()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage")
  endif ()
endif ()

################################################################################
#                                                                              #
# External programs/libraries: compilation does not require a molcas library,  #
# but needs to inherit the settings from our configuration.                    #
#                                                                              #
################################################################################

# global module directory
#========================

set (MAIN_MOD_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)
# The Intel compiler can generate automatic interfaces, in that case we want to keep them all in a single place
if (BIGOT AND CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set (SINGLE_MOD_DIR ON)
  set (MAIN_MOD_DIR ${MAIN_MOD_DIR}/_single)
endif ()
include_directories (${MAIN_MOD_DIR})

if (MSYM)
  message ("Configuring MSYM support:")

  set (proj ${extprojsdir}/libmsym)
  set (proj_bin ${PROJECT_BINARY_DIR}/External/libmsym)
  if (EXISTS ${proj}/CMakeLists.txt)
    set (MSYM_FOUND ON)

    # set and hide from the OpenMolcas CMake GUI
    make_internal (MSYM_BUILD_EXAMPLES OFF)
    make_internal (MSYM_BUILD_PYTHON   OFF)
    get_cmake_property (_in_vars VARIABLES)
    add_subdirectory (${proj} ${proj_bin} EXCLUDE_FROM_ALL)
    get_cmake_property (_out_vars VARIABLES)
    foreach (var ${_out_vars})
      get_property (in_cache CACHE ${var} PROPERTY TYPE)
      if (in_cache AND NOT var IN_LIST _in_vars)
        make_internal (${var})
      endif ()
    endforeach ()

    set (MSYM_INCLUDE ${proj}/src ${proj_bin})
    set (MSYM_LIBRARIES msym)

    message ("-- MSYM_INCLUDE_PATH: ${MSYM_INCLUDE}")
    message ("-- MSYM_LIBRARIES: ${MSYM_LIBRARIES}")

    list (APPEND EXTERNAL_LIBRARIES ${MSYM_LIBRARIES})
    list (APPEND EXTRA_INSTALL_TARGETS ${MSYM_LIBRARIES})
  else ()
    set (MSYM_FOUND OFF)
    message ("-- MSYM interface ignored (source not available, run \"${GIT} submodule update --init ${proj}\")")
  endif ()
else ()
  message ("MSYM support DISABLED")
endif ()

if (DMRG)
  message ("Configuring QCMaquis DMRG support:")

  if (DEFINED ENV{QCMaquis_ROOT})
    set (QCMaquis_ROOT $ENV{QCMaquis_ROOT})
  endif ()
  if (QCMaquis_ROOT)
    message ("-- External QCMaquis root specified at: ${QCMaquis_ROOT}")
    set (QCMaquis_ROOT ${QCMaquis_ROOT} CACHE STRING "Specify path to existing QCMaquis installation. This will allow to skip the installation of QCMaquis.")
    set (MAQUIS_DMRG_DIR ${QCMaquis_ROOT}/share/cmake)
    find_package (MAQUIS_DMRG PATHS ${QCMaquis_ROOT} NO_DEFAULT_PATH)
  endif ()

  foreach (prog
    dmrgscf
    last_energy
    mclr
    mcpdft
    mpssi
    nevpt2
    rasscf
    rassi
    wfn_util
  )
    list (APPEND ${prog}_defs "_DMRG_")
  endforeach ()

  set (custom_dir ${PROJECT_SOURCE_DIR}/cmake/custom)
  set (proj_bin ${PROJECT_BINARY_DIR}/External/qcmaquis)

  configure_file (${custom_dir}/qcmaquis.cmake ${proj_bin}/CMakeLists.txt COPYONLY)

  # Builds with build subdirectory outside the source tree require the second argument of add_subdirectory
  add_subdirectory (${proj_bin} ${proj_bin})

  # Make sure GlobalArrays is built first if we're building both
  if (TARGET GlobalArrays)
    add_dependencies (qcmaquis GlobalArrays)
  endif ()

  foreach (lib ${HDF5_QCM_LIBRARIES})
    get_filename_component (name ${lib} NAME_WE)
    if (NOT TARGET ${name})
      add_library (${name} UNKNOWN IMPORTED)
      set_target_properties (${name} PROPERTIES IMPORTED_LOCATION ${lib})
    endif ()
    target_link_libraries (${name} INTERFACE ${HDF5_C_LIBRARIES})

    list (APPEND EXTERNAL_LIBRARIES ${name})
  endforeach ()

  message ("Enabled QCMaquis DMRG support.")
  if (NOT MAQUIS_DMRG_FOUND)
    message ("-- QCMaquis libraries: ${MAQUIS_DMRG_LIBRARIES}")
  else ()
    message ("-- External QCMaquis libraries: ${MAQUIS_DMRG_LIBRARIES}")
  endif ()
  message ("-- HDF5 QCMaquis include path: ${HDF5_QCM_INCLUDE}")
  message ("-- HDF5 QCMaquis libraries: ${HDF5_QCM_LIBRARIES}")
else ()
  message ("QCMaquis DMRG support DISABLED")
endif ()

# NECI Full-CI Quantum Monte-Carlo
if (NECI)
  message(FATAL_ERROR "At the moment the embedded NECI does not work. "
    "Please use external NECI instead and contact the developers. "
    "If https://gitlab.com/Molcas/OpenMolcas/-/issues/346 is resolved we can enable it again."
  )
  if (NOT MPI)
    message (FATAL_ERROR "NECI requires MPI (-DMPI=ON).")
  endif ()
  message ("Configuring NECI support:")
  set (proj ${extprojsdir}/NECI)
  set (proj_bin ${PROJECT_BINARY_DIR}/External/NECI)
  if (EXISTS ${proj}/CMakeLists.txt)
    # set default NECI options when built from OpenMolcas
    # and hide them from the OpenMolcas CMake GUI
    make_internal (ENABLE_MOLCAS        ON)
    make_internal (ENABLE_HDF5          ON)
    make_internal (WARNINGS             OFF)
    make_internal (ENABLE_SHARED_MEMORY ON)
    make_internal (ENABLE_SHARED        ON)
    # make_internal (ENABLE_BUILD_HDF5    OFF)

    get_cmake_property (_in_vars VARIABLES)
    add_subdirectory (${proj} ${proj_bin} EXCLUDE_FROM_ALL)
    get_cmake_property (_out_vars VARIABLES)
    foreach (var ${_out_vars})
      get_property (in_cache CACHE ${var} PROPERTY TYPE)
      if (in_cache AND NOT var IN_LIST _in_vars)
        make_internal (${var})
      endif ()
    endforeach ()

    list (APPEND rasscf_defs "_NECI_")
    set (NECI_LIBRARIES libdneci)

    message ("-- NECI_INCLUDE_PATH: ${NECI_INCLUDE}")
    message ("-- NECI_LIBRARIES: ${NECI_LIBRARIES}")

    list (APPEND EXTERNAL_LIBRARIES ${NECI_LIBRARIES})
    list (APPEND EXTRA_INSTALL_TARGETS ${NECI_LIBRARIES})
  else ()
    message ("-- Internal NECI ignored (source not available, run \"${GIT} submodule update --init ${proj}\")")
    message ("   An externally configured NECI can still be used")
  endif ()
else ()
  message ("NECI interface DISABLED")
endif ()

# EFP support (effective fragment support)
if (EFPLIB)
  message ("Configuring EFP support:")
  set (proj ${extprojsdir}/efp)
  set (projname efplib)
  if (EXISTS ${proj}/CMakeLists.txt)
    if (ADDRMODE EQUAL 64)
      set (LIBEFP_FORTRAN_INT64 1)
    else ()
      set (LIBEFP_FORTRAN_INT64 0)
    endif ()

    include (ExternalProject)
    # Enabling source changes to keep ExternalProject happy
    set (CMAKE_DISABLE_SOURCE_CHANGES OFF)
    ExternalProject_Add (${projname}
                         PREFIX ${PROJECT_BINARY_DIR}/External
                         SOURCE_DIR ${proj}
                         CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
                                    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                    -DCMAKE_INSTALL_LIBDIR=lib
                                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
                                    -DENABLE_OPENMP=OFF  # OMP already in CMAKE_C_FLAGS
                                    -DFRAGLIB_UNDERSCORE_L=ON
                                    -DFRAGLIB_DEEP=ON
                                    -DLAPACK_LIBRARIES=${LINALG_LIBRARIES}
                         CMAKE_CACHE_ARGS -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                                          -DLIBEFP_FORTRAN_INT64:BOOL=${LIBEFP_FORTRAN_INT64}
                         INSTALL_COMMAND ""
    )
    set (CMAKE_DISABLE_SOURCE_CHANGES ON)

    ExternalProject_Get_Property (efplib BINARY_DIR)
    set (proj_bin ${BINARY_DIR})
    link_directories (${proj_bin})
    list (APPEND EXTERNAL_PROJECTS ${projname})

    list (APPEND integral_util_incs ${proj}/interface)

    set (EFP_LIBRARIES efp)

    message ("-- EFP_LIBRARIES: ${EFP_LIBRARIES}")

    list (APPEND EXTERNAL_LIBRARIES ${EFP_LIBRARIES})
    add_compile_definitions (_EFP_)
  else ()
    message ("-- EFP interface ignored (source not available, run git submodule update --init ${proj}\")")
  endif ()
else ()
  message ("EFP interface DISABLED")
endif ()

# Gen1Int interface
if (GEN1INT)
  message ("Configuring Gen1Int:")

  set (custom_dir ${PROJECT_SOURCE_DIR}/cmake/custom)
  set (proj_bin ${PROJECT_BINARY_DIR}/External/gen1int)

  configure_file (${custom_dir}/gen1int.cmake ${proj_bin}/CMakeLists.txt COPYONLY)
  add_subdirectory (${proj_bin} ${proj_bin})

  get_directory_property (${GEN1INT_INCLUDE} DIRECTORY ${proj_bin} DEFINITION ${GEN1INT_INCLUDE})
  get_directory_property (${GEN1INT_LIBRARIES} DIRECTORY ${proj_bin} DEFINITION ${GEN1INT_LIBRARIES})
  include_directories (${GEN1INT_INCLUDE})

  message ("-- GEN1INT_INCLUDE_PATH: ${GEN1INT_INCLUDE}")
  message ("-- GEN1INT_LIBRARIES:    ${GEN1INT_LIBRARIES}")

  list (APPEND EXTERNAL_LIBRARIES ${GEN1INT_LIBRARIES})
else ()
  message ("GEN1INT support DISABLED")
endif ()

# libwfa interface
if (WFA)
  if (NOT HDF5)
    message (FATAL_ERROR "WFA support requires enabling HDF5 too.")
  endif ()

  message ("Configuring libwfa support:")

  set (proj ${extprojsdir}/libwfa)
  set (proj_bin ${PROJECT_BINARY_DIR}/External/libwfa)
  if (EXISTS ${proj}/CMakeLists.txt)
    set (WFA_FOUND ON)
    mark_as_advanced (CLEAR ARMADILLO_INC)
    set (ARMADILLO_INC ${ARMADILLO_INC} CACHE PATH "Directory with armadillo header files")
    make_internal (MOLCAS_LIB ON)

    get_cmake_property (_in_vars VARIABLES)
    add_subdirectory (${proj}/libwfa)
    get_cmake_property (_out_vars VARIABLES)
    foreach (var ${_out_vars})
      get_property (in_cache CACHE ${var} PROPERTY TYPE)
      if (in_cache AND NOT var IN_LIST _in_vars)
        make_internal (${var})
      endif ()
    endforeach ()

    target_include_directories (wfa_molcas PRIVATE ${proj} ${ARMADILLO_INC} ${HDF5_INCLUDE_DIRS})
    target_compile_definitions (wfa_molcas PRIVATE "ARMA_DONT_USE_WRAPPER")
    target_link_libraries (wfa_molcas INTERFACE ${HDF5_CXX_LIBRARIES} ${LINALG_LIBRARIES})
    list (APPEND wfa_libs wfa_molcas)

    message ("-- wfa_libs: ${wfa_libs}")
    message ("-- ARMADILLO_INC: ${ARMADILLO_INC}")

    list (APPEND EXTERNAL_LIBRARIES ${wfa_libs})
  else ()
    set (WFA_FOUND OFF)
    mark_as_advanced (FORCE ARMADILLO_INC)
    message ("-- libwfa support ignored (source not available, run \"${GIT} submodule update --init ${proj}\")")
  endif ()
else ()
  mark_as_advanced (FORCE ARMADILLO_INC)
  message ("libwfa support DISABLED")
endif ()

# NEVPT2
if (NEVPT2)
  if (NOT DMRG)
    message (FATAL_ERROR "NEVPT2 support requires enabling DMRG too.")
  endif ()
  if (NOT (HDF5 AND HDF5_FOUND))
    message (FATAL_ERROR "NEVPT2 support requires enabling HDF5 too.")
  endif ()

  message ("Configuring NEVPT2:")

  set (custom_dir ${PROJECT_SOURCE_DIR}/cmake/custom)
  set (proj_bin ${PROJECT_BINARY_DIR}/External/nevpt2_ext)
  configure_file (${custom_dir}/nevpt2.cmake ${proj_bin}/CMakeLists.txt COPYONLY)

  get_cmake_property (_in_vars VARIABLES)
  add_subdirectory (${proj_bin} ${proj_bin})
  get_cmake_property (_out_vars VARIABLES)
  foreach (var ${_out_vars})
    get_property (in_cache CACHE ${var} PROPERTY TYPE)
    if (in_cache AND NOT var IN_LIST _in_vars)
      make_internal (${var})
    endif ()
  endforeach ()

  get_directory_property (${NEVPT2_INCLUDE} DIRECTORY ${proj_bin} DEFINITION ${NEVPT2_INCLUDE})
  get_directory_property (${NEVPT2_LIBRARIES} DIRECTORY ${proj_bin} DEFINITION ${NEVPT2_LIBRARIES})

  add_dependencies (nevpt2_ext qcmaquis-driver)

  message ("-- NEVPT2_LIBRARIES: ${NEVPT2_LIBRARIES}")
else ()
  message ("NEVPT2 support DISABLED")
endif ()

# MolGUI
if (MolGUI)
  message ("Configuring MolGUI:")

  set (projname MolGUI)
  set (projdir ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${projname})
  set (projinstall ${PROJECT_BINARY_DIR}/${projname})

  include (ExternalProject)
  ExternalProject_Add (${projname}
                       GIT_REPOSITORY https://gitlab.com/Molcas/MolGUI.git
                       GIT_TAG 599b75e79fab49501adb313f3e2dfab9d8211795
                       PREFIX ${projdir}
                       CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${projinstall}
                                  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                                  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                                  -DMOLCAS=${PROJECT_BINARY_DIR}
  )
  message ("-- MolGUI will be downloaded and compiled in ${projdir}")
  message ("-- MolGUI will be installed in ${projinstall}")

  list (APPEND EXTRA_INSTALL_DIRS ${projinstall})
else ()
  message ("MolGUI DISABLED")
endif ()

# Add the "bigot" flags after processing the external projects, so they are not affected
if (BIGOT)
  # Modify the "BUILD_TARGET" flags because we want them at the end
  string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TARGET)
  set (CMAKE_CXX_FLAGS_${BUILD_TARGET} "${CMAKE_CXX_FLAGS_${BUILD_TARGET}} ${CXXFLAGS_${CXXCID}_BIGOT}")
  set (CMAKE_C_FLAGS_${BUILD_TARGET} "${CMAKE_C_FLAGS_${BUILD_TARGET}} ${CFLAGS_${CCID}_BIGOT}")
  set (CMAKE_Fortran_FLAGS_${BUILD_TARGET} "${CMAKE_Fortran_FLAGS_${BUILD_TARGET}} ${FFLAGS_${FCID}_BIGOT}")
endif ()

################################################################################
#                                                                              #
# Configure settings for the OpenMolcas run-time environment                   #
#                                                                              #
################################################################################

message ("Configuring runtime environment settings:")

cmake_host_system_information (RESULT host_mem QUERY TOTAL_PHYSICAL_MEMORY)
cmake_host_system_information (RESULT host_name QUERY HOSTNAME)
cmake_host_system_information (RESULT host_system QUERY OS_NAME)

# memory/disk default sizes
#==========================

if (ADDRMODE EQUAL 64)
  set (mem "2048")
  set (disk "20000")
else ()
  set (mem "1024")
  set (disk "2047")
endif ()
math (EXPR max_mem "${host_mem}*3/4")
if (${max_mem} LESS ${mem})
  set (mem ${max_mem})
endif ()
set (DEFMOLCASMEM ${mem} CACHE STRING "Set default max memory (MiB).")
set (DEFMOLCASDISK ${disk} CACHE STRING "Set default max file size (MiB).")

message ("-- DEFMOLCASMEM:  ${DEFMOLCASMEM}")
message ("-- DEFMOLCASDISK: ${DEFMOLCASDISK}")

# launcher settings
#==================

# The runtime can be controlled by the user through two variables:
# SER_LAUNCHER: a command used to launch a serial program
# MPI_LAUNCHER: a command used to launch an MPI program
# Since these are advanced options, the user would typically modify
# them in the CMakeCache.txt. However, they can also be activated
# by other options for internal use, which is why their presence
# is detected instead of just caching a default string.

# Scripts are tricky, we don't want to set automatic wrappers
set (RUNSCRIPT "$program $input" CACHE STRING "Command to launch scripts.")

# Specify the launch command for a serial run
if (NOT SER_LAUNCHER)
  set (SER_LAUNCHER ""                CACHE STRING "Command to launch serial programs.")
else ()
  set (SER_LAUNCHER "${SER_LAUNCHER}" CACHE STRING "Command to launch serial programs.")
endif ()
set (RUNBINARYSER "${SER_LAUNCHER} $program")

# Specify the launch command for an MPI run
if (MPI)
  if (NOT MPI_LAUNCHER)
    set (MPI_LAUNCHER "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} $MOLCAS_NPROCS ${MPIEXEC_PREFLAGS}" CACHE STRING "Command to launch MPI programs.")
  endif ()
  set (RUNBINARY "${MPI_LAUNCHER} $program")
else ()
  set (RUNBINARY "${SER_LAUNCHER} $program")
endif ()

string (STRIP ${RUNSCRIPT} RUNSCRIPT)
string (STRIP ${RUNBINARY} RUNBINARY)
string (STRIP ${RUNBINARYSER} RUNBINARYSER)

message ("-- RUNSCRIPT:    ${RUNSCRIPT}")
message ("-- RUNBINARY:    ${RUNBINARY}")
message ("-- RUNBINARYSER: ${RUNBINARYSER}")

################################################################################
#                                                                              #
# Print summary                                                                #
#                                                                              #
################################################################################

message ("Build type: ${CMAKE_BUILD_TYPE}")
string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TARGET)
if (NOT BUILD_TARGET STREQUAL "NONE")
  message ("-- C compiler: ${CMAKE_C_COMPILER}")
  message ("-- C compiler flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${BUILD_TARGET}}")
  message ("-- Fortran compiler: ${CMAKE_Fortran_COMPILER}")
  message ("-- Fortran compiler flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${BUILD_TARGET}}")
else ()
  message ("-- C compiler: ${CMAKE_C_COMPILER}")
  message ("-- C compiler flags: ${CMAKE_C_FLAGS}")
  message ("-- Fortran compiler: ${CMAKE_Fortran_COMPILER}")
  message ("-- Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
endif ()
get_directory_property (defs COMPILE_DEFINITIONS)
message ("-- Definitions: ${defs}")
message ("-- Debug definitions: ${DEBUG_DEFS}")

################################################################################
#                                                                              #
# Utilities/Programs                                                           #
#                                                                              #
################################################################################

# Enable interprocedural optimization, if requested and supported, from now on
set (set_IPO OFF)
if (IPO AND ";RelWithDebInfo;Release;Fast;" MATCHES ";${CMAKE_BUILD_TYPE};")
  include (CheckIPOSupported)
  check_ipo_supported (RESULT set_IPO)
endif ()
set (CMAKE_INTERPROCEDURAL_OPTIMIZATION ${set_IPO})

# generate utility list
#======================
# FIXME: this should be split in several sublists, e.g. core, blas, etc.
# Then it will be easier to construct custom utility lists.

# autodiscover program and utility directories
foreach (dir ${source_roots})
  # find all files not starting with a dot
  file (GLOB dir_progs RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${dir}/[^.]*")
  file (GLOB dir_utils RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${dir}/[^.]*_util")
  # programs are those that do not end in "_util"
  if (dir_utils)
    list (REMOVE_ITEM dir_progs ${dir_utils})
  endif ()
  # and consider only directories not beginning with uppercase
  foreach (item ${dir_progs})
    if (NOT IS_DIRECTORY "${PROJECT_SOURCE_DIR}/${dir}/${item}")
      list (REMOVE_ITEM dir_progs ${item})
    endif ()
    if (item MATCHES "^[A-Z]")
      list (REMOVE_ITEM dir_progs ${item})
    endif ()
  endforeach ()
  foreach (item ${dir_utils})
    if (NOT IS_DIRECTORY "${PROJECT_SOURCE_DIR}/${dir}/${item}")
      list (REMOVE_ITEM dir_utils ${item})
    endif ()
    if (item MATCHES "^[A-Z]")
      list (REMOVE_ITEM dir_utils ${item})
    endif ()
  endforeach ()
  list (APPEND utils ${dir_utils})
  list (APPEND progs ${dir_progs})
endforeach ()
list (REMOVE_DUPLICATES utils)
list (REMOVE_DUPLICATES progs)
# some directories are special or managed manually, remove them from the automatic list
list (REMOVE_ITEM utils
  block_dmrg_util
  chemps2_util
  dice_util
  delayed_util
  dga_util
  embedding_util
  hdf5_util
  msym_util
)
list (REMOVE_ITEM progs
  parnell
  symmetrize
  wfa
  nevpt2
  dmrgscf
  mpssi
)

# set the location of the sources for the utilities
foreach (util ${utils})
  find_source (${util} ${source_roots})
  if (NOT DEFINED ${util}_src)
    list (REMOVE_ITEM utils ${util})
    message (WARNING "\"${util}\" not found in source directories, skipping.")
  endif ()
endforeach ()

# activate FDE interface
if (FDE)
  set (util embedding_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
  # directories using FDE
  foreach (prog
    gateway_util
    rasscf
    scf
    seward
  )
    list (APPEND ${prog}_defs "_FDE_")
  endforeach ()
endif ()

# If an external GA library is not present, we use the internal
# DGA code which resides in dga_util:
if (MPI AND NOT GA)
  set (util dga_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
    include_directories (${PROJECT_SOURCE_DIR}/${${util}_src})
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories. "
             "Enable Global Arrays (-DGA=ON) or disable MPI (-DMPI=OFF).")
  endif ()
  # workaround for OpenMPI
  if (MPI_IMPLEMENTATION STREQUAL "openmpi" AND OPENMP)
    list (APPEND ${util}_defs "_OPENMPI_")
    message (WARNING
             "MPI_THREAD_MULTIPLE support in OpenMPI is unreliable in "
             "some circumstances, using MPI_THREAD_SERIALIZED instead.")
  endif ()
endif ()

# Add optional interface to dynamic runtime BLAS/LAPACK
if (LINALG STREQUAL "Runtime")
  set (util delayed_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
  list (APPEND ${util}_incs ${internal_linalg})
endif ()

# activate HDF5 wrappers for I/O
if (HDF5 AND HDF5_FOUND)
  set (util hdf5_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
    list (APPEND ${util}_incs ${HDF5_INCLUDE_DIRS})
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
  # additional utilities doing internal HDF5 calls
  foreach (prog
    aniso_util
    gateway_util
    hyper_util
    molcas_ci_util
    msym_util
    oneint_util
    property_util
    restart_util
    runfile_util
    system_util
    wfn_util
  )
    list (APPEND ${prog}_defs "_HDF5_")
  endforeach ()
  if (DMRG)
    foreach (util
      transform_util
    )
      list (APPEND ${util}_defs "_HDF5_QCM_")
      list (APPEND ${util}_deps qcmaquis-hdf5-interface)
    endforeach ()
  endif ()
endif ()

if (NEVPT2)
  set (prog nevpt2)
  list (APPEND progs ${prog})
  list (APPEND ${prog}_incs ${NEVPT2_INCLUDE})
  list (APPEND ${prog}_deps nevpt2_ext)
  list (APPEND ${prog}_deplibs ${NEVPT2_LIBRARIES})
endif ()

# activate MSYM support
if (MSYM AND MSYM_FOUND)
  set (util msym_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
  list (APPEND ${util}_incs ${MSYM_INCLUDE})
  list (APPEND progs symmetrize)
endif ()

if (BLOCK)
  set (util block_dmrg_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
    list (APPEND ${util}_defs "_ENABLE_BLOCK_DMRG_")
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
endif ()

if (CHEMPS2)
  set (util chemps2_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
    list (APPEND ${util}_incs ${HDF5_INCLUDE_DIRS})
    list (APPEND ${util}_defs "_ENABLE_CHEMPS2_DMRG_")
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
endif ()

if (DICE)
  set (util dice_util)
  find_source (${util} ${source_roots})
  if (DEFINED ${util}_src)
    list (APPEND utils ${util})
    list (APPEND ${util}_defs "_ENABLE_DICE_SHCI_")
  else ()
    message (FATAL_ERROR "\"${util}\" not found in source directories.")
  endif ()
endif ()

# activate libwfa support
if (WFA AND WFA_FOUND)
  list (APPEND progs wfa)
endif ()

if (BLOCK OR DMRG)
  list (APPEND progs dmrgscf)
endif ()

if (DMRG)
  list (APPEND progs mpssi)
endif ()

# finally, sort the list
list (SORT utils)

# set the location of the sources for the programs
foreach (prog ${progs})
  find_source (${prog} ${source_roots})
  if (NOT DEFINED ${prog}_src)
    list (REMOVE_ITEM progs ${prog})
    message (WARNING "\"${prog}\" not found in source directories, skipping.")
  endif ()
endforeach ()

if (MSYM AND MSYM_FOUND)
  list (APPEND scf_defs "_MSYM_")
endif ()

if (WFA AND WFA_FOUND)
  list (APPEND wfa_defs "_WFA_")
endif ()

if (GEN1INT)
  foreach (prog
    seward
    integral_util
  )
    list (APPEND ${prog}_defs "_GEN1INT_")
  endforeach ()
endif ()

# programs doing HDF5 calls
if (HDF5 AND HDF5_FOUND)
  foreach (prog
    caspt2
    dynamix
    guessorb
    mcpdft
    nevpt2
    poly_aniso
    rasscf
    rassi
    rhodyn
    scf
    single_aniso
    slapaf
    surfacehop
  )
    list (APPEND ${prog}_defs "_HDF5_")
  endforeach ()
  if (DMRG)
    foreach (prog
      motra
      nevpt2
    )
      list (APPEND ${prog}_defs "_HDF5_QCM_")
      list (APPEND ${prog}_deplibs qcmaquis-hdf5-interface)
    endforeach ()
  endif ()
else ()
  message (WARNING "RHODYN program requires HDF5 interface, disabled")
  list (REMOVE_ITEM progs rhodyn)
endif ()

# Libraries that require QCMaquis-OpenMolcas interface
# must link with libstdc++ and against QCMaquis libraries
# (CXX linker does NOT work on some architectures)
# All other libraries do not require this (avoiding unnecessary linking saves compilation time)
if (DMRG)
  # Link QCMaquis libraries
  list (APPEND dmrg_progs
    casvb
    dmrgscf
    last_energy
    loprop
    mcpdft
    mpssi
    nevpt2
    numerical_gradient
    rasscf
    rassi
  )
  foreach (prog ${dmrg_progs})
    list (APPEND ${prog}_deplibs ${MAQUIS_DMRG_LIBRARIES} -lstdc++)
  endforeach ()
endif ()

if (GA AND MPI)
  foreach (prog
           caspt2
           cholesky_util
           cht3
           ga_util
           mclr
           rasscf
           ri_util
  )
    list (APPEND ${prog}_incs ${GA_INCLUDE_PATH})
  endforeach ()
endif ()

# finally, sort the list
list (SORT progs)

# set program dependencies of super modules
# (note dependencies are not recursive with object libraries, so child dependencies must be explicit)
set (casvb_deplibs
  rasscf
)
set (dmrgscf_deplibs
  rasscf
)
set (last_energy_deplibs
  caspt2
  ccsdt
  chcc
  cht3
  false
  mbpt2
  mcpdft
  motra
  rasscf
  scf
  seward
)
set (loprop_deplibs
  caspt2
  mbpt2
  rasscf
  scf
)
set (mpssi_deplibs
  rassi
)
set (numerical_gradient_deplibs
  caspt2
  ccsdt
  chcc
  cht3
  false
  mbpt2
  mcpdft
  motra
  rasscf
  scf
  seward
)
set (seward_deplibs
  guessorb
)

if (DMRG)
  list (APPEND last_energy_deplibs dmrgscf)
endif ()

################################################################################
#                                                                              #
# Disable some warnings that are buggy or on code we don't want to modify      #
#                                                                              #
################################################################################
# netlib blas/lapack produces warnings
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set (netlib_nowarn_flags " -w")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set (netlib_nowarn_flags " -w")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro")
  set (netlib_nowarn_flags " -erroff=%all")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  set (netlib_nowarn_flags " -Kieee -w")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set (netlib_nowarn_flags " -w -dcfuns")
endif ()
if (netlib_nowarn_flags)
  file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${LinAlg_internal_src}/*.f ${LinAlg_internal_src}/*.F90 ${delayed_util_src}/*.f ${delayed_util_src}/*.F90)
  set_property (SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS ${netlib_nowarn_flags})
endif ()

################################################################################
#                                                                              #
# Blacklisted sources                                                          #
#                                                                              #
################################################################################

if (DEFINED BLACKLIST)
  foreach (module ${BLACKLIST})
    file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR}
      ${${module}_src}/*.f
      ${${module}_src}/*.f90
    )
    set_property (
      SOURCE ${source}
      APPEND_STRING
      PROPERTY
      COMPILE_FLAGS " ${BLACKFLAGS}"
    )
  endforeach ()
  message ("Blacklisting enables, someone doesn't like certain modules...")
  message ("-- modules: ${BLACKLIST}")
  message ("-- *additional* flags: ${BLACKFLAGS}")
endif ()

################################################################################
#                                                                              #
# Build the utilities and the libmolcas library                                #
#                                                                              #
################################################################################

if (BUILD_SHARED_LIBS)
  set (CMAKE_POSITION_INDEPENDENT_CODE "TRUE")
endif ()

# first, build global modules separately, so other utils can add it as a
# dependency using the "mods_obj" target.
# note that mods_obj is empty at the beginning, source files will be added later
add_library (mods_obj OBJECT)
if (TARGET Libxc)
  add_dependencies (mods_obj Libxc)
endif ()
if (GA_BUILD AND MPI)
  add_dependencies (mods_obj GlobalArrays)
endif ()
if (GEN1INT)
  add_dependencies (mods_obj gen1int)
endif ()
if (BUILD_SHARED_LIBS)
  set_target_properties (mods_obj PROPERTIES POSITION_INDEPENDENT_CODE "TRUE")
endif ()
set_target_properties (mods_obj PROPERTIES Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR})

# set the module directory for internal linalg
if (LINALG STREQUAL "Internal")
  foreach (lib ${LINALG_LIBRARIES})
    set_module_directory (${lib})
  endforeach ()
endif ()

# now go through utility list
foreach (util ${utils})
  add_subdirectory (src/${util} ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${util})
  add_dependencies (${util}_obj mods_obj)
  if (${util}_mods)
    # move the module files to the mods_obj target
    # (note COMPILE_DEFINITIONS must be on target level to properly compute dependencies)
    pass_properties_to_target (${util}_obj mods_obj "COMPILE_DEFINITIONS")
    pass_properties_to_files (${util}_obj "${${util}_mods}" "INCLUDE_DIRECTORIES")
    target_sources (mods_obj PRIVATE ${${util}_mods})
  endif ()
  list (APPEND utils_obj ${util}_obj)
  list (APPEND ALL_SOURCES ${${util}_sources} ${${util}_mods})
  # get the list of "public" include directories, since these may be needed
  # when compiling programs, and not just for linking
  get_target_property(inc ${util}_obj INTERFACE_INCLUDE_DIRECTORIES)
  list (APPEND public_incs ${inc})
endforeach (util)

# build an actual libmolcas library
add_library (libmolcas "")
set_target_properties (libmolcas PROPERTIES PREFIX "")
set (MOLCAS_LIBRARIES "libmolcas")
target_link_libraries (libmolcas ${utils_obj})
target_link_libraries (libmolcas mods_obj)
# also create a static library!
if (BUILD_SHARED_LIBS AND BUILD_STATIC_LIBS)
  add_library (libmolcas_static STATIC "")
  set_target_properties (libmolcas_static PROPERTIES PREFIX "")
  set_target_properties (libmolcas_static PROPERTIES OUTPUT_NAME libmolcas)
  list (APPEND MOLCAS_LIBRARIES "libmolcas_static")
  target_link_libraries (libmolcas_static ${utils_obj})
  target_link_libraries (libmolcas_static mods_obj)
endif ()
# make sure these are installed
if (LINALG STREQUAL "Internal")
  list (APPEND MOLCAS_LIBRARIES ${INTERNAL_LINALG})
endif ()

################################################################################
#                                                                              #
# Build the program modules                                                    #
#                                                                              #
################################################################################

add_custom_target (only_objs
  DEPENDS ${utils_obj}
)

# set up proper RPATH for executables, use RUNPATH if on Linux
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH "TRUE")

# build the standard programs
foreach (prog ${progs})
  add_subdirectory (src/${prog} ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog})
  if (TARGET ${prog}_obj)
    add_dependencies (${prog}_obj ${utils_obj})
    add_dependencies (only_objs   ${prog}_obj)
  endif ()
  if (EXTERNAL_PROJECTS)
    add_dependencies (${prog}.exe ${EXTERNAL_PROJECTS})
  endif ()
  if (CUBLAS)
    set_target_properties (${prog}.exe PROPERTIES LINKER_LANGUAGE "CXX")
  endif ()
  list (APPEND progs_exe ${prog}.exe)
  list (APPEND ALL_SOURCES ${${prog}_sources})
  list (APPEND PROGRAM_EXECUTABLES ${prog}.exe)
  if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    if (${CMAKE_SYSTEM_VERSION} STREQUAL "15.0.0")
      if (LINALG STREQUAL "MKL")
        add_custom_target (dylibpath_${prog} ALL
                           DEPENDS $<TARGET_FILE:${prog}.exe>
                           COMMAND ${OPENMOLCAS_DIR}/sbin/fixdylib $<TARGET_FILE:${prog}.exe> libmkl ${MKLROOT}
                           COMMENT "Enforce absolute path to MKL in dynamically linked executables"
                           VERBATIM
        )
        add_dependencies (dylibpath_${prog} ${prog}.exe)
      endif ()
    endif ()
  endif ()
endforeach (prog)

# Set proper install RPATH if using internal QCMaquis
if (DMRG)
  foreach (prog ${dmrg_progs})
    if (NOT MAQUIS_DMRG_FOUND AND TARGET ${prog}.exe)
      set_target_properties (${prog}.exe PROPERTIES
        INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/qcmaquis/lib:${CMAKE_INSTALL_RPATH}"
      )
    endif ()
  endforeach ()
endif ()

# build parnell as a stand-alone executable
set (prog parnell)
if (MPI)
  list (APPEND ${prog}_deplibs ${MPI_C_LIBRARIES})
endif ()
add_subdirectory (src/${prog} ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${prog})
list (APPEND progs_exe ${prog}.exe)
list (APPEND ALL_SOURCES ${${prog}_sources})
list (APPEND PROGRAM_EXECUTABLES ${prog}.exe)

# Add custom debug flags per file:
#   For each file, tries to find out if each definition is used and modify the file's properties
# Limitation: it will rebuild the whole target (but not the whole OpenMolcas)
if (DEBUG_DEFS)
  get_directory_property (subdirs SUBDIRECTORIES)
  foreach (file ${ALL_SOURCES})
    file (READ ${file} contents)
    foreach (def ${DEBUG_DEFS})
      string (REGEX MATCH ${def} result "${contents}")
      if (NOT ${result} STREQUAL "")
        set_property (SOURCE ${file} DIRECTORY ${subdirs} APPEND PROPERTY COMPILE_DEFINITIONS ${def})
      endif ()
    endforeach ()
  endforeach ()
endif ()

################################################################################
#                                                                              #
# Tools, i.e. builds are dependent on libmolcas library                        #
#                                                                              #
################################################################################

if (TOOLS)
  message ("Configuring tools:")
  foreach (BASE_DIR ${basedirs})
    set (toolsdir ${BASE_DIR}/Tools)
    file (GLOB tools RELATIVE ${toolsdir} ${toolsdir}/*)
    list (SORT tools)
    foreach (tool ${tools})
      if (IS_DIRECTORY ${toolsdir}/${tool})
        if (EXISTS ${toolsdir}/${tool}/CMakeLists.txt)
          set (dir_error "")
          add_subdirectory (${toolsdir}/${tool} ${PROJECT_BINARY_DIR}/Tools/${tool})
          if ("${dir_error}" STREQUAL "")
            message ("-- ${tool}: added to targets")
          else ()
            message ("-- ${tool}: ${dir_error}")
            if ("${tool}" STREQUAL "pymolcas" AND (NOT DEFINED EXTRA_DIR))
              message (FATAL_ERROR "Failed to configure the pymolcas driver.")
            endif ()
          endif ()
        else ()
          message ("-- ${tool}: ignored (no CMake support)")
        endif ()
      endif ()
    endforeach ()
  endforeach ()
endif ()

################################################################################
#                                                                              #
# Post-build actions                                                           #
#                                                                              #
################################################################################

# create runtime environment file
file (WRITE ${PROJECT_BINARY_DIR}/molcas.rte
  "# molcas runtime environment\n"
  "OS='${OS}'\n"
  "PARALLEL='${MPI}'\n"
  "DEFMOLCASMEM='${DEFMOLCASMEM}'\n"
  "DEFMOLCASDISK='${DEFMOLCASDISK}'\n"
  "RUNSCRIPT='${RUNSCRIPT}'\n"
  "RUNBINARY='${RUNBINARY}'\n"
  "RUNBINARYSER='${RUNBINARYSER}'\n"
)

# create version file
file (WRITE ${PROJECT_BINARY_DIR}/.molcasversion
  ${OPENMOLCAS_VERSION}
)
if (DEFINED EXTRA_VERSION)
  file (APPEND ${PROJECT_BINARY_DIR}/.molcasversion
    "\n${EXTRA_VERSION}"
  )
endif ()
file (WRITE ${PROJECT_BINARY_DIR}/.molcashome)

# copy LICENSE and CONTRIBUTORS
file (COPY "${OPENMOLCAS_DIR}/LICENSE" DESTINATION ${PROJECT_BINARY_DIR})
file (COPY "${OPENMOLCAS_DIR}/CONTRIBUTORS.md" DESTINATION ${PROJECT_BINARY_DIR})

# copy the necessary files from data
file (REMOVE_RECURSE ${PROJECT_BINARY_DIR}/data)
configure_file (${OPENMOLCAS_DIR}/data/rysrw.ascii      data/rysrw            COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/abdata.ascii     data/abdata           COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/banner.src       data/banner.txt       COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/modalias.src     data/modalias.txt     COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/landing.src      data/landing.txt      COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/isotope_data.src data/isotope_data.txt COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/functionals.src  data/functionals.txt  COPYONLY)
file (COPY ${OPENMOLCAS_DIR}/data/inputs.tpl DESTINATION ${PROJECT_BINARY_DIR}/data)
if (DEFINED EXTRA_DIR)
  foreach (f alias.plx auto.prgm landing.txt)
    file (COPY ${EXTRA_DIR}/data/${f} DESTINATION ${PROJECT_BINARY_DIR}/data)
  endforeach ()
endif ()
# copy the necessary files from sbin
file (REMOVE_RECURSE ${PROJECT_BINARY_DIR}/sbin)
foreach (f setup verify version find_sources have_feature help_basis help_func help_doc)
  file (COPY ${OPENMOLCAS_DIR}/sbin/${f} DESTINATION ${PROJECT_BINARY_DIR}/sbin)
endforeach ()
if (DEFINED EXTRA_DIR)
  file (COPY ${EXTRA_DIR}/sbin DESTINATION ${PROJECT_BINARY_DIR})
endif ()
# copy basis sets
file (REMOVE_RECURSE ${PROJECT_BINARY_DIR}/basis_library)
foreach (BASE_DIR ${basedirs})
  if (EXISTS ${BASE_DIR}/basis_library)
    file (COPY ${BASE_DIR}/basis_library DESTINATION ${PROJECT_BINARY_DIR})
  endif ()
endforeach ()
# configure tests
file (REMOVE_RECURSE ${PROJECT_BINARY_DIR}/test)
file (APPEND ${PROJECT_BINARY_DIR}/test/testdirs
  ${OPENMOLCAS_DIR}/test\n
)
if (DEFINED EXTRA_DIR)
  file (APPEND ${PROJECT_BINARY_DIR}/test/testdirs
    ${EXTRA_DIR}/test\n
  )
endif ()
if (GEN1INT)
  file (APPEND ${PROJECT_BINARY_DIR}/test/testdirs
    ${extprojsdir}/gen1int-molcaslib/test\n
  )
endif ()

# find the molcas.driver script or use a dummy driver
if (DEFINED EXTRA_DIR)
  set (MOLCAS_DRIVER "${EXTRA_DIR}/sbin/molcas.driver")
  execute_process (
    COMMAND ${OPENMOLCAS_DIR}/sbin/install_driver.sh ${MOLCAS_DRIVER}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    RESULT_VARIABLE INSTALL_DRIVER_RC
  )
  if (INSTALL_DRIVER_RC)
    message (FATAL_ERROR "Failed to install the Molcas driver.")
  endif ()
else ()
  set (MOLCAS_DRIVER "${OPENMOLCAS_DIR}/sbin/dummy.driver")
  # if molcas-extra is not available, enforce building pymolcas
  if (NOT TOOLS)
    set (BASE_DIR "${OPENMOLCAS_DIR}")
    set (toolsdir "${OPENMOLCAS_DIR}/Tools")
    set (tool "pymolcas")
    add_subdirectory (${toolsdir}/${tool} ${PROJECT_BINARY_DIR}/Tools/${tool})
    if ("${dir_error}" STREQUAL "")
      message ("-- ${tool}: added to targets")
    else ()
      message ("-- ${tool}: ${dir_error}")
      message (FATAL_ERROR "Failed to configure the pymolcas driver.")
    endif ()
  endif ()
  if (MolGUI)
    add_dependencies (MolGUI pymolcas_target)
  endif ()
endif ()

get_directory_property (defs COMPILE_DEFINITIONS)
string (STRIP ${CMAKE_C_FLAGS} CMAKE_C_FLAGS)
string (STRIP ${CMAKE_Fortran_FLAGS} CMAKE_Fortran_FLAGS)
if (NOT BUILD_TARGET STREQUAL "NONE")
  file (WRITE "${PROJECT_BINARY_DIR}/data/info.txt"
    "\nconfiguration info\n"
    "------------------\n"
    "Host name: ${host_name} (${host_system})\n"
    "C Compiler ID: ${CMAKE_C_COMPILER_ID}\n"
    "C flags: ${CMAKE_C_FLAGS}\n"
    "Fortran Compiler ID: ${CMAKE_Fortran_COMPILER_ID}\n"
    "Fortran flags: ${CMAKE_Fortran_FLAGS}\n"
    "Definitions: ${defs}\n"
    "Parallel: ${MPI} (GA=${GA})\n"
  )
endif ()

# install git hooks
if (EXISTS ${GIT})
  foreach (BASE_DIR ${basedirs})
    execute_process (
      COMMAND ${OPENMOLCAS_DIR}/sbin/install_hooks.sh ${BASE_DIR}
      WORKING_DIRECTORY ${BASE_DIR}
    )
  endforeach ()
endif ()

# fetch the prebuilt executable (if using molcas-extra)
if (DEFINED EXTRA_DIR)
  file (WRITE ${PROJECT_BINARY_DIR}/xbin.cfg "PLATFORM = ${PLATFORM}\n")

  execute_process (
    COMMAND "${MOLCAS_DRIVER}" getemil
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    ERROR_VARIABLE GETEMIL_DESCRIBE_ERROR
    RESULT_VARIABLE GETEMIL_DESCRIBE_RC
  )
  if (GETEMIL_DESCRIBE_RC)
    message (FATAL_ERROR "Failed to fetch Molcas binary: ${GETEMIL_DESCRIBE_ERROR}")
  endif ()
  list (APPEND progs_exe "molcas.exe")
endif ()

# remove executables that should not be there
file (GLOB exe RELATIVE ${EXECUTABLE_OUTPUT_PATH} "${EXECUTABLE_OUTPUT_PATH}/*.exe")
foreach (f ${exe})
  if (NOT ";${progs_exe};" MATCHES ";${f};")
    file (REMOVE "${EXECUTABLE_OUTPUT_PATH}/${f}")
  endif ()
endforeach ()

# generate .prgm files from source
foreach (BASE_DIR ${basedirs})
  list (APPEND prgmincs "-I${BASE_DIR}/src/Driver")
  file (GLOB incs ${BASE_DIR}/src/Driver/*.inc)
  list (APPEND prgmincfiles ${incs})
endforeach (BASE_DIR)
string (REGEX REPLACE "([^;]+)" "-D\\1" defargs "${defs}")
message ("Definitions: ${defargs}")
foreach (BASE_DIR ${basedirs})
  file (GLOB source RELATIVE ${BASE_DIR}/src/Driver ${BASE_DIR}/src/Driver/*.prgm.src)
  if (DEFINED external_prgm)
    list (REVERSE external_prgm)
  endif ()
  foreach (prgmsrc ${external_prgm})
    file (GLOB ext_source RELATIVE ${BASE_DIR}/src/Driver ${prgmsrc})
    list (INSERT source 0 ${ext_source})
  endforeach ()
  foreach (prgmsrc ${source})
    string (REPLACE ".prgm.src" ".prgm" prgm ${prgmsrc})
    get_filename_component (prgm ${prgm} NAME)
    if (NOT ";${prgmtargets};" MATCHES ";data/${prgm};")
      add_custom_command (
        OUTPUT data/${prgm}
        DEPENDS ${BASE_DIR}/src/Driver/${prgmsrc} ${prgmincfiles}
        COMMAND cpp -P ${defargs} ${prgmincs} ${BASE_DIR}/src/Driver/${prgmsrc} data/${prgm}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )
      list (APPEND prgmtargets data/${prgm})
    endif ()
  endforeach (prgmsrc)
endforeach (BASE_DIR)
add_custom_target (prgms ALL
  DEPENDS ${prgmtargets}
)
add_dependencies (mods_obj prgms)

# generate help databases from doc
if (Python_FOUND)
  execute_process (
    COMMAND ${Python_EXECUTABLE} ${OPENMOLCAS_DIR}/sbin/extract_xmldoc.py ${OPENMOLCAS_DIR}/doc ${PROJECT_BINARY_DIR}/data
    ERROR_VARIABLE EXTRACT_ERROR
    RESULT_VARIABLE EXTRACT_RESULT
  )
  if (EXTRACT_ERROR)
    message (FATAL_ERROR "Running \"extract_xmldoc.py\" failed: ${EXTRACT_ERROR}")
  endif ()
  # check XML consistency, only if the program is found
  find_program (XMLLINT "xmllint")
  mark_as_advanced (FORCE XMLLINT)
  if (XMLLINT)
    execute_process (
      COMMAND ${XMLLINT} ${PROJECT_BINARY_DIR}/data/keyword.xml "-noout"
      ERROR_VARIABLE XMLLINT_ERROR
      RESULT_VARIABLE XMLLINT_RESULT
    )
    if (XMLLINT_ERROR)
      message (FATAL_ERROR "Extracted XML file is malformed: ${XMLLINT_ERROR}")
    endif ()
  endif ()
else ()
  message (WARNING "No Python found, help databases will not be created.")
endif ()

# generate rcodes.txt file
file (STRINGS "${OPENMOLCAS_DIR}/src/Include/warnings.h" lines)
file (WRITE "${PROJECT_BINARY_DIR}/data/rcodes.txt")
foreach (line ${lines})
  string (REGEX MATCH "^#define *(_RC_[^ ]*) *([0123456789]*)" match ${line})
  if (NOT "${match}" STREQUAL "")
    file (APPEND "${PROJECT_BINARY_DIR}/data/rcodes.txt"
      "${CMAKE_MATCH_1} = ${CMAKE_MATCH_2}\n"
    )
  endif ()
endforeach ()

################################################################################
#                                                                              #
# Testing                                                                      #
#                                                                              #
################################################################################

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  enable_testing ()
  add_subdirectory (unit_tests)
endif ()

################################################################################
#                                                                              #
# Documentation                                                                #
#                                                                              #
################################################################################

message ("Configuring documentation")

find_program (SPHINX_EXECUTABLE
  NAMES sphinx-build
  DOC "Sphinx Documentation Builder (sphinx-doc.org)"
)
mark_as_advanced (FORCE SPHINX_EXECUTABLE)

if (SPHINX_EXECUTABLE)
  message ("-- Sphinx compiler: ${SPHINX_EXECUTABLE}")

  # Get the short version (not using REGEX_REPLACE for backwards compatibility)
  string (REPLACE "v" "" tmp "${MOLCAS_VERSION}")
  string (FIND "${tmp}" "-" pos)
  string (SUBSTRING "${tmp}" 0 ${pos} SHORT_VERSION)

  set (SPHINX_BUILDDIR "${PROJECT_BINARY_DIR}/doc")
  set (SPHINX_OPTIONS -d ${SPHINX_BUILDDIR}/doctrees -D extract_dir=${SPHINX_BUILDDIR}/samples)

  add_custom_target (doc_all)

  add_custom_target (doc_html
    # workaround for sphinx bug #5966
    COMMAND touch ${PROJECT_SOURCE_DIR}/doc/source/index.rst
    COMMAND ${CMAKE_COMMAND} -E env OPENMOLCAS_RELEASE=${MOLCAS_VERSION} OPENMOLCAS_VERSION=${SHORT_VERSION}
            ${SPHINX_EXECUTABLE} -b html ${SPHINX_OPTIONS} ${PROJECT_SOURCE_DIR}/doc/source ${SPHINX_BUILDDIR}/html
  )
  add_dependencies (doc_all doc_html)

  find_program (LATEXMK "latexmk")
  mark_as_advanced (FORCE LATEXMK)
  if (LATEXMK)
    message ("-- LaTeX compiler: ${LATEXMK}")

    add_custom_target (doc_pdf
      COMMAND ${CMAKE_COMMAND} -E env OPENMOLCAS_RELEASE=${MOLCAS_VERSION} OPENMOLCAS_VERSION=${SHORT_VERSION}
              ${SPHINX_EXECUTABLE} -b latex ${SPHINX_OPTIONS} ${PROJECT_SOURCE_DIR}/doc/source ${SPHINX_BUILDDIR}/latex
      COMMAND make -C ${SPHINX_BUILDDIR}/latex all-pdf
    )
    add_dependencies (doc_all doc_pdf)
  else ()
    message ("-- LaTeX compiler: no latexmk available, pdf documentation disabled")
  endif ()

  add_custom_target (doc)
  add_dependencies (doc doc_html)
else ()
  message ("-- Sphinx compiler: no sphinx-build available, documentation disabled")
endif ()

################################################################################
#                                                                              #
# Installation                                                                 #
#                                                                              #
################################################################################

message ("Install directory: ${CMAKE_INSTALL_PREFIX}")

set (CMAKE_SKIP_INSTALL_ALL_DEPENDENCY "TRUE")

install (FILES
  ${PROJECT_BINARY_DIR}/.molcashome
  ${PROJECT_BINARY_DIR}/.molcasversion
  ${PROJECT_BINARY_DIR}/molcas.rte
  ${PROJECT_BINARY_DIR}/LICENSE
  ${PROJECT_BINARY_DIR}/CONTRIBUTORS.md
  DESTINATION ${CMAKE_INSTALL_PREFIX}
)

install (DIRECTORY
  ${PROJECT_BINARY_DIR}/basis_library
  ${PROJECT_BINARY_DIR}/data
  DESTINATION ${CMAKE_INSTALL_PREFIX}
)

if (BUILD_SHARED_LIBS AND BUILD_STATIC_LIBS)
  set (subdir "/static")
else ()
  set (subdir "")
endif ()
install (TARGETS
  ${PROGRAM_EXECUTABLES}
  ${MOLCAS_LIBRARIES}
  OPTIONAL
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${subdir}
)

install (PROGRAMS
  ${PROJECT_BINARY_DIR}/sbin/help_basis
  ${PROJECT_BINARY_DIR}/sbin/help_func
  ${PROJECT_BINARY_DIR}/sbin/help_doc
  ${PROJECT_BINARY_DIR}/sbin/setup
  ${PROJECT_BINARY_DIR}/sbin/version
  DESTINATION ${CMAKE_INSTALL_PREFIX}/sbin
)

if (INSTALL_TESTS)
  install (PROGRAMS
    ${PROJECT_BINARY_DIR}/sbin/verify
    ${PROJECT_BINARY_DIR}/sbin/have_feature
    DESTINATION ${CMAKE_INSTALL_PREFIX}/sbin
  )

  install (DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/test)
  install (CODE "execute_process (COMMAND ${OPENMOLCAS_DIR}/sbin/install_tests.sh ${CMAKE_INSTALL_PREFIX}/test)")
endif ()

if (TARGET pymolcas_target)
  install (PROGRAMS
    ${PYMOLCAS_SCRIPT}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/
    OPTIONAL
  )
endif ()

install (DIRECTORY
  ${PROJECT_BINARY_DIR}/doc/html
  DESTINATION ${CMAKE_INSTALL_PREFIX}/doc
  OPTIONAL
)

install (FILES
  ${PROJECT_BINARY_DIR}/doc/latex/Manual.pdf
  DESTINATION ${CMAKE_INSTALL_PREFIX}/doc
  OPTIONAL
)

if (DEFINED EXTRA_DIR)
  install (PROGRAMS
    ${PROJECT_BINARY_DIR}/sbin/getname.plx
    ${PROJECT_BINARY_DIR}/sbin/moclas
    ${PROJECT_BINARY_DIR}/sbin/molcas.driver
    ${PROJECT_BINARY_DIR}/sbin/sew2xyz
    DESTINATION ${CMAKE_INSTALL_PREFIX}/sbin
  )

  install (PROGRAMS
    ${PROJECT_BINARY_DIR}/bin/molcas.exe
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  )
endif ()

if (DMRG)
  install (DIRECTORY
    ${PROJECT_BINARY_DIR}/qcmaquis
    USE_SOURCE_PERMISSIONS
    DESTINATION ${CMAKE_INSTALL_PREFIX}
  )
  install (DIRECTORY
    ${PROJECT_BINARY_DIR}/Tools/qcmaquis
    USE_SOURCE_PERMISSIONS
    DESTINATION ${CMAKE_INSTALL_PREFIX}/Tools
  )
endif ()
if (NEVPT2)
  install (DIRECTORY
    ${PROJECT_BINARY_DIR}/Tools/distributed-4rdm
    USE_SOURCE_PERMISSIONS
    DESTINATION ${CMAKE_INSTALL_PREFIX}/Tools
  )
endif ()

# Additional targets and directories to install (from external projects)
if (EXTRA_INSTALL_TARGETS)
  install (TARGETS
    ${EXTRA_INSTALL_TARGETS}
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include
    OPTIONAL
  )
endif ()

foreach (dir ${EXTRA_INSTALL_DIRS})
  install (DIRECTORY
    ${dir}
    USE_SOURCE_PERMISSIONS
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    OPTIONAL
  )
endforeach ()

foreach (lib ${EXTRA_INSTALL_IMPORTED_LIBS})
  get_target_property (fname ${lib} IMPORTED_LOCATION)
  install (FILES
    ${fname}
    TYPE LIB
    OPTIONAL
  )
endforeach ()

if (IS_DIRECTORY "$ENV{HOME}/.Molcas")
  install (CODE "execute_process (COMMAND sh \"-c\" \"cd ${CMAKE_INSTALL_PREFIX}; echo `pwd` > $HOME/.Molcas/molcas 2> /dev/null\")")
endif ()

# Final warning(s)
if (GEN1INT)
  message (WARNING "OpenMolcas is LGPL code, but Gen1Int's original license is GPL.\nThe generated binaries *cannot* be redistributed.")
endif()
