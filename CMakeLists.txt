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
#               2014-2019, Ignacio Fdez. Galván                        *
#               2015, Stefan Knecht                                    *
#***********************************************************************
# Global CMakeLists.txt file for building Molcas with CMake.
#
# Steven Vancoillie - Spring 2014, on a Saturday afternoon, with coffee
# Ignacio Fdez. Galván made several improvements / additions
# Lasse Sørensen - Mac CMP0042 fix

cmake_minimum_required (VERSION 2.8.11)
#set (CMAKE_VERBOSE_MAKEFILE true)
if (POLICY CMP0042)
        cmake_policy(SET CMP0042 NEW)
endif ()
if (POLICY CMP0045)
        cmake_policy(SET CMP0045 NEW)
endif ()
if (POLICY CMP0074)
        cmake_policy(SET CMP0074 NEW)
endif ()
if (POLICY CMP0004)
        cmake_policy(SET CMP0004 NEW)
endif ()

# Prevent in-source builds
# If an in-source build is attempted, you will still need to clean up a few files manually.
set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)
get_filename_component (sourcedir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component (binarydir "${CMAKE_BINARY_DIR}" REALPATH)
if ("${sourcedir}" STREQUAL "${binarydir}")
        message(FATAL_ERROR "In-source builds in ${CMAKE_BINARY_DIR} are not "
                "allowed, please remove ./CMakeCache.txt and ./CMakeFiles/, create a "
                "separate build directory and run cmake from there.")
endif ()

# workaround for CMake bug 14874
# (i.e., fix for files with same name in different directories)
if (${CMAKE_VERSION} VERSION_LESS 3.1)
        if (NOT DEFINED CMAKE_Fortran_ARCHIVE_CREATE)
                set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> cq <TARGET> <LINK_FLAGS> <OBJECTS>")
        endif ()
        if (NOT DEFINED CMAKE_Fortran_ARCHIVE_APPEND)
                set(CMAKE_Fortran_ARCHIVE_APPEND "<CMAKE_AR> q  <TARGET> <LINK_FLAGS> <OBJECTS>")
        endif ()
        if (NOT DEFINED CMAKE_C_ARCHIVE_CREATE)
                set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> cq <TARGET> <LINK_FLAGS> <OBJECTS>")
        endif ()
        if (NOT DEFINED CMAKE_C_ARCHIVE_APPEND)
                set(CMAKE_C_ARCHIVE_APPEND "<CMAKE_AR> q  <TARGET> <LINK_FLAGS> <OBJECTS>")
        endif ()
        if (NOT DEFINED CMAKE_CXX_ARCHIVE_CREATE)
                set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> cq <TARGET> <LINK_FLAGS> <OBJECTS>")
        endif ()
        if (NOT DEFINED CMAKE_CXX_ARCHIVE_APPEND)
                set(CMAKE_CXX_ARCHIVE_APPEND "<CMAKE_AR> q  <TARGET> <LINK_FLAGS> <OBJECTS>")
        endif ()
endif ()

# discover some necessary tools
find_program (GIT "git")
find_program (PERL "perl")
mark_as_advanced(FORCE GIT PERL)

FIND_PACKAGE(PythonInterp 3.0)
if (NOT PYTHONINTERP_FOUND)
        FIND_PACKAGE(PythonInterp 2.7)
endif ()

################################################################################
#                                                                              #
# Options                                                                      #
#                                                                              #
################################################################################

# The build type, specifies the different optimization types.
if (NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
        set (CMAKE_BUILD_TYPE "Release" CACHE STRING
                "Type of build, options are: None (CFLAGS/FFLAGS can be used), Debug, Garble, RelWithDebInfo, Release, Fast."
                FORCE
                )
endif ()
set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS None Debug Garble RelWithDebInfo Release Fast)

# The place where Molcas is installed with make install.
if (NOT DEFINED CMAKE_INSTALL_PREFIX)
        set (CMAKE_INSTALL_PREFIX "/opt/molcas" CACHE STRING
                "Location where Molcas will be installed.")
endif ()

option (DEBUGPRINT "Enable excessive print statements for debugging." OFF)
option (DEBUGTRACE "Enable tracing print statements for debugging." OFF)

set (DEBUG_DEFS_DESC "Additional debugging definitions.")
if (DEFINED DEBUG_DEFS)
        separate_arguments(DEBUG_DEFS)
        set (DEBUG_DEFS "${DEBUG_DEFS}" CACHE STRING "${DEBUG_DEFS_DESC}" FORCE)
else ()
        set (DEBUG_DEFS "" CACHE STRING "${DEBUG_DEFS_DESC}")
endif ()

option (BUILD_SHARED_LIBS "Build dynamically-linked molcas library." OFF)
option (BUILD_STATIC_LIB "Build statically-linked molcas library too." OFF)
if (BUILD_SHARED_LIBS)
        mark_as_advanced (CLEAR BUILD_STATIC_LIB)
else ()
        mark_as_advanced (FORCE BUILD_STATIC_LIB)
endif ()

option (MPI "Enable MPI parallellization." OFF)
option (GA "Use Global Arrays library." OFF)

option (OPENMP "Enable multi-threading." OFF)

option (BOUNDS "Enable bounds checking (only gfortran >= 4.8)." OFF)

option (EXPERT "Show advanced CMake cache entries in GUI)." OFF)

option (BIGOT "Do not allow any compiler warning (treat them as errors)." OFF)

set (LINALG_OPTIONS "Choose linear algebra library, options: Internal (default), Runtime, MKL, OpenBLAS, ACML, Accelerate.")
if (DEFINED LINALG)
        set (LINALG "${LINALG}" CACHE STRING "${LINALG_OPTIONS}")
        if (NOT ";Internal;Runtime;MKL;OpenBLAS;ACML;Accelerate;" MATCHES ";${LINALG};")
                message (FATAL_ERROR "Unsuported LINALG: ${LINALG}")
        endif ()
else ()
        set (LINALG "Internal"  CACHE STRING "${LINALG_OPTIONS}")
endif ()
set_property (CACHE LINALG PROPERTY STRINGS Internal Runtime MKL OpenBLAS ACML Accelerate)

option (CUBLAS "Enable CUDA BLAS library." OFF)
option (NVBLAS "Enable NVidia BLAS library." OFF)

option (FDE "Enable Frozen-density-embedding (FDE) interface." OFF)
option (GROMACS "Compile Gromacs interface." OFF)
option (BLOCK "Activate BLOCK-DMRG support." OFF)
option (CHEMPS2 "Activate CheMPS2-DMRG support." OFF)

option (GPERFTOOLS "Activate gperftools CPU profiling." OFF)
option (GCOV "Activate code coverage profiling (use with Debug/-O0)." OFF)
option (HDF5 "Activate HDF5 support for wavefunction format." ON)

option (TOOLS "Compile supported tools." OFF)

option (CTEST "Enable CTest." OFF)
mark_as_advanced (CTEST)

# External projects
option (MSYM       "Activate MSYM support (requires External/libmsym submodule)."                               OFF)
option (DMRG       "Activate QCMaquis DMRG driver and library."                                                 OFF)
option (DMRG_DEBUG "Enable DEBUG print in QCMaquis driver."                                                     OFF)

if (DMRG)
        mark_as_advanced (CLEAR DMRG_DEBUG)
        # Name and e-mail required for downloading QCMaquis
        # Try to get the real name from the system account information
        if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
                set(QCMaquis_NAME_GET_CMD "id -F")
        else()
        # should we distinguish other platforms?
                set(QCMaquis_NAME_GET_CMD "sh" "-c" "getent passwd $USER | awk -F: '{print \$5}'")
        endif()
        execute_process (
                COMMAND ${QCMaquis_NAME_GET_CMD}
                OUTPUT_VARIABLE QCMaquis_NAME_
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set (QCMaquis_NAME ${QCMaquis_NAME_} CACHE STRING "Your name (required for downloading QCMaquis)")

        # Update QCMaquis flag. Off by default.
        # QCMaquis must be updated manually, as there might be QCMaquis upstream changes even when OpenMolcas is not changed,
        # or vice versa. Running cmake with QCMaquis_UPDATE=ON will force-update QCMaquis.
        if (NOT QCMaquis_UPDATE)
                set (QCMaquis_UPDATE OFF CACHE STRING "Trigger QCMaquis update" FORCE)
        endif()

        # try to get the e-mail address from git config (TODO: is there a better option?)
        execute_process (
                COMMAND ${GIT} "config" "user.email"
                OUTPUT_VARIABLE QCMaquis_EMAIL_
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        set (QCMaquis_EMAIL ${QCMaquis_EMAIL_} CACHE STRING "Your e-mail address (required for downloading QCMaquis)")

        # Prepare variables for QCMaquis: separate name and surname
        set (QCMaquis_NAME_LIST ${QCMaquis_NAME})

        separate_arguments (QCMaquis_NAME_LIST)
        list(LENGTH QCMaquis_NAME_LIST qcm_name_len)
        if (qcm_name_len LESS 2)
                message (FATAL_ERROR
                        "Your name has been specified as: " ${QCMaquis_NAME} "\n"
                        "Please specify both name and surname. If you have only one name, please specify it twice.")
        endif()

        # Separate first name from the last name
        # Currently, only the first word is treated by QCMaquis as the first name,
        # all the rest is considered surnames -- to accomodate multiple surnames common in the Spanish-speaking
        # areas or suffixes as Jr., III etc. Middle names are therefore incorrectly considered parts of
        # the last name. However, this has no influence on the QCMaquis download.
        list(GET QCMaquis_NAME_LIST 0 QCMaquis_GIVEN_NAME)
        list(REMOVE_AT QCMaquis_NAME_LIST 0)
        foreach(arg ${QCMaquis_NAME_LIST})
                      set(QCMaquis_SURNAME "${QCMaquis_SURNAME} ${arg}")
        endforeach()

        # Set QCMaquis download URL and path to the download script
        set (QCMaquis_URL "https://scine.ethz.ch/download/")
        set (QCMaquis_dl_SCRIPT "${CMAKE_SOURCE_DIR}/src/qcmaquis/get_qcmaquis.sh")
else ()
        mark_as_advanced (FORCE DMRG_DEBUG)
endif ()
option (GEN1INT    "Enable Gen1Int library for 1-electron integrals."                                           OFF)
option (NECI       "Activate NECI support (requires External/NECI submodule)."                                  OFF)
option (WFA        "Activate extended wavefunction analysis (requires External/libwfa submodule)."              OFF)
option (NEVPT2     "Activate (DMRG)-NEVPT2 support."                                                            OFF)
option (EFPLIB     "Enable EFPLib library for effective fragment potentials (requires External/efp submodule)." OFF)
option (GRID_IT    "Enable public GRID_IT code (requires External/grid_it submodule)."                          OFF)

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

if (EFPLIB)
        cmake_minimum_required (VERSION 3.0.2)
endif ()
if (DMRG)
        project (Molcas Fortran C CXX)
        cmake_minimum_required (VERSION 3.2.0)
        if (NOT DEFINED QCMaquis_ROOT)
                if (DEFINED ENV{QCMaquis_ROOT})
                        set (QCMaquis_ROOT $ENV{QCMaquis_ROOT})
                else ()
                        set (QCMaquis_ROOT "None")
                endif ()
                set (QCMaquis_ROOT ${QCMaquis_ROOT} CACHE STRING "Specify path to existing QCMaquis installation. This will allow to skip the installation of ALPS/BOOST and QCMaquis.")
        endif ()
elseif (WFA)
        project (Molcas Fortran C CXX)
else ()
        project (Molcas Fortran C)
endif ()

#avoid using CXX compiler for linking
set (CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)

add_definitions (-D_MOLCAS_)

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
                message (WARNING "The Molcas-extra directory ${EXTRA} does not exist.")
        else ()
                FILE (STRINGS "${EXTRA}/.molcashome" DIR_TYPE)
                if ("${DIR_TYPE}" STREQUAL "molcas-extra")
                        set (EXTRA_DIR ${EXTRA})
                        add_definitions (-D_HAVE_EXTRA_)
                else ()
                        message (FATAL_ERROR "The directory ${EXTRA} does not contain Molcas-extra.")
                endif ()
        endif ()
endif ()
set (EXTRA "${EXTRA}" CACHE PATH "Location of the Molcas-extra directory." FORCE)

# possible locations for source code
set (OPENMOLCAS_DIR ${CMAKE_CURRENT_LIST_DIR})
if (NOT ${PROJECT_SOURCE_DIR} STREQUAL ${OPENMOLCAS_DIR})
        set (basedirs ${PROJECT_SOURCE_DIR})
        set (EXTRA_DIR ${PROJECT_SOURCE_DIR})
        add_definitions (-D_HAVE_EXTRA_)
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

# set function to find source directories
#========================================

function (find_source name)
        foreach (src ${source_roots})
                if (EXISTS ${PROJECT_SOURCE_DIR}/${src}/${name})
                        file (GLOB tmplist RELATIVE ${PROJECT_SOURCE_DIR}/${src}/${name} ${PROJECT_SOURCE_DIR}/${src}/${name}/*)
                        # ignore possible leftover files from previous configure+make
                        list (REMOVE_ITEM tmplist Makefile 00dependencies 00sources)
                        if (tmplist)
                                set (found "${src}/${name}")
                                break ()
                        endif ()
                endif ()
        endforeach ()
        # do not overwrite existing definitions
        if (DEFINED ${name}_src)
                set (${name}_src "${${name}_src}" PARENT_SCOPE)
        else ()
                if (DEFINED found)
                        set (${name}_src "${found}" PARENT_SCOPE)
                #else ()
                #        message (WARNING "${name} not found, expect trouble")
                endif ()
        endif ()
endfunction ()

# set function to insert list
#============================

function (insert_before _list _item _new)
        list (GET ${_item} 0 _first)
        list (FIND ${_list} ${_first} _index)
        if (_index GREATER -1)
                list (INSERT ${_list} ${_index} ${_new})
        else ()
                list (APPEND ${_list} ${_new})
        endif ()
        set (${_list} ${${_list}} PARENT_SCOPE)
endfunction ()

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
                # workaround because "continue()" is not supported in CMake 2.8.11
                set (do ON)
                if ("${BASE_DIR}" STREQUAL "${OPENMOLCAS_DIR}")
                        if (EXISTS ${OPENMOLCAS_GIT_REPO})
                                set (MOLCAS_GIT_REPO ${OPENMOLCAS_GIT_REPO})
                        else ()
                                set (do OFF)
                        endif ()
                        set (label "")
                elseif ("${BASE_DIR}" STREQUAL "${EXTRA_DIR}")
                        if (EXISTS ${EXTRA_GIT_REPO})
                                set (MOLCAS_GIT_REPO ${EXTRA_GIT_REPO})
                        else ()
                                set (do OFF)
                        endif ()
                        set (label "-extra")
                endif ()
                if (do)
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
                string(STRIP ${GIT_BRANCH} GIT_BRANCH)
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
        add_definitions (-D_I8_)
else ()
        set (ADDRMODE 32)
endif ()

message ("-- ADDRMODE: ${ADDRMODE}")

# platform settings
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        add_definitions (-D_LINUX_)
        if ((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86") OR
            (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i686"))
                set (PLATFORM "LINUX")
        elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
                set (PLATFORM "LINUX64")
        elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ia64")
                set (PLATFORM "LINUX64_IA")
        elseif ((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc_64") OR
                (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64le"))
                set (PLATFORM "PPC64")
        endif ()
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        add_definitions (-D_LINUX_ -D_DARWIN_)
        set (CMAKE_MACOSX_RPATH 0)
        set (PLATFORM "MacOS")
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
        add_definitions (-D_LINUX_ -D_CYGWIN_)
        if (ADDRMODE EQUAL 32)
                set (PLATFORM "WIN")
        elseif (ADDRMODE EQUAL 64)
                set (PLATFORM "WIN64")
        endif ()
endif ()

if (NOT PLATFORM)
        message (FATAL_ERROR "unsupported platform")
else ()
        message ("-- PLATFORM: ${PLATFORM}")
endif ()

################################################################################
#                                                                              #
# Compiler-specific settings                                                   #
#                                                                              #
################################################################################

# Set the default compiler flags to pick from
#================================================
# For new compilers, check for the proper compiler ID on the following web page:
# http://www.cmake.org/cmake/help/v2.8.12/cmake.html#variable:CMAKE_LANG_COMPILER_ID
# The regular compiler options have to be compatible with other flags (so they
# might be added or removed depending on certain options), while the flags for
# build targets are mutually exclusive (only one of them will be active).

# Ideally, build targets should accomplish the features in the table below.
# legend: opt = optimizations; bt = backtrace (call stack); trap = catch
# floating point exceptions; warn = produce warnings during compilation
#                    opt  bt   trap  warn
#                    ===  ===  ====  ====
# - Debug            NO   YES  NO    YES
# - Garble           YES  YES  YES   YES
# - RelWithDebInfo   YES  YES  NO    YES
# - Release          YES  NO   NO    NO
# - Fast             YES  NO   NO    NO
#
# Note: trapping floating point exceptions is not compatible with LAPACK
# (which produces and handles the exceptions), so a better alternative
# has to be found

#======================================#
# GNU Compiler Collection (GCC)
#======================================#
# defines __GNUC__
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "GNU")
set (CXXFLAGS_GNU_BASIC "")
set (CXXFLAGS_GNU_OPENMP "-fopenmp")
# build targets
set (CXXFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELEASE "-O2")
set (CXXFLAGS_GNU_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "GNU")
set (CFLAGS_GNU_BASIC "-std=gnu99")
set (CFLAGS_GNU_OPENMP "-fopenmp")
# build targets
set (CFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_GNU_RELEASE "-O2")
set (CFLAGS_GNU_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "GNU")
set (FFLAGS_GNU_BASIC "-cpp")
set (FFLAGS_GNU_OPENMP "-fopenmp")
set (FFLAGS_GNU_ILP64 "-fdefault-integer-8")
# build targets
set (FFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (FFLAGS_GNU_GARBLE "-O2 -g -Wall -finit-real=snan -finit-integer=730432726 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
set (FFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (FFLAGS_GNU_RELEASE "-O2")
set (FFLAGS_GNU_FAST "-O3")

# Fix for GNU Fortran compiler version >= 4.8 to prevent wrong optimization
# version check for Fortran doesn't work yet, assume it is the same as the C compiler GCC version
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
    (NOT CMAKE_C_COMPILER_VERSION VERSION_LESS "4.8"))
        set (FFLAGS_GNU_BASIC "${FFLAGS_GNU_BASIC} -fno-aggressive-loop-optimizations")
        set (FFLAGS_GNU_BOUNDS "-fsanitize=address -fno-omit-frame-pointer")
endif ()

set (CXXFLAGS_GNU_BIGOT "-Wall -Wextra -Wpedantic -Werror")
set (CFLAGS_GNU_BIGOT "-Wall -Wextra -Wpedantic -Werror")
set (FFLAGS_GNU_BIGOT "-Wall -Werror")
set (GNU_TRAP "-ffpe-trap=invalid,zero,overflow")

#======================================#
# LLVM C frontend (Clang)
#======================================#
# defines __clang__
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Clang")
set (CXXFLAGS_Clang_BASIC "-std=gnu99")
set (CXXFLAGS_Clang_OPENMP "-fopenmp")
# build targets
set (CXXFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELEASE "-O2")
set (CXXFLAGS_Clang_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Clang")
set (CFLAGS_Clang_BASIC "-std=gnu99")
set (CFLAGS_Clang_OPENMP "-fopenmp")
# build targets
set (CFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_Clang_RELEASE "-O2")
set (CFLAGS_Clang_FAST "-O2")

set (CXXFLAGS_Clang_BIGOT "-Wall -Werror")
set (CFLAGS_Clang_BIGOT "-Wall -Werror")
#set (Clang_TRAP "???")

#======================================#
# Intel C++/C/Fortran Compilers
#======================================#
# defines __INTEL_COMPILER
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Intel")
set (CXXFLAGS_Intel_BASIC "")
set (CXXFLAGS_Intel_OPENMP "-qopenmp")
# build targets
set (CXXFLAGS_Intel_DEBUG "-debug -w3")
set (CXXFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELEASE "-O2")
set (CXXFLAGS_Intel_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Intel")
set (CFLAGS_Intel_BASIC "-std=gnu99")
set (CFLAGS_Intel_OPENMP "-qopenmp")
# build targets
set (CFLAGS_Intel_DEBUG "-debug -w3")
set (CFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CFLAGS_Intel_RELEASE "-O2")
set (CFLAGS_Intel_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "Intel")
set (FFLAGS_Intel_BASIC "-fpp")
set (FFLAGS_Intel_OPENMP "-qopenmp")
set (FFLAGS_Intel_ILP64 "-i8 -r8 -heap-arrays")
# build targets
set (FFLAGS_Intel_DEBUG "-debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_GARBLE "-O2 -debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_RELEASE "-O2 -traceback")
set (FFLAGS_Intel_FAST "-fast")

# Intel versions prior to 15 used -openmp
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
        set (CXXFLAGS_Intel_OPENMP "-openmp")
endif ()
if (CMAKE_C_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_C_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
        set (CFLAGS_Intel_OPENMP "-openmp")
endif ()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
        set (FFLAGS_Intel_OPENMP "-openmp")
endif ()

#fix for MAC OS X
if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
set (CXXFLAGS_Intel_BASIC "-headerpad_max_install_names ${CXXFLAGS_Intel_BASIC}")
set (CFLAGS_Intel_BASIC "-headerpad_max_install_names ${CFLAGS_Intel_BASIC}")
set (FFLAGS_Intel_BASIC "-headerpad_max_install_names ${FFLAGS_Intel_BASIC}")
endif ()

#note that "declarations" means "implicit none", which we don't want (yet)
set (CXXFLAGS_Intel_BIGOT "-w3 -Werror -diag-disable remark")
set (CFLAGS_Intel_BIGOT "-w3 -Werror -diag-disable remark")
set (FFLAGS_Intel_BIGOT "-warn all,nodeclarations -warn errors -diag-disable remark")
set (Intel_TRAP "-fpe0")

#======================================#
# Oracle Developer Studio (Sun Studio)
#======================================#
# defines __SUNPRO_F90
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "SunPro")
set (CFLAGS_SunPro_BASIC "-std=gnu99")
set (CFLAGS_SunPro_OPENMP "-xopenmp")
# build targets
set (CFLAGS_SunPro_DEBUG "-g -traceback -ftrap=%none")
set (CFLAGS_SunPro_GARBLE "-O -g -traceback -ftrap=%none -xcheck=init_local")
set (CFLAGS_SunPro_RELWITHDEBINFO "-O -g -traceback -ftrap=%none")
set (CFLAGS_SunPro_RELEASE "-O -ftrap=%none")
set (CFLAGS_SunPro_FAST "-O -ftrap=%none")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "SunPro")
set (FFLAGS_SunPro_BASIC "-fpp") # 32bit: -ftrap=common,no%division,%none
set (FFLAGS_SunPro_OPENMP "-xopenmp")
set (FFLAGS_SunPro_ILP64 "-dalign -xtypemap=real:64,double:64,integer:64")
# build targets
set (FFLAGS_SunPro_DEBUG "-g -traceback -ftrap=%none")
set (FFLAGS_SunPro_GARBLE "-O -g -traceback ftrap=%none -xcheck=init_local")
set (FFLAGS_SunPro_RELWITHDEBINFO "-O -g -traceback -ftrap=%none")
set (FFLAGS_SunPro_RELEASE "-O -ftrap=%none")
set (FFLAGS_SunPro_FAST "-fast -dbl_align_all=yes -ftrap=%none")

# older versions do not suport gnu99
if (CMAKE_C_COMPILER_ID STREQUAL "SunPro" AND
    CMAKE_C_COMPILER_VERSION VERSION_LESS "5.15.0")
        set (CFLAGS_SunPro_BASIC "-std=c99")
endif ()

set (CFLAGS_SunPro_BIGOT "-errwarn=%all")
set (FFLAGS_SunPro_BIGOT "-errwarn=%all")
set (SunPro_TRAP "-ftrap=common")

#======================================#
# Portland Group (PGI) Compilers
#======================================#
# defines __PGI
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "PGI")
set (CXXFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_PGI_OPENMP "-mp")
# build targets
set (CXXFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELEASE "-O2")
set (CXXFLAGS_PGI_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "PGI")
set (CFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CFLAGS_PGI_OPENMP "-mp")
# build targets
set (CFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELEASE "-O2")
set (CFLAGS_PGI_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "PGI")
set (FFLAGS_PGI_BASIC "-Mpreprocess -Mbackslash")
set (FFLAGS_PGI_OPENMP "-mp")
set (FFLAGS_PGI_ILP64 "-i8")
# build targets
set (FFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELEASE "-O2")
set (FFLAGS_PGI_FAST "-fast")

# The option does not exist, or I cannot find it
#set (CXXFLAGS_PGI_BIGOT "???")
#set (CFLAGS_PGI_BIGOT "???")
#set (FFLAGS_PGI_BIGOT "???")
#set (PGI_TRAP "???")

#======================================#
# Numerical Algorithms Group (NAG)
#======================================#
# defines NAGFOR
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "NAG")
set (FFLAGS_NAG_BASIC "-fpp -dusty -dcfuns -mismatch_all -w -kind=byte")
set (FFLAGS_NAG_OPENMP "-openmp")
set (FFLAGS_NAG_ILP64 "-i8")
# build targets
set (FFLAGS_NAG_DEBUG "-O0 -g -ieee=full")
set (FFLAGS_NAG_GARBLE "-O2 -g -ieee=full -nan")
set (FFLAGS_NAG_RELWITHDEBINFO "-O2 -g -ieee=full")
set (FFLAGS_NAG_RELEASE "-O2 -ieee=full")
set (FFLAGS_NAG_FAST "-O3 -ieee=full")

# The option does not exist, or I cannot find it
#set (FFLAGS_NAG_BIGOT "???")
set(NAG_TRAP "-ieee=stop")

# Generate the actual flags used in the build
#============================================

# Preprocessor
if (CMAKE_BUILD_TYPE STREQUAL "Garble")
        add_definitions (-D_GARBLE_ -D_NO_IEEECK_ -D_FPE_TRAP_)
        # Disable FPE traps for compilers that do not support IEEE Fortran modules
        # version check for Fortran doesn't work yet, assume it is the same as the C compiler GCC version
        if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND
            CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
            CMAKE_C_COMPILER_VERSION VERSION_LESS "5.0")
                remove_definitions (-D_FPE_TRAP_)
        elseif (NOT LINALG STREQUAL "Internal")
                message (WARNING
                        "Garble is active, but an external LINALG library is used! "
                        "The use of an external BLAS/LAPACK can result in floating point "
                        "exceptions handled by the library being trapped instead. "
                        "Please use -DLINALG=Internal for better chances for Garble to work correctly!")
        else ()
                message (WARNING
                        "Garble and floating point exception trapping are active! "
                        "It is possible that intended exceptions within BLAS/LAPACK crash the program. "
                        "Be warned!")
        endif ()
endif ()

if (DMRG OR WFA)
        # C++ compiler required for DMRG and WFA projects
        foreach (CXX_COMPILER ${SUPPORTED_CXX_COMPILERS})
                if (CXX_COMPILER STREQUAL "${CMAKE_CXX_COMPILER_ID}")
                        set (CXXCID ${CMAKE_CXX_COMPILER_ID})
                endif ()
        endforeach ()

        if (NOT CXXCID)
                message (FATAL_ERROR "unsupported CXX compiler: ${CMAKE_CXX_COMPILER_ID}")
        endif ()
endif ()

# C compiler
foreach (C_COMPILER ${SUPPORTED_C_COMPILERS})
        if (C_COMPILER STREQUAL "${CMAKE_C_COMPILER_ID}")
                set (CCID ${CMAKE_C_COMPILER_ID})
        endif ()
endforeach ()

if (NOT CCID)
        message (FATAL_ERROR "unsupported C compiler: ${CMAKE_C_COMPILER_ID}")
endif ()

# overwrite any CMake defaults (some are strange, so we set our own)
if (NOT CMAKE_CXX_FLAGS_USER)
        set (CMAKE_CXX_FLAGS_DEFAULT "$ENV{CXXFLAGS} ${CXXFLAGS_${CXXCID}_BASIC}"
                CACHE STRING "C++ compiler flags." FORCE)
        set (CMAKE_CXX_FLAGS_DEBUG ${CXXFLAGS_${CXXCID}_DEBUG}
                CACHE STRING "C++ compiler flags used for debug builds." FORCE)
        set (CMAKE_CXX_FLAGS_GARBLE ${CXXFLAGS_${CXXCID}_GARBLE}
                CACHE STRING "C++ compiler flags used for garble builds." FORCE)
        set (CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CXXFLAGS_${CXXCID}_RELWITHDEBINFO}
                CACHE STRING "C++ compiler flags used for release builds with debug info." FORCE)
        set (CMAKE_CXX_FLAGS_RELEASE ${CXXFLAGS_${CXXCID}_RELEASE}
                CACHE STRING "C++ compiler flags used for release builds." FORCE)
        set (CMAKE_CXX_FLAGS_FAST ${CXXFLAGS_${CXXCID}_FAST}
                CACHE STRING "C++ compiler flags used for fast release builds." FORCE)
        set (CMAKE_CXX_FLAGS_USER "TRUE"
                CACHE INTERNAL "set to FALSE to reset to default C++ flags" FORCE)
endif ()
if (EXPERT)
        mark_as_advanced(CLEAR
                CMAKE_CXX_FLAGS_DEFAULT
                CMAKE_CXX_FLAGS_DEBUG
                CMAKE_CXX_FLAGS_GARBLE
                CMAKE_CXX_FLAGS_RELWITHDEBINFO
                CMAKE_CXX_FLAGS_RELEASE
                CMAKE_CXX_FLAGS_FAST
                )
else ()
        mark_as_advanced(FORCE
                CMAKE_CXX_FLAGS_DEFAULT
                CMAKE_CXX_FLAGS_DEBUG
                CMAKE_CXX_FLAGS_GARBLE
                CMAKE_CXX_FLAGS_RELWITHDEBINFO
                CMAKE_CXX_FLAGS_RELEASE
                CMAKE_CXX_FLAGS_FAST
                )
endif ()
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEFAULT}")

# overwrite any CMake defaults (some are strange, so we set our own)
if (NOT CMAKE_C_FLAGS_USER)
        set (CMAKE_C_FLAGS_DEFAULT "$ENV{CFLAGS} ${CFLAGS_${CCID}_BASIC}"
                CACHE STRING "C compiler flags." FORCE)
        set (CMAKE_C_FLAGS_DEBUG ${CFLAGS_${CCID}_DEBUG}
                CACHE STRING "C compiler flags used for debug builds." FORCE)
        set (CMAKE_C_FLAGS_GARBLE ${CFLAGS_${CCID}_GARBLE}
                CACHE STRING "C compiler flags used for garble builds." FORCE)
        set (CMAKE_C_FLAGS_RELWITHDEBINFO ${CFLAGS_${CCID}_RELWITHDEBINFO}
                CACHE STRING "C compiler flags used for release builds with debug info." FORCE)
        set (CMAKE_C_FLAGS_RELEASE ${CFLAGS_${CCID}_RELEASE}
                CACHE STRING "C compiler flags used for release builds." FORCE)
        set (CMAKE_C_FLAGS_FAST ${CFLAGS_${CCID}_FAST}
                CACHE STRING "C compiler flags used for fast release builds." FORCE)
        set (CMAKE_C_FLAGS_USER "TRUE"
                CACHE INTERNAL "set to FALSE to reset to default C flags" FORCE)
endif ()
if (EXPERT)
        mark_as_advanced(CLEAR
                CMAKE_C_FLAGS_DEFAULT
                CMAKE_C_FLAGS_DEBUG
                CMAKE_C_FLAGS_GARBLE
                CMAKE_C_FLAGS_RELWITHDEBINFO
                CMAKE_C_FLAGS_RELEASE
                CMAKE_C_FLAGS_FAST
                )
else ()
        mark_as_advanced(FORCE
                CMAKE_C_FLAGS_DEFAULT
                CMAKE_C_FLAGS_DEBUG
                CMAKE_C_FLAGS_GARBLE
                CMAKE_C_FLAGS_RELWITHDEBINFO
                CMAKE_C_FLAGS_RELEASE
                CMAKE_C_FLAGS_FAST
                )
endif ()
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEFAULT}")

# Fortran compiler
foreach (Fortran_COMPILER ${SUPPORTED_Fortran_COMPILERS})
        if (Fortran_COMPILER STREQUAL "${CMAKE_Fortran_COMPILER_ID}")
                set (FCID ${CMAKE_Fortran_COMPILER_ID})
        endif ()
endforeach ()

if (NOT FCID)
        message (FATAL_ERROR "unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif ()

# overwrite any CMake defaults (some are strange, so we set our own)
if (NOT CMAKE_Fortran_FLAGS_USER)
        set (CMAKE_Fortran_FLAGS_DEFAULT "$ENV{FFLAGS} ${FFLAGS_${FCID}_BASIC}"
                CACHE STRING "Fortran compiler flags." FORCE)
        set (CMAKE_Fortran_FLAGS_DEBUG ${FFLAGS_${FCID}_DEBUG}
                CACHE STRING "Fortran compiler flags used for debug builds." FORCE)
        set (CMAKE_Fortran_FLAGS_GARBLE ${FFLAGS_${FCID}_GARBLE}
                CACHE STRING "Fortran compiler flags used for garble builds." FORCE)
        set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO ${FFLAGS_${FCID}_RELWITHDEBINFO}
                CACHE STRING "Fortran compiler flags used for release builds with debug info." FORCE)
        set (CMAKE_Fortran_FLAGS_RELEASE ${FFLAGS_${FCID}_RELEASE}
                CACHE STRING "Fortran compiler flags used for release builds." FORCE)
        set (CMAKE_Fortran_FLAGS_FAST ${FFLAGS_${FCID}_FAST}
                CACHE STRING "Fortran compiler flags used for fast release builds." FORCE)
        set (CMAKE_Fortran_FLAGS_USER "TRUE"
                CACHE INTERNAL "set to FALSE to reset to default Fortran flags" FORCE)
endif ()
if (EXPERT)
        mark_as_advanced(CLEAR
                CMAKE_Fortran_FLAGS_DEFAULT
                CMAKE_Fortran_FLAGS_DEBUG
                CMAKE_Fortran_FLAGS_GARBLE
                CMAKE_Fortran_FLAGS_RELWITHDEBINFO
                CMAKE_Fortran_FLAGS_RELEASE
                CMAKE_Fortran_FLAGS_FAST
                )
else ()
        mark_as_advanced(FORCE
                CMAKE_Fortran_FLAGS_DEFAULT
                CMAKE_Fortran_FLAGS_DEBUG
                CMAKE_Fortran_FLAGS_GARBLE
                CMAKE_Fortran_FLAGS_RELWITHDEBINFO
                CMAKE_Fortran_FLAGS_RELEASE
                CMAKE_Fortran_FLAGS_FAST
                )
endif ()
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEFAULT}")

if (ADDRMODE EQUAL 64)
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ILP64}")
endif ()

if (OPENMP)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_OPENMP}")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_OPENMP}")
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP}")
endif ()

if (BOUNDS)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_BOUNDS}")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_BOUNDS}")
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_BOUNDS}")
        add_definitions (-D_MALLOC_INTERCEPT_)
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

# Install directory
#==================
if (NOT IS_ABSOLUTE "${CMAKE_INSTALL_PREFIX}")
        set (CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}")
endif ()
get_filename_component (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)

# Special print statements
#=========================

if (DEBUGPRINT)
        add_definitions (-D_DEBUG_)
endif ()
if (DEBUGTRACE)
        add_definitions (-D_DEBUG_TRACE_)
endif ()

# Global Arrays settings
#=======================

# Global Arrays settings, non-portable!
if (GA AND MPI)
        message ("Configuring with external GA library:")

        find_path (GA_INCLUDE_PATH ga.h
                PATHS ENV GAROOT
                PATH_SUFFIXES include
                NO_DEFAULT_PATH)
        if (NOT GA_INCLUDE_PATH)
                message (FATAL_ERROR "GA_INCLUDE_PATH not found, make sure GAROOT environment variable is set.")
        endif ()
        find_library (LIBGA ga
                PATHS ENV GAROOT
                PATH_SUFFIXES lib
                NO_DEFAULT_PATH)
        if (NOT LIBGA)
                message (FATAL_ERROR "LIBGA not found, make sure GAROOT environment variable is set.")
        else ()
                list (APPEND GA_LIBRARIES ${LIBGA})
        endif ()
        find_library (LIBARMCI armci
                PATHS ENV GAROOT
                PATH_SUFFIXES lib
                NO_DEFAULT_PATH)
        if (NOT LIBARMCI)
                message (FATAL_ERROR "LIBARMCI not found, make sure GAROOT environment variable is set.")
        else ()
                list (APPEND GA_LIBRARIES ${LIBARMCI})
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
                        list (APPEND GA_LIBRARIES ${LIBIBVERBS})
                endif ()
        endif ()

        add_definitions (-D_GA_)
        include_directories (${GA_INCLUDE_PATH})

        message ("-- GA_INCLUDE_PATH: ${GA_INCLUDE_PATH}")
        message ("-- GA_LIBRARIES: ${GA_LIBRARIES}")

        list (APPEND EXTERNAL_LIBRARIES ${GA_LIBRARIES})
endif ()

# Parallel settings
#==================

if (MPI)
        message ("Configuring with MPI parallellization:")
        find_package (MPI REQUIRED)

        add_definitions (-D_MOLCAS_MPP_)
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
                message (FATAL_ERROR "Unknown or unsupported MPI implementation:\n${MPIEXEC_VERSION}")
        endif ()
        message ("-- MPI_IMPLEMENTATION: ${MPI_IMPLEMENTATION}")

        list (APPEND EXTERNAL_LIBRARIES ${MPI_Fortran_LIBRARIES})
endif ()

# Support for HDF5 I/O library
#=============================

if (HDF5)
        if (ADDRMODE EQUAL 32)
                message ("-- HDF5 interface not compatible with 32bit installations")
        else()
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
                        list (APPEND EXTERNAL_LIBRARIES ${HDF5_C_LIBRARIES})
                        if (HDF5_CXX_LIBRARIES)
                                message ("-- HDF5_CXX_LIBRARIES: ${HDF5_CXX_LIBRARIES}")
                                list (APPEND EXTERNAL_LIBRARIES ${HDF5_CXX_LIBRARIES})
                        endif ()
                else()
                        message ("-- HDF5 not found, it will be deactivated")
                endif()
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
                message (FATAL_ERROR "Use of MKL for 32bit installation not implemented")
        endif ()

        mark_as_advanced (CLEAR MKLROOT)
        if (MKLROOT STREQUAL "")
                set (MKLROOT $ENV{MKLROOT} CACHE PATH "MKL root directory." FORCE)
                if (NOT MKLROOT)
                        message (FATAL_ERROR
                                "You must set environment variable MKLROOT, "
                                "or specify -DMKLROOT=/path/to/mkl_root_dir "
                                "when running cmake."
                                )
                endif ()
        endif ()
        # at this point, MKLROOT should be defined and not empty
        message ("-- MKLROOT = ${MKLROOT}")

        # here we check if MKLROOT has changed after it was last used
        # succesfully, if so then we need to reprocess the checks here.
        if (NOT MKLROOT STREQUAL MKLROOT_LAST)
                # reset MKL paths
                set (MKL_INCLUDE_PATH "${MKLROOT}/include" CACHE PATH
                        "location of MKL include files." FORCE)
                if (ADDRMODE EQUAL 64)
                        set (libpath "${MKLROOT}/lib/intel64")
                elseif (ADDRMODE EQUAL 32)
                        set (libpath "${MKLROOT}/lib/ia32")
                endif ()
                set (MKL_LIBRARY_PATH ${libpath} CACHE PATH
                        "location of MKL libraries." FORCE)
                # uncache variables
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
        find_library (LIBMKL_CORE NAMES "mkl_core"
                PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        # compiler-specific library interface
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
                find_library (LIBMKL_INTERFACE NAMES "mkl_gf_ilp64"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
                find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI") # apparently PGI uses this too
                find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        endif ()
        # sequential/compiler-specific threading interface
        find_library (LIBMKL_SEQUENTIAL NAMES "mkl_sequential"
                PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
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
                find_library (LIBMKL_BLACS NAMES "mkl_blacs_ilp64"
                        PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
        endif ()

        # generate actual library list with paths
        if (MPI)
                add_definitions (-DSCALAPACK)
                list (APPEND MKL_LIBRARIES ${LIBMKL_SCALAPACK})
        endif ()
        list (APPEND MKL_LIBRARIES ${LIBMKL_INTERFACE})
        list (APPEND MKL_LIBRARIES ${LIBMKL_CORE})
        if (OPENMP)
                set (CMAKE_EXE_LINKER_FLAGS ${LIBMKL_OMP_LINK_FLAGS})
                list (APPEND MKL_LIBRARIES ${LIBMKL_THREADING})
        else ()
                list (APPEND MKL_LIBRARIES ${LIBMKL_SEQUENTIAL})
        endif ()
        if (MPI)
                list (APPEND MKL_LIBRARIES ${LIBMKL_BLACS})
        endif ()

        add_definitions (-D_MKL_)
        include_directories (${MKL_INCLUDE_PATH})
        set (LINALG_LIBRARIES ${MKL_LIBRARIES})
else ()
        mark_as_advanced(FORCE MKLROOT)
endif ()

# ACML settings, I don't think FindBLAS supports recent ACML either
set (ACMLROOT "" CACHE PATH "ACML path.")
if (LINALG STREQUAL "ACML")
        message ("-- Using AMD Core Math Library (ACML)")

        mark_as_advanced(CLEAR ACMLROOT)
        if (ACMLROOT STREQUAL "")
                set (ACMLROOT $ENV{ACMLROOT} CACHE PATH "ACML root directory." FORCE)
                if (NOT ACMLROOT)
                        message (FATAL_ERROR
                                "You must set environment variable ACMLROOT, "
                                "or specify -DACMLROOT=/path/to/acml_dir "
                                "when running cmake."
                                )
                endif ()
        endif ()

        unset (LIBACML CACHE)
        find_library (LIBACML
                NAMES acml
                PATHS ${ACMLROOT}
                PATH_SUFFIXES lib
                NO_DEFAULT_PATH)

        if (NOT LIBACML)
                message (FATAL_ERROR
                        "ACML library not found, please check that "
                        "the ACMLROOT variable is set and points to "
                        "a valid ACML installation directory."
                        )
        endif ()

        add_definitions (-D_ACML_)
        set (LINALG_LIBRARIES ${LIBACML})
else ()
        mark_as_advanced(FORCE ACMLROOT)
endif ()

# OpenBLAS
set (OPENBLASROOT "" CACHE PATH "OpenBLAS root directory.")
if (LINALG STREQUAL "OpenBLAS")
        message ("-- Using OpenBLAS LAPACK+BLAS libraries")

        # find out what OPENBLASROOT is
        mark_as_advanced(CLEAR OPENBLASROOT)
        if (OPENBLASROOT STREQUAL "")
                set (OPENBLASROOT $ENV{OPENBLASROOT} CACHE PATH "OpenBLAS root directory." FORCE)
                if (NOT OPENBLASROOT)
                        message (FATAL_ERROR
                                "You must set environment variable OPENBLASROOT, "
                                "or specify -DOPENBLASROOT=/path/to/openblas_dir "
                                "when running cmake."
                                )
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
                                "valid OpenBLAS installation directory."
                                )
                endif ()
                # check if the OpenBLAS installation was configured with 64bit integer support
                message ("-- Checking OpenBLAS for 64-bit integer interface...")
                include(CheckSymbolExists)
                unset (OPENBLAS_WITH_ILP64 CACHE)
                check_symbol_exists("OPENBLAS_USE64BITINT" ${OPENBLAS_CONFIG} OPENBLAS_WITH_ILP64)
                if (ADDRMODE EQUAL 64 AND NOT OPENBLAS_WITH_ILP64)
                        message (FATAL_ERROR
                                "OpenBLAS was not configured for 64-bit integer interface, "
                                "please build OpenBLAS with INTERFACE64=1 defined."
                                )
                endif ()
                # save the last location to check if it changed between configurations
                set (OPENBLASROOT_LAST ${OPENBLASROOT} CACHE INTERNAL "last value." FORCE)
        endif ()

        # search for the OpenBLAS library
        find_library (LIBOPENBLAS
                NAMES openblas
                PATHS ${OPENBLASROOT}
                PATH_SUFFIXES lib
                NO_DEFAULT_PATH)

        if (NOT LIBOPENBLAS)
                message (FATAL_ERROR
                        "OpenBLAS library not found, please check that "
                        "the OPENBLASROOT variable is set and points to "
                        "a valid OpenBLAS installation directory."
                        )
        endif ()

        # here we check if LIBOPENBLAS has changed after it was processed
        # succesfully, if not we do not need to rerun anything here.
        if (NOT LIBOPENBLAS STREQUAL LIBOPENBLAS_LAST)
                # check if the OpenBLAS library contains LAPACK functionality
                message ("-- Checking OpenBLAS for LAPACK functionality...")
                include(CheckFortranFunctionExists)
                set (CMAKE_REQUIRED_LIBRARIES ${LIBOPENBLAS})
                unset (OPENBLAS_WITH_LAPACK CACHE)
                check_fortran_function_exists("dsyev" OPENBLAS_WITH_LAPACK)
                if (NOT OPENBLAS_WITH_LAPACK)
                        unset (LIBOPENBLAS CACHE)
                        message (FATAL_ERROR
                                "LAPACK functionality missing from OpenBLAS library, "
                                "please build OpenBLAS with NO_LAPACK=0 defined."
                                )
                endif ()
                # save the last location to check if it changed between configurations
                set (LIBOPENBLAS_LAST ${LIBOPENBLAS} CACHE INTERNAL "last value." FORCE)
        endif ()

        add_definitions (-D_OPENBLAS_)
        set (LINALG_LIBRARIES ${LIBOPENBLAS})
else ()
        mark_as_advanced(FORCE OPENBLASROOT)
endif ()

if (LINALG STREQUAL "Accelerate")
        message ("-- Using XCode's veclib libraries")
        set (LINALG_LIBRARIES "-framework Accelerate")
endif ()

# Internal LAPACK+BLAS libraries shipped with Molcas source
if (LINALG STREQUAL "Internal")
        message ("-- Using Molcas' LAPACK+BLAS libraries (SLOW!)")
        set (LINALG_LIBRARIES "lapack" "blas")
endif ()

# Runtime loading of LAPACK+BLAS libraries (experimental)
if (LINALG STREQUAL "Runtime")
        message ("-- Using Runtime LAPACK+BLAS libraries (Internal as fallback)")
        set (LINALG_LIBRARIES "dl")

        add_definitions (-D_DELAYED_)
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
                        set_target_properties (cublas_fortran PROPERTIES COMPILE_DEFINITIONS "CUBLAS_GFORTRAN")
                endif ()
                target_include_directories (cublas_fortran PRIVATE ${CUDA_INCLUDE_PATH})
                list (APPEND objlibs "$<TARGET_OBJECTS:cublas_fortran>")

                # CUDA needs to link to the C++ standard library, so enable CXX here to
                # be able to use it as a linker later on for the program executables
                enable_language(CXX)

                add_definitions (-D_CUDA_BLAS_)

                find_library (LIBCUBLAS cublas
                        PATHS ${CUDAROOT}
                        PATH_SUFFIXES lib64
                        NO_DEFAULT_PATH)

                list (APPEND CUDA_LIBRARIES ${LIBCUBLAS})
        else ()
                if (NOT LINALG)
                        message (FATAL_ERROR "NVBLAS needs an external linear algebra library such as MKL or ACML")
                endif ()
                find_library (LIBNVBLAS nvblas
                        PATHS ${CUDAROOT}
                        PATH_SUFFIXES lib64
                        NO_DEFAULT_PATH)
                find_library (LIBCUBLAS cublas
                        PATHS ${CUDAROOT}
                        PATH_SUFFIXES lib64
                        NO_DEFAULT_PATH)

                list (APPEND CUDA_LIBRARIES ${LIBNVBLAS} ${LIBCUBLAS})
        endif ()

        message ("-- CUDA_LIBRARIES: ${CUDA_LIBRARIES}")

        list (APPEND EXTERNAL_LIBRARIES ${CUDA_LIBRARIES})
else ()
        mark_as_advanced(FORCE CUDAROOT)
endif ()

# Gromacs
#========

set (GROMACS_DIR "" CACHE PATH "Gromacs installation path.")
if (GROMACS)
        mark_as_advanced(CLEAR GROMACS_DIR)

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
        if (DEFINED ENV{GMXLDLIB})
                set (GMX_LIBDIR $ENV{GMXLDLIB})
        else ()
                set (GMX_LIBDIR "${GROMACS_DIR}/lib/*")
        endif ()

        find_library (LIBGMX
                NAMES gromacs_d gromacs
                PATHS ${GMX_LIBDIR}
                NO_DEFAULT_PATH)
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
                NO_DEFAULT_PATH)
        if (NOT INCGMX)
                message (FATAL_ERROR "INCGMX not found, make sure the GROMACS_DIR variable is set correctly.")
        else ()
                list (APPEND GMX_INCLUDE_PATH ${INCGMX})
        endif ()
        unset (INCGMX CACHE)

        foreach (util
                gateway_util
                espf_util
                )
                list (APPEND ${util}_incs ${GMX_INCLUDE_PATH})
                list (APPEND ${util}_defs "_GROMACS_")
        endforeach ()

        message ("-- GMX_LIBRARIES: ${GMX_LIBRARIES}")
        message ("-- GMX_INCLUDE_PATH: ${GMX_INCLUDE_PATH}")

        list (APPEND EXTERNAL_LIBRARIES ${GMX_LIBRARIES})
else ()
        mark_as_advanced(FORCE GROMACS_DIR)
        message ("Gromacs interface DISABLED")
endif ()

# Chan group's DMRG (Block code) support
#=======================================
if (BLOCK)
        mark_as_advanced(CLEAR BLOCK_DIR)

        message ("Configuring Block-DMRG support:")

        # Since CMake build is not available in current Block code,
        # suppose Block has been installed at ${BLOCK_DIR} for the time.
        set (BLOCK_DIR ${BLOCK_DIR} CACHE PATH "Block installation path." FORCE)
        if ("${BLOCK_DIR}" STREQUAL "")
                message (FATAL_ERROR "You must supply a path for BLOCK_DIR.")
        endif ()

        if (EXISTS ${BLOCK_DIR}/libqcdmrg.so)
                foreach (prog
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
        mark_as_advanced(FORCE BLOCK_DIR)
        message ("BLOCK interface DISABLED")
endif ()

# S.Wouters' DMRG (CheMPS2 code) support
#=======================================

if (CHEMPS2)
        mark_as_advanced(CLEAR CHEMPS2_DIR)

        message ("Configuring CheMPS2-DMRG support:")

        # suppose CheMPS2 has been installed at ${CHEMPS2_DIR} for the time.
        set (CHEMPS2_DIR ${CHEMPS2_DIR} CACHE PATH "CheMPS2 installation path." FORCE)
        if ("${CHEMPS2_DIR}" STREQUAL "")
                message (FATAL_ERROR "You must supply a path for CHEMPS2_DIR.")
        endif ()

        if (NOT HDF5)
                message (FATAL_ERROR "CheMPS2 interface requires enabling HDF5 with --enable-fortran2003.")
        endif ()

        find_library(HDF5_Fortran_LIBRARIES hdf5_fortran)
        if (HDF5_Fortran_LIBRARIES)
                message ("-- HDF5_Fortran_LIBRARIES: ${HDF5_Fortran_LIBRARIES}")
                list (APPEND EXTERNAL_LIBRARIES ${HDF5_Fortran_LIBRARIES})
        else ()
                message (FATAL_ERROR "CheMPS2 interface require HDF5 with fortran 2003 support")
        endif ()

        if (EXISTS ${CHEMPS2_DIR}/chemps2)
                foreach (prog
                        caspt2
                        rasscf
                        )
                        list (APPEND ${prog}_defs "_ENABLE_CHEMPS2_DMRG_")
                endforeach ()
                message ("-- CHEMPS2_DIR: ${CHEMPS2_DIR}")
        else ()
                message (FATAL_ERROR "No CheMPS2 installation found at ${CheMPS2_DIR}")
        endif ()

else ()
        mark_as_advanced(FORCE CHEMPS2_DIR)
        message ("CHEMPS2 interface DISABLED")
endif ()

if (BLOCK AND CHEMPS2)
        message (FATAL_ERROR "BLOCK and CHEMPS2 options are not compatible.")
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

# FDE interface
#==============

if (FDE)
        add_definitions (-D_FDE_)
endif ()

################################################################################
#                                                                              #
# External programs/libraries: compilation does not require a Molcas library,  #
# but needs to inherit the settings from our configuration.                    #
#                                                                              #
################################################################################

set (extprojsdir "${OPENMOLCAS_DIR}/External")

if (LINALG STREQUAL "Internal" OR LINALG STREQUAL "Runtime")
        set (proj ${extprojsdir}/lapack)
        if (EXISTS ${proj}/BLAS/SRC/ddot.f)
                list (APPEND internal_linalg ${proj}/SRC ${proj}/INSTALL ${proj}/BLAS/SRC)
        else ()
                message (FATAL_ERROR "LAPACK+BLAS sources not available, run \"${GIT} submodule update --init ${proj}\"")
        endif ()
endif ()

if (MSYM)
        message ("Configuring MSYM support:")

        set (proj ${extprojsdir}/libmsym)
        set (proj_bin ${PROJECT_BINARY_DIR}/External/libmsym)
        if (EXISTS ${proj}/CMakeLists.txt)
                set (MSYM_FOUND ON)
                add_subdirectory (${proj} ${proj_bin} EXCLUDE_FROM_ALL)
                # hide from the Molcas CMake GUI
                set (MSYM_BUILD_EXAMPLES OFF CACHE INTERNAL "Build example executables" FORCE)
                set (MSYM_BUILD_PYTHON   OFF CACHE INTERNAL "Build python binding"      FORCE)
                set (INSTALL_LIB_DIR     ${INSTALL_LIB_DIR}     CACHE INTERNAL "Installation directory for msym library." FORCE)
                set (INSTALL_INCLUDE_DIR ${INSTALL_INCLUDE_DIR} CACHE INTERNAL "Installation directory for header files." FORCE)
                set (INSTALL_CMAKE_DIR   ${INSTALL_CMAKE_DIR}   CACHE INTERNAL "Installation directory for CMake files."  FORCE)

                set (MSYM_INCLUDE ${proj}/src ${proj_bin})
                set (MSYM_LIBRARIES msym)

                message ("-- MSYM_INCLUDE_PATH: ${MSYM_INCLUDE}")
                message ("-- MSYM_LIBRARIES: ${MSYM_LIBRARIES}")

                list (APPEND EXTERNAL_LIBRARIES ${MSYM_LIBRARIES})
        else ()
                set (MSYM_FOUND OFF)
                message ("-- MSYM interface ignored (source not available, run \"${GIT} submodule update --init ${proj}\")")
        endif ()
else ()
        message ("MSYM support DISABLED")
endif ()

if (DMRG)
        message ("Configuring QCMaquis DMRG support:")

        if ( "${QCMaquis_NAME}" STREQUAL "" OR "${QCMaquis_EMAIL}" STREQUAL "")
          message (WARNING "Please set your name and e-mail with -DQCMaquis_NAME and -DQCMaquis_EMAIL to be able to download QCMaquis.")
        else ()
          message ("-- Your name (used to download QCMaquis): " ${QCMaquis_NAME})
          message ("-- Your email (used to download QCMaquis): " ${QCMaquis_EMAIL})
        endif ()

        add_definitions (-D_DMRG_)
        if (DMRG_DEBUG)
                add_definitions (-D_DMRG_DEBUG_)
        endif ()

        set (proj ${PROJECT_SOURCE_DIR}/src/qcmaquis)
        set (proj_bin ${PROJECT_BINARY_DIR}/External/qcmaquis-full)

        configure_file(${proj}/CMakeLists.txt.in ${proj_bin}/CMakeLists.txt COPYONLY)
        ## Builds with build subdirectory outside the source tree require the second argument of add_subdirectory
        add_subdirectory(${proj_bin} ${proj_bin})

        get_directory_property(${DMRG_INCLUDE} DIRECTORY ${proj_bin} DEFINITION ${DMRG_INCLUDE})
        get_directory_property(${DMRG_LIBRARIES} DIRECTORY ${proj_bin} DEFINITION ${DMRG_LIBRARIES})
        include_directories(${DMRG_INCLUDE})

        message ("-- DMRG_INCLUDE_PATH: ${DMRG_INCLUDE}")
        message ("-- DMRG_LIBRARIES: ${DMRG_LIBRARIES}")

        list (APPEND EXTERNAL_LIBRARIES ${DMRG_LIBRARIES})
else ()
        message ("DMRG support DISABLED")
endif ()

# HDF5 QCMaquis interface
# QCMaquis uses its own HDF5 interface, which is used here
# TODO: Leon: Actually, it's mostly NEVPT2 code that's using this interface, therefore it should be independent from DMRG/QCMaquis once we have a CASSCF-NEVPT2 interface
# and in that case we should add another variable here (or another variable for NEVPT2 and link the HDF5_QCM support to it)
if(HDF5 AND HDF5_FOUND AND DMRG)
        message ("Configuring HDF5 QCMaquis support:")

        set (proj ${PROJECT_SOURCE_DIR}/src/hdf5_qcm_util)
        set (proj_bin ${PROJECT_BINARY_DIR}/External/hdf5_qcm_util)

        configure_file(${proj}/CMakeLists.txt.in ${proj_bin}/CMakeLists.txt COPYONLY)
        add_subdirectory(${proj_bin} ${proj_bin})

        get_directory_property(${HDF5_QCM_INCLUDE} DIRECTORY ${proj_bin} DEFINITION ${HDF5_QCM_INCLUDE})
        get_directory_property(${HDF5_QCM_LIBRARIES} DIRECTORY ${proj_bin} DEFINITION ${HDF5_QCM_LIBRARIES})
        include_directories(${HDF5_QCM_INCLUDE})

        message ("-- HDF5_QCM_INCLUDE_PATH: ${HDF5_QCM_INCLUDE}")
        message ("-- HDF5_QCM_LIBRARIES: ${HDF5_QCM_LIBRARIES}")

        insert_before(EXTERNAL_LIBRARIES HDF5_C_LIBRARIES "${HDF5_QCM_LIBRARIES}")
else ()
        message ("HDF5 QCMaquis interface DISABLED")
endif ()

# NECI Full-CI Quantum Monte-Carlo
if (NECI)
        if (NOT MPI)
            message (FATAL_ERROR "NECI requires MPI (-DMPI=ON).")
        endif()
        message ("Configuring NECI support:")
        set (proj ${extprojsdir}/NECI)
        set (proj_bin ${PROJECT_BINARY_DIR}/External/NECI)
        if (EXISTS ${proj}/CMakeLists.txt)
                # set default NECI options when built from Molcas
                # and hide them from the Molcas CMake GUI
                set (ENABLE_MOLCAS ON CACHE INTERNAL "Enable Molcas interface.")
                set (ENABLE_HDF5 ON CACHE INTERNAL "Enable HDF5 popsfiles.")
                set (WARNINGS OFF CACHE INTERNAL "Enable lots of compile warnings.")
                set (ENABLE_SHARED_MEMORY ON CACHE INTERNAL "Enable shared memory.")
                set (ENABLE_SHARED ON CACHE INTERNAL "Enable shared libraries.")
                # set (ENABLE_BUILD_HDF5 OFF CACHE INTERNAL "Provide precompiled HDF5.")

                list (APPEND rasscf_defs "_NECI_")
                add_subdirectory (${proj} ${proj_bin} EXCLUDE_FROM_ALL)
                set (NECI_LIBRARIES libdneci)

                message ("-- NECI_INCLUDE_PATH: ${NECI_INCLUDE}")
                message ("-- NECI_LIBRARIES: ${NECI_LIBRARIES}")

                list (APPEND EXTERNAL_LIBRARIES ${NECI_LIBRARIES})
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
                        set(LIBEFP_FORTRAN_INT64 1)
                else ()
                        set(LIBEFP_FORTRAN_INT64 0)
                endif ()

                include(ExternalProject)
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
                link_directories(${proj_bin})
                list (APPEND EXTERNAL_PROJECTS ${projname})

                set (EFP_LIBRARIES efp)

                message ("-- EFP_LIBRARIES: ${EFP_LIBRARIES}")

                list (APPEND EXTERNAL_LIBRARIES ${EFP_LIBRARIES})
                add_definitions (-D_EFP_)
        else()
                message ("-- EFP interface ignored (source not available, run git submodule update --init ${proj}\")")
        endif ()
else ()
        message ("EFP interface DISABLED")
endif ()

# GEN1INT interface
if (GEN1INT)
        message ("Configuring GEN1INT")

        set (proj ${PROJECT_SOURCE_DIR}/src/gen1int_util)
        set (proj_bin ${PROJECT_BINARY_DIR}/External/gen1int_util)

        configure_file(${proj}/CMakeLists.txt.in ${proj_bin}/CMakeLists.txt COPYONLY)
        add_subdirectory(${proj_bin} ${proj_bin})

        get_directory_property(${GEN1INT_INCLUDE} DIRECTORY ${proj_bin} DEFINITION ${GEN1INT_INCLUDE})
        get_directory_property(${GEN1INT_LIBRARIES} DIRECTORY ${proj_bin} DEFINITION ${GEN1INT_LIBRARIES})
        include_directories(${GEN1INT_INCLUDE})

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

        message("Configuring libwfa support:")

        set (proj ${extprojsdir}/libwfa)
        set (proj_bin ${PROJECT_BINARY_DIR}/External/libwfa)
        if (EXISTS ${proj}/CMakeLists.txt)
                set (WFA_FOUND ON)
                mark_as_advanced(CLEAR ARMADILLO_INC)
                set (ARMADILLO_INC ${ARMADILLO_INC} CACHE PATH "Directory with armadillo header files")
                set (MOLCAS_LIB ON CACHE INTERNAL "Compile Molcas support")

                add_subdirectory (${proj}/libwfa)
                target_include_directories (wfa_molcas PRIVATE ${proj} ${ARMADILLO_INC} ${HDF5_INCLUDE_DIRS})
                target_compile_definitions (wfa_molcas PRIVATE "ARMA_DONT_USE_WRAPPER")
                list (APPEND wfa_libs wfa_molcas)

                message ("-- wfa_libs: ${wfa_libs}")
                message ("-- ARMADILLO_INC: ${ARMADILLO_INC}")

                insert_before(EXTERNAL_LIBRARIES HDF5_C_LIBRARIES "${wfa_libs}")
        else ()
                set (WFA_FOUND OFF)
                mark_as_advanced(FORCE ARMADILLO_INC)
                message ("-- libwfa support ignored (source not available, run \"${GIT} submodule update --init ${proj}\")")
        endif ()
else ()
        mark_as_advanced(FORCE ARMADILLO_INC)
        message("libwfa support DISABLED")
endif ()

# NEVPT2
if (NEVPT2)
        if (NOT DMRG)
                message (FATAL_ERROR "NEVPT2 support requires enabling DMRG too.")
        endif ()
        if (NOT (HDF5 AND HDF5_FOUND))
                message (FATAL_ERROR "NEVPT2 support requires enabling HDF5 too.")
        endif ()

        set (proj ${PROJECT_SOURCE_DIR}/src/nevpt2)
        set (proj_bin ${PROJECT_BINARY_DIR}/External/nevpt2_ext)

        configure_file(${proj}/CMakeLists.txt.in ${proj_bin}/CMakeLists.txt COPYONLY)

        add_subdirectory(${proj_bin} ${proj_bin})

        get_directory_property(${NEVPT2_INCLUDE} DIRECTORY ${proj_bin} DEFINITION ${NEVPT2_INCLUDE})
        get_directory_property(${NEVPT2_LIBRARIES} DIRECTORY ${proj_bin} DEFINITION ${NEVPT2_LIBRARIES})

        add_dependencies(nevpt2_ext hdf5-qcm qcmaquis-driver)

        message ("-- NEVPT2_LIBRARIES: ${NEVPT2_LIBRARIES}")

        insert_before(EXTERNAL_LIBRARIES HDF5_QCM_LIBRARIES "${NEVPT2_LIBRARIES}")
else ()
        message("NEVPT2 support DISABLED")
endif ()

# external GRID_IT
if (GRID_IT)
        message("Configuring external GRID_IT:")
        set (proj ${extprojsdir}/grid_it)
        if (EXISTS ${proj}/main.f)
                set (GRID_IT_FOUND ON)
                set (prog "grid_it")
                list (APPEND progs ${prog})
                file (RELATIVE_PATH ${prog}_src ${PROJECT_SOURCE_DIR} ${proj})
                # some routines in the external code are referenced by libmolcas,
                # so they must be compiled with the molcas library and not with the "program" library
                list (APPEND objlibs "$<TARGET_OBJECTS:${prog}>")
                list (APPEND proglib ${prog})
                message ("-- Source found, external GRID_IT enabled")
        else()
                set (GRID_IT_FOUND OFF)
                message ("-- External GRID_IT ignored (source not available, run \"${GIT} submodule update --init ${proj}\")")
        endif ()
else ()
        message("External GRID_IT DISABLED")
endif ()

# Add the "bigot" flags after processing the external projects, so they are not affected
if (BIGOT)
        # Modify the "BUILD_TARGET" flags because we want them at the end
        string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TARGET)
        set (CMAKE_CXX_FLAGS_${BUILD_TARGET} "${CMAKE_CXX_FLAGS_${BUILD_TARGET}} ${CXXFLAGS_${CXXCID}_BIGOT}")
        set (CMAKE_C_FLAGS_${BUILD_TARGET} "${CMAKE_C_FLAGS_${BUILD_TARGET}} ${CFLAGS_${CCID}_BIGOT}")
        set (CMAKE_Fortran_FLAGS_${BUILD_TARGET} "${CMAKE_Fortran_FLAGS_${BUILD_TARGET}} ${FFLAGS_${FCID}_BIGOT}")
        # The Intel compiler can generate automatic interfaces, in that case we want to keep them all in a single place
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
                set(SINGLE_MOD_DIR "1")
        endif ()
endif ()

# global module directory
#========================

set (MAIN_MOD_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)
if (SINGLE_MOD_DIR)
        include_directories (${MAIN_MOD_DIR}/_single)
else ()
        include_directories (${MAIN_MOD_DIR})
endif ()

################################################################################
#                                                                              #
# Configure settings for the Molcas run-time environment                       #
#                                                                              #
################################################################################

message ("Configuring runtime environment settings:")

# memory/disk default sizes
#==========================

if (ADDRMODE EQUAL 64)
        set (DEFMOLCASMEM "2048" CACHE STRING "Set default max memory (MB).")
        set (DEFMOLCASDISK "20000" CACHE STRING "Set default max file size (MB).")
else ()
        set (DEFMOLCASMEM "1024" CACHE STRING "Set default max memory (MB).")
        set (DEFMOLCASDISK "2047"  CACHE STRING "Set default max file size (MB).")
endif ()

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
get_directory_property(defs COMPILE_DEFINITIONS)
message ("-- Definitions: ${defs}")
message ("-- Debug definitions: ${DEBUG_DEFS}")

################################################################################
#                                                                              #
# Utilities/Programs                                                           #
#                                                                              #
################################################################################

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
        blas_m_util
        block_dmrg_util
        chemps2_util
        delayed_util
        dga_util
        embedding_util
        hdf5_util
        lapack_m_util
        libmsym_util # remove when deleted
        msym_util
        )
list (REMOVE_ITEM progs
        parnell
        symmetrize
        wfa
        nevpt2 # although not necessarily
        dmrgscf
        mpssi
        )

# set the location of the sources for the utilities
foreach (util ${utils})
        find_source(${util})
        if (NOT DEFINED ${util}_src)
                list (REMOVE_ITEM utils ${util})
                message (WARNING "\"${util}\" not found in source directories, skipping")
        endif ()
endforeach ()

if (FDE)
        set (util embedding_util)
        find_source(${util})
        if (DEFINED ${util}_src)
                list (APPEND utils ${util})
                include_directories (${${util}_src})
        else ()
                message (FATAL_ERROR "\"${util}\" not found in source directories")
        endif ()
        foreach (util
                gateway_util
                )
                list (APPEND ${util}_incs ${${util}_src})
        endforeach ()
        foreach (prog
                seward
                scf
                rasscf
                grid_it
                )
                list (APPEND ${prog}_incs ${${util}_src})
        endforeach ()
endif ()

# If an external GA library is not present, we use the internal
# DGA code which resides in dga_util:
if (MPI AND NOT GA)
        set (util dga_util)
        find_source(${util})
        if (DEFINED ${util}_src)
                list (APPEND utils ${util})
                include_directories (${${util}_src})
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
        find_source(${util})
        if (DEFINED ${util}_src)
                list (APPEND utils ${util})
        else ()
                message (FATAL_ERROR "\"${util}\" not found in source directories")
        endif ()
        list (APPEND ${util}_incs ${internal_linalg})
endif ()

if (GEN1INT)
        list (APPEND integral_util_defs "_GEN1INT_")
endif()

# activate HDF5 wrappers for I/O
if (HDF5 AND HDF5_FOUND)
        set (util hdf5_util)
        find_source(${util})
        if (DEFINED ${util}_src)
                list (APPEND utils ${util})
                list (APPEND ${util}_incs ${HDF5_INCLUDE_DIRS})
        else ()
                message (FATAL_ERROR "\"${util}\" not found in source directories")
        endif ()
        # additional utilities doing internal HDF5 calls
        foreach (util
                Modules
                aniso_util
                guessorb_util
                hyper_util
                molcas_ci_util
                msym_util
                oneint_util
                property_util
                restart_util
                runfile_util
                )
                list (APPEND ${util}_defs "_HDF5_")
        endforeach ()
        if (DMRG)
                foreach (util
                        transform_util
                        )
                        list (APPEND ${util}_defs "_HDF5_QCM_")
                endforeach ()
        endif ()
endif ()

if (NEVPT2)
       set (prog nevpt2)
       list (APPEND progs ${prog})
       list (APPEND ${prog}_incs ${NEVPT2_INCLUDE})
#       list (APPEND nevpt2_defs "_NEVPT2_") # if we have conditional compilation, we theoretically don't need this flag
#       set (nevpt2_deps $<TARGET_OBJECTS:caspt2>) # TODO: and somehow hdf5_qcm and qcmaquis driver
endif ()

# activate MSYM support
if (MSYM AND MSYM_FOUND)
        set (util msym_util)
        find_source(${util})
        if (DEFINED ${util}_src)
                list (APPEND utils ${util})
        else ()
                message (FATAL_ERROR "\"${util}\" not found in source directories")
        endif ()
        list (APPEND ${util}_incs ${MSYM_INCLUDE})
        list (APPEND progs symmetrize)
endif ()

if (BLOCK)
        set (util block_dmrg_util)
        find_source(${util})
        if (DEFINED ${util}_src)
                list (APPEND utils ${util})
                list (APPEND ${util}_defs "_ENABLE_BLOCK_DMRG_")
        else ()
                message (FATAL_ERROR "\"${util}\" not found in source directories")
        endif ()
endif ()

if (CHEMPS2)
        set (util chemps2_util)
        find_source(${util})
        if (DEFINED ${util}_src)
                list (APPEND utils ${util})
                list (APPEND ${util}_incs ${HDF5_INCLUDE_DIRS})
                list (APPEND ${util}_defs "_ENABLE_CHEMPS2_DMRG_")
        else ()
                message (FATAL_ERROR "\"${util}\" not found in source directories")
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
        find_source(${prog})
        if (NOT DEFINED ${prog}_src)
                list (REMOVE_ITEM progs ${prog})
                message (WARNING "\"${prog}\" not found in source directories, skipping")
        endif ()
endforeach ()

if (MSYM AND MSYM_FOUND)
        list (APPEND scf_defs "_MSYM_")
endif ()

if (WFA AND WFA_FOUND)
        list (APPEND wfa_defs "_WFA_")
endif ()

if (GEN1INT)
        list (APPEND seward_defs "_GEN1INT_")
endif()

# programs doing HDF5 calls
if (HDF5 AND HDF5_FOUND)
        foreach (prog
                caspt2
                dynamix
                mcpdft
                nevpt2
                poly_aniso
                rasscf
                rassi
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
                endforeach ()
        endif ()
endif ()

# activate grid_it
if (GRID_IT AND GRID_IT_FOUND)
        # Some additional files and configuration needed if using external grid_it
        configure_file (${grid_util_src}/grid.fh.in ${MAIN_INCLUDE_DIR}/grid.fh)
        list (APPEND external_prgm "${grid_util_src}/grid_it.prgm.src")
        list (APPEND grid_util_defs "_EXTERNAL_GRID_IT_")
else ()
        file (REMOVE ${MAIN_INCLUDE_DIR}/grid.fh)
endif ()
# grid_it could be external or from molcas_extra
if (DEFINED grid_it_src)
        list (APPEND grid_util_defs "_HAVE_GRID_IT_")
endif ()

# finally, sort the list
list (SORT progs)

# set program dependencies of super modules
set (casvb_deps
        $<TARGET_OBJECTS:rasscf>
        )
set (dmrgscf_deps
        $<TARGET_OBJECTS:rasscf>
        )
set (mpssi_deps
        $<TARGET_OBJECTS:rassi>
        )
set (espf_deps
        $<TARGET_OBJECTS:seward>
        $<TARGET_OBJECTS:alaska>
        )
set (loprop_deps
        $<TARGET_OBJECTS:scf>
        $<TARGET_OBJECTS:rasscf>
        $<TARGET_OBJECTS:caspt2>
        $<TARGET_OBJECTS:mbpt2>
        )
set (numerical_gradient_deps
        $<TARGET_OBJECTS:seward>
        $<TARGET_OBJECTS:scf>
        $<TARGET_OBJECTS:rasscf>
        $<TARGET_OBJECTS:caspt2>
        $<TARGET_OBJECTS:mbpt2>
        $<TARGET_OBJECTS:ccsdt>
        $<TARGET_OBJECTS:motra>
        $<TARGET_OBJECTS:chcc>
        $<TARGET_OBJECTS:cht3>
        $<TARGET_OBJECTS:mcpdft>
        )
set (last_energy_deps
        $<TARGET_OBJECTS:seward>
        $<TARGET_OBJECTS:scf>
        $<TARGET_OBJECTS:rasscf>
        $<TARGET_OBJECTS:caspt2>
        $<TARGET_OBJECTS:mbpt2>
        $<TARGET_OBJECTS:ccsdt>
        $<TARGET_OBJECTS:motra>
        $<TARGET_OBJECTS:chcc>
        $<TARGET_OBJECTS:cht3>
        $<TARGET_OBJECTS:mcpdft>
        )
if (BLOCK OR DMRG)
list (APPEND last_energy_deps
        $<TARGET_OBJECTS:dmrgscf>
        )
endif ()

################################################################################
#                                                                              #
# Disable some warnings that are buggy or on code we don't want to modify      #
#                                                                              #
#################################################################################
# netlib blas/lapack produces warnings
# some external codes (grid_it) too
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set (netlib_nowarn_flags " -w")
        set (ext_nowarn_flags " -w")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        set (netlib_nowarn_flags " -w")
        set (ext_nowarn_flags " -w")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro")
        set (netlib_nowarn_flags " -erroff=%all")
        set (ext_nowarn_flags " -erroff=%all")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
        set (netlib_nowarn_flags " -Kieee -w")
        set (ext_nowarn_flags " -w")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
        set (netlib_nowarn_flags " -w")
        set (ext_nowarn_flags " -w")
endif ()
find_source(LinAlg_internal)
if (netlib_nowarn_flags)
        file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${LinAlg_internal_src}/*.f ${delayed_util_src}/*.f)
        set_property (SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS ${netlib_nowarn_flags})
endif ()
if (ext_nowarn_flags)
        file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${grid_it_src}/*.f)
        set_property (SOURCE ${source} APPEND_STRING PROPERTY COMPILE_FLAGS ${ext_nowarn_flags})
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
# Build the utilities and the Molcas library                                   #
#                                                                              #
################################################################################

# first, build global modules separately, so other utils can add it as a
# dependency using the "mods_obj" target.
set (util "Modules")
file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${OPENMOLCAS_DIR}/src/${util}/*.f ${OPENMOLCAS_DIR}/src/${util}/*.f90)
add_library (mods_obj OBJECT ${source})
list (APPEND ALL_SOURCES ${source})
if (DEFINED ${util}_defs)
        set_target_properties (mods_obj PROPERTIES COMPILE_DEFINITIONS "${${util}_defs}")
endif ()
if (DMRG)
        if (HDF5 AND HDF5_FOUND)
                add_dependencies(mods_obj hdf5-qcm)
        else ()
                message (FATAL_ERROR "QCMaquis DMRG support requires the HDF5 library.")
        endif ()

        add_dependencies(mods_obj qcmaquis-driver)
endif ()
if (GEN1INT)
        add_dependencies(mods_obj gen1int)
endif ()
if (NEVPT2)
        add_dependencies(mods_obj nevpt2_ext)
endif ()
if (BUILD_SHARED_LIBS)
        set_target_properties (mods_obj PROPERTIES
                POSITION_INDEPENDENT_CODE "TRUE"
        )
endif ()
if (SINGLE_MOD_DIR)
        set_target_properties (mods_obj PROPERTIES
                Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR}/_single)
else ()
        set_target_properties (mods_obj PROPERTIES
                Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR})
endif ()
target_include_directories (mods_obj PRIVATE
        ${OPENMOLCAS_DIR}/src/${util}
        )
list (APPEND utils_obj mods_obj)
list (APPEND objlibs "$<TARGET_OBJECTS:mods_obj>")

# now go through utility list
foreach (util ${utils})
        file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${${util}_src}/*.f ${${util}_src}/*.f90 ${${util}_src}/*.c)
        list (APPEND ALL_SOURCES ${source})
        # ignore empty directories
        if (NOT "${source}" STREQUAL "")
                add_library (${util}_obj OBJECT ${source})
                if (SINGLE_MOD_DIR)
                        set (mod_dir ${MAIN_MOD_DIR}/_single)
                else ()
                        set (mod_dir ${MAIN_MOD_DIR}/${util})
                endif ()
                set_target_properties (${util}_obj PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})

                if (BUILD_SHARED_LIBS)
                        set_target_properties (${util}_obj PROPERTIES
                                POSITION_INDEPENDENT_CODE "TRUE")
                endif ()

                # utility-specific compile definitions
                if (DEFINED ${util}_defs)
                        set_target_properties (${util}_obj PROPERTIES COMPILE_DEFINITIONS "${${util}_defs}")
                endif ()
                # utility-specific include directories
                list (APPEND ${util}_incs ${PROJECT_SOURCE_DIR}/${${util}_src} ${mod_dir})
                foreach (inc ${${util}_incs})
                        target_include_directories (${util}_obj PRIVATE ${inc})
                endforeach ()

                add_dependencies (${util}_obj mods_obj)

                list (APPEND utils_obj ${util}_obj)
                list (APPEND objlibs "$<TARGET_OBJECTS:${util}_obj>")
                #add_library (${util} $<TARGET_OBJECTS:${util}_obj>)

                # if (EXISTS ${PROJECT_SOURCE_DIR}/${${util}_src}/custom.cmake)
                #         include(${PROJECT_SOURCE_DIR}/${${util}_src}/custom.cmake)
                # endif ()
        endif ()
endforeach (util)

# build an actual Molcas library
add_library (molcas ${objlibs})
set (MOLCAS_LIBRARIES "molcas")
if (BUILD_SHARED_LIBS)
        # flatten namespace on Mac OS X
        if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
                set_target_properties (molcas PROPERTIES
                        LINK_FLAGS "-Wl,-flat_namespace")
        endif ()
        # also create a static library!
        if (BUILD_STATIC_LIB)
                add_library (molcas_static STATIC ${objlibs})
                list (APPEND MOLCAS_LIBRARIES "molcas_static")
        endif ()
endif ()

# if necessary, build BLAS/LAPACK
if (LINALG STREQUAL "Internal")
        set (util "LinAlg_internal")
        if (NOT DEFINED ${util}_src)
                message (FATAL_ERROR "\"${util}\" not found in source directories")
        endif ()
        foreach (lib "lapack" "blas")
                set (source ${${util}_src}/${lib}.f)
                add_library (${lib} ${source})
                list (APPEND ALL_SOURCES ${source})
                list (APPEND ${util}_incs ${internal_linalg})
                foreach (inc ${${util}_incs})
                        target_include_directories (${lib} PRIVATE ${inc})
                endforeach ()
                if (BUILD_SHARED_LIBS)
                        set_target_properties (${lib} PROPERTIES
                                POSITION_INDEPENDENT_CODE "TRUE")
                endif ()
        endforeach ()
endif ()

################################################################################
#                                                                              #
# Build the program modules                                                    #
#                                                                              #
################################################################################

add_custom_target(only_objs
        DEPENDS ${utils_obj}
        )

foreach (prog ${progs})
        file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${${prog}_src}/*.f ${${prog}_src}/*.f90 ${${prog}_src}/*.F90 ${${prog}_src}/*.c)
        list (APPEND ALL_SOURCES ${source})
        # ignore empty programs
        if ("${source}" STREQUAL "")
                list (REMOVE_ITEM progs ${prog})
        else ()
                # FIXME: the checks below remove any main.* file, which is needed in
                # order to be able to make a library out of the program sources. This
                # is required for supermodules.
                #------------------------------------------------------------------
                foreach (lext "f" "f90" "c")
                        if (EXISTS ${PROJECT_SOURCE_DIR}/${${prog}_src}/main.${lext})
                                list (REMOVE_ITEM source ${${prog}_src}/main.${lext})
                                set (${prog}_main ${${prog}_src}/main.${lext})
                        endif ()
                endforeach ()

                if (NOT "${source}" STREQUAL "")
                        add_library (${prog} OBJECT ${source})
                        if (SINGLE_MOD_DIR)
                                set (mod_dir ${MAIN_MOD_DIR}/_single)
                        else ()
                                set (mod_dir ${MAIN_MOD_DIR}/${prog})
                        endif ()
                        set_target_properties (${prog} PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})

                        # programs in "proglib" have their libraries compiled as utils
                        if (BUILD_SHARED_LIBS AND ";${proglib};" MATCHES ";${prog};")
                                set_target_properties (${prog} PROPERTIES
                                        POSITION_INDEPENDENT_CODE "TRUE")
                        endif ()

                        # program-specific compile definitions
                        if (DEFINED ${prog}_defs)
                                set_target_properties (${prog} PROPERTIES COMPILE_DEFINITIONS "${${prog}_defs}")
                        endif ()
                        # program-specific include directories
                        list (APPEND ${prog}_incs ${PROJECT_SOURCE_DIR}/${${prog}_src} ${mod_dir})
                        foreach (inc ${${prog}_incs})
                                target_include_directories (${prog} PRIVATE ${inc})
                        endforeach ()

                        add_dependencies (${prog} ${utils_obj})
                        add_dependencies (only_objs ${prog})
                endif ()
        endif ()
endforeach ()

# set file-specific flags
if (BOUNDS AND (${FFLAGS_GNU_BOUNDS} MATCHES "-fsanitize=address"))
        set_property(SOURCE ${rassi_src}/prprop.f APPEND_STRING PROPERTY COMPILE_FLAGS " --param=max-vartrack-size=120000000")
endif ()

# set up proper RPATH for executables, use RUNPATH if on Linux
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

foreach (prog ${progs})
        # FIXME: Actual program target only depends on main file now, but if
        # supermodules are gone, the dependency should be on ${source} instead.
        #----------------------------------------------------------------------
        # programs in "proglib" have their library included in libmolcas
        if (TARGET ${prog} AND NOT ";${proglib};" MATCHES ";${prog};")
                list (APPEND ${prog}_deps $<TARGET_OBJECTS:${prog}>)
        endif ()
        add_executable (${prog}.exe ${${prog}_main} ${${prog}_deps})
        list (APPEND progs_exe ${prog}.exe)
        if (SINGLE_MOD_DIR)
                set_target_properties (${prog}.exe PROPERTIES
                        Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR}/_single)
        else ()
                set_target_properties (${prog}.exe PROPERTIES
                        Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR}/${prog})
        endif ()
        foreach (inc ${${prog}_incs})
                target_include_directories (${prog}.exe PRIVATE ${inc})
        endforeach ()

        target_link_libraries (${prog}.exe molcas)
        target_link_libraries (${prog}.exe ${EXTERNAL_LIBRARIES})
        if (EXTERNAL_PROJECTS)
                add_dependencies (${prog}.exe ${EXTERNAL_PROJECTS})
        endif ()

        if (CUBLAS)
                set_target_properties (${prog}.exe PROPERTIES LINKER_LANGUAGE "CXX")
        endif ()

        list (APPEND PROGRAM_EXECUTABLES ${prog}.exe)
endforeach (prog)

# build parnell as a stand-alone executable
set (prog parnell)
find_source(${prog})
if (NOT DEFINED ${prog}_src)
        message (FATAL_ERROR "\"${prog}\" not found in source directories")
endif ()
file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${${prog}_src}/*.f ${${prog}_src}/*.f90 ${${prog}_src}/*.c)
add_executable (${prog}.exe ${source})
list (APPEND progs_exe ${prog}.exe)
list (APPEND ALL_SOURCES ${source})
list (APPEND PROGRAM_EXECUTABLES ${prog}.exe)
if (MPI)
        target_link_libraries (${prog}.exe ${MPI_C_LIBRARIES})
endif ()

if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        if (${CMAKE_SYSTEM_VERSION} STREQUAL "15.0.0")
                if (LINALG STREQUAL "MKL")
                        foreach (prog ${PROGRAM_EXECUTABLES})
                                # workaround because "continue()" is not supported in CMake 2.8.11
                                if (NOT ${prog} STREQUAL "parnell.exe")
                                        add_custom_target(dylibpath_${prog} ALL
                                                          DEPENDS ${prog}
                                                          VERBATIM
                                                         )
                                        add_dependencies(dylibpath_${prog} ${prog})
                                        add_custom_command(TARGET dylibpath_${prog}
                                                           POST_BUILD
                                                           DEPENDS dylibpath_${prog}
                                                           COMMAND ${OPENMOLCAS_DIR}/sbin/fixdylib  ${EXECUTABLE_OUTPUT_PATH}/${prog} libmkl ${MKLROOT}
                                                           COMMENT "Enforce absolute path to MKL in dynamically linked executables"
                                                           VERBATIM
                                                          )
                                endif ()
                        endforeach (prog)
                endif ()
        endif ()
endif ()

# Add custom debug flags per file:
#   For each file, tries to find out if each definition is used and modify the file's properties
# Limitation: it will rebuild the whole target (but not the whole molcas)
if (DEBUG_DEFS)
        foreach (file ${ALL_SOURCES})
                file (READ ${file} contents)
                foreach (def ${DEBUG_DEFS})
                        string (REGEX MATCH ${def} result "${contents}")
                        if (NOT ${result} STREQUAL "")
                                set_property (SOURCE ${file} APPEND PROPERTY COMPILE_DEFINITIONS ${def})
                        endif ()
                endforeach ()
        endforeach ()
endif ()

################################################################################
#                                                                              #
# Tools, i.e. builds are dependent on Molcas library                           #
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
                                        add_subdirectory (${toolsdir}/${tool} ${PROJECT_BINARY_DIR}/Tools/${tool})
                                        message ("-- ${tool}: added to targets")
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
configure_file (${OPENMOLCAS_DIR}/data/rysrw.ascii  data/rysrw        COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/abdata.ascii data/abdata       COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/banner.src   data/banner.txt   COPYONLY)
configure_file (${OPENMOLCAS_DIR}/data/modalias.src data/modalias.txt COPYONLY)
file (COPY ${OPENMOLCAS_DIR}/data/inputs.tpl DESTINATION ${PROJECT_BINARY_DIR}/data)
if (DEFINED EXTRA_DIR)
        foreach (f alias.plx auto.prgm landing.txt)
                file (COPY ${EXTRA_DIR}/data/${f} DESTINATION ${PROJECT_BINARY_DIR}/data)
        endforeach ()
endif ()
# copy the necessary files from sbin
file (REMOVE_RECURSE ${PROJECT_BINARY_DIR}/sbin)
foreach (f verify version chkunprint.plx updatetest.plx find_sources have_feature help_basis help_doc)
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
                message (FATAL_ERROR "Failed to install the Molcas driver")
        endif ()
else ()
        set (MOLCAS_DRIVER "${OPENMOLCAS_DIR}/sbin/dummy.driver")
        # if molcas-extra is not available, enforce building pymolcas
        if (NOT TOOLS)
                set (BASE_DIR "${OPENMOLCAS_DIR}")
                set (toolsdir "${OPENMOLCAS_DIR}/Tools")
                set (tool "pymolcas")
                add_subdirectory (${toolsdir}/${tool} ${PROJECT_BINARY_DIR}/Tools/${tool})
                message ("-- ${tool}: added to targets")
        endif ()
endif ()

get_directory_property(defs COMPILE_DEFINITIONS)
string(STRIP ${CMAKE_C_FLAGS} CMAKE_C_FLAGS)
string(STRIP ${CMAKE_Fortran_FLAGS} CMAKE_Fortran_FLAGS)
if (NOT BUILD_TARGET STREQUAL "NONE")
        file (WRITE "${PROJECT_BINARY_DIR}/data/info.txt"
                "\nconfiguration info\n"
                "------------------\n"
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
                        COMMAND "${OPENMOLCAS_DIR}/sbin/install_hooks.sh" "${BASE_DIR}"
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
string(REGEX REPLACE "([^;]+)" "-D\\1" defargs "${defs}")
message("Definitions: ${defargs}")
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
                get_filename_component(prgm ${prgm} NAME)
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
add_custom_target(prgms ALL
        DEPENDS ${prgmtargets}
        )
add_dependencies(mods_obj prgms)

# generate help databases from doc
if (PYTHONINTERP_FOUND)
        execute_process (
                COMMAND ${PYTHON_EXECUTABLE} ${OPENMOLCAS_DIR}/sbin/extract_xmldoc.py ${OPENMOLCAS_DIR}/doc ${PROJECT_BINARY_DIR}/data
                ERROR_VARIABLE EXTRACT_ERROR
                RESULT_VARIABLE EXTRACT_RESULT
                )
        if (EXTRACT_ERROR)
                message(FATAL_ERROR "Running \"extract_xmldoc.py\" failed: ${EXTRACT_ERROR}")
        endif ()
        # check XML consistency, only if the program is found
        find_program (XMLLINT "xmllint")
        mark_as_advanced(FORCE XMLLINT)
        if (XMLLINT)
                execute_process (
                        COMMAND ${XMLLINT} "${PROJECT_BINARY_DIR}/data/keyword.xml" "-noout"
                        ERROR_VARIABLE XMLLINT_ERROR
                        RESULT_VARIABLE XMLLINT_RESULT
                        )
                if (XMLLINT_ERROR)
                        message(FATAL_ERROR "Extracted XML file is malformed: ${XMLLINT_ERROR}")
                endif ()
        endif ()
else ()
        message (WARNING "No python found, help databases will not be created.")
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

if (CTEST)
        enable_testing()

        set (MODULENAMES ${PROJECT_SOURCE_DIR}/sbin/modulenames.plx)
        file (STRINGS ${PROJECT_BINARY_DIR}/test/testdirs testdirs)
        set (known_groups "")
        foreach (dir ${testdirs})
                file (GLOB groups RELATIVE ${dir} ${dir}/*)
                foreach (group ${groups})
                        list (FIND known_groups ${group} _index)
                        if (_index GREATER -1)
                                break ()
                        else ()
                                list (APPEND known_groups ${group})
                        endif ()
                        file (GLOB inputs RELATIVE ${dir}/${group} ${dir}/${group}/*)
                        list (SORT inputs)
                        foreach (input ${inputs})
                                string (REGEX MATCH "^([^.]+).input$" match ${input})
                                if (NOT ${match} STREQUAL "")
                                        set (number ${CMAKE_MATCH_1})
                                        execute_process (
                                                COMMAND ${MODULENAMES} ${dir}/${group}/${input}
                                                WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                                                OUTPUT_VARIABLE module_names
                                                )
                                        add_test (
                                                NAME "${group}-${CMAKE_MATCH_1}__${module_names}"
                                                COMMAND ${MOLCAS_DRIVER} verify -q "${group}:${number}"
                                                )
                                        set_property (
                                                TEST "${group}-${CMAKE_MATCH_1}__${module_names}"
                                                PROPERTY ENVIRONMENT MOLCAS=${PROJECT_BINARY_DIR}
                                                )
                                endif ()
                        endforeach ()
                endforeach ()
        endforeach ()
endif ()

################################################################################
#                                                                              #
# Documentation                                                                #
#                                                                              #
################################################################################

message ("Configuring documentation")

find_program(SPHINX_EXECUTABLE
        NAMES sphinx-build
        DOC "Sphinx Documentation Builder (sphinx-doc.org)"
        )
mark_as_advanced(FORCE SPHINX_EXECUTABLE)

if (SPHINX_EXECUTABLE)
        message ("-- Sphinx compiler: ${SPHINX_EXECUTABLE}")

        # Get the short version (not using REGEX_REPLACE for backwards compatibility)
        string (REPLACE "v" "" tmp ${MOLCAS_VERSION})
        string (FIND "${tmp}" "-" pos)
        string (SUBSTRING "${tmp}" 0 ${pos} SHORT_VERSION)

        set (SPHINX_BUILDDIR "${PROJECT_BINARY_DIR}/doc")
        set (SPHINX_OPTIONS -d ${SPHINX_BUILDDIR}/doctrees -D extract_dir=${SPHINX_BUILDDIR}/samples)

        add_custom_target (doc_all)

        add_custom_target (
                doc_html
                # workaround for sphinx bug #5966
                COMMAND touch ${PROJECT_SOURCE_DIR}/doc/source/index.rst
                COMMAND OPENMOLCAS_RELEASE=${MOLCAS_VERSION} OPENMOLCAS_VERSION=${SHORT_VERSION}
                        ${SPHINX_EXECUTABLE} -b html ${SPHINX_OPTIONS} ${PROJECT_SOURCE_DIR}/doc/source ${SPHINX_BUILDDIR}/html
                )
        add_dependencies (doc_all doc_html)

        find_program (LATEXMK "latexmk")
        mark_as_advanced(FORCE LATEXMK)
        if (LATEXMK)
                message ("-- LaTeX compiler: ${LATEXMK}")

                add_custom_target (
                        doc_pdf
                        COMMAND OPENMOLCAS_RELEASE=${MOLCAS_VERSION} OPENMOLCAS_VERSION=${SHORT_VERSION}
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

install (TARGETS ${PROGRAM_EXECUTABLES} ${MOLCAS_LIBRARIES}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        )

install (PROGRAMS
        ${PROJECT_BINARY_DIR}/sbin/help_basis
        ${PROJECT_BINARY_DIR}/sbin/help_doc
        DESTINATION ${CMAKE_INSTALL_PREFIX}/sbin
        )

if (DEFINED PYMOLCAS_TARGET)
        install (PROGRAMS
                ${PYMOLCAS_TARGET}
                DESTINATION ${CMAKE_INSTALL_PREFIX}/sbin
                RENAME pymolcas
                )

        install (CODE "execute_process(COMMAND ${OPENMOLCAS_DIR}/sbin/install_pymolcas.sh ${PYMOLCAS_TARGET})")
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

        install (CODE "execute_process(COMMAND ${OPENMOLCAS_DIR}/sbin/install_driver.sh ${MOLCAS_DRIVER})")
endif ()

if (IS_DIRECTORY "$ENV{HOME}/.Molcas")
        install (CODE "execute_process(COMMAND sh \"-c\" \"cd ${CMAKE_INSTALL_PREFIX}; echo `pwd` > $HOME/.Molcas/molcas 2> /dev/null\")")
endif ()
