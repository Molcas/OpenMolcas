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
# Copyright (C) 2017, Stefan Knecht                                    *
#               2020, Leon Freitag                                     *
#***********************************************************************
#                                                                      *
#***********************************************************************
# CMakeLists.txt for QCMaquis                                          *
#***********************************************************************

project (qcmaquis_ext)
# load External Project macro
include (ExternalProject)

set (extprojpath "${CMAKE_BINARY_DIR}/External/qcmaquis")

if (NOT QCMaquis_ROOT)

  if (NOT LOCAL_QCM_INSTALL_PATH)
    set (LOCAL_QCM_INSTALL_PATH ${CMAKE_BINARY_DIR}/qcmaquis)
  endif ()
  # Set up compilation of QCMaquis components

  # QCMaquis does not know profile
  if (CMAKE_BUILD_TYPE MATCHES "profile")
    set (QCM_BUILD_TYPE "release")
  else ()
    set (QCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
  endif ()
  list (APPEND CMAKE_MODULE_PATH ${CMAKE_ROOT})

#else ()
#  set (MAQUIS_DMRG_DIR ${QCMaquis_ROOT}/share/cmake)
#  find_package (MAQUIS_DMRG PATHS ${QCMaquis_ROOT} NO_DEFAULT_PATH)
endif ()

# save openmolcas Tools subdirectory
set (OPENMOLCAS_TOOLS_DIR ${CMAKE_BINARY_DIR}/Tools/qcmaquis)
if (SINGLE_MOD_DIR)
  set (mod_dir ${MAIN_MOD_DIR}/_single)
else ()
  set (mod_dir ${MAIN_MOD_DIR}) # making a subdirectory for qcmaquis dirs doesn't work
endif ()
set (CMAKE_DISABLE_SOURCE_CHANGES ON)

if (MPI)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D_MOLCAS_MPP_")

  # Add GA include path to the Fortran script, required to find "mafdecls.fh"
  if (GA_INCLUDE_PATH)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${GA_INCLUDE_PATH}")
  else ()
    message (FATAL_ERROR "Could not find GA include path")
  endif ()

  # Workaround for gfortran 10
  if (GA AND "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif ()
endif ()

list (APPEND QCMaquisCMakeArgs
      -DCMAKE_BUILD_TYPE=${QCM_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/External
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS=${CMake_C_FLAGS}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS=${QCM_CMake_CXX_FLAGS}
      -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
)
if (HDF5_ROOT)
  list (APPEND QCMaquisCMakeArgs
        -DHDF5_ROOT=${HDF5_ROOT}
        -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
  )
endif ()

if (NOT MAQUIS_DMRG_FOUND) # Does the opposite work?

  # Pass BOOST_ROOT on if needed
  if (BOOST_ROOT)
    list (APPEND QCMaquisCMakeArgs
          -DBOOST_ROOT=${BOOST_ROOT}
    )
  endif (BOOST_ROOT)

  if (LINALG STREQUAL "Manual")
    target_files (LINALG_LIBRARIES_FILES ${LINALG_LIBRARIES})
    list (APPEND LINALG_LIBRARIES_FILES ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
    string (REPLACE ";" '\' LINALG_LIBRARIES_FILES "${LINALG_LIBRARIES_FILES}")
    list (APPEND QCMaquisCMakeArgs
          "-DBLAS_LAPACK_SELECTOR=manual"
          "-DMAQUISLapack_LIBRARIES=${LINALG_LIBRARIES_FILES}"
    )
  elseif (LINALG STREQUAL "MKL")
    list (APPEND QCMaquisCMakeArgs
          "-DBLAS_LAPACK_SELECTOR=mkl_sequential"
    )
  elseif (LINALG STREQUAL "OpenBLAS")
    list (APPEND QCMaquisCMakeArgs
          "-DBLAS_LAPACK_SELECTOR=openblas"
          "-DOPENBLASROOT=${OPENBLASROOT}"
    )
  elseif (LINALG STREQUAL "Accelerate")
    list (APPEND QCMaquisCMakeArgs
          "-DBLAS_LAPACK_SELECTOR:STRING=veclib"
    )
  elseif (LINALG STREQUAL "Internal")

    # To link QCMaquis with Fortran static libraries, we
    # need to add -lgfortran for gfortran
    # It seems that ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
    # is not suited for this because it contains also other unnecessary libraries

    # for some reason, the list does not work if the generator expression -lgfortran is not first
    # but for correct linking it needs to be last AND with a prepended "-l"
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set (Fortran_RUNTIME_LIBRARY "gfortran")
    endif ()

    list (APPEND QCMaquisCMakeArgs
          "-DBLAS_LAPACK_SELECTOR=manual"
          "-DMAQUISLapack_LIBRARIES=$<$<BOOL:Fortran_RUNTIME_LIBRARY>:${Fortran_RUNTIME_LIBRARY}\ >$<TARGET_FILE:blas>\ $<TARGET_FILE:lapack>\ $<TARGET_FILE:blas>\ -l$<$<BOOL:Fortran_RUNTIME_LIBRARY>:${Fortran_RUNTIME_LIBRARY}>"
    )
  else ()
    message (FATAL_ERROR "LINALG=${LINALG} is not supported by QCMaquis")
  endif ()

  # Enabling source changes to keep ExternalProject happy
  set (CMAKE_DISABLE_SOURCE_CHANGES OFF)
  ########
  ### New QCMaquis installation
  ########
  # HDF5 #
  ########
  message ("Check for cmake-configured HDF5")
  find_package (HDF5 NAMES hdf5)
  if (HDF5_FOUND)
    include_directories (${HDF5_INCLUDE_DIR})
    set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
  else (HDF5_FOUND)
    message ("Check for non-cmake-configured HDF5")
    find_package (HDF5) # Find non-cmake built HDF5
    if (HDF5_FOUND)
      include_directories (${HDF5_INCLUDE_DIR})
      set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
    else (HDF5_FOUND)
      set (LINK_LIBS "")
    endif (HDF5_FOUND)
  endif (HDF5_FOUND)

  if (NOT HDF5_FOUND)
    message (FATAL_ERROR "Could not find HDF5 support which is required for QCMaquis")
  endif (NOT HDF5_FOUND)

  ##########################
  # GNU scientific library #
  ##########################
  find_package (GSL REQUIRED)
  if (GSL_FOUND)
    list (APPEND MAQUIS_DMRG_LIBRARIES ${GSL_LIBRARIES})
  endif (GSL_FOUND)

  ##########
  # OpenMP #
  ##########
  if (OPENMP)
    find_package (OpenMP)
    if (NOT OPENMP_FOUND)
      message (FATAL_ERROR "Could not configure OpenMP.")
    endif (NOT OPENMP_FOUND)
  endif (OPENMP)

  ################################################
  # set CXX FLAGS for ALPS/BOOST and QCMaquis    #
  ################################################
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-attributes -Wno-deprecated-declarations")

  # fix for Intel compiler
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    if (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
      set (CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    endif ()
  endif ()

  # default
  set (QCMaquis_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

  # OpenMP flags
  if (OPENMP_FOUND)
    set (QCMaquis_CXX_FLAGS "${QCMaquis_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (QCMaquis_OPENMP "-DENABLE_OMP:BOOL=ON")
  else ()
    set (QCMaquis_OPENMP "-DENABLE_OMP:BOOL=OFF")
  endif ()
  unset (OPENMP_FOUND)

  list (APPEND CMAKE_MODULE_PATH ${CMAKE_ROOT})
  #list (APPEND CMAKE_MODULE_PATH ${extprojpath}/scripts/common/cmake)

  set (EP_PROJECT "qcmaquis")

  if (ADDRMODE EQUAL 64)
    set (EP_CMAKE_ARGS "${QCMaquisCMakeArgs}" "-DLAPACK_64_BIT:BOOL=ON")
  endif ()

  set (EP_CMAKE_CACHE_ARGS
       "-DBUILD_SYMMETRIES:STRING=TwoU1PG;SU2U1PG"
       "${QCMaquis_OPENMP}"
       "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
       "-DDMRG_NUMSYMM:STRING=6"
       "-DBUILD_DMRG:BOOL=ON"
       "-DBUILD_MPS_TRANSFORM:BOOL=ON"
       "-DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}"
       "-DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}"
       "-DCMAKE_CXX_FLAGS:STRING=${QCMaquis_CXX_FLAGS}"
       "-DCMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT:BOOL=OFF"
       "-DBUILD_OPENMOLCAS_INTERFACE:BOOL=ON"
  )

  if (LINALG STREQUAL "MKL")
    set (EP_CMAKE_CACHE_ARGS ${EP_CMAKE_CACHE_ARGS} "-DMKLROOT:STRING=${MKLROOT}")
  endif ()

  if (MPI AND GA)
    target_files (GA_LIBRARIES_FILES ${GA_LIBRARIES})
    set (EP_CMAKE_CACHE_ARGS
         ${EP_CMAKE_CACHE_ARGS}
         "-DBUILD_OPENMOLCAS_MPI:BOOL=ON"
         "-DGA_INCLUDE_DIR:STRING=${GA_INCLUDE_PATH}"
         "-DGA_LIBRARIES:STRING=${GA_LIBRARIES_FILES}"
    )
  endif ()

  # Boost for QCMaquis is required already here
  set (Boost_requirements program_options filesystem system serialization thread)
  set (Boost_NO_BOOST_CMAKE ON)

  find_package (Boost 1.56 COMPONENTS ${Boost_requirements})
  if (Boost_FOUND)
    list (APPEND MAQUIS_DMRG_LIBRARIES ${Boost_LIBRARIES})
  else (Boost_FOUND)
    message (FATAL_ERROR "Boost >= 1.56 is required for QCMaquis")
  endif (Boost_FOUND)

  ###############################
  # git references for QCMaquis #
  ###############################

  set (reference_git_repo https://github.com/qcscine/qcmaquis.git)
  set (reference_git_commit release-3.1.4)

  set (last_hash "None")
  set (hash_file ${extprojpath}/${EP_PROJECT}.hash)
  if (EXISTS ${hash_file})
    file (READ ${hash_file} last_hash)
    string (REGEX REPLACE "\n$" "" last_hash "${last_hash}")
  endif ()
  if (last_hash STREQUAL ${reference_git_commit})
    set (EP_SkipUpdate ON)
  else ()
    set (EP_SkipUpdate OFF)
  endif ()

  ExternalProject_Add (${EP_PROJECT}
                       PREFIX ${extprojpath}
                       GIT_REPOSITORY ${reference_git_repo}
                       GIT_TAG ${reference_git_commit}
                       UPDATE_DISCONNECTED ${EP_SkipUpdate}

                       SOURCE_SUBDIR dmrg
                       CMAKE_ARGS ${EP_CMAKE_ARGS}
                       CMAKE_CACHE_ARGS ${EP_CMAKE_CACHE_ARGS}
                       LIST_SEPARATOR '\'
                       INSTALL_DIR ${LOCAL_QCM_INSTALL_PATH}
                       UPDATE_COMMAND ""
  )

  ExternalProject_Add_Step (${EP_PROJECT} update_hash
                            COMMAND echo ${reference_git_commit} > ${hash_file}
                            DEPENDEES build
  )

  # Retrieve information about linking to the shared library
  # Unfortunately with external project we need to hard-code the library paths and all the rpath
  # because at the configure time this information is unknown
  # NOTE: the library paths are platform-specific
  # If we use external QCMaquis installation, all the paths are handled by find_package correctly

  #add_library (alps SHARED IMPORTED GLOBAL)
  add_library (maquis_dmrg SHARED IMPORTED GLOBAL)
  add_library (qcmaquis-driver SHARED IMPORTED GLOBAL)
  add_library (qcmaquis-hdf5-interface SHARED IMPORTED GLOBAL)

  add_dependencies (maquis_dmrg qcmaquis)
  add_dependencies (qcmaquis-driver qcmaquis)
  add_dependencies (qcmaquis-hdf5-interface qcmaquis)

  if (APPLE)
    set_target_properties (maquis_dmrg PROPERTIES
                           IMPORTED_LOCATION "${LOCAL_QCM_INSTALL_PATH}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}maquis_dmrg.dylib"
                           IMPORTED_SONAME "${CMAKE_FIND_LIBRARY_PREFIXES}maquis_dmrg.dylib"
    )
    set (ALPS_LIBRARY "${LOCAL_QCM_INSTALL_PATH}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}alps.dylib")
  else ()
    set_target_properties (maquis_dmrg PROPERTIES
                           IMPORTED_LOCATION "${LOCAL_QCM_INSTALL_PATH}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}maquis_dmrg.so"
                           IMPORTED_SONAME "${CMAKE_FIND_LIBRARY_PREFIXES}maquis_dmrg.so"
    )
    set (ALPS_LIBRARY "${LOCAL_QCM_INSTALL_PATH}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}alps.so")
  endif ()

  set_target_properties (qcmaquis-driver PROPERTIES
                         IMPORTED_LOCATION "${LOCAL_QCM_INSTALL_PATH}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}qcmaquis-driver.a"
                         IMPORTED_SONAME "${CMAKE_FIND_LIBRARY_PREFIXES}qcmaquis-driver.a"
  )

  set_target_properties (qcmaquis-hdf5-interface PROPERTIES
                         IMPORTED_LOCATION "${LOCAL_QCM_INSTALL_PATH}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}qcmaquis-hdf5-interface.a"
                         IMPORTED_SONAME "${CMAKE_FIND_LIBRARY_PREFIXES}qcmaquis-hdf5-interface.a"
  )

  ExternalProject_Get_Property (${EP_PROJECT} install_dir)
  ExternalProject_Add_Step (${EP_PROJECT} install_tools
                            DEPENDEES install
                            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${install_dir}/bin/qcm_checkpoint_rename.py ${OPENMOLCAS_TOOLS_DIR}/qcm_checkpoint_rename.py)

  set (QCM_MOD_SUBDIR share/qcmaquis/fortran_mod)
  file (MAKE_DIRECTORY ${mod_dir})
  # copy QCMaquis interface module files to a place where OpenMolcas can find them
  ExternalProject_Add_Step (${EP_PROJECT} install_modules
                            DEPENDEES install
                            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface_mpssi.mod
                                                                          ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface_utility_routines.mod
                                                                          ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface.mod
                                                                          ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface_cfg.mod
                                                                          ${install_dir}/${QCM_MOD_SUBDIR}/hdf5_utils.mod
                                                                          ${mod_dir}
  )

  # Fix shebang line
  add_custom_command (TARGET qcmaquis POST_BUILD
                      COMMAND sed -i '1 s@ python$$@ ${Python_EXECUTABLE}@' ${OPENMOLCAS_TOOLS_DIR}/qcm_checkpoint_rename.py
  )

  # set MAQUIS_DMRG_DIR so that future find_package calls can find QCMaquis
  #if (NOT MAQUIS_DMRG_DIR)
  #  set (MAQUIS_DMRG_DIR ${install_dir}/share/cmake)
  #endif ()
else ()
  set (install_dir ${MAQUIS_DMRG_DIR}/../../)
  set (QCM_MOD_SUBDIR share/qcmaquis/fortran_mod)
  # copy the above files manually
  file (COPY ${install_dir}/bin/qcm_checkpoint_rename.py DESTINATION ${OPENMOLCAS_TOOLS_DIR})
  file (COPY ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface_mpssi.mod
             ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface_utility_routines.mod
             ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface.mod
             ${install_dir}/${QCM_MOD_SUBDIR}/qcmaquis_interface_cfg.mod
             ${install_dir}/${QCM_MOD_SUBDIR}/hdf5_utils.mod
        DESTINATION ${mod_dir}
  )
endif (NOT MAQUIS_DMRG_FOUND)

set (CMAKE_DISABLE_SOURCE_CHANGES ON)

set (DMRG_INCLUDE ${mod_dir} PARENT_SCOPE)

# kszenes: only include gsl library and not other imported targets from GSL
# normally, it also includes cblas and blas which would overwrite the OpenMolcas BLAS
get_target_property(GSL_LOCATION GSL::gsl IMPORTED_LOCATION)

# set library paths
if (MAQUIS_DMRG_FOUND)
  set (MAQUIS_DMRG_LIBRARIES qcmaquis-driver ${Boost_LIBRARIES} ${GSL_LOCATION} PARENT_SCOPE)
else ()
  # add static QCMaquis libraries
  set (MAQUIS_DMRG_LIBRARIES
       qcmaquis-driver
       maquis_dmrg
       ${ALPS_LIBRARY}
       ${CMAKE_BINARY_DIR}/qcmaquis/lib/${CMAKE_FIND_LIBRARY_PREFIXES}dmrg_models.a
       ${CMAKE_BINARY_DIR}/qcmaquis/lib/${CMAKE_FIND_LIBRARY_PREFIXES}dmrg_utils.a
       ${Boost_LIBRARIES}
       ${GSL_LOCATION}
       PARENT_SCOPE
  )
endif ()

# add HDF5 QCMaquis interface libraries
set (HDF5_QCM_INCLUDE ${mod_dir} PARENT_SCOPE)
set (HDF5_QCM_LIBRARIES qcmaquis-hdf5-interface PARENT_SCOPE)
