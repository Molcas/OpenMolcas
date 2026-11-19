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
# Tinker for the QM/MM interface (ESPF/Gateway)
#
# Tinker is not linked, it is a separate program run through $TINKER
# (default $MOLCAS/tinker/bin). There are no cross-dependencies, so one
# patched Tinker can serve any number of OpenMolcas installations.
#
# This directory is always configured by the main CMakeLists.txt (it is
# NOT part of TOOLS):
#  * always: "molcas get_tinker" is made available (sbin/get_tinker)
#  * TINKER=ON: provide Tinker in ${PROJECT_BINARY_DIR}/tinker, either
#    - an existing patched Tinker (TINKER_DIR), which is only linked, or
#    - downloaded, patched and compiled with get_tinker (TINKER_DIR empty)
#  * TINKER=ON, make install: ${CMAKE_INSTALL_PREFIX}/tinker is provided
#    the same way and it is checked that the interface program works

set (tinker_script ${CMAKE_CURRENT_LIST_DIR}/get_tinker)
file (STRINGS ${tinker_script} tinker_version_line REGEX "^TINKERVERSION=")
string (REGEX REPLACE "^TINKERVERSION='(.*)'.*" "\\1" tinker_version "${tinker_version_line}")
set (tinker_patch ${CMAKE_CURRENT_LIST_DIR}/patch_${tinker_version}.diff)
if (NOT EXISTS ${tinker_patch})
  message (FATAL_ERROR "Tinker patch not found: ${tinker_patch}")
endif ()

# "molcas get_tinker" (the script finds the patch in ../data/tinker)
file (COPY ${tinker_script} DESTINATION ${PROJECT_BINARY_DIR}/sbin)
configure_file (${tinker_patch} ${PROJECT_BINARY_DIR}/data/tinker/patch_${tinker_version}.diff COPYONLY)
install (PROGRAMS ${PROJECT_BINARY_DIR}/sbin/get_tinker
         DESTINATION ${CMAKE_INSTALL_PREFIX}/sbin
)

set (TINKER_DIR "" CACHE PATH "Existing patched Tinker (directory containing bin/tkr2qm_s) to use instead of compiling it.")
set (TINKER_SOURCE "" CACHE STRING "Local Tinker tarball (or directory containing it, relative to the build directory) or URL (empty: default upstream URL).")
set (TINKER_MAXATM "" CACHE STRING "Maximum number of atoms in Tinker (empty: Tinker default, 100000).")
set (TINKER_MAKE_JOBS "4" CACHE STRING "Parallel make jobs when compiling Tinker.")
set (TINKER_MAKE_RETRIES "1" CACHE STRING "Attempts at compiling Tinker (works around intermittent compiler crashes).")
set (TINKER_Fortran_COMPILER "" CACHE FILEPATH "GNU-compatible Fortran compiler for Tinker (empty: autodetect).")
set (tinker_vars TINKER_DIR TINKER_SOURCE TINKER_MAXATM TINKER_MAKE_JOBS TINKER_MAKE_RETRIES TINKER_Fortran_COMPILER)

if (NOT TINKER)
  mark_as_advanced (FORCE ${tinker_vars})
  message ("Tinker DISABLED")
  return ()
endif ()

mark_as_advanced (CLEAR ${tinker_vars})
message ("Configuring Tinker:")
set (tinker_build_dir ${PROJECT_BINARY_DIR}/tinker)

if (TINKER_DIR)

  # Existing Tinker: check it, and link it into the build directory
  get_filename_component (tinker_dir "${TINKER_DIR}" ABSOLUTE BASE_DIR ${PROJECT_BINARY_DIR})
  if (NOT EXISTS ${tinker_dir}/bin/tkr2qm_s)
    message (FATAL_ERROR "TINKER_DIR does not contain a patched Tinker (no bin/tkr2qm_s): ${tinker_dir}")
  endif ()
  message ("-- Tinker: using ${tinker_dir}")
  if (IS_SYMLINK ${tinker_build_dir} OR NOT EXISTS ${tinker_build_dir})
    file (REMOVE ${tinker_build_dir})
    execute_process (COMMAND ${CMAKE_COMMAND} -E create_symlink ${tinker_dir} ${tinker_build_dir})
  else ()
    message (FATAL_ERROR "${tinker_build_dir} exists and is not a link, remove it to use TINKER_DIR.")
  endif ()

else ()

  # The Tinker Makefile only supports GNU Fortran flags
  set (tinker_fc "${TINKER_Fortran_COMPILER}")
  if (NOT tinker_fc)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set (tinker_fc "${CMAKE_Fortran_COMPILER}")
    else ()
      find_program (tinker_gfortran gfortran)
      mark_as_advanced (tinker_gfortran)
      if (NOT tinker_gfortran)
        message (FATAL_ERROR "Tinker needs GNU Fortran, set TINKER_Fortran_COMPILER.")
      endif ()
      set (tinker_fc "${tinker_gfortran}")
    endif ()
  endif ()

  # get_tinker makes tinker a link to ${tinker_version}, any other link
  # is left over from a previous TINKER_DIR configuration
  if (IS_SYMLINK ${tinker_build_dir})
    get_filename_component (tinker_target ${tinker_build_dir} REALPATH)
    get_filename_component (tinker_own ${PROJECT_BINARY_DIR}/${tinker_version} REALPATH)
    if (NOT tinker_target STREQUAL tinker_own)
      file (REMOVE ${tinker_build_dir})
    endif ()
  endif ()
  set (tinker_stamp ${PROJECT_BINARY_DIR}/External/tinker/tinker_args.txt)
  if (NOT EXISTS ${tinker_build_dir}/bin/tkr2qm_s)
    # force (re)compilation
    file (REMOVE ${tinker_stamp})
  endif ()

  set (tinker_args -y -p ${PROJECT_BINARY_DIR} -P ${tinker_patch} -c ${tinker_fc}
                   -j ${TINKER_MAKE_JOBS} -r ${TINKER_MAKE_RETRIES})
  if (TINKER_MAXATM)
    list (APPEND tinker_args -a ${TINKER_MAXATM})
  endif ()
  if (TINKER_SOURCE MATCHES "^[a-zA-Z]+://")
    list (APPEND tinker_args -u ${TINKER_SOURCE})
  elseif (TINKER_SOURCE)
    # relative paths are relative to the build directory
    get_filename_component (tinker_tarball "${TINKER_SOURCE}" ABSOLUTE BASE_DIR ${PROJECT_BINARY_DIR})
    if (NOT EXISTS ${tinker_tarball})
      message (FATAL_ERROR "TINKER_SOURCE not found: ${tinker_tarball}")
    endif ()
    list (APPEND tinker_args -t ${tinker_tarball})
  endif ()
  message ("-- Tinker version: ${tinker_version}")
  message ("-- Tinker compiler: ${tinker_fc}")

  # Rebuild when the script, the patch or the options change
  file (WRITE ${tinker_stamp}.tmp "${tinker_args}\n")
  configure_file (${tinker_stamp}.tmp ${tinker_stamp} COPYONLY)

  include (ExternalProject)
  set (CMAKE_DISABLE_SOURCE_CHANGES OFF)
  ExternalProject_Add (tinker
                       PREFIX ${PROJECT_BINARY_DIR}/External/tinker
                       DOWNLOAD_COMMAND ""
                       CONFIGURE_COMMAND ""
                       BUILD_COMMAND sh ${tinker_script} ${tinker_args}
                       BUILD_IN_SOURCE ON
                       INSTALL_COMMAND ""
                       BUILD_BYPRODUCTS ${tinker_build_dir}/bin/tkr2qm_s
  )
  ExternalProject_Add_StepDependencies (tinker build ${tinker_script} ${tinker_patch} ${tinker_stamp})
  set (CMAKE_DISABLE_SOURCE_CHANGES ON)

endif ()

# Installation: an existing Tinker is linked, a compiled one is copied
# (only bin and params are needed at run time)
set (tinker_install_dir ${CMAKE_INSTALL_PREFIX}/tinker)
if (TINKER_DIR)
  install (CODE "
    if (IS_SYMLINK \"\$ENV{DESTDIR}${tinker_install_dir}\")
      file (REMOVE \"\$ENV{DESTDIR}${tinker_install_dir}\")
    endif ()
    if (EXISTS \"\$ENV{DESTDIR}${tinker_install_dir}\")
      message (FATAL_ERROR \"${tinker_install_dir} exists and is not a link, cannot link it to ${tinker_dir}\")
    endif ()
    message (STATUS \"Linking: \$ENV{DESTDIR}${tinker_install_dir} -> ${tinker_dir}\")
    execute_process (COMMAND \"${CMAKE_COMMAND}\" -E create_symlink \"${tinker_dir}\" \"\$ENV{DESTDIR}${tinker_install_dir}\")
  ")
else ()
  install (CODE "
    if (IS_SYMLINK \"\$ENV{DESTDIR}${tinker_install_dir}\")
      file (REMOVE \"\$ENV{DESTDIR}${tinker_install_dir}\")
    endif ()
  ")
  install (DIRECTORY
    ${tinker_build_dir}/bin
    ${tinker_build_dir}/params
    DESTINATION ${tinker_install_dir}
    USE_SOURCE_PERMISSIONS
  )
endif ()

# Make sure the installed OpenMolcas can run Tinker: $MOLCAS/tinker/bin
# is where ESPF/Gateway look for it when $TINKER is not set
install (CODE "
  set (tkr \"\$ENV{DESTDIR}${tinker_install_dir}/bin/tkr2qm_s\")
  if (NOT EXISTS \"\${tkr}\" OR IS_DIRECTORY \"\${tkr}\")
    message (FATAL_ERROR \"Tinker interface program not installed: \${tkr}\")
  endif ()
  execute_process (COMMAND test -x \"\${tkr}\" RESULT_VARIABLE rc)
  if (rc)
    message (FATAL_ERROR \"Tinker interface program not executable: \${tkr}\")
  endif ()
  if (NOT IS_DIRECTORY \"\$ENV{DESTDIR}${tinker_install_dir}/params\")
    message (FATAL_ERROR \"Tinker parameter files not installed: \$ENV{DESTDIR}${tinker_install_dir}/params\")
  endif ()
  message (STATUS \"Tinker interface OK: \${tkr}\")
")
