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
# Copyright (C) 2016,2018,2020, Ignacio Fdez. Galván                   *
#               2020, Oskar Weser                                      *
#***********************************************************************

# Find the location of a non-empty directory called Name as a subdirectory
# of one of the elements of Source_Roots.
# The result is returned as ${Name}_src
function (find_source Name Source_Roots)
  foreach (src ${Source_Roots})
    if (EXISTS ${PROJECT_SOURCE_DIR}/${src}/${Name})
      file (GLOB tmplist RELATIVE ${PROJECT_SOURCE_DIR}/${src}/${Name} ${PROJECT_SOURCE_DIR}/${src}/${Name}/*)
      # ignore possible leftover files from previous configure+make
      list (REMOVE_ITEM tmplist Makefile 00dependencies 00sources)
      if (tmplist)
        set (found "${src}/${Name}")
        break ()
      endif ()
    endif ()
  endforeach ()
  # do not overwrite existing definitions
  if (DEFINED ${Name}_src)
    set (${Name}_src "${${Name}_src}" PARENT_SCOPE)
  else ()
    if (DEFINED found)
      set (${Name}_src "${found}" PARENT_SCOPE)
    endif ()
  endif ()
endfunction ()

# Insert element New before Item in List, if Item does
# not exist, append New at the end
function (insert_before List Item New)
  list (GET ${Item} 0 first)
  list (FIND ${List} ${first} index)
  if (index GREATER -1)
    list (INSERT ${List} ${index} ${New})
  else ()
    list (APPEND ${List} ${New})
  endif ()
  set (${List} ${${List}} PARENT_SCOPE)
endfunction ()

# Has the same signature as add_library, but adds the
# Fortran_MODULE_DIRECTORY property in addition.
function (add_Fortran_library Target)
  add_library (${ARGV})
  get_target_property (LIB_DIR ${Target} BINARY_DIR)
  set_target_properties (${Target} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/mod)
  target_include_directories (${Target} INTERFACE ${LIB_DIR}/mod)
endfunction ()

function (add_directory Dir)
  add_subdirectory (${Dir} bin)
endfunction ()

# Same as add_library, but sets the Fortran_MODULE_DIRECTORY too.
function (add_Molcas_library Target)
  add_library (${ARGV})
  set_module_directory (${Target})
endfunction ()

# Same as add_executable, but sets the Fortran_MODULE_DIRECTORY too.
function (add_Molcas_executable Target)
  add_executable (${ARGV})
  set_module_directory (${Target})
endfunction ()

# Helper function for add_Molcas_{library,executable}
function (set_module_directory Target)
  if (SINGLE_MOD_DIR)
    set (mod_dir ${MAIN_MOD_DIR})
  else ()
    set (mod_dir ${MAIN_MOD_DIR}/${Target})
  endif ()
  set_target_properties (${Target} PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
  target_include_directories (${Target} INTERFACE ${mod_dir})
endfunction ()

# Get the specified Properties from Target, and pass them to Files
# This may be useful if a file is moved to a different virtual target
function (pass_properties_to_files Target Files Properties)
  foreach (prop ${Properties})
    set_source_files_properties (${Files} PROPERTIES ${prop} $<TARGET_PROPERTY:${Target},${prop}>)
  endforeach ()
endfunction ()

# Get the specified Properties from Target_Source, and pass them to Target
function (pass_properties_to_target Target_Source Target Properties)
  foreach (prop ${Properties})
    set_target_properties (${Target} PROPERTIES ${prop} $<TARGET_PROPERTY:${Target_Source},${prop}>)
  endforeach ()
endfunction ()

# Remove specified source files from the Target
# This may be useful to move some files to a different virtual target
function (target_remove_sources Target)
  get_target_property (sources ${Target} SOURCES)
  list (REMOVE_ITEM sources ${ARGN})
  set_target_properties (${Target} PROPERTIES SOURCES "${sources}")
endfunction ()

# Convert a list of relative paths to absolute paths
function (set_absolute_paths Output Base_Dir)
  foreach (file ${ARGN})
    get_filename_component (file_path "${file}" REALPATH BASE_DIR "${Base_Dir}")
    if (NOT EXISTS ${file_path})
      message (FATAL_ERROR "${file_path} does not exist!")
    endif ()
    list (APPEND absolute_paths ${file_path})
  endforeach ()
  set (${Output} ${absolute_paths} PARENT_SCOPE)
endfunction ()
