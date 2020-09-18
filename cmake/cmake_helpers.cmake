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
# Copyright (C) 2020, Oskar Weser                                      *
#***********************************************************************

function(add_Fortran_library Target)
# Has the same signature as add_library, but adds the
# Fortran_MODULE_DIRECTORY property in addition.
    add_library(${ARGV})
    get_target_property(LIB_DIR ${Target} BINARY_DIR)
    set_target_properties(${Target} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/mod)
    target_include_directories(${Target} INTERFACE ${LIB_DIR}/mod)
endfunction()

function(add_directory Dir)
    add_subdirectory(${Dir} bin)
endfunction()
