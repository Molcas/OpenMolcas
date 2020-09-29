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

include(${PROJECT_SOURCE_DIR}/cmake/cmake_helpers.cmake)

set(not_supported ";SunPro;PGI;")

if (NOT ((CMAKE_Fortran_COMPILER_ID IN_LIST not_supported) OR (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "14")))
    add_unit_test(test_fortran_strings test_fortran_strings.F90)
endif()


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_source_files_properties(test_fortran_strings.F90 COMPILE_FLAGS -Wno-maybe-uninitialized)
endif()
