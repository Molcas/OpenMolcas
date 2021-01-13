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

include_guard(GLOBAL)

include(${PROJECT_SOURCE_DIR}/cmake/cmake_helpers.cmake)

# unsupported compilers
set(not_supported ";SunPro;PGI;")

set(UNIT_TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

function(add_unit_test Target)
    add_directory(../fruit_molcas)

    add_executable(${Target} ${ARGN})

    target_link_libraries(${Target} fruit_molcas mods_obj molcas ${EXTERNAL_LIBRARIES})

    set_target_properties(${Target}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${UNIT_TEST_BUILD_DIR}/bin")

    add_test(
        NAME
            ${Target}
        COMMAND
            $<TARGET_FILE:${Target}>)

    set_tests_properties(
        ${Target}
        PROPERTIES
            ENVIRONMENT "MOLCAS_MEM=1000")
endfunction()


# start adding tests, but don't bother with shitty compilers
if (NOT ((CMAKE_Fortran_COMPILER_ID IN_LIST not_supported) OR (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "14")))
    add_subdirectory(linalg_mod/)
    add_subdirectory(fortran_strings/)
    add_subdirectory(sorting/)
endif()
