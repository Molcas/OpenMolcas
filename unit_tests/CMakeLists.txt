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

set(all_tests "")

function(add_unit_test Target)
    add_directory(../fruit_molcas)

    add_executable(${Target} ${ARGN})

    target_link_libraries(${Target} fruit_molcas libmolcas ${EXTERNAL_LIBRARIES})

    set_target_properties(${Target}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
            "${UNIT_TEST_BUILD_DIR}/bin")

    add_test(
        NAME
            ${Target}
        COMMAND
            $<TARGET_FILE:${Target}>)

    list(APPEND all_tests ${Target})
    set(all_tests ${all_tests} PARENT_SCOPE)

    set_tests_properties(
        ${Target}
        PROPERTIES
            ENVIRONMENT "MOLCAS_MEM=256"
            FIXTURES_REQUIRED "test_fixture")
endfunction()

macro(add_test_directory SubDir)
    add_subdirectory(${SubDir}/)
    get_directory_property(all_tests DIRECTORY ${SubDir} DEFINITION all_tests)
endmacro()

# start adding tests, but don't bother with shitty compilers
if (NOT ((CMAKE_Fortran_COMPILER_ID IN_LIST not_supported) OR (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "14")))
    add_test_directory(linalg_mod/)
    add_test_directory(fortran_strings/)
    add_test_directory(sorting/)
    add_test_directory(filesystem/)

    add_custom_target(build_tests)
    add_dependencies(build_tests ${all_tests})

    add_test(test_build
        "${CMAKE_COMMAND}"
        --build "${CMAKE_BINARY_DIR}"
        --config "$<CONFIG>"
        --target build_tests
    )

    set_tests_properties(
        test_build
        PROPERTIES
            FIXTURES_REQUIRED "test_fixture")
endif()
