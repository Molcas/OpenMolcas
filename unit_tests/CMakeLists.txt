include(${PROJECT_SOURCE_DIR}/cmake/cmake_helpers.cmake)

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


# And start adding tests
add_subdirectory(linalg_mod/)

