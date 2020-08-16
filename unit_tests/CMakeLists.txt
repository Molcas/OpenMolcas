set(UNIT_TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

function(add_unit_test Target Source)
    add_subdirectory(${PROJECT_SOURCE_DIR}/unit_tests/fruit fruit)

    add_executable(${Target} ${Source})

    target_link_libraries(${Target} fruit mods_obj molcas ${EXTERNAL_LIBRARIES})

    set_target_properties(${Target}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${UNIT_TEST_BUILD_DIR}/bin")
endfunction()


# And start adding tests
add_subdirectory(linalg_mod/)

