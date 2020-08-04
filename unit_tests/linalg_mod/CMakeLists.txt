include(${PROJECT_SOURCE_DIR}/cmake_helpers.cmake)

add_subdirectory(${PROJECT_SOURCE_DIR}/unit_tests/fruit fruit)

add_executable(test_linalg test_linalg.F90)

target_link_libraries(test_linalg fruit mods_obj molcas ${EXTERNAL_LIBRARIES})

set_target_properties(test_linalg
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${UNIT_TEST_BUILD_DIR}/bin"
)
