# macro (add_unit_test test_name)
#
#           neci_add_test(
#             TARGET test_neci_${test_name}
#         SOURCES test_${test_name}.F90
#         LIBS libneci fruit
#         LINKER_LANGUAGE Fortran
#         META_TARGET test_neci
#         DEFINITIONS
#                 "${${PROJECT_NAME}_neci_DEFINITIONS}"
#     )
# endmacro()

# Make some stuff available to all the tests
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Build FRUIT for unit testing
add_subdirectory( fruit_src )

# And start adding tests
add_subdirectory( linalg_mod )

