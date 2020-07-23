
# Libfruit provides the machinery for unit testing in Fortran

add_library(fruit)

get_target_property(LIB_DIR fruit BINARY_DIR)

set_target_properties(fruit PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/mod)

target_include_directories(fruit INTERFACE ${LIB_DIR}/mod)
