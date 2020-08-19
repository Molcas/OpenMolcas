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
